/*
 * Copyright 2016-2018 Leo McCormack
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

/**
 * @file saf_sh.c
 * @ingroup SH
 * @brief Public source for the Spherical Harmonic Transform and Spherical Array
 *        Processing module (#SAF_SH_MODULE)
 *
 * A collection of spherical harmonic related functions. Many of which have been
 * derived from the MATLAB libraries found in [1-3].
 *
 * @see [1] https://github.com/polarch/Spherical-Harmonic-Transform
 *          Copyright (c) 2015, Archontis Politis, BSD-3-Clause License
 * @see [2] https://github.com/polarch/Array-Response-Simulator
 *          Copyright (c) 2015, Archontis Politis, BSD-3-Clause License
 * @see [3] https://github.com/polarch/Spherical-Array-Processing
 *          Copyright (c) 2016, Archontis Politis, BSD-3-Clause License
 *
 * @author Leo McCormack
 * @date 22.05.2016
 * @license ISC
 */

#include "saf_sh.h"
#include "saf_sh_internal.h"

/** First-order ACN/N3D to FuMa [without sqrt(4pi) term] conversion matrix */
const float wxyzCoeffs[4][4] = {
    {3.544907701811032f, 0.0f, 0.0f, 0.0f},
    {0.0f, 0.0f, 0.0f, 2.046653415892977f},
    {0.0f, 2.046653415892977f, 0.0f, 0.0f},
    {0.0f, 0.0f, 2.046653415892977f, 0.0f} };


/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

void unnorm_legendreP
(
    int n,
    double* x,
    int lenX,
    double* y /* FLAT: (n+1) x lenX  */
)
{
    int i, m;
    double s, norm, scale;
    double* P, *s_n, *tc, *sqrt_n;
    
    if(n==0){
        for(i=0; i<lenX; i++)
            y[i] = 1.0;
        return;
    }
    
    /* alloc */
    P = calloc1d((n+3)*lenX,sizeof(double));
    s_n = malloc1d(lenX*sizeof(double));
    tc = malloc1d(lenX*sizeof(double));
    sqrt_n = malloc1d((2*n+1)*sizeof(double));
    
    /* init */
    for(i=0; i<lenX; i++){
        s = sqrt(1.0-pow(x[i],2.0)) + 2.23e-20;
        s_n[i] = pow(-s, (double)n);
        tc[i] = -2.0 * x[i]/s;
    }
    for(i=0; i<2*n+1; i++)
        sqrt_n[i] = sqrt((double)i);
    norm = 1.0;
    for(i=1; i<=n; i++)
        norm *= 1.0 - 1.0/(2.0*(double)i);
    
    /* Starting values for downwards recursion */
    for(i=0; i<lenX; i++){
        P[(n)*lenX+i] = sqrt(norm)*s_n[i];
        P[(n-1)*lenX+i] = P[(n)*lenX+i] * tc[i] * (double)n/sqrt_n[2*n];
    }
    
    /* 3-step downwards recursion to m == 0 */
    for(m=n-2; m>=0; m--)
        for(i=0; i<lenX; i++)
            P[(m)*lenX+i] = (P[(m+1)*lenX+i]*tc[i]*((double)m+1.0) - P[(m+2)*lenX+i] * sqrt_n[n+m+2] * sqrt_n[n-m-1])/(sqrt_n[n+m+1]*sqrt_n[n-m]);
    
    /* keep up to the last 3 elements in P */
    for(i=0; i<n+1; i++)
        memcpy(&(y[i*lenX]), &(P[i*lenX]), lenX*sizeof(double));
    
    /* Account for polarity when x == -/+1 for first value of P */
    for(i=0; i<lenX; i++)
        if(sqrt(1.0-pow(x[i],2.0))==0)
            y[i] = pow(x[i],(double)n);
    
    /* scale each row by: sqrt((n+m)!/(n-m)!) */
    for(m=1; m<n; m++){
        scale = 1.0;
        for(i=n-m+1; i<n+m+1; i++)
            scale*=sqrt_n[i];
        for(i=0; i<lenX; i++)
            y[m*lenX+i] *= scale;
    }
    scale = 1.0;
    for(i=1; i<2*n+1; i++)
        scale*=sqrt_n[i];
    for(i=0; i<lenX; i++)
        y[n*lenX+i] *= scale;
    
    free(P);
    free(s_n);
    free(tc);
    free(sqrt_n);
}
 
void unnorm_legendreP_recur
(
    int n,
    float* x,
    int lenX,
    float* Pnm_minus1,
    float* Pnm_minus2,
    float* Pnm
)
{
    int i, m, k, kk;
    float x2, one_min_x2, dfact_k;
    
    if(n==0){
        for(i=0; i<lenX; i++)
            Pnm[i] = 1.0f;
        return;
    }
    
    for(i=0; i<lenX; i++){
        x2 = (x[i])*(x[i]);
        switch(n) {
            case 1:
                Pnm[0*lenX+i] = x[i];
                Pnm[1*lenX+i] = sqrtf(1.0f-x2);
                break;
            case 2:
                Pnm[0*lenX+i] = (3.0f*x2-1.0f)/2.0f;
                Pnm[1*lenX+i] = (x[i])*3.0f*sqrtf(1.0f-x2);
                Pnm[2*lenX+i] = 3.0f*(1.0f-x2);
                break;
            default:
                one_min_x2 = 1.0f-x2;
                /* last term m=n */
                k = 2*n-1;
                dfact_k = 1.0f;
                if ((k % 2) == 0)
                    for (kk=1; kk<k/2+1; kk++)
                        dfact_k *= 2.0f*(float)kk;
                else
                    for (kk=1; kk<(k+1)/2+1; kk++)
                        dfact_k *= (2.0f*(float)kk-1.0f);
                
                Pnm[n*lenX+i] = dfact_k * powf(one_min_x2, (float)n/2.0f);
                /* before last term */
                /* P_{n(n-1)} = (2*n-1)*x*P_{(n-1)(n-1)} */
                Pnm[(n-1)*lenX+i] = (float)k * (x[i]) *Pnm_minus1[(n-1)*lenX+i];
                /* three term recurence for the rest */
                for (m=0; m<n-1; m++)
                    /* P_l = ( (2l-1)xP_(l-1) - (l+m-1)P_(l-2) )/(l-m) */
                    Pnm[m*lenX+i] = ( ((float)k * (x[i]) *Pnm_minus1[m*lenX+i]) - ((float)(n+m-1) * Pnm_minus2[m*lenX+i])) / (float)(n-m);
                break;
        }
    }
}


/* ========================================================================== */
/*                    SH and Beamforming related Functions                    */
/* ========================================================================== */

void getSHreal
(
    int order,
    float* dirs_rad,
    int nDirs,
    float* Y  /* the SH weights: (order+1)^2 x nDirs */
)
{
    int dir, j, n, m, idx_Y;
    double* Lnm, *CosSin;
    double *p_nm, *cos_incl;
    double *norm_real;

    if(nDirs<1)
        return;
    
    Lnm = malloc1d((2*order+1)*nDirs*sizeof(double));
    norm_real = malloc1d((2*order+1)*sizeof(double));
    CosSin = malloc1d((2*order+1)*sizeof(double));
    cos_incl = malloc1d(nDirs*sizeof(double));
    p_nm = malloc1d((order+1)*nDirs * sizeof(double));
    for (dir = 0; dir<nDirs; dir++)
        cos_incl[dir] = cos((double)dirs_rad[dir*2+1]);
    
    idx_Y = 0;
    for(n=0; n<=order; n++){
        /* vector of unnormalised associated Legendre functions of current order */
        unnorm_legendreP(n, cos_incl, nDirs, p_nm); /* includes Condon-Shortley phase term */
        
        for(dir=0; dir<nDirs; dir++){
            /* cancel the Condon-Shortley phase from the definition of the Legendre functions to result in signless real SH */
            if (n != 0)
                for(m=-n, j=0; m<=n; m++, j++)
                    Lnm[j*nDirs+dir] = pow(-1.0, (double)abs(m)) * p_nm[abs(m)*nDirs+dir];
            else
                Lnm[dir] = p_nm[dir];
        }
        
        /* normalisation */
        for(m=-n, j=0; m<=n; m++, j++)
            norm_real[j] = sqrt( (2.0*(double)n+1.0) * (double)factorial(n-abs(m)) / (4.0*SAF_PId*(double)factorial(n+abs(m))) );
        
        /* norm_real * Lnm_real .* CosSin; */
        for(dir=0; dir<nDirs; dir++){
            for(m=-n, j=0; m<=n; m++, j++){
                if(j<n)
                    Y[(j+idx_Y)*nDirs+dir] = (float)(norm_real[j] * Lnm[j*nDirs+dir] * sqrt(2.0)*sin((double)(n-j)*(double)dirs_rad[dir*2]));
                else if(j==n)
                    Y[(j+idx_Y)*nDirs+dir] = (float)(norm_real[j] * Lnm[j*nDirs+dir]);
                else /* (j>n) */
                    Y[(j+idx_Y)*nDirs+dir] = (float)(norm_real[j] * Lnm[j*nDirs+dir] * sqrt(2.0)*cos((double)(abs(m))*(double)dirs_rad[dir*2]));
            }
        }
        
        /* increment */
        idx_Y = idx_Y + (2*n+1);
    }
    
    free(p_nm);
    free(Lnm);
    free(norm_real);
    free(CosSin);
    free(cos_incl);
}

void getSHreal_recur
(
    int N,
    float* dirs_rad,
    int nDirs,
    float* Y
)
{
    int n, m, i, dir, index_n;
    float Nn0, Nnm;
    float sleg_n[8], sleg_n_1[8], sleg_n_2[8], scos_incl, sfactorials_n[15];
    float* leg_n, *leg_n_1, *leg_n_2, *cos_incl, *factorials_n;

    if(nDirs<1)
        return;

    if(N<=7 && nDirs == 1){
        /* Single direction optimisation for up to 7th order */
        leg_n = sleg_n;
        leg_n_1 = sleg_n_1;
        leg_n_2 = sleg_n_2;
        cos_incl = &scos_incl;
        factorials_n = sfactorials_n;
    }
    else{
        factorials_n = malloc1d((2*N+1)*sizeof(float));
        leg_n = malloc1d((N+1)*nDirs * sizeof(float));
        leg_n_1 = malloc1d((N+1)*nDirs * sizeof(float));
        leg_n_2 = malloc1d((N+1)*nDirs * sizeof(float));
        cos_incl = malloc1d(nDirs * sizeof(float));
    }
    index_n = 0;

    /* precompute factorials */
    for (i = 0; i < 2*N+1; i++)
        factorials_n[i] = (float)factorial(i); 

    /* cos(inclination) = sin(elevation) */
    for (dir = 0; dir<nDirs; dir++)
        cos_incl[dir] = cosf(dirs_rad[dir*2+1]);

    /* compute SHs with the recursive Legendre function */
    for (n = 0; n<N+1; n++) {
        if (n==0) {
            for (dir = 0; dir<nDirs; dir++)
                Y[n*nDirs+dir] = 1.0f/SQRT4PI;
            index_n = 1;
        }
        else {
            unnorm_legendreP_recur(n, cos_incl, nDirs, leg_n_1, leg_n_2, leg_n); /* does NOT include Condon-Shortley phase term */

            Nn0 = sqrtf(2.0f*(float)n+1.0f);
            for (dir = 0; dir<nDirs; dir++){
                for (m = 0; m<n+1; m++) {
                    if (m==0)
                        Y[(index_n+n)*nDirs+dir] = Nn0/SQRT4PI * leg_n[m*nDirs+dir];
                    else {
                        Nnm = Nn0* sqrtf( 2.0f * factorials_n[n-m]/factorials_n[n+m] );
                        Y[(index_n+n-m)*nDirs+dir] = Nnm/SQRT4PI * leg_n[m*nDirs+dir] * sinf((float)m * (dirs_rad[dir*2]));
                        Y[(index_n+n+m)*nDirs+dir] = Nnm/SQRT4PI * leg_n[m*nDirs+dir] * cosf((float)m * (dirs_rad[dir*2]));
                    }
                }
            }
            index_n += 2*n+1;
        }
        utility_svvcopy(leg_n_1, (N+1)*nDirs, leg_n_2);
        utility_svvcopy(leg_n,   (N+1)*nDirs, leg_n_1);
    }

    if(N>7 || nDirs > 1){
        free(factorials_n);
        free(leg_n);
        free(leg_n_1);
        free(leg_n_2);
        free(cos_incl);
    }
}

void getSHcomplex
(
    int order,
    float* dirs_rad,
    int nDirs,
    float_complex* Y
)
{
    int dir, j, n, m, idx_Y;
    double *norm_real;
    double *Lnm, *cos_incl;
    double_complex Ynm;
    
    Lnm = malloc1d((order+1)*nDirs*sizeof(double));
    norm_real = malloc1d((order+1)*sizeof(double));
    cos_incl = malloc1d(nDirs*sizeof(double));
    for (dir = 0; dir<nDirs; dir++)
        cos_incl[dir] = cos((double)dirs_rad[dir*2+1]);
    
    idx_Y = 0;
    for(n=0; n<=order; n++){
        /* vector of unnormalised associated Legendre functions of current order */
        unnorm_legendreP(n, cos_incl, nDirs, Lnm); /* includes Condon-Shortley phase term */
        
        /* normalisation */
        for(m=0; m<=n; m++)
            norm_real[m] = sqrt( (2.0*(double)n+1.0)*(double)factorial(n-m) / (4.0*SAF_PId*(double)factorial(n+m)) );
        
        /* norm_real .* Lnm_real .* CosSin; */
        for(dir=0; dir<nDirs; dir++){
            for(m=-n, j=0; m<=n; m++, j++){
                if(m<0){
                    Ynm = crmul(conj(crmul(cexp(cmplx(0.0, (double)abs(m)*(double)dirs_rad[dir*2])), norm_real[abs(m)] * Lnm[abs(m)*nDirs+dir])), pow(-1.0, (double)abs(m)));
                    Y[(j+idx_Y)*nDirs+dir] = cmplxf((float)creal(Ynm), (float)cimag(Ynm));
                }
                else {/* (m>=0) */
                    Ynm = crmul(cexp(cmplx(0.0, (double)abs(m)*(double)dirs_rad[dir*2])), norm_real[m] * Lnm[m*nDirs+dir]);
                    Y[(j+idx_Y)*nDirs+dir] = cmplxf((float)creal(Ynm), (float)cimag(Ynm));
                }
            }
        }
        
        /* increment */
        idx_Y = idx_Y + (2*n+1);
    }
    
    free(Lnm);
    free(norm_real);
    free(cos_incl);
}

void complex2realSHMtx
(
    int order,
    float_complex* T_c2r
)
{
    int n, m, q, p, idx, nSH;
    
    nSH = ORDER2NSH(order);
    memset(T_c2r, 0, nSH*nSH*sizeof(float_complex));
    T_c2r[0] = cmplxf(1.0f, 0.0f);
    if(order == 0)
        return;
    
    idx = 1;
    for(n=1, q = 1; n<=order; n++){
        idx += (2*n+1);
        for(m=-n, p=0; m<=n; m++, q++, p++){
            if(m<0){
                T_c2r[(q)*nSH+(q)] = cmplxf(0.0f, 1.0f/sqrtf(2.0f));
                T_c2r[(idx-p-1)*nSH+(q)] = cmplxf(1.0f/sqrtf(2.0f), 0.0f);
            }
            else if(m==0)
                T_c2r[(q)*nSH+(q)] = cmplxf(1.0f, 0.0f);
            else{
                T_c2r[(q)*nSH+(q)] = cmplxf(powf(-1.0f,(float)m)/sqrtf(2.0f), 0.0f);
                T_c2r[(idx-p-1)*nSH+(q)] = cmplxf(0.0f, -powf(-1.0f, (float)abs(m))/sqrtf(2.0f));
            } 
        }
    }
}

void real2complexSHMtx
(
    int order,
    float_complex* T_r2c
)
{
    int n, m, q, p, idx, nSH;
    
    nSH = ORDER2NSH(order);
    memset(T_r2c, 0, nSH*nSH*sizeof(float_complex));
    T_r2c[0] = cmplxf(1.0f, 0.0f);
    if(order == 0)
        return;
    
    idx = 1;
    for(n=1, q = 1; n<=order; n++){
        idx += (2*n+1);
        for(m=-n, p=0; m<=n; m++, q++, p++){
            if(m<0){
                T_r2c[(q)*nSH+(q)] = cmplxf(0.0f, -1.0f/sqrtf(2.0f));
                T_r2c[(idx-p-1)*nSH+(q)] = cmplxf(0.0f, powf(-1.0f, (float)abs(m))/sqrtf(2.0f)); //cmplxf(1.0f/sqrtf(2.0f), 0.0f);
            }
            else if(m==0)
                T_r2c[(q)*nSH+(q)] = cmplxf(1.0f, 0.0f);
            else{
                T_r2c[(q)*nSH+(q)] = cmplxf(powf(-1.0f,(float)m)/sqrtf(2.0f), 0.0f);
                T_r2c[(idx-p-1)*nSH+(q)] = cmplxf(1.0f/sqrtf(2.0f), 0.0f); //cmplxf(0.0f, -powf(-1.0f, (float)abs(m))/sqrtf(2.0f));
            }
        }
    }
}

void complex2realCoeffs
(
    int order,
    float_complex* C_N,
    int K,
    float* R_N
)
{
    int i, nSH;
    float_complex* T_c2r, *R_N_c;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH = ORDER2NSH(order);
    T_c2r = malloc1d(nSH*nSH*sizeof(float_complex));
    R_N_c = malloc1d(nSH*K*sizeof(float_complex));
    complex2realSHMtx(order, T_c2r);
    for(i=0; i<nSH*nSH; i++)
        T_c2r[i] = conjf(T_c2r[i]);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, K, nSH, &calpha,
                T_c2r, nSH,
                C_N, K, &cbeta,
                R_N_c, K);
    for(i=0; i<nSH*K; i++)
        R_N[i] = crealf(R_N_c[i]);
    
    free(T_c2r);
    free(R_N_c);
}

/* Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 * by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
void getSHrotMtxReal
(
    float Rxyz[3][3],
    float* RotMtx/*(L+1)^2 x (L+1)^2 */,
    int L
)
{
    int i, j, M, l, m, n, d, bandIdx, denom;
    float u, v, w;
    float R_1[3][3], _R_lm1[64*64], _R_l[64*64];
    float* R_lm1, *R_l;

    /* Prep */
    M = (L+1) * (L+1);
    if(L<=7){
        R_lm1 = _R_lm1;
        R_l = _R_l;
    }
    else{
        R_lm1 = malloc1d(M*M*sizeof(float));
        R_l = malloc1d(M*M*sizeof(float));
    }
    memset(RotMtx, 0, M*M*sizeof(float));
    
    /* zeroth-band (l=0) is invariant to rotation */
    RotMtx[0] = 1;
    
    /* the first band (l=1) is directly related to the rotation matrix */
    R_1[0][0] = Rxyz[1][1];
    R_1[0][1] = Rxyz[1][2];
    R_1[0][2] = Rxyz[1][0];
    R_1[1][0] = Rxyz[2][1];
    R_1[1][1] = Rxyz[2][2];
    R_1[1][2] = Rxyz[2][0];
    R_1[2][0] = Rxyz[0][1];
    R_1[2][1] = Rxyz[0][2];
    R_1[2][2] = Rxyz[0][0];
    for (i=1; i<4; i++){
        R_lm1[(i-1)*M+0] = R_1[i-1][0];
        R_lm1[(i-1)*M+1] = R_1[i-1][1];
        R_lm1[(i-1)*M+2] = R_1[i-1][2];
        for (j=1; j<4; j++)
            RotMtx[i*M+j] = R_1[i-1][j-1];
    }
    
    /* compute rotation matrix of each subsequent band recursively */
    bandIdx = 4;
    for(l = 2; l<=L; l++){
        for(i=0; i<2*l+1; i++)
            memset(R_l + i*M, 0, (2*l+1) * sizeof(float));
        for(m=-l; m<=l; m++){
            for(n=-l; n<=l; n++){
                /* compute u,v,w terms of Eq.8.1 (Table I) */
                d = m == 0 ? 1 : 0; /* the delta function d_m0 */
                denom = abs(n) == l ? (2*l)*(2*l-1) : (l*l-n*n);
                u = sqrtf( (float)((l*l-m*m)) /  (float)denom);
                v = sqrtf( (float)((1+d)*(l+abs(m)-1)*(l+abs(m))) /  (float)denom) * (float)(1-2*d)*0.5f;
                w = sqrtf( (float)((l-abs(m)-1)*(l-abs(m))) / (float)denom) * (float)(1-d)*(-0.5f);
                
                /* computes Eq.8.1 */
                if (u!=0)
                    u = u* getU(M,l,m,n,R_1,R_lm1);
                if (v!=0)
                    v = v* getV(M,l,m,n,R_1,R_lm1);
                if (w!=0)
                    w = w* getW(M,l,m,n,R_1,R_lm1);
                
                R_l[(m+l)*M+(n+l)] = u+v+w;
            }
        }
        
        for(i=0; i<2*l+1; i++)
            for(j=0; j<2*l+1; j++)
                RotMtx[(bandIdx + i)*M + (bandIdx + j)] = R_l[i*M+j];
        for(i=0; i<2*l+1; i++)
            memcpy(R_lm1+i*M, R_l + i*M, (2*l+1) * sizeof(float));
        bandIdx += 2*l+1;
    }

    /* clean-up */
    if(L>7){
        free(R_lm1);
        free(R_l);
    }
}

void computeVelCoeffsMtx
(
    int sectorOrder,
    float_complex* A_xyz
)
{
    int i, j, Nxyz, Ns, nC_xyz, nC_s;
    float x1, x3, z2, y1, y3;
    float* G_mtx;
    
    Ns = sectorOrder;
    Nxyz = Ns+1;
    nC_xyz = (Nxyz+1)*(Nxyz+1);
    nC_s = (Ns+1)*(Ns+1);
    x1 = sqrtf(2.0f*SAF_PI/3.0f);
    x3 = -x1;
    y1 = y3 = sqrtf(2.0f*SAF_PI/3.0f);
    z2 = sqrtf(4.0f*SAF_PI/3.0f);
    G_mtx = malloc1d(nC_s*4*nC_xyz*sizeof(float));
    gaunt_mtx(Ns, 1, Nxyz, G_mtx);
    
    for (i=0; i<nC_xyz; i++){
        for (j=0; j<nC_s; j++){
            A_xyz[i*nC_s*3 + j*3 + 0] = cmplxf(x1*G_mtx[j*4*nC_xyz + 1*nC_xyz + i] + x3*G_mtx[j*4*nC_xyz + 3*nC_xyz + i], 0.0f);
            A_xyz[i*nC_s*3 + j*3 + 1] = cmplxf(0.0f, y1*G_mtx[j*4*nC_xyz + 1*nC_xyz + i] + y3*G_mtx[j*4*nC_xyz + 3*nC_xyz + i]);
            A_xyz[i*nC_s*3 + j*3 + 2] = cmplxf(z2*G_mtx[j*4*nC_xyz + 2*nC_xyz + i], 0.0f);
        }
    }
     
    free(G_mtx);
}

float computeSectorCoeffsEP
(
    int orderSec,
    float_complex* A_xyz, /* FLAT: (sectorOrder+2)^2 x (sectorOrder+1)^2 x 3 */
    SECTOR_PATTERNS pattern,
    float* sec_dirs_deg,
    int nSecDirs,
    float* sectorCoeffs /* FLAT: (nSecDirs*4) x (orderSec+2)^2 */
)
{
    int i, j, ns, orderVel, nSH;
    float normSec, azi_sec, elev_sec, Q;
    float* b_n, *c_nm, *xyz_nm;
    
    if(orderSec==0){
        memcpy(sectorCoeffs, wxyzCoeffs, 16*sizeof(float)); /* ACN/N3D to WXYZ */
        normSec = 1.0f;
    }
    else{
        orderVel = orderSec+1;
        nSH = (orderSec+2)*(orderSec+2);
        b_n = malloc1d((orderSec+1)*sizeof(float));
        c_nm = calloc1d((orderVel+1)*(orderVel+1), sizeof(float)); /* pad with zeros */
        xyz_nm = malloc1d((orderVel+1)*(orderVel+1)*3*sizeof(float));
        switch(pattern){
            case SECTOR_PATTERN_PWD:
                beamWeightsHypercardioid2Spherical(orderSec, b_n);
                Q = (float)((orderSec+1) * (orderSec+1));
                break;
            case SECTOR_PATTERN_MAXRE:
                beamWeightsMaxEV(orderSec, b_n);
                cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 1, 1, orderSec+1, 1.0f,
                            b_n, 1,
                            b_n, 1, 0.0f,
                            &Q, 1);
                Q = 4.0f*SAF_PI/(Q);
                break;
            case SECTOR_PATTERN_CARDIOID:
                beamWeightsCardioid2Spherical(orderSec, b_n);
                Q = 2.0f*(float)orderSec + 1.0f;
                break;
        }
        normSec = Q/(float)nSecDirs; /* directivity factor / number of sectors */
        
        for(ns=0; ns<nSecDirs; ns++){
            /* rotate the pattern by rotating the coefficients */
            azi_sec = sec_dirs_deg[ns*2] * SAF_PI/180.0f;
            elev_sec = sec_dirs_deg[ns*2+1] * SAF_PI/180.0f; /* from elevation to inclination */
            rotateAxisCoeffsReal(orderSec, b_n, SAF_PI/2.0f-elev_sec, azi_sec, c_nm);
            beamWeightsVelocityPatternsReal(orderSec, b_n, azi_sec, elev_sec, A_xyz, xyz_nm);
     
            /* store coefficients */
            for(j=0; j<nSH; j++){
                sectorCoeffs[ns*4*nSH + 0*nSH +j] = sqrtf(normSec) * c_nm[j];
                for(i=0; i<3; i++)
                    sectorCoeffs[ns*4*nSH + (i+1)*nSH +j] = sqrtf(normSec) * xyz_nm[j*3+i];
            }
        }
        
        free(b_n);
        free(c_nm);
        free(xyz_nm);
    }
    return normSec;
}

float computeSectorCoeffsAP
(
    int orderSec,
    float_complex* A_xyz,
    SECTOR_PATTERNS pattern,
    float* sec_dirs_deg,
    int nSecDirs,
    float* sectorCoeffs
)
{
    int i, j, ns, orderVel, nSH;
    float normSec, azi_sec, elev_sec;
    float* b_n, *c_nm, *xyz_nm;
    
    if(orderSec==0){
        memcpy(sectorCoeffs, wxyzCoeffs, 16*sizeof(float)); /* ACN/N3D to WXYZ */
        normSec = 1.0f;
    }
    else{
        orderVel = orderSec+1;
        nSH = (orderSec+2)*(orderSec+2);
        b_n = malloc1d((orderSec+1)*sizeof(float));
        c_nm = calloc1d((orderVel+1)*(orderVel+1), sizeof(float)); /* pad with zeros */
        xyz_nm = malloc1d((orderVel+1)*(orderVel+1)*3*sizeof(float));
        switch(pattern){
            case SECTOR_PATTERN_PWD: beamWeightsHypercardioid2Spherical(orderSec, b_n); break;
            case SECTOR_PATTERN_MAXRE: beamWeightsMaxEV(orderSec, b_n); break;
            case SECTOR_PATTERN_CARDIOID: beamWeightsCardioid2Spherical(orderSec, b_n); break;
        }
        normSec = (float)(orderSec+1)/(float)nSecDirs;
        
        for(ns=0; ns<nSecDirs; ns++){
            /* rotate the pattern by rotating the coefficients */
            azi_sec = sec_dirs_deg[ns*2] * SAF_PI/180.0f;
            elev_sec = sec_dirs_deg[ns*2+1] * SAF_PI/180.0f;
            rotateAxisCoeffsReal(orderSec, b_n, SAF_PI/2.0f-elev_sec, azi_sec, c_nm);
            beamWeightsVelocityPatternsReal(orderSec, b_n, azi_sec, elev_sec, A_xyz, xyz_nm);
            
            /* store coefficients */
            for(j=0; j<nSH; j++){
                sectorCoeffs[ns*4*nSH + 0*nSH +j] = normSec * c_nm[j];
                for(i=0; i<3; i++)
                    sectorCoeffs[ns*4*nSH + (i+1)*nSH +j] = normSec * xyz_nm[j*3+i];
            }
        }
        
        free(b_n);
        free(c_nm);
        free(xyz_nm);
    }
    return normSec;
}

void beamWeightsCardioid2Spherical
(
    int N,
    float* b_n
)
{
    int n;
    
    /* The coefficients can be derived by the binomial expansion of the cardioid function */
    for(n=0; n<N+1; n++) {
        b_n[n] = sqrtf(4.0f*SAF_PI*(2.0f*(float)n+1.0f)) *
                 (float)factorial(N)* (float)factorial(N+1)/
                 ((float)factorial(N+n+1)*(float)factorial(N-n))/
                 ((float)N+1.0f);
    }
}

void beamWeightsHypercardioid2Spherical
(
    int N,
    float* b_n
)
{
    int n;
    float* c_n;
    float dirs_rad[2] = {0.0f, 0.0f};
    
    c_n = malloc1d((N+1)*(N+1)*sizeof(float));
    getSHreal(N, dirs_rad, 1, c_n);
    for(n=0; n<N+1; n++)
        b_n[n] = c_n[(n+1)*(n+1)-n-1] * 4.0f * SAF_PI/(powf((float)N+1.0f, 2.0f));
    
    free(c_n);
}

void beamWeightsMaxEV
(
    int N,
    float* b_n
)
{
    int n;
    float norm;
    double temp_i;
    double* temp_o;
    
    temp_o = malloc1d( (N+1)*sizeof(double));
    norm = 0.0f;
    for (n=0; n<=N; n++) {
        temp_i = cos(2.4068f/((double)N+1.51));
        unnorm_legendreP(n, &temp_i, 1, temp_o);
        b_n[n] = sqrtf((2.0f*(float)n+1.0f)/(4.0f*SAF_PI))*(float)temp_o[0];
        norm +=  sqrtf((2.0f*(float)n+1.0f)/(4.0f*SAF_PI))*b_n[n];
    }
    
    /* normalise to unity response on look-direction */
    for (n=0; n<=N; n++)
        b_n[n] /= norm;
    
    free(temp_o);
}

void beamWeightsVelocityPatternsReal
(
    int order,
    float* b_n,
    float azi_rad,
    float elev_rad,
    float_complex* A_xyz,
    float* velCoeffs
)
{
    int nSH;
    float_complex* velCoeffs_c;
    
    nSH = ORDER2NSH(order+1);
    velCoeffs_c = malloc1d(nSH*3*sizeof(float_complex));
    beamWeightsVelocityPatternsComplex(order, b_n, azi_rad, elev_rad, A_xyz, velCoeffs_c);
    complex2realCoeffs(order+1, velCoeffs_c, 3, velCoeffs);
      
    free(velCoeffs_c);
}

void beamWeightsVelocityPatternsComplex
(
    int order,
    float* b_n,
    float azi_rad,
    float elev_rad,
    float_complex* A_xyz,
    float_complex* velCoeffs
)
{
    int i, j, d3, nSH_l, nSH;
    float_complex* c_nm, *A_1, *velCoeffs_T;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH_l = ORDER2NSH(order);
    nSH = ORDER2NSH(order+1);
    c_nm = malloc1d(nSH_l*sizeof(float_complex));
    A_1 = malloc1d(nSH*nSH_l*sizeof(float_complex));
    velCoeffs_T = malloc1d(3*nSH*sizeof(float_complex));
    rotateAxisCoeffsComplex(order, b_n, SAF_PI/2.0f-elev_rad, azi_rad, c_nm);
    
    /* x_nm, y_nm, z_nm */
    for(d3 = 0; d3<3; d3++){
        for(i=0; i<nSH; i++)
            for(j=0; j<nSH_l; j++)
                A_1[i*nSH_l+j] = A_xyz[i*nSH_l*3 + j*3 + d3];
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, 1, nSH_l, &calpha,
                    A_1, nSH_l,
                    c_nm, 1, &cbeta,
                    &velCoeffs_T[d3*nSH], 1);
    }
    for(d3 = 0; d3<3; d3++)
        for(i=0; i<nSH; i++)
            velCoeffs[i*3+d3] = velCoeffs_T[d3*nSH+i]; /* transpose */

    free(c_nm);
    free(A_1);
    free(velCoeffs_T);
}

void rotateAxisCoeffsReal
(
    int order,
    float* c_n,
    float theta_0, /* inclination*/
    float phi_0,  /* azimuth */
    float* c_nm
)
{
    int nSH;
    float_complex* c_nm_c;

    nSH = ORDER2NSH(order);
    c_nm_c = malloc1d(nSH*sizeof(float_complex));
    rotateAxisCoeffsComplex(order, c_n, theta_0, phi_0, c_nm_c);
    complex2realCoeffs(order, c_nm_c, 1, c_nm);
    
    free(c_nm_c);
}

void rotateAxisCoeffsComplex
(
    int order,
    float* c_n,
    float theta_0,  /* inclination*/
    float phi_0,    /* azimuth */
    float_complex* c_nm
)
{
    int n, m, q, nSH;
    float phi_theta[2];
    float_complex* Y_N;
    
    phi_theta[0] = phi_0;
    phi_theta[1] = theta_0;
    nSH = ORDER2NSH(order);
    Y_N = malloc1d(nSH*sizeof(float_complex));
    getSHcomplex(order, (float*)phi_theta, 1, Y_N);
    for(n=0, q = 0; n<=order; n++)
        for(m=-n; m<=n; m++, q++)
            c_nm[q] = crmulf(conjf(Y_N[q]), sqrtf(4.0f*SAF_PI/(2.0f*(float)n+1.0f)) * c_n[n]);
    
    free(Y_N);
}

void checkCondNumberSHTReal
(
    int order,
    float* dirs_rad,
    int nDirs,
    float* w,
    float* cond_N
)
{
    int n, i, j, nSH, nSH_n, ind;
    float minVal, maxVal;
    float *YY_n, *W, *W_Yn, *s;
    float** Y_N, **Y_n;
    
    /* get SH */
    nSH = ORDER2NSH(order);
    Y_N = (float**)malloc2d(nSH, nDirs, sizeof(float));
    Y_n = (float**)malloc2d(nDirs, nSH, sizeof(float));
    YY_n = malloc1d(nSH*nSH*sizeof(float));
    getSHreal(order, dirs_rad, nDirs, FLATTEN2D(Y_N));
    
    /* diagonalise the integration weights, if available */
    if(w!=NULL){
        W = calloc1d(nDirs*nDirs, sizeof(float));
        W_Yn = malloc1d(nDirs*nSH*sizeof(float));
        for(i=0; i<nDirs; i++)
            W[i*nDirs+i] = w[i]; 
    }
    else{
        W = NULL;
        W_Yn = NULL;
    }
    
    /* compute the condition number for each order up to N */
    s = malloc1d(nSH*sizeof(float));
    for(n=0; n<=order; n++){
        nSH_n = (n+1)*(n+1);
        for(i=0; i<nDirs; i++)
            for(j=0; j<nSH_n; j++)
                Y_n[i][j] = Y_N[j][i]; /* truncate to current order and transpose */
        if(w==NULL){
            cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nSH_n, nSH_n, nDirs, 1.0f,
                        FLATTEN2D(Y_n), nSH,
                        FLATTEN2D(Y_n), nSH, 0.0f,
                        YY_n, nSH_n);
        }
        else{
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nDirs, nSH_n, nDirs, 1.0f,
                        W, nDirs,
                        FLATTEN2D(Y_n), nSH, 0.0f,
                        W_Yn, nSH_n);
            cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nSH_n, nSH_n, nDirs, 1.0f,
                        FLATTEN2D(Y_n), nSH,
                        W_Yn, nSH_n, 0.0f,
                        YY_n, nSH_n);
        }
        
        /* condition number = max(singularValues)/min(singularValues) */
        utility_ssvd(NULL, YY_n, nSH_n, nSH_n, NULL, NULL, NULL, s);
        utility_simaxv(s, nSH_n, &ind);
        maxVal = s[ind];
        utility_siminv(s, nSH_n, &ind);
        minVal = s[ind];
        cond_N[n] = maxVal/(minVal+2.23e-7f);
    }
    
    free(Y_N);
    free(Y_n);
    free(YY_n);
    free(W);
    free(W_Yn); 
    free(s);
}


int calculateGridWeights
(
    float* dirs_rad,
    int nDirs,
    int order,
    float* w
)
{
    int i, j, nSH;
    float sumW;
    float** Y_N, **Y_N_T, **Y_leftinv;

    if(order<0){
        int nSH, ind;
        float minVal, maxVal, cond_N;

        float *YY_N, *s;

        /* get SH */
        s = NULL;
        Y_N = NULL;
        YY_N = NULL;

        for(int n=1; n<100; n++){
            /* compute the condition number for order N */
            nSH = ORDER2NSH(n);
            Y_N = (float**)realloc2d((void**)Y_N, nSH, nDirs, sizeof(float));
            YY_N = (float*)realloc1d(YY_N, nSH*nSH*sizeof(float));
            s = (float*) realloc1d(s, nSH*sizeof(float));
            getSHreal(n, dirs_rad, nDirs, FLATTEN2D(Y_N));

            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, nSH, nDirs, 1.0f,
                        FLATTEN2D(Y_N), nDirs,
                        FLATTEN2D(Y_N), nDirs, 0.0f,
                        YY_N, nSH);

            /* condition number = max(singularValues)/min(singularValues) */
            utility_ssvd(NULL, YY_N, nSH, nSH, NULL, NULL, NULL, s);
            utility_simaxv(s, nSH, &ind);
            maxVal = s[ind];
            utility_siminv(s, nSH, &ind);
            minVal = s[ind];
            cond_N = maxVal/(minVal+2.23e-7f);

            if(cond_N > 1.5f * n){
                order = n-1;
                break;
            }
        }
    }
    assert(order>0);
        
    nSH = ORDER2NSH(order);
    Y_N = (float**)malloc2d(nSH, nDirs, sizeof(float));
    Y_N_T = (float**)malloc2d(nDirs, nSH, sizeof(float));
    Y_leftinv = (float**)malloc2d(nSH, nDirs, sizeof(float));

    getSHreal(order, dirs_rad, nDirs, FLATTEN2D(Y_N));

    for(i=0; i<nDirs; i++)
        for(j=0; j<nSH; j++)
            Y_N_T[i][j] = Y_N[j][i]; /* truncate to current order and transpose */

    utility_spinv(NULL, FLATTEN2D(Y_N_T), nDirs, nSH, FLATTEN2D(Y_leftinv));

    sumW=0.f;
    for(int idx=0; idx<nDirs; idx++){
        w[idx] = sqrtf(FOURPI)*Y_leftinv[0][idx];
        sumW += w[idx];
    }

    if(fabs(sumW - FOURPI) > 0.001){
        saf_print_warning("Grid weights no bueno!");
    }
    return order;
}

/* ========================================================================== */
/*                     Localisation Functions in the  SHD                     */
/* ========================================================================== */

void sphPWD_create
(
    void ** const phPWD,
    int order,
    float* grid_dirs_deg,
    int nDirs
)
{
    *phPWD = malloc1d(sizeof(sphPWD_data));
    sphPWD_data *h = (sphPWD_data*)(*phPWD);
    int i, j;
    float** grid_dirs_rad, **grid_svecs_tmp;

    h->order = order;
    h->nSH = ORDER2NSH(h->order);
    h->nDirs = nDirs;

    /* steering vectors for each grid direction  */
    h->grid_svecs = malloc1d(h->nSH * (h->nDirs) * sizeof(float_complex));
    grid_dirs_rad  = (float**)malloc2d(h->nDirs, 2, sizeof(float));
    grid_svecs_tmp = (float**)malloc2d(h->nSH, h->nDirs, sizeof(float));
    for(i=0; i<h->nDirs; i++){
        grid_dirs_rad[i][0] = grid_dirs_deg[i*2]*SAF_PI/180.0f;
        grid_dirs_rad[i][1] = SAF_PI/2.0f - grid_dirs_deg[i*2+1]*SAF_PI/180.0f;
    }
    getSHreal(h->order, FLATTEN2D(grid_dirs_rad), h->nDirs, FLATTEN2D(grid_svecs_tmp));
    for(i=0; i<h->nSH; i++)
        for(j=0; j<h->nDirs; j++)
            h->grid_svecs[j*(h->nSH)+i] = cmplxf(grid_svecs_tmp[i][j], 0.0f);

    /* store cartesian coords of scanning directions (for optional peak finding) */
    h->grid_dirs_xyz = malloc1d(h->nDirs * 3 * sizeof(float));
    unitSph2cart(grid_dirs_deg, h->nDirs, 1, h->grid_dirs_xyz);

    /* for run-time */
    h->A_Cx = malloc1d((h->nSH) * sizeof(float_complex));
    h->pSpec = malloc1d(h->nDirs*sizeof(float));
    h->P_minus_peak = malloc1d(h->nDirs*sizeof(float));
    h->VM_mask = malloc1d(h->nDirs*sizeof(float));
    h->P_tmp = malloc1d(h->nDirs*sizeof(float));

    /* clean-up */
    free(grid_dirs_rad);
    free(grid_svecs_tmp);
}

void sphPWD_destroy
(
    void ** const phPWD
)
{
    sphPWD_data *h = (sphPWD_data*)(*phPWD);

    if (h != NULL) {
        free(h->grid_dirs_xyz);
        free(h->grid_svecs);
        free(h->A_Cx);
        free(h->pSpec);
        free(h->P_minus_peak);
        free(h->P_tmp);
        free(h->VM_mask);
        free(h);
        h = NULL;
        *phPWD = NULL;
    }
}

void sphPWD_compute
(
    void* const hPWD,
    float_complex* Cx,
    int nSrcs,
    float* P_map,
    int* peak_inds
)
{
    sphPWD_data *h = (sphPWD_data*)(hPWD);
    int i, k, peak_idx;
    float kappa, scale;
    float VM_mean[3];
    float_complex A_Cx_A;
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta = cmplxf(0.0f, 0.0f);

    /* derive the power-map value for each grid direction */ 
    for (i = 0; i < (h->nDirs); i++){ 
        cblas_cgemv(CblasRowMajor, CblasNoTrans, h->nSH, h->nSH, &calpha,
                    Cx, h->nSH,
                    &(h->grid_svecs[i*(h->nSH)]), 1, &cbeta,
                    h->A_Cx, 1); 
        cblas_cdotu_sub(h->nSH, h->A_Cx, 1, &(h->grid_svecs[i*(h->nSH)]), 1, &A_Cx_A);
        h->pSpec[i] = crealf(A_Cx_A);
    }

    /* Output power-map */
    if(P_map!=NULL)
        cblas_scopy(h->nDirs, h->pSpec, 1, P_map, 1);

    /* Peak-finding */
    if(peak_inds!=NULL){
        kappa = 50.0f;
        scale = kappa/(2.0f*SAF_PI*expf(kappa)-expf(-kappa));
        cblas_scopy(h->nDirs, h->pSpec, 1, h->P_minus_peak, 1);

        /* Loop over the number of sources */
        for(k=0; k<nSrcs; k++){
            utility_simaxv(h->P_minus_peak, h->nDirs, &peak_idx);
            peak_inds[k] = peak_idx;
            if(k==nSrcs-1)
                break;
            VM_mean[0] = h->grid_dirs_xyz[peak_idx*3];
            VM_mean[1] = h->grid_dirs_xyz[peak_idx*3+1];
            VM_mean[2] = h->grid_dirs_xyz[peak_idx*3+2];

            /* Apply mask for next iteration */
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, h->nDirs, 1, 3, 1.0f,
                        h->grid_dirs_xyz, 3,
                        (const float*)VM_mean, 3, 0.0f,
                        h->VM_mask, 1);
            cblas_sscal(h->nDirs, kappa, h->VM_mask, 1);
            for(i=0; i<h->nDirs; i++)
                h->VM_mask[i] = expf(h->VM_mask[i]);             /* VM distribution */
            cblas_sscal(h->nDirs, scale, h->VM_mask, 1);
            for(i=0; i<h->nDirs; i++)
                h->VM_mask[i] = 1.0f/(0.00001f+(h->VM_mask[i])); /* inverse VM distribution */
            utility_svvmul(h->P_minus_peak, h->VM_mask, h->nDirs, h->P_tmp);
            cblas_scopy(h->nDirs, h->P_tmp, 1, h->P_minus_peak, 1);
        }
    }
}

void sphMUSIC_create
(
    void ** const phMUSIC,
    int order,
    float* grid_dirs_deg,
    int nDirs
)
{
    *phMUSIC = malloc1d(sizeof(sphMUSIC_data));
    sphMUSIC_data *h = (sphMUSIC_data*)(*phMUSIC);
    int i, j;
    float** grid_dirs_rad, **grid_svecs_tmp;

    h->order = order;
    h->nSH = ORDER2NSH(h->order);
    h->nDirs = nDirs;

    /* steering vectors for each grid direction  */
    h->grid_svecs = malloc1d(h->nSH * (h->nDirs) * sizeof(float_complex));
    grid_dirs_rad  = (float**)malloc2d(h->nDirs, 2, sizeof(float));
    grid_svecs_tmp = (float**)malloc2d(h->nSH, h->nDirs, sizeof(float));
    for(i=0; i<h->nDirs; i++){
        grid_dirs_rad[i][0] = grid_dirs_deg[i*2]*SAF_PI/180.0f;
        grid_dirs_rad[i][1] = SAF_PI/2.0f - grid_dirs_deg[i*2+1]*SAF_PI/180.0f;
    }
    getSHreal(h->order, FLATTEN2D(grid_dirs_rad), h->nDirs, FLATTEN2D(grid_svecs_tmp));
    for(i=0; i<h->nSH; i++)
        for(j=0; j<h->nDirs; j++)
            h->grid_svecs[i*(h->nDirs)+j] = cmplxf(grid_svecs_tmp[i][j], 0.0f);

    /* store cartesian coords of scanning directions (for optional peak finding) */
    h->grid_dirs_xyz = malloc1d(h->nDirs * 3 * sizeof(float));
    unitSph2cart(grid_dirs_deg, h->nDirs, 1, h->grid_dirs_xyz);

    /* for run-time */
    h->VnA = malloc1d(h->nSH * (h->nDirs) * sizeof(float_complex));
    h->abs_VnA = malloc1d(h->nSH * (h->nDirs) * sizeof(float));
    h->pSpec = malloc1d(h->nDirs*sizeof(float));
    h->pSpecInv = malloc1d(h->nDirs*sizeof(float));
    h->P_minus_peak = malloc1d(h->nDirs*sizeof(float));
    h->P_tmp = malloc1d(h->nDirs*sizeof(float));
    h->VM_mask = malloc1d(h->nDirs*sizeof(float));

    /* clean-up */
    free(grid_dirs_rad);
    free(grid_svecs_tmp);
}

void sphMUSIC_destroy
(
    void ** const phMUSIC
)
{
    sphMUSIC_data *h = (sphMUSIC_data*)(*phMUSIC);

    if (h != NULL) {
        free(h->grid_dirs_xyz);
        free(h->grid_svecs);
        free(h->VnA);
        free(h->abs_VnA);
        free(h->pSpec);
        free(h->pSpecInv);
        free(h->P_minus_peak);
        free(h->P_tmp);
        free(h->VM_mask);
        free(h);
        h = NULL;
        *phMUSIC = NULL;
    }
}

void sphMUSIC_compute
(
    void* const hMUSIC,
    float_complex *Vn, /* nSH x (nSH - nSrcs) */
    int nSrcs,
    float* P_music,
    int* peak_inds
)
{
    sphMUSIC_data *h = (sphMUSIC_data*)(hMUSIC);
    int i, k, VnD2, peak_idx;
    float kappa, scale;
    float VM_mean[3];
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta = cmplxf(0.0f, 0.0f);

    VnD2 = h->nSH - nSrcs; /* noise subspace second dimension length */

    /* derive the pseudo-spectrum value for each grid direction */
    cblas_cgemm(CblasRowMajor, CblasTrans, CblasNoTrans, h->nDirs, VnD2, h->nSH, &calpha,
                h->grid_svecs, h->nDirs,
                Vn, VnD2, &cbeta,
                h->VnA, VnD2);
    utility_cvabs(h->VnA, (h->nDirs)*VnD2, h->abs_VnA);
    for (i = 0; i < (h->nDirs); i++)
        h->pSpecInv[i] = cblas_sdot(VnD2, &(h->abs_VnA[i*VnD2]), 1, &(h->abs_VnA[i*VnD2]), 1);
    //h->pSpec[i] = 1.0f / (h->pSpecInv[i] + 2.23e-10f);
    utility_svrecip(h->pSpecInv, h->nDirs, h->pSpec);

    /* Output pseudo-spectrum */
    if(P_music!=NULL)
        cblas_scopy(h->nDirs, h->pSpec, 1, P_music, 1);

    /* Peak-finding */
    if(peak_inds!=NULL){
        kappa = 50.0f;
        scale = kappa/(2.0f*SAF_PI*expf(kappa)-expf(-kappa));
        cblas_scopy(h->nDirs, h->pSpec, 1, h->P_minus_peak, 1);

        /* Loop over the number of sources */
        for(k=0; k<nSrcs; k++){
            utility_simaxv(h->P_minus_peak, h->nDirs, &peak_idx);
            peak_inds[k] = peak_idx;
            if(k==nSrcs-1)
                break;
            VM_mean[0] = h->grid_dirs_xyz[peak_idx*3];
            VM_mean[1] = h->grid_dirs_xyz[peak_idx*3+1];
            VM_mean[2] = h->grid_dirs_xyz[peak_idx*3+2];

            /* Apply mask for next iteration */
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, h->nDirs, 1, 3, 1.0f,
                        h->grid_dirs_xyz, 3,
                        (const float*)VM_mean, 3, 0.0f,
                        h->VM_mask, 1);
            cblas_sscal(h->nDirs, kappa, h->VM_mask, 1);
            for(i=0; i<h->nDirs; i++)
                h->VM_mask[i] = expf(h->VM_mask[i]);             /* VM distribution */
            cblas_sscal(h->nDirs, scale, h->VM_mask, 1);
            for(i=0; i<h->nDirs; i++)
                h->VM_mask[i] = 1.0f/(0.00001f+(h->VM_mask[i])); /* inverse VM distribution */
            utility_svvmul(h->P_minus_peak, h->VM_mask, h->nDirs, h->P_tmp);
            cblas_scopy(h->nDirs, h->P_tmp, 1, h->P_minus_peak, 1);
        }
    }
}

void sphESPRIT_create
(
    void ** const phESPRIT,
    int order
)
{
    *phESPRIT = malloc1d(sizeof(sphESPRIT_data));
    sphESPRIT_data *h = (sphESPRIT_data*)(*phESPRIT);
    int i, j;

    h->N = order;
    h->NN = order*order;
    h->maxK = h->NN;

    /* pre-compute indices and matrices */
    for(i=0; i<6; i++){
        h->rWVnimu[i] = malloc1d((order*order) * (order*order) * sizeof(double));
        h->WVnimu[i]  = malloc1d((order*order) * (order*order) * sizeof(double_complex));
    }
    h->nIdx[0] = h->nIdx[1] = h->nIdx[4] = h->nIdx[5] = h->nIdx[10] = h->nIdx[11] = (order*order);
    h->nIdx[2] = h->nIdx[3] = h->nIdx[6] = h->nIdx[7] = h->nIdx[8]  = h->nIdx[9]  = ((order-1)*(order-1));
    for(i=0; i<12; i++){
        if(h->nIdx[i] == 0)
            h->idx_from_Ynm2Ynimu[i] = NULL;
        else
            h->idx_from_Ynm2Ynimu[i] = calloc1d(h->nIdx[i], sizeof(int));
    }
    getWnimu(order, 1, 1,-1, h->rWVnimu[0]);
    getWnimu(order,-1, 0, 0, h->rWVnimu[1]);
    getWnimu(order,-1, 1,-1, h->rWVnimu[2]);
    getWnimu(order, 1, 0, 0, h->rWVnimu[3]);
    getVnimu(order, 0, 0,    h->rWVnimu[4]);
    getVnimu(order, 1, 0,    h->rWVnimu[5]);
    for(i=0; i<6; i++){
        for(j=0; j<(order*order) * (order*order); j++)
            h->WVnimu[i][j] = cmplx(h->rWVnimu[i][j], 0.0);
    }
    muni2q(order, 1,-1, h->idx_from_Ynm2Ynimu[0],  h->idx_from_Ynm2Ynimu[1]);
    muni2q(order,-1,-1, h->idx_from_Ynm2Ynimu[2],  h->idx_from_Ynm2Ynimu[3]);
    muni2q(order, 1, 1, h->idx_from_Ynm2Ynimu[4],  h->idx_from_Ynm2Ynimu[5]);
    muni2q(order,-1, 1, h->idx_from_Ynm2Ynimu[6],  h->idx_from_Ynm2Ynimu[7]);
    muni2q(order,-1, 0, h->idx_from_Ynm2Ynimu[8],  h->idx_from_Ynm2Ynimu[9]);
    muni2q(order, 1, 0, h->idx_from_Ynm2Ynimu[10], h->idx_from_Ynm2Ynimu[11]);

    /* memory allocations for run-time matrices */
    utility_zpinv_create(&(h->hZpinv), h->maxK, h->maxK);
    utility_zeigmp_create(&(h->hZeigmp), h->maxK);
    utility_zglslv_create(&(h->hZglslv), h->maxK, h->maxK);
    h->Us_1m1  = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->Us_m1m1 = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->Us_11   = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->Us_m11  = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->Us_m10  = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->Us_10   = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->Us_00   = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->WVnimu0_Us1m1  = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->WVnimu1_Usm1m1 = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->WVnimu2_Us11   = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->WVnimu3_Usm11  = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->WVnimu4_Usm10  = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->WVnimu5_Us10   = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->LambdaXYp      = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->LambdaXYm      = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->LambdaZ        = malloc1d((h->NN) * (h->maxK) * sizeof(double_complex));
    h->pinvUs = malloc1d((h->maxK) * (h->NN)   * sizeof(double_complex));
    h->PsiXYp = malloc1d((h->maxK) * (h->maxK) * sizeof(double_complex));
    h->PsiXYm = malloc1d((h->maxK) * (h->maxK) * sizeof(double_complex));
    h->PsiZ   = malloc1d((h->maxK) * (h->maxK) * sizeof(double_complex));
    h->tmp_KK = malloc1d((h->maxK) * (h->maxK) * sizeof(double_complex));
    h->V      = malloc1d((h->maxK) * (h->maxK) * sizeof(double_complex));
    h->PhiXYp = malloc1d((h->maxK) * (h->maxK) * sizeof(double_complex));
    h->PhiXYm = malloc1d((h->maxK) * (h->maxK) * sizeof(double_complex));
    h->PhiZ   = malloc1d((h->maxK) * (h->maxK) * sizeof(double_complex));
}

void sphESPRIT_destroy
(
    void ** const phESPRIT
)
{
    sphESPRIT_data *h = (sphESPRIT_data*)(*phESPRIT);
    int i;

    if (h != NULL) {
        for(i=0; i<6; i++){
            free(h->rWVnimu[i]);
            free(h->WVnimu[i]);
        }
        for(i=0; i<12; i++)
            free(h->idx_from_Ynm2Ynimu[i]);
        utility_zpinv_destroy(&(h->hZpinv));
        utility_zeigmp_destroy(&(h->hZeigmp));
        utility_zglslv_destroy(&(h->hZglslv));
        free(h->Us_1m1);
        free(h->Us_m1m1);
        free(h->Us_11);
        free(h->Us_m11);
        free(h->Us_m10);
        free(h->Us_10);
        free(h->Us_00);
        free(h->WVnimu0_Us1m1);
        free(h->WVnimu1_Usm1m1);
        free(h->WVnimu2_Us11);
        free(h->WVnimu3_Usm11);
        free(h->WVnimu4_Usm10);
        free(h->WVnimu5_Us10);
        free(h->LambdaXYp);
        free(h->LambdaXYm);
        free(h->LambdaZ);
        free(h->pinvUs);
        free(h->PsiXYp);
        free(h->PsiXYm);
        free(h->PsiZ);
        free(h->tmp_KK);
        free(h->V);
        free(h->PhiXYp);
        free(h->PhiXYm);
        free(h->PhiZ);

        free(h);
        h = NULL;
        *phESPRIT = NULL;
    }
}

void sphESPRIT_estimateDirs
(
    void * const hESPRIT,
    float_complex* Us, /* nSH * K */
    int K,
    float* src_dirs_rad /*  K x 2 */
)
{
    sphESPRIT_data *h = (sphESPRIT_data*)(hESPRIT);
    int i, j;
    const double_complex i2_ = cmplx(0.0, 2.0);
    const double_complex calpha = cmplx(1.0, 0.0); const double_complex cbeta = cmplx(0.0, 0.0); /* blas */
    double phiX, phiY;

    /* Fill matrices */
    memset(h->Us_1m1, 0, (h->NN) * K * sizeof(double_complex));
    memset(h->Us_m1m1, 0, (h->NN) * K * sizeof(double_complex));
    memset(h->Us_11, 0, (h->NN) * K * sizeof(double_complex));
    memset(h->Us_m11, 0, (h->NN) * K * sizeof(double_complex));
    memset(h->Us_m10, 0, (h->NN) * K * sizeof(double_complex));
    memset(h->Us_10, 0, (h->NN) * K * sizeof(double_complex));
    memset(h->Us_00, 0, (h->NN) * K * sizeof(double_complex));
    for (i = 0; i < K; i++) {
        for (j = 0; j < h->nIdx[0]; j++)
            h->Us_1m1[h->idx_from_Ynm2Ynimu[0][j] * K + i] = cmplx(crealf(Us[h->idx_from_Ynm2Ynimu[1][j] * K + i]), cimagf(Us[h->idx_from_Ynm2Ynimu[1][j] * K + i]));
        for (j = 0; j < h->nIdx[2]; j++)
            h->Us_m1m1[h->idx_from_Ynm2Ynimu[2][j] * K + i] = cmplx(crealf(Us[h->idx_from_Ynm2Ynimu[3][j] * K + i]), cimagf(Us[h->idx_from_Ynm2Ynimu[3][j] * K + i]));
        for (j = 0; j < h->nIdx[4]; j++)
            h->Us_11[h->idx_from_Ynm2Ynimu[4][j] * K + i] = cmplx(crealf(Us[h->idx_from_Ynm2Ynimu[5][j] * K + i]), cimagf(Us[h->idx_from_Ynm2Ynimu[5][j] * K + i]));
        for (j = 0; j < h->nIdx[6]; j++)
            h->Us_m11[h->idx_from_Ynm2Ynimu[6][j] * K + i] = cmplx(crealf(Us[h->idx_from_Ynm2Ynimu[7][j] * K + i]), cimagf(Us[h->idx_from_Ynm2Ynimu[7][j] * K + i]));
        for (j = 0; j < h->nIdx[8]; j++)
            h->Us_m10[h->idx_from_Ynm2Ynimu[8][j] * K + i] = cmplx(crealf(Us[h->idx_from_Ynm2Ynimu[9][j] * K + i]), cimagf(Us[h->idx_from_Ynm2Ynimu[9][j] * K + i]));
        for (j = 0; j < h->nIdx[10]; j++)
            h->Us_10[h->idx_from_Ynm2Ynimu[10][j] * K + i] = cmplx(crealf(Us[h->idx_from_Ynm2Ynimu[11][j] * K + i]), cimagf(Us[h->idx_from_Ynm2Ynimu[11][j] * K + i]));
        for (j = 0; j < (h->NN); j++)
            h->Us_00[j*K + i] = cmplx(crealf(Us[j*K + i]), cimagf(Us[j*K + i]));
    }

    /*  */
    cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, (h->NN), K, (h->NN), &calpha,
                h->WVnimu[0], (h->NN),
                h->Us_1m1, K, &cbeta,
                h->WVnimu0_Us1m1, K);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (h->NN), K, (h->NN), &calpha,
                h->WVnimu[1], (h->NN),
                h->Us_m1m1, K, &cbeta,
                h->WVnimu1_Usm1m1, K);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (h->NN), K, (h->NN), &calpha,
                h->WVnimu[2], (h->NN),
                h->Us_11, K, &cbeta,
                h->WVnimu2_Us11, K);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (h->NN), K, (h->NN), &calpha,
                h->WVnimu[3], (h->NN),
                h->Us_m11, K, &cbeta,
                h->WVnimu3_Usm11, K);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (h->NN), K, (h->NN), &calpha,
                h->WVnimu[4], (h->NN),
                h->Us_m10, K, &cbeta,
                h->WVnimu4_Usm10, K);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (h->NN), K, (h->NN), &calpha,
                h->WVnimu[5], (h->NN),
                h->Us_10, K, &cbeta,
                h->WVnimu5_Us10, K);
    utility_zvvsub(h->WVnimu0_Us1m1, h->WVnimu1_Usm1m1, h->NN*K, h->LambdaXYp);
    cblas_dscal(/*re+im*/2*h->NN*K, -1.0, (double*)h->WVnimu2_Us11, 1);
    utility_zvvadd(h->WVnimu2_Us11, h->WVnimu3_Usm11, h->NN*K, h->LambdaXYm);
    utility_zvvadd(h->WVnimu4_Usm10, h->WVnimu5_Us10, h->NN*K, h->LambdaZ);

    /*  */
    utility_zpinv(h->hZpinv, h->Us_00, h->NN, K, h->pinvUs);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, K, K, (h->NN), &calpha,
                h->pinvUs, (h->NN),
                h->LambdaXYp, K, &cbeta,
                h->PsiXYp, K);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, K, K, (h->NN), &calpha,
                h->pinvUs, (h->NN),
                h->LambdaXYm, K, &cbeta,
                h->PsiXYm, K);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, K, K, (h->NN), &calpha,
                h->pinvUs, (h->NN),
                h->LambdaZ, K, &cbeta,
                h->PsiZ, K);

    /*  */
    utility_zeigmp(h->hZeigmp, h->PsiXYp, h->PsiZ, K,  NULL, h->V, NULL);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, K, K, K, &calpha,
                h->PsiXYp, K,
                h->V, K, &cbeta,
                h->tmp_KK, K);
    utility_zglslv(h->hZglslv, h->V, K, h->tmp_KK, K, h->PhiXYp);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, K, K, K, &calpha,
                h->PsiXYm, K,
                h->V, K, &cbeta,
                h->tmp_KK, K);
    utility_zglslv(h->hZglslv, h->V, K, h->tmp_KK, K, h->PhiXYm);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, K, K, K, &calpha,
                h->PsiZ, K,
                h->V, K, &cbeta,
                h->tmp_KK, K);
    utility_zglslv(h->hZglslv, h->V, K, h->tmp_KK, K, h->PhiZ);

    /* extract DoAs */
    for(i=0; i<K; i++){
        phiX = (creal(h->PhiXYp[i*K+i])+creal(h->PhiXYm[i*K+i]))/2.0;
        phiY = creal(ccdiv(ccsub(h->PhiXYp[i*K+i], h->PhiXYm[i*K+i]), i2_));
        src_dirs_rad[i*2] = (float)atan2(phiY, phiX);
        src_dirs_rad[i*2+1] = (float)SAF_MIN(atan2(creal(h->PhiZ[i*K+i]), sqrt(phiX*phiX+phiY*phiY)), SAF_PI/2.0f);
    }
}

void generatePWDmap
(
    int order,
    float_complex* Cx,
    float_complex* Y_grid,
    int nGrid_dirs,
    float* pmap
)
{
    int i, j, nSH;
    float_complex* Cx_Y, *Y_Cx_Y, *Cx_Y_s, *Y_grid_s;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH = ORDER2NSH(order);
    Cx_Y = malloc1d(nSH * nGrid_dirs * sizeof(float_complex));
    Y_Cx_Y = malloc1d(nGrid_dirs*sizeof(float_complex));
    Cx_Y_s = malloc1d(nSH*sizeof(float_complex));
    Y_grid_s = malloc1d(nSH*sizeof(float_complex));
    
    /* Calculate PWD powermap: real(diag(Y_grid.'*C_x*Y_grid)) */
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, nGrid_dirs, nSH, &calpha,
                Cx, nSH,
                Y_grid, nGrid_dirs, &cbeta,
                Cx_Y, nGrid_dirs);
    for(i=0; i<nGrid_dirs; i++){
        for(j=0; j<nSH; j++){
            Cx_Y_s[j] = Cx_Y[j*nGrid_dirs+i];
            Y_grid_s[j] = Y_grid[j*nGrid_dirs+i];
        }
        /* faster to perform the dot-product for each vector seperately */
        utility_cvvdot(Y_grid_s, Cx_Y_s, nSH, NO_CONJ, &Y_Cx_Y[i]);
    }
    
    for(i=0; i<nGrid_dirs; i++)
        pmap[i] = crealf(Y_Cx_Y[i]);
    
    free(Cx_Y);
    free(Y_Cx_Y);
    free(Cx_Y_s);
    free(Y_grid_s);
}

void generateMVDRmap
(
    int order,
    float_complex* Cx,
    float_complex* Y_grid,
    int nGrid_dirs,
    float regPar,
    float* pmap,
    float_complex* w_MVDR_out
)
{
    int i, j, nSH;
    float Cx_trace;
    float_complex *Cx_d, *invCx_Ygrid, *w_MVDR, *invCx_Ygrid_s, *Y_grid_s;
    float_complex denum;
    
    nSH = ORDER2NSH(order);
    w_MVDR = malloc1d(nSH * nGrid_dirs*sizeof(float_complex));
    Cx_d = malloc1d(nSH*nSH*sizeof(float_complex));
    invCx_Ygrid = malloc1d(nSH*nGrid_dirs*sizeof(float_complex));
    invCx_Ygrid_s = malloc1d(nSH*sizeof(float_complex));
    Y_grid_s = malloc1d(nSH*sizeof(float_complex));
    
    /* apply diagonal loading */
    Cx_trace = 0.0f;
    for(i=0; i<nSH; i++)
        Cx_trace += crealf(Cx[i*nSH+i]);
    Cx_trace /= (float)nSH;
    memcpy(Cx_d, Cx, nSH*nSH*sizeof(float_complex));
    for(i=0; i<nSH; i++)
        Cx_d[i*nSH+i] = craddf(Cx_d[i*nSH+i], regPar*Cx_trace);
    
    /* solve the numerator part of the MVDR weights for all grid directions: Cx^-1 * Y */
    utility_cslslv(NULL, Cx_d, nSH, Y_grid, nGrid_dirs, invCx_Ygrid);
    for(i=0; i<nGrid_dirs; i++){
        /* solve the denumerator part of the MVDR weights for each grid direction: Y^T * Cx^-1 * Y */
        for(j=0; j<nSH; j++){
            invCx_Ygrid_s[j] = conjf(invCx_Ygrid[j*nGrid_dirs+i]);
            Y_grid_s[j] = Y_grid[j*nGrid_dirs+i];
        }
        /* faster to perform the dot-product for each vector seperately */
        utility_cvvdot(Y_grid_s, invCx_Ygrid_s, nSH, NO_CONJ, &denum);
        
        /* calculate the MVDR weights per grid direction: (Cx^-1 * Y) * (Y^T * Cx^-1 * Y)^-1 */
        for(j=0; j<nSH; j++)
            w_MVDR[j*nGrid_dirs +i] = ccdivf(invCx_Ygrid[j*nGrid_dirs +i], denum);
    }
    
    /* generate MVDR powermap, by using the generatePWDmap function with the MVDR weights instead */
    generatePWDmap(order, Cx, w_MVDR, nGrid_dirs, pmap);
    
    /* optional output of the beamforming weights */
    if (w_MVDR_out!=NULL)
        memcpy(w_MVDR_out, w_MVDR, nSH * nGrid_dirs*sizeof(float_complex));
    
    free(w_MVDR);
    free(Cx_d);
    free(invCx_Ygrid);
    free(invCx_Ygrid_s);
    free(Y_grid_s);
}

/* EXPERIMENTAL
 * Delikaris-Manias, S., Vilkamo, J., & Pulkki, V. (2016). Signal-dependent spatial filtering based on
 * weighted-orthogonal beamformers in the spherical harmonic domain. IEEE/ACM Transactions on Audio,
 * Speech and Language Processing (TASLP), 24(9), 1507-1519. */
void generateCroPaCLCMVmap
(
    int order,
    float_complex* Cx,
    float_complex* Y_grid,
    int nGrid_dirs,
    float regPar,
    float lambda,
    float* pmap  
)
{
    int i, j, k, nSH;
    float Cx_trace, S, G;
    float* mvdr_map;
    float_complex* Cx_d, *A, *invCxd_A, *invCxd_A_tmp, *w_LCMV_s, *w_CroPaC, *wo, *Cx_Y, *Cx_Y_s;
    float_complex b[2]; 
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    float_complex A_invCxd_A[2][2];
    float_complex Y_wo_xspec;
    
    b[0] = cmplxf(1.0f, 0.0f);
    b[1] = cmplxf(0.0f, 0.0f);
    nSH = ORDER2NSH(order);
    Cx_Y = malloc1d(nSH * nGrid_dirs * sizeof(float_complex));
    Cx_d = malloc1d(nSH*nSH*sizeof(float_complex));
    A = malloc1d(nSH*2*sizeof(float_complex));
    invCxd_A = malloc1d(nSH*2*sizeof(float_complex));
    invCxd_A_tmp = malloc1d(nSH*2*sizeof(float_complex));
    w_LCMV_s = malloc1d(2*nGrid_dirs*sizeof(float_complex));
    w_CroPaC = malloc1d(nSH*nGrid_dirs*sizeof(float_complex));
    wo = malloc1d(nSH*sizeof(float_complex));
    mvdr_map = malloc1d(nGrid_dirs*sizeof(float));
    Cx_Y_s = malloc1d(nSH*sizeof(float_complex));
    
    /* generate MVDR map and weights to use as a basis */
    generateMVDRmap(order, Cx, Y_grid, nGrid_dirs, regPar, mvdr_map, w_CroPaC);
    
    /* first half of the cross-spectrum */
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, nGrid_dirs, nSH, &calpha,
                Cx, nSH,
                Y_grid, nGrid_dirs, &cbeta,
                Cx_Y, nGrid_dirs);
    
    /* apply diagonal loading to cov matrix */
    Cx_trace = 0.0f;
    for(i=0; i<nSH; i++)
        Cx_trace += crealf(Cx[i*nSH+i]);
    Cx_trace /= (float)nSH;
    memcpy(Cx_d, Cx, nSH*nSH*sizeof(float_complex));
    for(i=0; i<nSH; i++)
        Cx_d[i*nSH+i] = craddf(Cx_d[i*nSH+i], regPar*Cx_trace);
    
    /* calculate CroPaC beamforming weights for each grid direction */
    for(i=0; i<nGrid_dirs; i++){
        for(j=0; j<nSH; j++){
            A[j*2] = Y_grid[j*nGrid_dirs+i];
            A[j*2+1] = ccmulf(A[j*2], Cx[j*nSH+j]);
        }
        
        /* solve for minimisation problem for LCMV weights: (Cx^-1 * A) * (A^H * Cx^-1 * A)^-1 * b */
        utility_cslslv(NULL, Cx_d, nSH, A, 2, invCxd_A);
        for(j=0; j<nSH*2; j++)
            invCxd_A_tmp[j] = conjf(invCxd_A[j]);
        cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, 2, 2, nSH, &calpha,
                    A, 2,
                    invCxd_A_tmp, 2, &cbeta,
                    A_invCxd_A, 2);
        for(j=0; j<nSH; j++)
            for(k=0; k<2; k++)
                invCxd_A_tmp[k*nSH+j] = invCxd_A[j*2+k];
        utility_cglslv(NULL, (float_complex*)A_invCxd_A, 2, invCxd_A_tmp, nSH, w_LCMV_s);
        cblas_cgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nSH, 1, 2, &calpha,
                    w_LCMV_s, nSH,
                    b, 1, &cbeta,
                    wo, 1);
        
        /* calculate the cross-spectrum between static beam Y, and adaptive beam wo (LCMV) */
        for(j=0; j<nSH; j++)
            Cx_Y_s[j] = Cx_Y[j*nGrid_dirs+i];
        utility_cvvdot(wo, Cx_Y_s, nSH, NO_CONJ, &Y_wo_xspec);
        
        /* derive CroPaC weights  */
        S = SAF_MIN(cabsf(Y_wo_xspec), mvdr_map[i]); /* ensures distortionless response */
        G = sqrtf(S/(mvdr_map[i]+2.23e-10f));
        G = SAF_MAX(lambda, G); /* optional spectral floor parameter, to control harshness of attenuation (good for demos) */
        for(j=0; j<nSH; j++)
            w_CroPaC[j*nGrid_dirs + i] = crmulf(w_CroPaC[j*nGrid_dirs + i], G);
    }
    
    /* generate CroPaC powermap, by using the generatePWDmap function with the CroPaC weights instead */
    generatePWDmap(order, Cx, w_CroPaC, nGrid_dirs, pmap);
    
    free(mvdr_map);
    free(Cx_d);
    free(A);
    free(invCxd_A);
    free(invCxd_A_tmp);
    free(w_LCMV_s);
    free(w_CroPaC);
    free(wo);
    free(Cx_Y);
    free(Cx_Y_s);
}

void generateMUSICmap
(
    int order,
    float_complex* Cx,
    float_complex* Y_grid,
    int nSources,
    int nGrid_dirs,
    int logScaleFlag,
    float* pmap
)
{
    int i, j, nSH;
    float_complex* V, *Vn, *Vn_Y;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    float_complex tmp;
    
    nSH = ORDER2NSH(order);
    nSources = SAF_MIN(nSources, nSH/2);
    V = malloc1d(nSH*nSH*sizeof(float_complex));
    Vn = malloc1d(nSH*(nSH-nSources)*sizeof(float_complex));
    Vn_Y = malloc1d((nSH-nSources)*nGrid_dirs*sizeof(float_complex));
    
    /* obtain eigenvectors */
    //utility_ceig(Cx, nSH, 1, NULL, V, NULL, NULL);
    utility_cseig(NULL, Cx, nSH, 1, V, NULL, NULL);
    
    /* truncate, to obtain noise sub-space */
    for (i = 0; i < nSH; i++)
        for (j = 0; j < nSH - nSources; j++)
            Vn[i*(nSH - nSources) + j] = V[i*nSH + j + nSources];
    
    /* derive the pseudo-spectrum value for each grid direction */
    cblas_cgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nSH-nSources, nGrid_dirs, nSH, &calpha,
                Vn, nSH-nSources,
                Y_grid, nGrid_dirs, &cbeta,
                Vn_Y, nGrid_dirs);
    for(i=0; i<nGrid_dirs; i++){
        tmp = cmplxf(0.0f,0.0f);
        for(j=0; j<nSH-nSources; j++)
            tmp = ccaddf(tmp, ccmulf(conjf(Vn_Y[j*nGrid_dirs+i]),Vn_Y[j*nGrid_dirs+i]));
        pmap[i] = logScaleFlag ? logf(1.0f/(crealf(tmp)+2.23e-10f)) : 1.0f/(crealf(tmp)+2.23e-10f);
    }
    
    free(V);
    free(Vn);
    free(Vn_Y);
}

void generateMinNormMap
(
    int order,
    float_complex* Cx,
    float_complex* Y_grid,
    int nSources,
    int nGrid_dirs,
    int logScaleFlag,
    float* pmap
)
{
    int i, j, nSH;
    float_complex* V, *Vn, *Vn1, *Un, *Un_Y;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    float_complex Vn1_Vn1H;
    
    nSH = ORDER2NSH(order);
    nSources = SAF_MIN(nSources, nSH/2);
    V = malloc1d(nSH*nSH*sizeof(float_complex));
    Vn = malloc1d(nSH*(nSH-nSources)*sizeof(float_complex));
    Vn1 = malloc1d((nSH-nSources)*sizeof(float_complex));
    Un = malloc1d(nSH*sizeof(float_complex));
    Un_Y = malloc1d(nGrid_dirs*sizeof(float_complex));
    
    /* obtain eigenvectors */
    utility_ceig(NULL, Cx, nSH, NULL, V, NULL, NULL);
    
    /* truncate, to obtain noise sub-space */
    for(i=0; i<nSH; i++)
        for(j=0; j<nSH-nSources; j++)
            Vn[i*(nSH-nSources)+j] = V[i*nSH + j + nSources];
    for(j=0; j<nSH-nSources; j++)
        Vn1[j] = V[j + nSources];
    
    /* derive the pseudo-spectrum value for each grid direction */
    utility_cvvdot(Vn1, Vn1, nSH-nSources, NO_CONJ, &Vn1_Vn1H);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, 1, nSH-nSources, &calpha,
                Vn, nSH-nSources,
                Vn1, nSH-nSources, &cbeta,
                Un, 1);
    for(i=0; i<nSH; i++)
        Un[i] = ccdivf(Un[i], craddf(Vn1_Vn1H, 2.23e-9f));
    cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, 1, nGrid_dirs, nSH, &calpha,
                Un, 1,
                Y_grid, nGrid_dirs, &cbeta,
                Un_Y, nGrid_dirs);
    for(i=0; i<nGrid_dirs; i++)
        pmap[i] = logScaleFlag ? logf(1.0f/(powf(cabsf(Un_Y[i]),2.0f) + 2.23e-9f)) : 1.0f/(powf(cabsf(Un_Y[i]),2.0f) + 2.23e-9f);
    
    free(V);
    free(Vn);
    free(Vn1);
    free(Un);
    free(Un_Y);
}


/* ========================================================================== */
/*              Microphone/Hydrophone array processing functions              */
/* ========================================================================== */

void cylModalCoeffs
(
    int order,
    double* kr,
    int nBands,
    ARRAY_CONSTRUCTION_TYPES arrayType,
    double_complex* b_N
)
{
    int i, n;
    double* Jn, *Jnprime;
    double_complex* Hn2, *Hn2prime;
    
    memset(b_N, 0, nBands*(order+1)*sizeof(double_complex));
    switch(arrayType){
        default:
        case ARRAY_CONSTRUCTION_OPEN:
            /* compute spherical Bessels of the first kind */
            Jn = malloc1d(nBands*(order+1)*sizeof(double));
            bessel_Jn_ALL(order, kr, nBands, Jn, NULL);
            
            /* modal coefficients for open spherical array (omni sensors): 1i^n * jn; */
            for(n=0; n<order+1; n++)
                for(i=0; i<nBands; i++)
                    b_N[i*(order+1)+n] = crmul(cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)), Jn[i*(order+1)+n]);
            
            free(Jn);
            break;
            
        case ARRAY_CONSTRUCTION_RIGID:
            /* compute spherical Bessels/Hankels and their derivatives */
            Jn = malloc1d(nBands*(order+1)*sizeof(double));
            Jnprime = malloc1d(nBands*(order+1)*sizeof(double));
            Hn2 = malloc1d(nBands*(order+1)*sizeof(double_complex));
            Hn2prime = malloc1d(nBands*(order+1)*sizeof(double_complex));
            bessel_Jn_ALL(order, kr, nBands, Jn, Jnprime);
            hankel_Hn2_ALL(order, kr, nBands, Hn2, Hn2prime);
            
            /* modal coefficients for rigid spherical array: 1i^n * (jn-(jnprime./hn2prime).*hn2); */
            for(i=0; i<nBands; i++){
                for(n=0; n<order+1; n++){
                    if(n==0 && kr[i]<=1e-20)
                        b_N[i*(order+1)+n] = cmplx(1.0, 0.0);
                    else if(kr[i] <= 1e-20)
                        b_N[i*(order+1)+n] = cmplx(0.0, 0.0);
                    else{
                        b_N[i*(order+1)+n] = ccmul(cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)), ( ccsub(cmplx(Jn[i*(order+1)+n], 0.0),
                                             ccmul(ccdiv(cmplx(Jnprime[i*(order+1)+n],0.0), Hn2prime[i*(order+1)+n]), Hn2[i*(order+1)+n]))));
                    }
                }
            }
            
            free(Jn);
            free(Jnprime);
            free(Hn2);
            free(Hn2prime);
            break;
            
        case ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL:
        case ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL:
            saf_print_error("Unsupported array type");
            break;
    }
}

float sphArrayAliasLim
(
    float r,
    float c,
    int maxN
)
{
   return c*(float)maxN/(2.0f*SAF_PI*r);
}

void sphArrayNoiseThreshold
(
    int maxN,
    int Nsensors,
    float r,
    float c,
    ARRAY_CONSTRUCTION_TYPES arrayType,
    double dirCoeff,
    float maxG_db,
    float* f_lim
)
{
    int n;
    float kR_lim, maxG;
    double kr;
    double_complex* b_N;
    
    maxG = powf(10.0f, maxG_db/10.0f);
    kr = 1.0f;
    for (n=1; n<maxN+1; n++){
        b_N = malloc1d((n+1) * sizeof(double_complex));
        sphModalCoeffs(n, &kr, 1, arrayType, dirCoeff, b_N);
        kR_lim = powf(maxG*(float)Nsensors* powf((float)cabs(b_N[n])/(4.0f*SAF_PI), 2.0f), (-10.0f*log10f(2.0f)/(6.0f*n)));
        f_lim[n-1] = kR_lim*c/(2.0f*SAF_PI*r);
        free(b_N);
    }
}

void sphModalCoeffs
(
    int order,
    double* kr,
    int nBands,
    ARRAY_CONSTRUCTION_TYPES arrayType,
    double dirCoeff,
    double_complex* b_N
)
{
    int i, n, maxN, maxN_tmp;
    double* jn, *jnprime;
    double_complex* hn2, *hn2prime;
    
    memset(b_N, 0, nBands*(order+1)*sizeof(double_complex));
    switch(arrayType){
        default:
        case ARRAY_CONSTRUCTION_OPEN:
            /* compute spherical Bessels of the first kind */
            jn = malloc1d(nBands*(order+1)*sizeof(double));
            bessel_jn_ALL(order, kr, nBands, &maxN, jn, NULL);
            
            /* modal coefficients for open spherical array (omni sensors): 4*pi*1i^n * jn; */
            for(n=0; n<maxN+1; n++)
                for(i=0; i<nBands; i++)
                    b_N[i*(order+1)+n] = crmul(crmul(cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)), 4.0*SAF_PId), jn[i*(order+1)+n]);
            
            free(jn);
            break;
            
        case ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL:
            /* compute spherical Bessels of the first kind + derivatives */
            jn = malloc1d(nBands*(order+1)*sizeof(double));
            jnprime = malloc1d(nBands*(order+1)*sizeof(double));
            bessel_jn_ALL(order, kr, nBands, &maxN, jn, jnprime);

            /* modal coefficients for open spherical array (directional sensors): 4*pi*1i^n * (dirCoeff*jn - 1i*(1-dirCoeff)*jnprime); */
            for(n=0; n<maxN+1; n++)
                for(i=0; i<nBands; i++)
                    b_N[i*(order+1)+n] = ccmul(crmul(cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)), 4.0*SAF_PId), ccsub(cmplx(dirCoeff*jn[i*(order+1)+n], 0.0),
                                         cmplx(0.0, (1.0-dirCoeff)*jnprime[i*(order+1)+n]))  );
            
            free(jn);
            free(jnprime);
            break;
            
        case ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL:
            /* equivalent to "ARRAY_CONSTRUCTION_RIGID", as long as the sensor radius is the same as the scatterer radius. Call
               "sphScattererModalCoeffs" or "sphScattererDirModalCoeffs", if sensors protrude from the rigid baffle. */
            
        case ARRAY_CONSTRUCTION_RIGID:
            /* compute spherical Bessels/Hankels and their derivatives */
            jn = malloc1d(nBands*(order+1)*sizeof(double));
            jnprime = malloc1d(nBands*(order+1)*sizeof(double));
            hn2 = malloc1d(nBands*(order+1)*sizeof(double_complex));
            hn2prime = malloc1d(nBands*(order+1)*sizeof(double_complex));
            maxN = 1000000000;
            bessel_jn_ALL(order, kr, nBands, &maxN_tmp, jn, jnprime);
            maxN = SAF_MIN(maxN_tmp, maxN);
            hankel_hn2_ALL(order, kr, nBands, &maxN_tmp, hn2, hn2prime);
            maxN = SAF_MIN(maxN_tmp, maxN); /* maxN being the minimum highest order that was computed for all values in kr */

            /* modal coefficients for rigid spherical array: 4*pi*1i^n * (jn-(jnprime./hn2prime).*hn2); */
            for(i=0; i<nBands; i++){
                for(n=0; n<maxN+1; n++){
                    if(n==0 && kr[i]<=1e-20)
                        b_N[i*(order+1)+n] = cmplx(4.0*SAF_PId, 0.0);
                    else if(kr[i] <= 1e-20)
                        b_N[i*(order+1)+n] = cmplx(0.0, 0.0);
                    else{
                        b_N[i*(order+1)+n] = ccmul(crmul(cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)), 4.0*SAF_PId), ( ccsub(cmplx(jn[i*(order+1)+n], 0.0),
                                             ccmul(ccdiv(cmplx(jnprime[i*(order+1)+n],0.0), hn2prime[i*(order+1)+n]), hn2[i*(order+1)+n]))));
                    }
                }
            }
            
            free(jn);
            free(jnprime);
            free(hn2);
            free(hn2prime);
            break;
    }
}

void sphScattererModalCoeffs
(
    int order,
    double* kr,
    double* kR,
    int nBands,
    double_complex* b_N
)
{
    int i, n, maxN, maxN_tmp;
    double* jn, *jnprime;
    double_complex* hn2, *hn2prime;
    
    /* compute spherical Bessels/Hankels and their derivatives */
    jn = malloc1d(nBands*(order+1)*sizeof(double));
    jnprime = malloc1d(nBands*(order+1)*sizeof(double));
    hn2 = malloc1d(nBands*(order+1)*sizeof(double_complex));
    hn2prime = malloc1d(nBands*(order+1)*sizeof(double_complex));
    maxN = 1000000000;
    bessel_jn_ALL(order, kr, nBands, &maxN_tmp, jn, NULL);
    maxN = SAF_MIN(maxN_tmp, maxN);
    bessel_jn_ALL(order, kR, nBands, &maxN_tmp, NULL, jnprime);
    maxN = SAF_MIN(maxN_tmp, maxN);
    hankel_hn2_ALL(order, kr, nBands, &maxN_tmp, hn2, NULL);
    maxN = SAF_MIN(maxN_tmp, maxN);
    hankel_hn2_ALL(order, kR, nBands, &maxN_tmp, NULL, hn2prime);
    maxN = SAF_MIN(maxN_tmp, maxN); /* maxN being the minimum highest order that was computed for all values in kr */
    
    /* modal coefficients for rigid spherical array (OMNI): 4*pi*1i^n * (jn_kr-(jnprime_kr./hn2prime_kr).*hn2_kr); */
    /* modal coefficients for rigid spherical scatterer (OMNI): 4*pi*1i^n * (jn_kr-(jnprime_kR./hn2prime_kR).*hn2_kr); */
    for(i=0; i<nBands; i++){
        for(n=0; n<maxN+1; n++){
            if(n==0 && kr[i]<=1e-20)
                b_N[i*(order+1)+n] = cmplx(4.0*SAF_PId, 0.0);
            else if(kr[i] <= 1e-20)
                b_N[i*(order+1)+n] = cmplx(0.0, 0.0);
            else{
                b_N[i*(order+1)+n] = ccmul(crmul(cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)), 4.0*SAF_PId), ( ccsub(cmplx(jn[i*(order+1)+n], 0.0),
                                     ccmul(ccdiv(cmplx(jnprime[i*(order+1)+n],0.0), hn2prime[i*(order+1)+n]), hn2[i*(order+1)+n]))));
            }
        }
    }
    
    free(jn);
    free(jnprime);
    free(hn2);
    free(hn2prime);
}

void sphScattererDirModalCoeffs
(
    int order,
    double* kr,
    double* kR,
    int nBands,
    double dirCoeff, /* 0.0 gives NaNs */
    double_complex* b_N
)
{
    int i, n, maxN, maxN_tmp;
    double* jn_kr, *jnprime_kr, *jnprime_kR;
    double_complex* hn2_kr, *hn2prime_kr, *hn2prime_kR;
    
    /* compute spherical Bessels/Hankels and their derivatives */
    jn_kr = malloc1d(nBands*(order+1)*sizeof(double));
    jnprime_kr = malloc1d(nBands*(order+1)*sizeof(double));
    jnprime_kR = malloc1d(nBands*(order+1)*sizeof(double));
    hn2_kr = malloc1d(nBands*(order+1)*sizeof(double_complex));
    hn2prime_kr = malloc1d(nBands*(order+1)*sizeof(double_complex));
    hn2prime_kR = malloc1d(nBands*(order+1)*sizeof(double_complex));
    maxN = 1000000000;
    bessel_jn_ALL(order, kr, nBands, &maxN_tmp, jn_kr, jnprime_kr);
    maxN = SAF_MIN(maxN_tmp, maxN);
    bessel_jn_ALL(order, kR, nBands, &maxN_tmp, NULL, jnprime_kR);
    maxN = SAF_MIN(maxN_tmp, maxN);
    hankel_hn2_ALL(order, kr, nBands, &maxN_tmp, hn2_kr, hn2prime_kr);
    maxN = SAF_MIN(maxN_tmp, maxN);
    hankel_hn2_ALL(order, kR, nBands, &maxN_tmp, NULL, hn2prime_kR);
    maxN = SAF_MIN(maxN_tmp, maxN); /* maxN being the minimum highest order that was computed for all values in kr */
    
    /* modal coefficients for rigid spherical array (OMNI): 4*pi*1i^n * (jn_kr-(jnprime_kr./hn2prime_kr).*hn2_kr); */
    /* modal coefficients for rigid spherical scatterer (OMNI): 4*pi*1i^n * (jn_kr-(jnprime_kR./hn2prime_kR).*hn2_kr); */
    /* modal coefficients for rigid spherical scatterer (DIRECTIONAL):
           4*pi*1i^n * [ (beta*jn_kr - i(1-beta)*jnprime_kr) - (jnprime_kR/hn2prime_kR) * (beta*hn2_kr - i(1-beta)hn2prime_kr) ] */
    for(i=0; i<nBands; i++){
        for(n=0; n<maxN+1; n++){
            if(n==0 && kr[i]<=1e-20)
                b_N[i*(order+1)+n] = cmplx(4.0*SAF_PId, 0.0);
            else if(kr[i] <= 1e-20)
                b_N[i*(order+1)+n] = cmplx(0.0, 0.0);
            else{ // Dear god, what happened here?!...
//#if __STDC_VERSION__ >= 199901L
//                b_N[i*(order+1)+n] = 4.0f*PI*cpowf(I,(float)n) * ( (dirCoeff*jn_kr[i*(order+1)+n] - I*(1.0f-dirCoeff)*jnprime_kr[i*(order+1)+n]) -
//                                                                   (jnprime_kR[i*(order+1)+n]/hn2prime_kR[i*(order+1)+n]) * (dirCoeff*hn2_kr[i*(order+1)+n] -
//                                                                     I*(1.0f-dirCoeff)*hn2prime_kr[i*(order+1)+n]) );
//#else
                b_N[i*(order+1)+n] = cmplx(dirCoeff * jn_kr[i*(order+1)+n], -(1.0-dirCoeff)* jnprime_kr[i*(order+1)+n]);
                b_N[i*(order+1)+n] = ccsub(b_N[i*(order+1)+n], ccmul(ccdiv(cmplx(jnprime_kR[i*(order+1)+n], 0.0), hn2prime_kR[i*(order+1)+n]),
                                    (ccsub(crmul(hn2_kr[i*(order+1)+n], dirCoeff), ccmul(cmplx(0.0f,1.0-dirCoeff), hn2prime_kr[i*(order+1)+n])))));
                b_N[i*(order+1)+n] = crmul(ccmul(cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)), b_N[i*(order+1)+n]), 4.0*SAF_PId/dirCoeff); /* had to scale by directivity to preserve amplitude? */ 
//                b_N[i*(order+1)+n] = dirCoeff * jn_kr[i*(order+1)+n] - I*(1.0-dirCoeff)* jnprime_kr[i*(order+1)+n];
//                b_N[i*(order+1)+n] = b_N[i*(order+1)+n] - (jnprime_kR[i*(order+1)+n]/hn2prime_kR[i*(order+1)+n])*(dirCoeff*hn2_kr[i*(order+1)+n] - I*(1.0-dirCoeff)*hn2prime_kr[i*(order+1)+n]);
//                b_N[i*(order+1)+n] = cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)) * b_N[i*(order+1)+n] * 4.0*SAF_PId/dirCoeff; /* had to scale by directivity to preserve amplitude? */
//#endif
            }
        }
    }
    
    free(jn_kr);
    free(jnprime_kr);
    free(jnprime_kR);
    free(hn2_kr);
    free(hn2prime_kr);
    free(hn2prime_kR);
}

void sphDiffCohMtxTheory
(
    int order,
    float* sensor_dirs_rad,
    int N_sensors,
    ARRAY_CONSTRUCTION_TYPES arrayType,
    double dirCoeff,
    double* kr,
    int nBands,
    double* M_diffcoh
)
{
    int i, j, k, n;
    float cosangle;
    float* sensor_dirs_xyz, *ppm, *ppm_z1, *ppm_z2;
    double* b_N2, *Pn;
    double_complex* b_N;
    
    /* sph->cart */
    sensor_dirs_xyz = malloc1d(N_sensors*3*sizeof(float));
    for(i=0; i<N_sensors; i++){
        sensor_dirs_xyz[i*3] = cosf(sensor_dirs_rad[i*2+1]) * cosf(sensor_dirs_rad[i*2]);
        sensor_dirs_xyz[i*3+1] = cosf(sensor_dirs_rad[i*2+1]) * sinf(sensor_dirs_rad[i*2]);
        sensor_dirs_xyz[i*3+2] = sinf(sensor_dirs_rad[i*2+1]);
    }
    
    /* calculate modal coefficients */
    b_N = malloc1d(nBands * (order+1) * sizeof(double_complex));
    b_N2 = malloc1d(nBands * (order+1) * sizeof(double));
    switch (arrayType){
        case ARRAY_CONSTRUCTION_OPEN:
            sphModalCoeffs(order, kr, nBands, ARRAY_CONSTRUCTION_OPEN, 1.0, b_N); break;
        case ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL:
            sphModalCoeffs(order, kr, nBands, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, dirCoeff, b_N); break;
        case ARRAY_CONSTRUCTION_RIGID: /* fall through */
        case ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL:
            sphModalCoeffs(order, kr, nBands, ARRAY_CONSTRUCTION_RIGID, 1.0, b_N);
            break;
    }
    for(i=0; i<nBands * (order+1); i++)
        b_N2[i] = pow(cabs(ccdiv(b_N[i], cmplx(4.0*SAF_PId, 0.0))), 2.0);
    
    /* determine theoretical diffuse-coherence matrix for sensor array */
    ppm = malloc1d((order+1)*sizeof(float));
    ppm_z1 = malloc1d((order+1)*sizeof(float));
    ppm_z2 = malloc1d((order+1)*sizeof(float));
    Pn = malloc1d((order+1)*sizeof(double));
    for(i=0; i<N_sensors; i++){
        for(j=i; j<N_sensors; j++){
            cosangle = 0.0f;
            for(k=0; k<3; k++)
                cosangle += sensor_dirs_xyz[j*3+k] * sensor_dirs_xyz[i*3+k];
            cosangle = cosangle>1.0f ? 1.0f : (cosangle<-1.0f ? -1.0f : cosangle);
            for(n=0; n<order+1; n++){
                unnorm_legendreP_recur(n, &cosangle, 1, ppm_z1, ppm_z2, ppm);  
                Pn[n] =  (2.0*(double)n+1.0) * 4.0f*SAF_PI * (double)ppm[0];
                memcpy(ppm_z2, ppm_z1, (order+1)*sizeof(float));
                memcpy(ppm_z1, ppm, (order+1)*sizeof(float));
            }
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nBands, 1, order+1, 1.0,
                        b_N2, order+1,
                        Pn, 1, 0.0,
                        &M_diffcoh[j*N_sensors*nBands + i*nBands], 1);

            memcpy(&M_diffcoh[i*N_sensors*nBands + j*nBands], &M_diffcoh[j*N_sensors*nBands + i*nBands], nBands*sizeof(double));
        }
    }
    
    free(b_N);
    free(b_N2);
    free(sensor_dirs_xyz);
    free(ppm);
    free(ppm_z1);
    free(ppm_z2);
    free(Pn);
}

void simulateCylArray /*untested*/
(
    int order,
    double* kr,
    int nBands,
    float* sensor_dirs_rad,
    int N_sensors,
    float* src_dirs_deg,
    int N_srcs,
    ARRAY_CONSTRUCTION_TYPES arrayType,
    float_complex* H_array
)
{
    int i, j, n, band;
    double angle;
    double_complex* b_N, *C, *b_NC;
    const double_complex calpha = cmplx(1.0, 0.0), cbeta = cmplx(0.0, 0.0);
    
    /* calculate modal coefficients */
    b_N = malloc1d(nBands * (order+1) * sizeof(double_complex));
    cylModalCoeffs(order, kr, nBands, arrayType, b_N);
    
    /* Compute angular-dependent part of the array responses */
    C = malloc1d((order+1)*N_sensors*sizeof(double_complex));
    b_NC = malloc1d(nBands*N_sensors*sizeof(double_complex));
    for(i=0; i<N_srcs; i++){
        for(j=0; j<N_sensors; j++){
            angle = sensor_dirs_rad[i*2] - src_dirs_deg[i*2]*SAF_PId/180.0;
            for(n=0; n<order+1; n++){
                /* Jacobi-Anger expansion */
                if(n==0)
                    C[n*N_sensors+j] = cmplx(1.0, 0.0);
                else
                    C[n*N_sensors+j] = cmplx(2.0*cos((double)n*angle), 0.0);
            }
        }
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nBands, N_sensors, order+1, &calpha,
                    b_N, order+1,
                    C, N_sensors, &cbeta,
                    b_NC, N_sensors);
        
        /* store array response per frequency, sensors and plane wave dirs */
        for(band=0; band<nBands; band++)
            for(j=0; j<N_sensors; j++)
                H_array[band*N_sensors*N_srcs + j*N_srcs + i] = cmplxf((float)creal(b_NC[band*N_sensors+j]), (float)cimag(b_NC[band*N_sensors+j]));
    }
    
    free(b_N);
    free(C);
    free(b_NC);
}

void simulateSphArray
(
    int order,
    double* kr,
    double* kR,
    int nBands,
    float* sensor_dirs_rad,
    int N_sensors,
    float* src_dirs_deg,
    int N_srcs,
    ARRAY_CONSTRUCTION_TYPES arrayType,
    double dirCoeff,
    float_complex* H_array
)
{
    int i, j, n, band;
    double dcosangle;
    double *ppm;
    float cosangle;
    float* U_sensors, *U_srcs;
    double_complex* b_N, *P, *b_NP;
    const double_complex calpha = cmplx(1.0, 0.0), cbeta = cmplx(0.0, 0.0);
    
    /* calculate modal coefficients */
    b_N = malloc1d(nBands * (order+1) * sizeof(double_complex));
    switch (arrayType){
        case ARRAY_CONSTRUCTION_OPEN:
            sphModalCoeffs(order, kr, nBands, ARRAY_CONSTRUCTION_OPEN, 1.0, b_N); break;
        case ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL:
            sphModalCoeffs(order, kr, nBands, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, dirCoeff, b_N); break;
        case ARRAY_CONSTRUCTION_RIGID: /* fall through */
        case ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL:
            if(kR==NULL)
                sphModalCoeffs(order, kr, nBands, ARRAY_CONSTRUCTION_RIGID, 1.0, b_N); /* if kr==kR, dirCoeff is irrelevant */
            else
                sphScattererDirModalCoeffs(order, kr, kR, nBands, dirCoeff, b_N);
            break;
    }
    
    /* calculate (unit) cartesian coords for sensors and plane waves */
    U_sensors = malloc1d(N_sensors*3*sizeof(float));
    U_srcs = malloc1d(N_srcs*3*sizeof(float));
    unitSph2cart(sensor_dirs_rad, N_sensors, 0, U_sensors);
    unitSph2cart(src_dirs_deg, N_srcs, 1, U_srcs); 
    
    /* Compute angular-dependent part of the array responses */
    ppm = malloc1d((order+1)*sizeof(double));
    P = malloc1d((order+1)*N_sensors*sizeof(double_complex));
    b_NP = malloc1d(nBands*N_sensors*sizeof(double_complex));
    for(i=0; i<N_srcs; i++){
        for(j=0; j<N_sensors; j++){
            utility_svvdot((const float*)&U_sensors[j*3], (const float*)&U_srcs[i*3], 3, &cosangle);
            for(n=0; n<order+1; n++){
                /* Legendre polynomials correspond to the angular dependency */
                dcosangle = (double)cosangle;
                unnorm_legendreP(n, &dcosangle, 1, ppm);
                P[n*N_sensors+j] = cmplx((2.0*(double)n+1.0)/(4.0*SAF_PId) * ppm[0], 0.0);
            }
        }
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nBands, N_sensors, order+1, &calpha,
                    b_N, order+1,
                    P, N_sensors, &cbeta,
                    b_NP, N_sensors);
        
        /* store array response per frequency, sensors and plane wave dirs */
        for(band=0; band<nBands; band++)
            for(j=0; j<N_sensors; j++)
                H_array[band*N_sensors*N_srcs + j*N_srcs + i] = cmplxf((float)creal(b_NP[band*N_sensors+j]), (float)cimag(b_NP[band*N_sensors+j]));
    }
    
    free(U_sensors);
    free(U_srcs);
    free(b_N);
    free(ppm);
    free(P);
    free(b_NP);
}

void evaluateSHTfilters
(
    int order,
    float_complex* M_array2SH,
    int nSensors,
    int nBands,
    float_complex* H_array,
    int nDirs,
    float_complex* Y_grid,
    float* cSH,
    float* lSH
#if 0
    , float* WNG
#endif
)
{
    int band, i, n, m, nSH, q;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    float w_uni_grid, lSH_n, lSH_nm;
    float_complex cSH_n, cSH_nm, yre_yre_dot, yre_yid_dot;
    float_complex *y_recon_kk, *y_recon_nm, *w_y_recon_nm, *y_ideal_nm, *MH_M, *EigV;
    
    nSH = ORDER2NSH(order);
    w_uni_grid = 1.0f/(float)nDirs;
    y_recon_kk = malloc1d(nSH*nDirs*sizeof(float_complex));
    y_recon_nm = malloc1d(nDirs*sizeof(float_complex));
    w_y_recon_nm = malloc1d(nDirs*sizeof(float_complex));
    y_ideal_nm = malloc1d(nDirs*sizeof(float_complex));
    MH_M = malloc1d(nSensors*nSensors*sizeof(float_complex));
    EigV = malloc1d(nSensors*nSensors*sizeof(float_complex));
    for(band=0; band<nBands; band++){
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, nDirs, nSensors, &calpha,
                    &M_array2SH[band*nSH*nSensors], nSensors,
                    &H_array[band*nSensors*nDirs], nDirs, &cbeta,
                    y_recon_kk, nDirs);
        for(n=0; n<order+1; n++){
            cSH_n = cmplxf(0.0f, 0.0f);
            lSH_n = 0.0f;
            for(m=-n; m<=n; m++){
                q = n*n+n+m;
                for(i=0; i<nDirs; i++){
                    y_recon_nm[i] = y_recon_kk[q*nDirs+i];
                    w_y_recon_nm[i] = crmulf(y_recon_nm[i], w_uni_grid);
                    y_ideal_nm[i] = Y_grid[q*nDirs+i];
                }
                utility_cvvdot(w_y_recon_nm, y_recon_nm, nDirs, CONJ, &yre_yre_dot);
                utility_cvvdot(w_y_recon_nm, y_ideal_nm, nDirs, CONJ, &yre_yid_dot);
                cSH_nm = ccdivf(yre_yid_dot, ccaddf(csqrtf(yre_yre_dot), cmplxf(2.23e-9f, 0.0f)));
                cSH_n = ccaddf(cSH_n, cSH_nm);
                lSH_nm = crealf(yre_yre_dot);
                lSH_n += lSH_nm;
            }
            cSH[band*(order+1)+n] = SAF_MAX(SAF_MIN(cabsf(cSH_n)/(2.0f*(float)n+1.0f),1.0f),0.0f);
            lSH[band*(order+1)+n] = 10.0f*log10f(lSH_n/(2.0f*(float)n+1.0f)+2.23e-9f);
        } 
    }
    
#if 0
    /* find maximum noise amplification of all filters in matrix */
    for(band=0; band<nBands; band++){
        cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, nSensors, nSensors, nSH, &calpha,
                    &M_array2SH[band*nSH*nSensors], nSensors,
                    &M_array2SH[band*nSH*nSensors], nSensors, &cbeta,
                    MH_M, nSensors);
        utility_ceig(NULL, MH_M, nSensors, 1, NULL, NULL, EigV); /* eigenvalues in decending order */
        WNG[band] = 10.0f*log10f(crealf(EigV[0])+2.23e-9f);
    }
#endif
    
    free(y_recon_kk);
    free(y_recon_nm);
    free(w_y_recon_nm);
    free(y_ideal_nm);
    free(MH_M);
    free(EigV);
}

