/*
 * Copyright 2018 Leo McCormack
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
 * @file saf_hoa_internal.c
 * @ingroup HOA
 * @brief Internal source for the higher-order Ambisonics module
 *        (#SAF_HOA_MODULE)
 *
 * A collection of Ambisonics related functions. Many of which are derived from
 * the MATLAB library found in [1].
 *
 * @see [1] https://github.com/polarch/Higher-Order-Ambisonics
 *          Copyright (c) 2015, Archontis Politis, BSD-3-Clause License
 *
 * @author Leo McCormack
 * @date 19.03.2018
 * @license ISC
 */
 
#include "../saf_hoa/saf_hoa.h"
#include "../saf_hoa/saf_hoa_internal.h"

/* ========================================================================== */
/*                       Loudspeaker Ambisonic Decoders                       */
/* ========================================================================== */

void getEPAD
(
    int order,
    float* ls_dirs_deg,
    int nLS,
    float* decMtx
)
{
    int i, j, nSH;
    float scale;
    float* Y_ls, *U, *V, *U_tr, *V_tr;
    
    nSH = ORDER2NSH(order);
    scale = 1.0f/SQRT4PI;

    /* Prep */
    Y_ls = malloc1d(nSH*nLS*sizeof(float));
    U = malloc1d(nSH*nSH*sizeof(float));
    V = malloc1d(nLS*nLS*sizeof(float));
    getRSH(order, ls_dirs_deg, nLS, Y_ls);
    cblas_sscal(nLS*nSH, scale, Y_ls, 1);
    utility_ssvd(NULL, Y_ls, nSH, nLS, U, NULL, V, NULL);

    /* Apply truncation */
    if(nSH>nLS){
        /* truncate the U matrix */
        U_tr = malloc1d(nSH*nLS*sizeof(float));
        for(i=0; i<nSH; i++)
            for(j=0; j<nLS; j++)
                U_tr[i*nLS+j] = U[i*nSH+j];
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nLS, nSH, nLS, 1.0f,
                    V, nLS,
                    U_tr, nLS, 0.0f,
                    decMtx, nSH);
        free(U_tr);
    }
    else{
        /* truncate the V matrix (NOT V^T!) */
        V_tr = malloc1d(nLS*nSH*sizeof(float));
        for(i=0; i<nLS; i++)
            for(j=0; j<nSH; j++)
                V_tr[i*nSH+j] = V[i*nLS+j];
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nLS, nSH, nSH, 1.0f,
                    V_tr, nSH,
                    U, nSH, 0.0f,
                    decMtx, nSH);
        free(V_tr);
    }

    /* Apply normalisation, and scale by number of loudspeakers */
    scale = sqrtf(4.0f*SAF_PI/(float)nLS);
    utility_svsmul(decMtx, &scale, nLS*nSH, decMtx);

    /* clean-up */
    free(U);
    free(V);
    free(Y_ls);
}
 
void getAllRAD
(
    int order,
    float* ls_dirs_deg,
    int nLS,
    float* decMtx
)
{
    int nDirs_td, N_gtable, nGroups, nSH;
    float scale;
    float* Y_td, *G_td, *t_dirs;
    
    nSH = ORDER2NSH(order);
    scale = 1.0f/SQRT4PI;

#if 0
    /* define a sufficiently dense t-design for this decoding order, as to
     * conserve omni energy */
    t = 4*order;
    G_td = NULL;
    if(t<=21){
        /* suitable for up to 5th order */
        nDirs_td = __Tdesign_nPoints_per_degree[t-1];
        t_dirs = (float*)__HANDLES_Tdesign_dirs_deg[t-1];
    }
    else if (order > 7){
        /* suitable for >7th order */
        nDirs_td = 5100; /* Minimum t-design of degree 100 has 5100 points */
        t_dirs = (float*)__Tdesign_degree_100_dirs_deg;
    }
    else{
        /* suitable for 6th & 7th order */
        nDirs_td = 480; /* Minimum t-design of degree 30 has 480 points (sufficient for up to 7th order) */
        t_dirs = (float*)__Tdesign_degree_30_dirs_deg;
    }
#else
    nDirs_td = 5100; /* Minimum t-design of degree 100 has 5100 points */
    t_dirs = (float*)__Tdesign_degree_100_dirs_deg;
#endif
    
    /* calculate vbap gains and SH matrix for this t-design */ 
    generateVBAPgainTable3D_srcs(t_dirs, nDirs_td, ls_dirs_deg, nLS, 0, 0, 0.0f, &G_td, &N_gtable, &nGroups);
    Y_td = malloc1d(nSH*nDirs_td*sizeof(float));
    getRSH(order, t_dirs, nDirs_td, Y_td);
    cblas_sscal(nDirs_td*nSH, scale, Y_td, 1);
    
    /* AllRAD decoder is simply (G_td * T_td * 1/nDirs_td) */
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasTrans, nLS, nSH, nDirs_td, 1.0f,
                G_td, nLS,
                Y_td, nDirs_td, 0.0f,
                decMtx, nSH);
    cblas_sscal(nLS*nSH, (4.0f*SAF_PI)/(float)nDirs_td, decMtx, 1);

    free(Y_td);
    free(G_td);
}


/* ========================================================================== */
/*                         Binaural Ambisonic Decoders                        */
/* ========================================================================== */

void getBinDecoder_LS
(
    float_complex* hrtfs,  /* the HRTFs; FLAT: N_bands x 2 x N_dirs */
    float* hrtf_dirs_deg,
    int N_dirs,
    int N_bands,
    int order,
    float* weights,
    float_complex* decMtx /* N_bands x 2 x (order+1)^2  */
)
{
    int i, j, nSH, band;
    float* Y_tmp;
    float_complex* W, *Y_na, *Yna_W, *Yna_W_Yna, *Yna_W_H, *B;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH = ORDER2NSH(order);
    
    /* SH */
    Y_tmp = malloc1d(nSH*N_dirs*sizeof(float));
    Y_na = malloc1d(nSH*N_dirs*sizeof(float_complex));
    B = malloc1d(nSH * 2 * sizeof(float_complex));
    getRSH(order, hrtf_dirs_deg, N_dirs, Y_tmp);
    for(i=0; i<nSH*N_dirs; i++)
        Y_na[i] = cmplxf(Y_tmp[i], 0.0f);
    free(Y_tmp);
    
    /* compute decoding matrix, incorporating integration weights */
    W = calloc1d(N_dirs*N_dirs, sizeof(float_complex));
    if(weights!=NULL)
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(weights[i], 0.0f);
    else
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(1.0f/(float)N_dirs, 0.0f);

    /* calculate decoding matrix per band */
    Yna_W = malloc1d(nSH * N_dirs*sizeof(float_complex));
    Yna_W_Yna = malloc1d(nSH * nSH * sizeof(float_complex));
    Yna_W_H = malloc1d(nSH * 2 * sizeof(float_complex));
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, N_dirs, N_dirs, &calpha,
                Y_na, N_dirs,
                W, N_dirs, &cbeta,
                Yna_W, N_dirs);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, nSH, N_dirs, &calpha,
                Yna_W, N_dirs,
                Y_na, N_dirs, &cbeta,
                Yna_W_Yna, nSH);
    for(band=0; band<N_bands; band++){
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, 2, N_dirs, &calpha,
                    Yna_W, N_dirs,
                    &hrtfs[band*2*N_dirs], N_dirs, &cbeta,
                    Yna_W_H, 2);
        utility_cglslv(NULL, Yna_W_Yna, nSH, Yna_W_H, 2, B);
        for(i=0; i<nSH; i++)
            for(j=0; j<2; j++)
                decMtx[band*2*nSH + j*nSH + i] = conjf(B[i*2+j]); /* ^H */
    }
    
    /* clean-up */
    free(W);
    free(Yna_W);
    free(Yna_W_Yna);
    free(Yna_W_H);
    free(Y_na);
    free(B);
}

void getBinDecoder_LSDIFFEQ
(
    float_complex* hrtfs,  /* the HRTFs; FLAT: N_bands x 2 x N_dirs */
    float* hrtf_dirs_deg,
    int N_dirs,
    int N_bands,
    int order,
    float* weights,
    float_complex* decMtx /* N_bands x 2 x (order+1)^2  */
)
{
    int i, j, nSH, band;
    float Gh;
    float* Y_tmp;
    float_complex* W, *Y_na, *Yna_W, *Yna_W_Yna, *Yna_W_H, *B_ls, *hrtfs_ls, *H_W;
    float_complex C_ref[2][2], C_ls[2][2];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH = ORDER2NSH(order);
    
    /* integration weights */
    W = calloc1d(N_dirs*N_dirs, sizeof(float_complex));
    if(weights!=NULL)
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(weights[i], 0.0f);
    else
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(1.0f/(float)N_dirs, 0.0f);
    
    /* SH */
    Y_tmp = malloc1d(nSH*N_dirs*sizeof(float));
    Y_na = malloc1d(nSH*N_dirs*sizeof(float_complex));
    getRSH(order, hrtf_dirs_deg, N_dirs, Y_tmp);
    for(i=0; i<nSH*N_dirs; i++)
        Y_na[i] = cmplxf(Y_tmp[i], 0.0f);
    free(Y_tmp);
    
    /* calculate decoding matrix per band */
    Yna_W = malloc1d(nSH * N_dirs*sizeof(float_complex));
    Yna_W_Yna = malloc1d(nSH * nSH * sizeof(float_complex));
    Yna_W_H = malloc1d(nSH * 2 * sizeof(float_complex));
    B_ls = malloc1d(nSH * 2 * sizeof(float_complex));
    hrtfs_ls = malloc1d(2*N_dirs*sizeof(float_complex));
    H_W = malloc1d(2*N_dirs*sizeof(float_complex));
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, N_dirs, N_dirs, &calpha,
                Y_na, N_dirs,
                W, N_dirs, &cbeta,
                Yna_W, N_dirs);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, nSH, N_dirs, &calpha,
                Yna_W, N_dirs,
                Y_na, N_dirs, &cbeta,
                Yna_W_Yna, nSH);
    for(band=0; band<N_bands; band++){
        /* find least-squares decoding matrix */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, 2, N_dirs, &calpha,
                    Yna_W, N_dirs,
                    &hrtfs[band*2*N_dirs], N_dirs, &cbeta,
                    Yna_W_H, 2);
        utility_cglslv(NULL, Yna_W_Yna, nSH, Yna_W_H, 2, B_ls);
        cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, 2, N_dirs, nSH, &calpha,
                    B_ls, 2,
                    Y_na, N_dirs, &cbeta,
                    hrtfs_ls, N_dirs);
        
        /* Diffuse-field responses */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, N_dirs, N_dirs, &calpha,
                    &hrtfs[band*2*N_dirs], N_dirs,
                    W, N_dirs, &cbeta,
                    H_W, N_dirs);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, 2, 2, N_dirs, &calpha,
                    H_W, N_dirs,
                    &hrtfs[band*2*N_dirs], N_dirs, &cbeta,
                    (float_complex*)C_ref, 2);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, N_dirs, N_dirs, &calpha,
                    hrtfs_ls, N_dirs,
                    W, N_dirs, &cbeta,
                    H_W, N_dirs);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, 2, 2, N_dirs, &calpha,
                    H_W, N_dirs,
                    hrtfs_ls, N_dirs, &cbeta,
                    (float_complex*)C_ls, 2);
        
        /* Diffuse-Equalisation factor */
        Gh = (sqrtf(crealf(C_ref[0][0])/(crealf(C_ls[0][0])+2.23e-7f)) +
              sqrtf(crealf(C_ref[1][1])/(crealf(C_ls[1][1])+2.23e-7f))) /2.0f;
        
        /* apply diff-EQ */
        for(i=0; i<nSH; i++)
            for(j=0; j<2; j++)
                decMtx[band*2*nSH + j*nSH + i] = crmulf(conjf(B_ls[i*2+j]), Gh); /* ^H */
    }
    
    free(W);
    free(Y_na);
    free(Yna_W);
    free(Yna_W_Yna);
    free(Yna_W_H);
    free(B_ls);
    free(hrtfs_ls);
    free(H_W);
}

void getBinDecoder_SPR
(
    float_complex* hrtfs, /* the HRTFs; FLAT: N_bands x 2 x N_dirs */
    float* hrtf_dirs_deg,
    int N_dirs,
    int N_bands,
    int order,
    float* weights,
    float_complex* decMtx /* N_bands x 2 x (order+1)^2  */
)
{
    int band, i, j, nSH, nSH_nh, Nh_max, Nh, K_td;
    float* hrtf_dirs_rad, *W, *cnd_num, *Y_nh, *Y_na, *tdirs_deg, *Y_td, *Ynh_Ytd, *tmp;
    float_complex* Y_td_cmplx, *W_Ynh_Ytd, *hrtfs_td, *B;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH = ORDER2NSH(order);
    
    /* integration weights */
    W = calloc1d(N_dirs*N_dirs, sizeof(float));
    if(weights!=NULL)
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = weights[i]/(4.0f*SAF_PI);
    else
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = 1.0f/(float)N_dirs;
    
    /* find SH-order for interpolation of the HRTF set */
    Nh_max = SAF_MIN((int)(sqrtf((float)N_dirs)-1.0f), 20); /* Cap to something more sensible if needed... */
    hrtf_dirs_rad = malloc1d(N_dirs*2*sizeof(float));
    cnd_num = malloc1d((Nh_max+1)*sizeof(float));
    for(i=0; i<N_dirs; i++){ /* [azi, elev] degrees, to: [azi, inclination] radians */
        hrtf_dirs_rad[i*2]   = hrtf_dirs_deg[i*2]*(SAF_PI/180.0f);
        hrtf_dirs_rad[i*2+1] = SAF_PI/2.0f - hrtf_dirs_deg[i*2+1]*(SAF_PI/180.0f);
    }
    checkCondNumberSHTReal(Nh_max, hrtf_dirs_rad, N_dirs, weights, cnd_num);
    for(i=0, Nh=0; i<Nh_max+1; i++)
        Nh = cnd_num[i] < 100.0f ? i : Nh;
    saf_assert(Nh>=order, "Input order exceeds the modal order of the spatial grid");
    nSH_nh = (Nh+1)*(Nh+1);
    Y_nh = malloc1d(nSH_nh*N_dirs*sizeof(float));
    getRSH(Nh, hrtf_dirs_deg, N_dirs, Y_nh);
    Y_na = malloc1d(nSH * N_dirs *sizeof(float));
    for(i=0; i<nSH; i++)
        for(j=0; j<N_dirs; j++)
            Y_na[i*N_dirs+j] = Y_nh[i*N_dirs+j]; 
    
    /* Get t-design SH for ambisonic signals */
    tdirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[2*order-1];
    K_td = __Tdesign_nPoints_per_degree[2*order-1];
    Y_td = malloc1d(nSH_nh*K_td*sizeof(float));
    getRSH(Nh, tdirs_deg, K_td, Y_td);
    Y_td_cmplx = malloc1d(nSH_nh*K_td*sizeof(float_complex));
    for(i=0; i<nSH_nh*K_td; i++)
        Y_td_cmplx[i] = cmplxf(Y_td[i], 0.0f);
    
    /* calculate decoding matrix per band */
    Ynh_Ytd = malloc1d(N_dirs * K_td * sizeof(float));
    tmp = malloc1d(N_dirs * K_td * sizeof(float));
    W_Ynh_Ytd = malloc1d(N_dirs * K_td * sizeof(float_complex));
    hrtfs_td = malloc1d(2*K_td*sizeof(float_complex));
    B = malloc1d(nSH * 2 * sizeof(float_complex));
    for(band=0; band<N_bands; band++){
        cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N_dirs, K_td, nSH_nh, 1.0f,
                    Y_nh, N_dirs,
                    Y_td, K_td, 0.0f,
                    Ynh_Ytd, K_td);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N_dirs, K_td, N_dirs, 1.0f,
                    W, N_dirs,
                    Ynh_Ytd, K_td, 0.0f,
                    tmp, K_td);
        for(i=0; i<N_dirs*K_td; i++)
            W_Ynh_Ytd[i] = cmplxf(tmp[i], 0.0f);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, K_td, N_dirs, &calpha,
                    &hrtfs[band*2*N_dirs], N_dirs,
                    W_Ynh_Ytd, K_td, &cbeta,
                    hrtfs_td, K_td);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, 2, K_td, &calpha,
                    Y_td_cmplx, K_td,
                    hrtfs_td, K_td, &cbeta,
                    B, 2);
        for(i=0; i<nSH; i++)
            for(j=0; j<2; j++)
                decMtx[band*2*nSH + j*nSH + i] = crmulf(conjf(B[i*2+j]), 1.0f/(float)K_td); /* ^H */
    }
    
    free(hrtf_dirs_rad);
    free(W);
    free(cnd_num);
    free(Y_nh);
    free(Y_na);
    free(Y_td);
    free(Y_td_cmplx);
    free(Ynh_Ytd);
    free(tmp);
    free(W_Ynh_Ytd);
    free(hrtfs_td);
    free(B);
}

void getBinDecoder_TA
(
    float_complex* hrtfs, /* the HRTFs; FLAT: N_bands x 2 x N_dirs */
    float* hrtf_dirs_deg,
    int N_dirs,
    int N_bands,
    int order,
    float* freqVector,
    float* itd_s,
    float* weights,
    float_complex* decMtx /* N_bands x 2 x (order+1)^2  */
)
{
    int i, j, nSH, band, band_cutoff;
    float cutoff, minVal;
    float* Y_tmp;
    float_complex* W, *Y_na, *hrtfs_mod, *Yna_W, *Yna_W_Yna, *Yna_W_H, *B;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH = ORDER2NSH(order);
    
    /* integration weights */
    W = calloc1d(N_dirs*N_dirs, sizeof(float_complex));
    if(weights!=NULL)
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(weights[i], 0.0f);
    else
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(1.0f/(float)N_dirs, 0.0f);
    
    /* SH */
    Y_tmp = malloc1d(nSH*N_dirs*sizeof(float));
    Y_na = malloc1d(nSH*N_dirs*sizeof(float_complex));
    getRSH(order, hrtf_dirs_deg, N_dirs, Y_tmp);
    for(i=0; i<nSH*N_dirs; i++)
        Y_na[i] = cmplxf(Y_tmp[i], 0.0f);
    free(Y_tmp);
    
    /* find band index for cutoff frequency */
    cutoff = 1.5e3f;
    minVal = 2.23e10f;
    for(band=0, band_cutoff=0; band<N_bands; band++){
        if(minVal>fabsf(freqVector[band]-cutoff)){
            minVal = fabsf(freqVector[band]-cutoff);
            band_cutoff = band;
        }
    }
    
    /* calculate decoding matrix per band */
    Yna_W = malloc1d(nSH * N_dirs*sizeof(float_complex));
    Yna_W_Yna = malloc1d(nSH * nSH * sizeof(float_complex));
    Yna_W_H = malloc1d(nSH * 2 * sizeof(float_complex));
    B = malloc1d(nSH * 2 * sizeof(float_complex));
    hrtfs_mod = malloc1d(2*N_dirs*sizeof(float_complex));
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, N_dirs, N_dirs, &calpha,
                Y_na, N_dirs,
                W, N_dirs, &cbeta,
                Yna_W, N_dirs);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, nSH, N_dirs, &calpha,
                Yna_W, N_dirs,
                Y_na, N_dirs, &cbeta,
                Yna_W_Yna, nSH);
    for(band=0; band < N_bands; band++){
        /* Remove itd from high frequency HRTFs */
        if(band>=band_cutoff){
            for(j=0; j<N_dirs; j++){
                hrtfs_mod[0*N_dirs+j] = ccmulf(hrtfs[band_cutoff*2*N_dirs + 0*N_dirs + j],
                                               cexpf( crmulf(cmplxf(0.0f, 0.0f), (itd_s[j]/2.0f))));
                hrtfs_mod[1*N_dirs+j] = ccmulf(hrtfs[band_cutoff*2*N_dirs + 1*N_dirs + j],
                                               cexpf( crmulf(cmplxf(0.0f, 0.0f), (-itd_s[j]/2.0f))));
            }
        }
        else
            memcpy(hrtfs_mod, &hrtfs[band*2*N_dirs], 2*N_dirs*sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, 2, N_dirs, &calpha,
                    Yna_W, N_dirs,
                    hrtfs_mod, N_dirs, &cbeta,
                    Yna_W_H, 2);
        utility_cglslv(NULL, Yna_W_Yna, nSH, Yna_W_H, 2, B);
        for(i=0; i<nSH; i++)
            for(j=0; j<2; j++)
                decMtx[band*2*nSH + j*nSH + i] = conjf(B[i*2+j]); /* ^H */
    }
    
    free(Y_na);
    free(W);
    free(Yna_W);
    free(Yna_W_Yna);
    free(Yna_W_H);
    free(B);
    free(hrtfs_mod);
}

void getBinDecoder_MAGLS
(
    float_complex* hrtfs, /* the HRTFs; FLAT: N_bands x 2 x N_dirs */
    float* hrtf_dirs_deg,
    int N_dirs,
    int N_bands,
    int order,
    float* freqVector,
    float* weights,
    float_complex* decMtx /* N_bands x 2 x (order+1)^2  */
)
{
    int i, j, nSH, band, band_cutoff;
    float cutoff, minVal;
    float* Y_tmp;
    float_complex* W, *Y_na, *hrtfs_ls, *Yna_W, *Yna_W_Yna, *Yna_W_H, *H_mod, *B_magls;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    
    nSH = ORDER2NSH(order);

    /* integration weights */
    W = calloc1d(N_dirs*N_dirs, sizeof(float_complex));
    if(weights!=NULL)
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(weights[i], 0.0f);
    else
        for(i=0; i<N_dirs; i++)
            W[i*N_dirs+i] = cmplxf(1.0f/(float)N_dirs, 0.0f);
    
    /* SH */
    Y_tmp = malloc1d(nSH*N_dirs*sizeof(float));
    Y_na = malloc1d(nSH*N_dirs*sizeof(float_complex));
    getRSH(order, hrtf_dirs_deg, N_dirs, Y_tmp);
    for(i=0; i<nSH*N_dirs; i++)
        Y_na[i] = cmplxf(Y_tmp[i], 0.0f);
    free(Y_tmp);
    
    /* find band index for cutoff frequency */
    cutoff = 1.5e3f;
    minVal = 2.23e10f;
    for(band=0, band_cutoff=0; band<N_bands; band++){
        if(minVal>fabsf(freqVector[band]-cutoff)){
            minVal = fabsf(freqVector[band]-cutoff);
            band_cutoff = band;
        }
    }
    
    /* calculate decoding matrix per band */
    Yna_W = malloc1d(nSH * N_dirs*sizeof(float_complex));
    Yna_W_Yna = malloc1d(nSH * nSH * sizeof(float_complex));
    Yna_W_H = malloc1d(nSH * 2 * sizeof(float_complex));
    B_magls = malloc1d(nSH * 2 * sizeof(float_complex));
    hrtfs_ls = malloc1d(2*N_dirs*sizeof(float_complex));
    H_mod = malloc1d(2*N_dirs*sizeof(float_complex));
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, N_dirs, N_dirs, &calpha,
                Y_na, N_dirs,
                W, N_dirs, &cbeta,
                Yna_W, N_dirs);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, nSH, N_dirs, &calpha,
                Yna_W, N_dirs,
                Y_na, N_dirs, &cbeta,
                Yna_W_Yna, nSH);
    for (band=0; band<N_bands; band++){
        if(band<=band_cutoff){
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, 2, N_dirs, &calpha,
                        Yna_W, N_dirs,
                        &hrtfs[band*2*N_dirs], N_dirs, &cbeta,
                        Yna_W_H, 2);
            utility_cglslv(NULL, Yna_W_Yna, nSH, Yna_W_H, 2, B_magls);
        }
        else{
            /* Remove itd from high frequency HRTFs */
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, N_dirs, nSH, &calpha,
                        &decMtx[(band-1)*2*nSH] , nSH,
                        Y_na, N_dirs, &cbeta,
                        H_mod, N_dirs);
            for(i=0; i<2*N_dirs; i++)
                H_mod[i] = ccmulf(cmplxf(cabsf(hrtfs[band*2*N_dirs + i]), 0.0f), cexpf(cmplxf(0.0f, atan2f(cimagf(H_mod[i]), crealf(H_mod[i])))));
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, 2, N_dirs, &calpha,
                        Yna_W, N_dirs,
                        H_mod, N_dirs, &cbeta,
                        Yna_W_H, 2);
            utility_cglslv(NULL, Yna_W_Yna, nSH, Yna_W_H, 2, B_magls);
        }
        
        for(i=0; i<nSH; i++)
            for(j=0; j<2; j++)
                decMtx[band*2*nSH + j*nSH + i] = conjf(B_magls[i*2+j]); /* ^H */
    }
    
    free(W);
    free(Y_na);
    free(Yna_W);
    free(Yna_W_Yna);
    free(Yna_W_H);
    free(B_magls);
    free(hrtfs_ls);
    free(H_mod);
}
