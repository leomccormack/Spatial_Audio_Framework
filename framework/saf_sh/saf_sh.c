    /*
 Copyright 2016-2018 Leo McCormack
 
 Permission to use, copy, modify, and/or distribute this software for any purpose with or
 without fee is hereby granted, provided that the above copyright notice and this permission
 notice appear in all copies.
 
 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 OR PERFORMANCE OF THIS SOFTWARE.
*/
/*
 * Filename:
 *     saf_sh.c  
 * Description:
 *     A collection of spherical harmonic related functions. Some of which have been
 *     derived from the Matlab library by Archontis Politis; found here:
 *     https://github.com/polarch/Spherical-Harmonic-Transform
 *     and MATLAB code by Symeon Delikaris-Manias
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 22.05.2016
 */

#include "saf_sh.h"
#include "saf_sh_internal.h"

static double Jn(int n, double z)
{
#ifdef __APPLE__
    return jn(n,z);
#else
    return _jn(n,z);
#endif
}

static double Yn(int n, double z)
{
#ifdef __APPLE__
    return yn(n,z);
#else
    return _yn(n,z);
#endif
}

static int MSTA1(double, int);
static int MSTA2(double,int,int);
static double ENVJ(int N, double X);

static unsigned long factorial(unsigned long f)
{
    if ( f == 0 )
        return 1;
    else
        return(f * factorial(f - 1));
}

/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
static void SPHI(int N, double X, int *NM, double *SI, double *DI)
{
    int K, M;
    double CS, F, F0, F1, SI0;
    
    *NM=N;
    if (fabs(X) < 1e-100) {
        for (K=0; K<=N; K++) {
            SI[K]=0.0;
            DI[K]=0.0;
        }
        SI[0]=1.0;
        DI[1]=0.333333333333333;
        return;
    }
    SI[0]=sinh(X)/X;
    SI[1]=-(sinh(X)/X-cosh(X))/X;
    SI0=SI[0];
    if (N >= 2) {
        M=MSTA1(X,200);
        if (M < N)
            *NM=M;
        else
            M=MSTA2(X,N,15);
        F0=0.0;
        F1=1.0-100;
        for (K=M; K>-1; K--) {
            F=(2.0*K+3.0)*F1/X+F0;
            if (K <= *NM) SI[K]=F;
            F0=F1;
            F1=F;
        }
        CS=SI0/F;
        for (K=0; K<=*NM; K++)  SI[K] *= CS;
    }
    DI[0]=SI[1];
    for (K=1; K<=*NM; K++)
        DI[K]=SI[K-1]-(K+1.0)/X*SI[K];
}

/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
static void SPHK(int N, double X, int *NM, double *SK, double *DK)
{
    int K;
    double F, F0, F1;
    
    *NM=N;
    if (X < 1e-60) {
        for (K=0; K<=N; K++) {
            SK[K]=1.0e+300;
            DK[K]=-1.0e+300;
        }
        return;
    }
    SK[0]=0.5*PI/X*exp(-X);
    SK[1]=SK[0]*(1.0+1.0/X);
    F0=SK[0];
    F1=SK[1];
    for (K=2; K<=N; K++) {
        F=(2.0*K-1.0)*F1/X+F0;
        SK[K]=F;
        if (fabs(F) > 1.0e+300) goto e20;
        F0=F1;
        F1=F;
    }
e20:    *NM=K-1;
    DK[0]=-SK[1];
    for (K=1; K<=*NM; K++)
        DK[K]=-SK[K-1]-(K+1.0)/X*SK[K];
}

/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
static void SPHJ(int N, double X, int *NM, double *SJ, double *DJ)
{
    int K, M;
    double CS, F, F0, F1, SA, SB;
    
    *NM=N;
    if (fabs(X) < 1e-100) {
        for (K=0; K<=N; K++) {
            SJ[K]=0.0;
            DJ[K]=0.0;
        }
        SJ[0]=1.0;
        DJ[1]=0.333333333333333;
        return;
    }
    SJ[0]=sin(X)/X;
    SJ[1]=(SJ[0]-cos(X))/X;
    if (N >= 2) {
        SA=SJ[0];
        SB=SJ[1];
        M=MSTA1(X,200);
        if (M < N)
            *NM=M;
        else
            M=MSTA2(X,N,15);
        F0=0.0;
        F1=1.0-100;
        for (K=M; K>-1; K--) {
            F=(2.0*K+3.0)*F1/X-F0;
            if (K <= *NM)  SJ[K]=F;
            F0=F1;
            F1=F;
        }
        if (fabs(SA) > fabs(SB))  CS=SA/F;
        if (fabs(SA) <= fabs(SB)) CS=SB/F0;
        for (K=0; K<=*NM; K++)  SJ[K] *= CS;
    }
    DJ[0]=(cos(X)-sin(X)/X)/X;
    for (K=1; K<=*NM; K++)
        DJ[K]=SJ[K-1]-(K+1.0)*SJ[K]/X;
}

/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
static int MSTA1(double X, int MP)
{
    double A0,F,F0,F1;
    int IT,NN,N0,N1;
    
    A0=fabs(X);
    N0=floor(1.1*A0)+1;
    F0=ENVJ(N0,A0)-MP;
    N1=N0+5;
    F1=ENVJ(N1,A0)-MP;
    for (IT=1; IT<=20; IT++) {
        NN=N1-(N1-N0)/(1.0-F0/F1);
        F=ENVJ(NN,A0)-MP;
        if (abs(NN-N1) < 1) goto e20;
        N0=N1;
        F0=F1;
        N1=NN;
        F1=F;
    }
e20:    return NN;
}

/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
static int MSTA2(double X, int N, int MP)
{
    double A0,EJN,F,F0,F1,HMP,OBJ;
    int IT,N0,N1,NN;
    
    A0=fabs(X);
    HMP=0.5*MP;
    EJN=ENVJ(N,A0);
    if (EJN <= HMP) {
        OBJ=MP;
        N0=floor(1.1*A0);
    }
    else {
        OBJ=HMP+EJN;
        N0=N;
    }
    F0=ENVJ(N0,A0)-OBJ;
    N1=N0+5;
    F1=ENVJ(N1,A0)-OBJ;
    for (IT=1; IT<=20; IT++) {
        NN=N1-(N1-N0)/(1.0-F0/F1);
        F=ENVJ(NN,A0)-OBJ;
        if (abs(NN-N1) < 1) goto e20;
        N0=N1;
        F0=F1;
        N1=NN;
        F1=F;
    }
e20:    return NN+10;
}

/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
static double ENVJ(int N, double X) {
    return (0.5*log(6.28*N)-N*log(1.36*X/N));
}

/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
static void SPHY(int N, double X, int *NM, double *SY, double *DY)
{
    int K;
    double F, F0, F1;
    
    *NM=N;
    if (X < 1e-60) {
        for (K=0; K<=N; K++) {
            SY[K]=-1.0e+300;
            DY[K]=1e+300;
        }
        return;
    }
    SY[0]=-cos(X)/X;
    SY[1]=(SY[0]-sin(X))/X;
    F0=SY[0];
    F1=SY[1];
    for (K=2; K<=N; K++) {
        F=(2.0*K-1.0)*F1/X-F0;
        SY[K]=F;
        if (fabs(F) >= 1e+300) goto e20;
        F0=F1;
        F1=F;
    }
e20:    *NM=K-1;
    DY[0]=(sin(X)+cos(X)/X)/X;
    for (K=1; K<=*NM; K++)
        DY[K]=SY[K-1]-(K+1.0)*SY[K]/X;
}


/* Used in the calculation of spherical harmonic rotation matrices
 * Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 * by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
static float getP(int i, int l, int a, int b, float** R_1, float** R_lm1)
{
    float ret, ri1, rim1, ri0;
    //ret = 0.0f;

    ri1 = R_1[i + 1][1 + 1];
    rim1 = R_1[i + 1][-1 + 1];
    ri0 = R_1[i + 1][0 + 1];

    if (b == -l)
        ret = ri1 * R_lm1[a + l - 1][0] + rim1 * R_lm1[a + l - 1][2 * l - 2];
    else {
        if (b == l)
            ret = ri1*R_lm1[a + l - 1][2 * l - 2] - rim1 * R_lm1[a + l - 1][0];
        else
            ret = ri0 * R_lm1[a + l - 1][b + l - 1];
    }

    return ret;
}

/* Used in the calculation of spherical harmonic rotation matrices
 * Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 * by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
static float getU(int l, int m, int n, float** R_1, float** R_lm1)
{
    return getP(0, l, m, n, R_1, R_lm1);
}

/* Used in the calculation of spherical harmonic rotation matrices
 * Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 * by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
static float getV(int l, int m, int n, float** R_1, float** R_lm1)
{
    int d;
    float ret, p0, p1;
    ret = 0.0f;

    if (m == 0) {
        p0 = getP(1, l, 1, n, R_1, R_lm1);
        p1 = getP(-1, l, -1, n, R_1, R_lm1);
        ret = p0 + p1;
    }
    else {
        if (m>0) {
            d = m == 1 ? 1 : 0;
            p0 = getP(1, l, m - 1, n, R_1, R_lm1);
            p1 = getP(-1, l, -m + 1, n, R_1, R_lm1);
            ret = p0*sqrtf(1.0f + d) - p1*(1.0f - d);
        }
        else {
            d = m == -1 ? 1 : 0;
            p0 = getP(1, l, m + 1, n, R_1, R_lm1);
            p1 = getP(-1, l, -m - 1, n, R_1, R_lm1);
            ret = p0*(1.0f - (float)d) + p1*sqrtf(1.0f + (float)d);
        }
    }

    return ret;
}

/* Used in the calculation of spherical harmonic rotation matrices
 * Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 * by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
static float getW(int l, int m, int n, float** R_1, float** R_lm1)
{
    float ret, p0, p1;
    ret = 0.0f;

    if (m != 0) {
        if (m>0) {
            p0 = getP(1, l, m + 1, n, R_1, R_lm1);
            p1 = getP(-1, l, -m - 1, n, R_1, R_lm1);
            ret = p0 + p1;
        }
        else {
            p0 = getP(1, l, m - 1, n, R_1, R_lm1);
            p1 = getP(-1, l, -m + 1, n, R_1, R_lm1);
            ret = p0 - p1;
        }
    }
    return ret;
}

/* Will be removed in a later version: use "unnorm_legendreP" */
void legendreP
(
    int l,
    float x,
    float* ppm
)
{
    int m, j, p_x;
    float c_f, c_l, maxc_f, xx, p;

    memset(ppm, 0, (l+1)*sizeof(float));
    if (l>0){
        for (m=0;m<=l;m++){
            c_f = 1.0f;
            c_l = powf((-1.0f), (float)m) * c_f * (float)factorial(2*l) /
                ( powf(2.0f,(float)l) * (float)factorial(l)* (float)factorial(l-m));
            maxc_f = fabsf(c_l);
            p_x = l-m;
            xx = x*x;

            /* Calculate P_l^m (x)/sqrt(1-x^2)^(m/2) */
            p=c_l;
            for (j=l-1; j>=0; j--){
                if(p_x>=2){
                    c_l=-(2.0f*(float)j+2.0f-(float)l-(float)m)*(2.0f*(float)j+1.0f-(float)l-(float)m)
                        /(2.0f*(2.f*(float)j+1.0f)*((float)l-(float)j))*c_l;
                    if (maxc_f < fabsf(c_l))
                        maxc_f = fabsf(c_l);
                    p=p*xx + c_l;
                    p_x=p_x-2;
                }
            }
            if(p_x==1)
                p*=x;
            if(m!=0){
                xx=1.0f-xx;
                for (j=1; j<=(int)((float)m/2.0f); j++)
                    p=p*xx;
                if(m != 2*(int)((float)m/2.0f))
                    p = p*sqrtf(xx);
            }
            ppm[m] = p;
        }
    }
    else
        ppm[0] = 1.0f;
}

void unnorm_legendreP(int n, double* x, int lenX, double* y)
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
    P = calloc((n+3)*lenX,sizeof(double));
    s_n = malloc(lenX*sizeof(double));
    tc = malloc(lenX*sizeof(double));
    sqrt_n = malloc((2*n+1)*sizeof(double));
    
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

void getRSH
(
    int N,
    float* dirs_deg,
    int nDirs,
    float** Y
)
{
    int i, j, nSH;
    float scale;
    float* Y_dir;
    
    scale = sqrtf(4.0f*M_PI);
    nSH = (N+1)*(N+1);
    if((*Y)!=NULL)
        free(*Y);
    (*Y) = malloc(nSH*nDirs*sizeof(float));
    Y_dir = malloc(nSH*sizeof(float));
    for(i=0; i<nDirs; i++){
        /* compute spherical harmonics for each direction */
        getSHreal(N, dirs_deg[i*2]*M_PI/180.0f, M_PI/2.0f- dirs_deg[i*2+1]*M_PI/180.0f, Y_dir);
        for( j=0; j<nSH; j++)
            (*Y)[j*nDirs + i] = Y_dir[j]*scale;
    }
    free(Y_dir);
}

void getSHreal
(
    int N,
    float azi,
    float incl,
    float* Y
)
{
    int i, j, n, m, idx_Y;
    float *p_mm, *Lnm_real, *condon, *norm_real, *Nnm_real, *CosSin, *Ynm;
    //Nharm = (N+1)*(N+1);

    p_mm = (float*)malloc((N+1) * sizeof(float));
    norm_real = (float*)malloc((N+1) * sizeof(float));
    CosSin = (float*)malloc((N*2+1) * sizeof(float));
    Ynm = (float*)malloc((N*2+1) * sizeof(float));
    Nnm_real = (float*)malloc((N*2+1) * sizeof(float));
    Lnm_real = (float*)malloc((N*2+1) * sizeof(float));
    condon = (float*)malloc((N*2+1) * sizeof(float));

    idx_Y = 0;
    for(n=0; n<=N; n++){
        legendreP(n, cosf(incl), p_mm);
        if (n != 0){
            for(i=-n, j=0; i<=n; i++, j++){
                condon[j] = powf(-1.0f, fabsf((float)i));
                Lnm_real[j] = condon[j] * p_mm[abs(i)];
            }
        }
        else
            Lnm_real[0] = p_mm[0];
        for(m=0; m <= n; m++)
            norm_real[m] = sqrtf( (2.0f*(float)n+1.0f)* (float)factorial(n-m) / (4.0f*M_PI*(float)factorial(n+m)) );
        if (n != 0){
            for(i=-n, j=0; i<=n; i++, j++){
                Nnm_real[j] = norm_real[abs(i)];
            }
        }
        else
            Nnm_real[0] = norm_real[0];
        memset(CosSin, 0, (2*n+1)*sizeof(float));
        CosSin[n] = 1.0f;
        if (n != 0){
            for(j=0; j<2*n+1; j++){
                if (j>n)
                    CosSin[j] = sqrtf(2.0f)*cosf((float)(j-n)*azi);
                else if (j<n)
                    CosSin[j] = sqrtf(2.0f)*sinf((float)(n-j)*azi);
            }
        }
        for(j=0; j<2*n+1; j++){
            Ynm[j] = Nnm_real[j] * Lnm_real[j] * CosSin[j];
            Y[idx_Y+j] = Ynm[j];
        }
        idx_Y = idx_Y + (2*n+1);
    }

    free(p_mm);
    free(norm_real);
    free(CosSin);
    free(Ynm);
    free(Nnm_real);
    free(Lnm_real);
    free(condon);
}

void yawPitchRoll2Rzyx
(
    float yaw,
    float pitch,
    float roll,
    int rollPitchYawFLAG, /* use Rxyz, i.e. apply roll, pitch and then yaw */
    float R[3][3]
)
{
    int m,n,k;
    float Rtmp[3][3] = {{0.0f}};
    float Rx[3][3] = { {1.0f ,0.0f ,0.0f }, { 0.0f ,1.0f ,0.0f }, { 0.0f ,0.0f ,1.0f } };
    float Ry[3][3] = { {1.0f ,0.0f ,0.0f }, { 0.0f ,1.0f ,0.0f }, { 0.0f ,0.0f ,1.0f } };
    float Rz[3][3] = { {1.0f ,0.0f ,0.0f }, { 0.0f ,1.0f ,0.0f }, { 0.0f ,0.0f ,1.0f } };
    
    /* var Rx, Ry, Rz; */
    if (roll != 0) {
        Rx[1][1] =  cosf(roll); Rx[1][2] = sinf(roll);
        Rx[2][1] = -sinf(roll); Rx[2][2] = cosf(roll);
    }
    if (pitch != 0){
        Ry[0][0] = cosf(pitch); Ry[0][2] = -sinf(pitch);
        Ry[2][0] = sinf(pitch); Ry[2][2] =  cosf(pitch);
    }
    if (yaw != 0){
        Rz[0][0] =  cosf(yaw); Rz[0][1] = sinf(yaw);
        Rz[1][0] = -sinf(yaw); Rz[1][1] = cosf(yaw);
    }
    if(rollPitchYawFLAG){
        /* rotation order: roll-pitch-yaw */
        for (m=0;m<3; m++){
            memset(R[m], 0, 3*sizeof(float));
            for(n=0;n<3; n++)
                for(k=0; k<3; k++)
                    Rtmp[m][n] += Ry[m][k] * Rx[k][n];
        }
        for (m=0;m<3; m++)
            for(n=0;n<3; n++)
                for(k=0; k<3; k++)
                    R[m][n] += Rz[m][k] * Rtmp[k][n];
    }
    else{
        /* rotation order: yaw-pitch-roll */
        for (m=0;m<3; m++){
            memset(R[m], 0, 3*sizeof(float));
            for(n=0;n<3; n++)
                for(k=0; k<3; k++)
                    Rtmp[m][n] += Ry[m][k] * Rz[k][n];
        }
        for (m=0;m<3; m++)
            for(n=0;n<3; n++)
                for(k=0; k<3; k++)
                    R[m][n] += Rx[m][k] * Rtmp[k][n];
    }
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
    float** R_1,  **R_lm1, **R_l;
    
    M = (L+1) * (L+1);
    R_1 = (float**)calloc2d(3, 3, sizeof(float));
    R_lm1 = (float**)calloc2d(M, M, sizeof(float));
    R_l = (float**)calloc2d(M, M, sizeof(float));
    memset(RotMtx, 0, M*M*sizeof(float));
    
    /* zeroth-band (l=0) is invariant to rotation */
    RotMtx[0] = 1;
    
    /* the first band (l=1) is directly related to the rotation matrix */
    R_1[-1+1][-1+1] = Rxyz[1][1];
    R_1[-1+1][0+1] = Rxyz[1][2];
    R_1[-1+1][1+1] = Rxyz[1][0];
    R_1[ 0+1][-1+1] = Rxyz[2][1];
    R_1[ 0+1][0+1] = Rxyz[2][2];
    R_1[ 0+1][1+1] = Rxyz[2][0];
    R_1[ 1+1][-1+1] = Rxyz[0][1];
    R_1[ 1+1][0+1] = Rxyz[0][2];
    R_1[ 1+1][1+1] = Rxyz[0][0];
    for (i=1; i<4; i++){
        memcpy(R_lm1[i-1], R_1[i-1], 3*sizeof(float));
        for (j=1; j<4; j++)
            RotMtx[i*M+j] = R_1[i-1][j-1];
    }
    
    /* compute rotation matrix of each subsequent band recursively */
    bandIdx = 4;
    for(l = 2; l<=L; l++){
        for(i=0; i<2*l+1; i++)
            memset(R_l[i], 0, (2*l+1) * sizeof(float));
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
                    u = u* getU(l,m,n,R_1,R_lm1);
                if (v!=0)
                    v = v* getV(l,m,n,R_1,R_lm1);
                if (w!=0)
                    w = w* getW(l,m,n,R_1,R_lm1);
                
                R_l[m+l][n+l] = u+v+w;
            }
        }
        
        for(i=0; i<2*l+1; i++)
            for(j=0; j<2*l+1; j++)
                RotMtx[(bandIdx + i)*M + (bandIdx + j)] = R_l[i][j];
                //RotMtx[(bandIdx+i)*(2*l+1) +(bandIdx+j)] = R_l[i][j];
        for(i=0; i<2*l+1; i++)
            memcpy(R_lm1[i], R_l[i], (2*l+1) * sizeof(float));
        bandIdx += 2*l+1;
    }
    
    free2d((void**)R_1, 3);
    free2d((void**)R_lm1, M);
    free2d((void**)R_l, M);
}

void calcBFweights
(
    BEAMFORMING_WEIGHT_TYPES BFW_type,
    int order,
    float azi,
    float elev,
    float* weights
)
{
    int i, j, nSH;
    int o[9] = { 0,1,4,9,16,25,36,49,64 };
    float *d, *Y;
    nSH = (order + 1)*(order + 1);

    /* compute real spherical hamonics */
    Y = (float*)malloc(nSH * sizeof(float));
    getSHreal(order, azi, M_PI / 2.0f - elev, Y);

    /* calculate beamforming weights */
    switch (BFW_type) {
        case BFW_BASIC:
            for (i = 0; i<nSH; i++)
                weights[i] = Y[i];
            break;

        case BFW_MAX_RE:
            d = (float*)calloc((order + 1), sizeof(float));
            maxre3d(order, d);
            for (i = 0; i< (order + 1); i++)
                for (j = o[i]; j< o[i + 1]; j++)
                    weights[j] = Y[j] * d[i];
            free(d);
            break;

        case BFW_DOLPH_CHEBY_MAIN:
            d = (float*)calloc((order + 1), sizeof(float));
            dolph_chebyshev(order, d, 0);
            for (i = 0; i< (order + 1); i++)
                for (j = o[i]; j< o[i + 1]; j++)
                    weights[j] = Y[j] * d[i];
            free(d);
            break;

        case BFW_DOLPH_CHEBY_DESIRED:
            d = (float*)calloc((order + 1), sizeof(float));
            dolph_chebyshev(order, d, 1);
            for (i = 0; i< (order + 1); i++) 
                for (j = o[i]; j< o[i + 1]; j++)
                    weights[j] = Y[j] * d[i];
            free(d);
            break;

        default:
            break;
    }
    free(Y);
}

void unitSph2Cart
(
    float azi_rad,
    float elev_rad,
    float xyz[3]
)
{
    xyz[0] = cosf(elev_rad) * cosf(azi_rad);
    xyz[1] = cosf(elev_rad) * sinf(azi_rad);
    xyz[2] = sinf(elev_rad);
}

void unitCart2Sph
(
    float xyz[3],
    float AziElev_rad[2]
)
{
    float hypotxy = sqrtf(powf(xyz[0], 2.0f) + powf(xyz[1], 2.0f));
    AziElev_rad[0] = atan2f(xyz[1], xyz[0]);
    AziElev_rad[1] = atan2f(xyz[2], hypotxy);
}

void unitCart2Sph_aziElev
(
    float xyz[3],
    float* azi_rad,
    float* elev_rad
)
{
    float hypotxy = sqrtf(powf(xyz[0], 2.0f) + powf(xyz[1], 2.0f));
    (*azi_rad) = atan2f(xyz[1], xyz[0]);
    (*elev_rad) = atan2f(xyz[2], hypotxy);
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
    
    nSH = (order+1)*(order+1);
    Cx_Y = malloc(nSH * nGrid_dirs * sizeof(float_complex));
    Y_Cx_Y = malloc(nGrid_dirs*sizeof(float_complex));
    Cx_Y_s = malloc(nSH*sizeof(float_complex));
    Y_grid_s = malloc(nSH*sizeof(float_complex));
    
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
    
    nSH = (order+1)*(order+1);
    w_MVDR = malloc(nSH * nGrid_dirs*sizeof(float_complex));
    Cx_d = malloc(nSH*nSH*sizeof(float_complex));
    invCx_Ygrid = malloc(nSH*nGrid_dirs*sizeof(float_complex));
    invCx_Ygrid_s = malloc(nSH*sizeof(float_complex));
    Y_grid_s = malloc(nSH*sizeof(float_complex));
    
    /* apply diagonal loading */
    Cx_trace = 0.0f;
    for(i=0; i<nSH; i++)
        Cx_trace += crealf(Cx[i*nSH+i]);
    Cx_trace /= (float)nSH;
    memcpy(Cx_d, Cx, nSH*nSH*sizeof(float_complex));
    for(i=0; i<nSH; i++)
        Cx_d[i*nSH+i] = craddf(Cx_d[i*nSH+i], regPar*Cx_trace);
    
    /* solve the numerator part of the MVDR weights for all grid directions: Cx^-1 * Y */
    utility_cslslv(Cx_d, nSH, Y_grid, nGrid_dirs, invCx_Ygrid);
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
    nSH = (order+1)*(order+1);
    Cx_Y = malloc(nSH * nGrid_dirs * sizeof(float_complex));
    Cx_d = malloc(nSH*nSH*sizeof(float_complex));
    A = malloc(nSH*2*sizeof(float_complex));
    invCxd_A = malloc(nSH*2*sizeof(float_complex));
    invCxd_A_tmp = malloc(nSH*2*sizeof(float_complex));
    w_LCMV_s = malloc(2*nGrid_dirs*sizeof(float_complex));
    w_CroPaC = malloc(nSH*nGrid_dirs*sizeof(float_complex));
    wo = malloc(nSH*sizeof(float_complex));
    mvdr_map = malloc(nGrid_dirs*sizeof(float));
    Cx_Y_s = malloc(nSH*sizeof(float_complex));
    
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
        utility_cslslv(Cx_d, nSH, A, 2, invCxd_A);
        for(j=0; j<nSH*2; j++)
            invCxd_A_tmp[j] = conjf(invCxd_A[j]);
        cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, 2, 2, nSH, &calpha,
                    A, 2,
                    invCxd_A_tmp, 2, &cbeta,
                    A_invCxd_A, 2);
        for(j=0; j<nSH; j++)
            for(k=0; k<2; k++)
                invCxd_A_tmp[k*nSH+j] = invCxd_A[j*2+k];
        utility_cglslv((float_complex*)A_invCxd_A, 2, invCxd_A_tmp, nSH, w_LCMV_s);
        cblas_cgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nSH, 1, 2, &calpha,
                    w_LCMV_s, nSH,
                    b, 1, &cbeta,
                    wo, 1);
        
        /* calculate the cross-spectrum between static beam Y, and adaptive beam wo (LCMV) */
        for(j=0; j<nSH; j++)
            Cx_Y_s[j] = Cx_Y[j*nGrid_dirs+i];
        utility_cvvdot(wo, Cx_Y_s, nSH, NO_CONJ, &Y_wo_xspec);
        
        /* derive CroPaC weights  */
        S = MIN(cabsf(Y_wo_xspec), mvdr_map[i]); /* ensures distortionless response */
        G = sqrtf(S/(mvdr_map[i]+2.23e-10f));
        G = MAX(lambda, G); /* optional spectral floor parameter, to control harshness of attenuation (good for demos) */
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
    
    nSH = (order+1)*(order+1);
    nSources = MIN(nSources, nSH/2);
    V = malloc(nSH*nSH*sizeof(float_complex));
    Vn = malloc(nSH*(nSH-nSources)*sizeof(float_complex));
    Vn_Y = malloc((nSH-nSources)*nGrid_dirs*sizeof(float_complex));
    
    /* obtain eigenvectors */
    utility_ceig(Cx, nSH, 1, NULL, V, NULL, NULL);
    
    /* truncate, to obtain noise sub-space */
    for(i=0; i<nSH; i++)
        for(j=0; j<nSH-nSources; j++)
            Vn[i*(nSH-nSources)+j] = V[i*nSH + j + nSources];
    
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
    
    nSH = (order+1)*(order+1);
    nSources = MIN(nSources, nSH/2);
    V = malloc(nSH*nSH*sizeof(float_complex));
    Vn = malloc(nSH*(nSH-nSources)*sizeof(float_complex));
    Vn1 = malloc((nSH-nSources)*sizeof(float_complex));
    Un = malloc(nSH*sizeof(float_complex));
    Un_Y = malloc(nGrid_dirs*sizeof(float_complex));
    
    /* obtain eigenvectors */
    utility_ceig(Cx, nSH, 1, NULL, V, NULL, NULL);
    
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

void bessel_Jn /* untested */
(
    int N,
    double* z,
    int nZ,
    double* J_n,
    double* dJ_n
)
{
    int n, i;
    
    for(i=0; i<nZ; i++){
        if(z[i] <= 1e-15){
            if(J_n!=NULL)
                memset(&J_n[0], 0, (N+1)*sizeof(double));
            if(dJ_n!=NULL)
                memset(&dJ_n[0], 0, (N+1)*sizeof(double));
        }
        else{
            for(n=0; n<N; n++){
                if(J_n!=NULL)
                    J_n[i*(N+1)+n] = Jn(n, z[i]);
                if(n==0 && dJ_n!=NULL)
                    dJ_n[i*(N+1)+n] = -Jn(1, z[i]);
                else if(dJ_n!=NULL)
                    dJ_n[i*(N+1)+n] = (Jn(n-1, z[i])-Jn(n+1, z[i]))/2.0;
            }
        }
    }
}

void bessel_Yn /* untested */
(
    int N,
    double* z,
    int nZ,
    double* Y_n,
    double* dY_n
)
{
    int n, i;
    
    for(i=0; i<nZ; i++){
        if(z[i] <= 1e-15){
            if(Y_n!=NULL)
                memset(&Y_n[0], 0, (N+1)*sizeof(double));
            if(dY_n!=NULL)
                memset(&dY_n[0], 0, (N+1)*sizeof(double));
        }
        else{
            for(n=0; n<N; n++){
                if(Y_n!=NULL)
                    Y_n[i*(N+1)+n] = Yn(n, z[i]);
                if(n==0 && dY_n!=NULL)
                    dY_n[i*(N+1)+n] = -Yn(1, z[i]);
                else if(dY_n!=NULL)
                    dY_n[i*(N+1)+n] = (Yn(n-1, z[i])-Yn(n+1, z[i]))/2.0;
            }
        }
    }
}

void hankel_Hn1 /* untested */
(
    int N,
    double* z,
    int nZ,
    double_complex* H_n1,
    double_complex* dH_n1
)
{
    int n, i;
    
    for(i=0; i<nZ; i++){
        if(z[i] <= 1e-15){
            if(H_n1!=NULL)
                memset(&H_n1[0], 0, (N+1)*sizeof(double_complex));
            if(dH_n1!=NULL)
                memset(&dH_n1[0], 0, (N+1)*sizeof(double_complex));
        }
        else{
            for(n=0; n<N; n++){
                if(H_n1!=NULL)
                    H_n1[i*(N+1)+n] = cmplx(Jn(n, z[i]), Yn(n, z[i]));
                if(dH_n1!=NULL)
                    dH_n1[i*(N+1)+n] = ccsub(crmul(cmplx(Jn(n, z[i]), Yn(n, z[i])), (double)n/MAX(z[i],2.23e-13f)), cmplx(Jn(n+1, z[i]), Yn(n+1, z[i])));
            }
        }
    }
}

void hankel_Hn2 /* untested */
(
    int N,
    double* z,
    int nZ,
    double_complex* H_n2,
    double_complex* dH_n2
)
{
    int n, i;
    
    for(i=0; i<nZ; i++){
        if(z[i] <= 1e-15){
            if(H_n2!=NULL)
                memset(&H_n2[0], 0, (N+1)*sizeof(double_complex));
            if(dH_n2!=NULL)
                memset(&dH_n2[0], 0, (N+1)*sizeof(double_complex));
        }
        else{
            for(n=0; n<N; n++){
                if(H_n2!=NULL)
                    H_n2[i*(N+1)+n] = cmplx(Jn(n, z[i]), -Yn(n, z[i]));
                if(n==0 && dH_n2!=NULL)
                    dH_n2[i*(N+1)+n] = crmul(ccsub(ccmul(cmplx(Jn(1, z[i]), Yn(1, z[i])), cexp(cmplx(0.0, -M_PI))), cmplx(Jn(1, z[i]), -Yn(1, z[i]))), 0.5);
                else if(dH_n2!=NULL)
                    dH_n2[i*(N+1)+n] = crmul(ccsub(cmplx(Jn(n-1, z[i]), -Yn(n-1, z[i])), cmplx(Jn(n+1, z[i]), -Yn(n+1, z[i]))), 0.5);
            }
        }
    }
}

void bessel_jn
(
    int N,
    double* z,
    int nZ,
    int* maxN,
    double* j_n,
    double* dj_n
)
{
    int n, i, NM;
    double* j_n_tmp, *dj_n_tmp;
    
    j_n_tmp = malloc((N+1)*sizeof(double));
    dj_n_tmp = malloc((N+1)*sizeof(double));
    *maxN = 1e8;
    for(i=0; i<nZ; i++){
        if(z[i] <= 1e-15){
            if(j_n!=NULL){
                memset(&j_n[0], 0, (N+1)*sizeof(double));
                j_n[0] = 1.0f;
            }
            if(dj_n!=NULL){
                memset(&dj_n[0], 0, (N+1)*sizeof(double));
                if(N>0)
                    dj_n[1] = 1.0/3.0;
            }
        }
        else{
            SPHJ(N, z[i], &NM, j_n_tmp, dj_n_tmp);
            *maxN = MIN(NM, *maxN );
            for(n=0; n<NM; n++){
                if(j_n!=NULL)
                    j_n [i*(N+1)+n] = j_n_tmp[n];
                if(dj_n!=NULL)
                    dj_n[i*(N+1)+n] = dj_n_tmp[n];
            }
            for(; n<N+1; n++){
                if(j_n!=NULL)
                    j_n [i*(N+1)+n] = 0.0;
                if(dj_n!=NULL)
                    dj_n [i*(N+1)+n] = 0.0;
            }
        }
    }
    *maxN = *maxN==1e8 ? 0 : *maxN; /* maximum order that could be computed */
    
    free(j_n_tmp);
    free(dj_n_tmp);
}

void bessel_in /* untested */
(
    int N,
    double* z,
    int nZ,
    int* maxN,
    double* i_n,
    double* di_n
)
{
    int n, i, NM;
    double* i_n_tmp, *di_n_tmp;
    
    i_n_tmp = malloc((N+1)*sizeof(double));
    di_n_tmp = malloc((N+1)*sizeof(double));
    *maxN = 1e8;
    for(i=0; i<nZ; i++){
        if(z[i] <= 1e-15){
            if(i_n!=NULL){
                memset(&i_n[0], 0, (N+1)*sizeof(double));
                i_n[0] = 1.0f;
            }
            if(di_n!=NULL){
                memset(&di_n[0], 0, (N+1)*sizeof(double));
                if(N>0)
                    di_n[1] = 1.0/3.0;
            }
        }
        else{
            SPHI(N, z[i], &NM, i_n_tmp, di_n_tmp);
            *maxN = MIN(NM, *maxN );
            for(n=0; n<NM; n++){
                if(i_n!=NULL)
                    i_n [i*(N+1)+n] = i_n_tmp[n];
                if(di_n!=NULL)
                    di_n[i*(N+1)+n] = di_n_tmp[n];
            }
            for(; n<N+1; n++){
                if(i_n!=NULL)
                    i_n [i*(N+1)+n] = 0.0;
                if(di_n!=NULL)
                    di_n [i*(N+1)+n] = 0.0;
            }
        }
    }
    *maxN = *maxN==1e8 ? 0 : *maxN; /* maximum order that could be computed */
    
    free(i_n_tmp);
    free(di_n_tmp);
}

void bessel_yn
(
    int N,
    double* z,
    int nZ,
    int* maxN,
    double* y_n,
    double* dy_n
)
{
    int n, i, NM;
    double* y_n_tmp, *dy_n_tmp;
    
    y_n_tmp = malloc((N+1)*sizeof(double));
    dy_n_tmp = malloc((N+1)*sizeof(double));
    *maxN = 1e8;
    for(i=0; i<nZ; i++){
        if(z[i] <= 1e-15){
            if(y_n!=NULL)
                memset(&y_n[0], 0, (N+1)*sizeof(double));
            if(dy_n!=NULL)
                memset(&dy_n[0], 0, (N+1)*sizeof(double));
        }
        else{
            SPHY(N, z[i], &NM, y_n_tmp, dy_n_tmp);
            *maxN = MIN(NM, *maxN );
            for(n=0; n<NM; n++){
                if(y_n!=NULL)
                    y_n [i*(N+1)+n] = y_n_tmp[n];
                if(dy_n!=NULL)
                    dy_n[i*(N+1)+n] = dy_n_tmp[n];
            }
            for(; n<N+1; n++){
                if(y_n!=NULL)
                    y_n [i*(N+1)+n] = 0.0;
                if(dy_n!=NULL)
                    dy_n [i*(N+1)+n] = 0.0;
            }
        }
    }
    *maxN = *maxN==1e8 ? 0 : *maxN; /* maximum order that could be computed */
    
    free(y_n_tmp);
    free(dy_n_tmp);
}

void bessel_kn /* untested */
(
    int N,
    double* z,
    int nZ,
    int* maxN,
    double* k_n,
    double* dk_n
)
{
    int n, i, NM;
    double* k_n_tmp, *dk_n_tmp;
    
    k_n_tmp = malloc((N+1)*sizeof(double));
    dk_n_tmp = malloc((N+1)*sizeof(double));
    *maxN = 1e8;
    for(i=0; i<nZ; i++){
        if(z[i] <= 1e-15){
            if(k_n!=NULL)
                memset(&k_n[0], 0, (N+1)*sizeof(double));
            if(dk_n!=NULL)
                memset(&dk_n[0], 0, (N+1)*sizeof(double));
        }
        else{
            SPHK(N, z[i], &NM, k_n_tmp, dk_n_tmp);
            *maxN = MIN(NM, *maxN );
            for(n=0; n<NM; n++){
                if(k_n!=NULL)
                    k_n [i*(N+1)+n] = k_n_tmp[n];
                if(dk_n!=NULL)
                    dk_n[i*(N+1)+n] = dk_n_tmp[n];
            }
            for(; n<N+1; n++){
                if(k_n!=NULL)
                    k_n [i*(N+1)+n] = 0.0;
                if(dk_n!=NULL)
                    dk_n [i*(N+1)+n] = 0.0;
            }
        }
    }
    *maxN = *maxN==1e8 ? 0 : *maxN; /* maximum order that could be computed */
    
    free(k_n_tmp);
    free(dk_n_tmp);
}

void hankel_hn1
(
    int N,
    double* z,
    int nZ,
    int* maxN,
    double_complex* h_n1,
    double_complex* dh_n1
)
{
    int n, i, NM1, NM2;
    double* j_n_tmp, *dj_n_tmp, *y_n_tmp, *dy_n_tmp;
    
    j_n_tmp = calloc((N+1),sizeof(double));
    dj_n_tmp = calloc((N+1),sizeof(double));
    y_n_tmp = calloc((N+1),sizeof(double));
    dy_n_tmp = calloc((N+1),sizeof(double));
    *maxN = 1e8;
    for(i=0; i<nZ; i++){
        if(z[i] <= 1e-15){
            if(h_n1!=NULL){
                memset(&h_n1[0], 0, (N+1)*sizeof(double_complex));
                h_n1[0] = cmplx(1.0, 0.0);
            }
            if(dh_n1!=NULL)
                memset(&dh_n1[0], 0, (N+1)*sizeof(double_complex));
        }
        else{
            SPHJ(N, z[i], &NM1, j_n_tmp, dj_n_tmp);
            *maxN = MIN(NM1, *maxN );
            SPHY(N, z[i], &NM2, y_n_tmp, dy_n_tmp);
            *maxN = MIN(NM2, *maxN );
            for(n=0; n<MIN(NM1,NM2); n++){
                if(h_n1!=NULL)
                    h_n1 [i*(N+1)+n] = cmplx(j_n_tmp[n], y_n_tmp[n]);
                if(dh_n1!=NULL)
                    dh_n1[i*(N+1)+n] = cmplx(dj_n_tmp[n], dy_n_tmp[n]);
            }
            for(; n<N+1; n++){
                if(h_n1!=NULL)
                    h_n1 [i*(N+1)+n] = cmplx(0.0,0.0);
                if(dh_n1!=NULL)
                    dh_n1 [i*(N+1)+n] = cmplx(0.0,0.0);
            }
        }
    }
    *maxN = *maxN==1e8 ? 0 : *maxN; /* maximum order that could be computed */
    
    free(j_n_tmp);
    free(dj_n_tmp);
    free(y_n_tmp);
    free(dy_n_tmp);
}

void hankel_hn2
(
    int N,
    double* z,
    int nZ,
    int* maxN,
    double_complex* h_n2,
    double_complex* dh_n2
)
{
    int n, i, NM1, NM2;
    double* j_n_tmp, *dj_n_tmp, *y_n_tmp, *dy_n_tmp;
    
    j_n_tmp = calloc((N+1),sizeof(double));
    dj_n_tmp = calloc((N+1),sizeof(double));
    y_n_tmp = calloc((N+1),sizeof(double));
    dy_n_tmp = calloc((N+1),sizeof(double));
    *maxN = 1e8;
    for(i=0; i<nZ; i++){
        if(z[i] <= 1e-15){
            if(h_n2!=NULL){
                memset(&h_n2[0], 0, (N+1)*sizeof(double_complex));
                h_n2[0] = cmplx(1.0, 0.0);
            }
            if(dh_n2!=NULL)
                memset(&dh_n2[0], 0, (N+1)*sizeof(double_complex));
        }
        else{
            SPHJ(N, z[i], &NM1, j_n_tmp, dj_n_tmp);
            *maxN = MIN(NM1, *maxN );
            SPHY(N, z[i], &NM2, y_n_tmp, dy_n_tmp);
            *maxN = MIN(NM2, *maxN );
            for(n=0; n<MIN(NM1,NM2); n++){
                if(h_n2!=NULL)
                    h_n2 [i*(N+1)+n] = cmplx(j_n_tmp[n], -y_n_tmp[n]);
                if(dh_n2!=NULL)
                    dh_n2[i*(N+1)+n] = cmplx(dj_n_tmp[n], -dy_n_tmp[n]);
            }
            for(; n<N+1; n++){
                if(h_n2!=NULL)
                    h_n2 [i*(N+1)+n] = cmplx(0.0,0.0);
                if(dh_n2!=NULL)
                    dh_n2 [i*(N+1)+n] = cmplx(0.0,0.0);
            }
        }
    }
    *maxN = *maxN==1e8 ? 0 : *maxN; /* maximum order that could be computed */
    
    free(j_n_tmp);
    free(dj_n_tmp);
    free(y_n_tmp);
    free(dy_n_tmp);
}

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
    
    assert(arrayType!=ARRAY_CONSTRUCTION_DIRECTIONAL);
    
    memset(b_N, 0, nBands*(order+1)*sizeof(double_complex));
    switch(arrayType){
        default:
        case ARRAY_CONSTRUCTION_OPEN:
            /* compute spherical Bessels of the first kind */
            Jn = malloc(nBands*(order+1)*sizeof(double));
            bessel_Jn(order, kr, nBands, Jn, NULL);
            
            /* modal coefficients for open spherical array (omni sensors): 1i^n * jn; */
            for(n=0; n<order; n++)
                for(i=0; i<nBands; i++)
                    b_N[i*(order+1)+n] = crmul(cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)), Jn[i*(order+1)+n]);
            
            free(Jn);
            break;
            
        case ARRAY_CONSTRUCTION_RIGID:
            /* compute spherical Bessels/Hankels and their derivatives */
            Jn = malloc(nBands*(order+1)*sizeof(double));
            Jnprime = malloc(nBands*(order+1)*sizeof(double));
            Hn2 = malloc(nBands*(order+1)*sizeof(double_complex));
            Hn2prime = malloc(nBands*(order+1)*sizeof(double_complex));
            bessel_Jn(order, kr, nBands, Jn, Jnprime);
            hankel_Hn2(order, kr, nBands, Hn2, Hn2prime);
            
            /* modal coefficients for rigid spherical array: 1i^n * (jn-(jnprime./hn2prime).*hn2); */
            for(i=0; i<nBands; i++){
                for(n=0; n<order; n++){
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
            
        case ARRAY_CONSTRUCTION_DIRECTIONAL:
            /* not supported */
            break;
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
            jn = malloc(nBands*(order+1)*sizeof(double));
            bessel_jn(order, kr, nBands, &maxN, jn, NULL);
            
            /* modal coefficients for open spherical array (omni sensors): 4*pi*1i^n * jn; */
            for(n=0; n<maxN; n++)
                for(i=0; i<nBands; i++)
                    b_N[i*(order+1)+n] = crmul(crmul(cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)), 4.0*M_PI), jn[i*(order+1)+n]);
            
            free(jn);
            break;
            
        case ARRAY_CONSTRUCTION_RIGID:
            /* compute spherical Bessels/Hankels and their derivatives */
            jn = malloc(nBands*(order+1)*sizeof(double));
            jnprime = malloc(nBands*(order+1)*sizeof(double));
            hn2 = malloc(nBands*(order+1)*sizeof(double_complex));
            hn2prime = malloc(nBands*(order+1)*sizeof(double_complex));
            maxN = 1e8;
            bessel_jn(order, kr, nBands, &maxN_tmp, jn, jnprime);
            maxN = MIN(maxN_tmp, maxN);
            hankel_hn2(order, kr, nBands, &maxN_tmp, hn2, hn2prime);
            maxN = MIN(maxN_tmp, maxN); /* maxN being the minimum highest order that was computed for all values in kr */
            
            /* modal coefficients for rigid spherical array: 4*pi*1i^n * (jn-(jnprime./hn2prime).*hn2); */
            for(i=0; i<nBands; i++){
                for(n=0; n<maxN; n++){
                    if(n==0 && kr[i]<=1e-20)
                        b_N[i*(order+1)+n] = cmplx(4.0*M_PI, 0.0);
                    else if(kr[i] <= 1e-20)
                        b_N[i*(order+1)+n] = cmplx(0.0, 0.0);
                    else{
                        b_N[i*(order+1)+n] = ccmul(crmul(cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)), 4.0*M_PI), ( ccsub(cmplx(jn[i*(order+1)+n], 0.0),
                                             ccmul(ccdiv(cmplx(jnprime[i*(order+1)+n],0.0), hn2prime[i*(order+1)+n]), hn2[i*(order+1)+n]))));
                    }
                }
            }
            
            free(jn);
            free(jnprime);
            free(hn2);
            free(hn2prime);
            break;
            
        case ARRAY_CONSTRUCTION_DIRECTIONAL:
            /* compute spherical Bessels of the first kind + derivatives */
            jn = malloc(nBands*(order+1)*sizeof(double));
            jnprime = malloc(nBands*(order+1)*sizeof(double));
            bessel_jn(order, kr, nBands, &maxN, jn, jnprime);
            
            /* modal coefficients for open spherical array (directional sensors): 4*pi*1i^n * (dirCoeff*jn - 1i*(1-dirCoeff)*jnprime); */
            for(n=0; n<maxN; n++)
                for(i=0; i<nBands; i++)
                    b_N[i*(order+1)+n] = ccmul(crmul(cpow(cmplx(0.0,1.0), cmplx((double)n,0.0)), 4.0*M_PI), ccsub(cmplx(dirCoeff*jn[i*(order+1)+n], 0.0),
                                         cmplx(0.0, (1.0-dirCoeff)*jnprime[i*(order+1)+n]))  );
            
            free(jn);
            free(jnprime);
            break;
    }
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
    b_N = malloc(nBands * (order+1) * sizeof(double_complex));
    cylModalCoeffs(order, kr, nBands, arrayType, b_N);
    
    /* Compute angular-dependent part of the array responses */
    C = malloc((order+1)*N_sensors*sizeof(double_complex));
    b_NC = malloc(nBands*N_sensors*sizeof(double_complex));
    for(i=0; i<N_srcs; i++){
        for(j=0; j<N_sensors; j++){
            angle = sensor_dirs_rad[i*2] - src_dirs_deg[i*2]*M_PI/180.0;
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
    b_N = malloc(nBands * (order+1) * sizeof(double_complex));
    sphModalCoeffs(order, kr, nBands, arrayType, dirCoeff, b_N); /* double precision recommended for high orders of sph Bessels  */
    
    /* calculate (unit) cartesian coords for sensors and plane waves */
    U_sensors = malloc(N_sensors*3*sizeof(float));
    U_srcs = malloc(N_srcs*3*sizeof(float));
    for(i=0; i<N_sensors; i++)
        unitSph2Cart(sensor_dirs_rad[i*2], sensor_dirs_rad[i*2+1], (float*)&U_sensors[i*3]);
    for(i=0; i<N_srcs; i++)
        unitSph2Cart(src_dirs_deg[i*2]*M_PI/180.0f, src_dirs_deg[i*2+1]*M_PI/180.0f, (float*)&U_srcs[i*3]);
    
    /* Compute angular-dependent part of the array responses */
    ppm = malloc((order+1)*sizeof(double));
    P = malloc((order+1)*N_sensors*sizeof(double_complex));
    b_NP = malloc(nBands*N_sensors*sizeof(double_complex));
    for(i=0; i<N_srcs; i++){
        for(j=0; j<N_sensors; j++){
            utility_svvdot((const float*)&U_sensors[j*3], (const float*)&U_srcs[i*3], 3, &cosangle);
            for(n=0; n<order+1; n++){
                /* Legendre polynomials correspond to the angular dependency */
                dcosangle = (double)cosangle;
                unnorm_legendreP(n, &dcosangle, 1, ppm);
                P[n*N_sensors+j] = cmplx((2.0*(double)n+1.0)/(4.0*M_PI) * ppm[0], 0.0);
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
    
    nSH = (order+1)*(order+1);
    w_uni_grid = 1.0f/(float)nDirs;
    y_recon_kk = malloc(nSH*nDirs*sizeof(float_complex));
    y_recon_nm = malloc(nDirs*sizeof(float_complex));
    w_y_recon_nm = malloc(nDirs*sizeof(float_complex));
    y_ideal_nm = malloc(nDirs*sizeof(float_complex));
    MH_M = malloc(nSensors*nSensors*sizeof(float_complex));
    EigV = malloc(nSensors*nSensors*sizeof(float_complex));
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
            cSH[band*(order+1)+n] = MAX(MIN(cabsf(cSH_n)/(2.0f*(float)n+1.0f),1.0f),0.0f);
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
        utility_ceig(MH_M, nSensors, 1, NULL, NULL, EigV); /* eigenvalues in decending order */
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


