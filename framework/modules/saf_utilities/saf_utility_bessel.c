/*
 * Copyright 2020 Leo McCormack
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
 * @file saf_bessel.c
 * @brief Collection of routines for computing spherical and cylindrical Bessel
 *        and Hankel functions, including their derivatives
 * @author Leo McCormack
 * @date 26.05.2020
 */

#include "saf_utility_bessel.h"
#include "saf_utilities.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ========================================================================== */
/*                            Internal Functions                              */
/* ========================================================================== */

/**
 * Helper function, used when computing spherical bessel function values
 *
 * Original Fortran code: "Fortran Routines for Computation of Special
 * Functions": jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
static double ENVJ
(
    int N,
    double X
)
{
    return (0.5*log(6.28*N)-N*log(1.36*X/N));
}

/**
 * Helper function, used when computing spherical bessel function values
 *
 * Original Fortran code: "Fortran Routines for Computation of Special
 * Functions": jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
static int MSTA1
(
    double X,
    int MP
)
{
    double A0,F,F0,F1;
    int IT,NN,N0,N1;

    NN = 0;
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

/**
 * Helper function, used when computing spherical bessel function values
 *
 * Original Fortran code: "Fortran Routines for Computation of Special
 * Functions": jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
static int MSTA2
(
    double X,
    int N,
    int MP
)
{
    double A0,EJN,F,F0,F1,HMP,OBJ;
    int IT,N0,N1,NN;

    NN = 0;
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

/**
 * Helper function for bessel_in
 *
 * Original Fortran code: "Fortran Routines for Computation of Special
 * Functions": jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr)
 *
 * @note This function has been modified to avoid numerical instability, at the
 *       cost of slightly higher numerical inprecision (but only when such cases
 *       arise)
 */
static void SPHI
(
    int N,
    double X,
    int *NM,
    double *SI,
    double *DI
)
{
    int K, M, i;
    double CS, F, F0, F1, SI0;

    *NM=N;
    if (fabs(X) < 1e-20) {
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
        /* I had to add this while loop to avoid NaNs and sacrifice some precision, but only when needed */
        i=0;
        while (M < 0) {
            M=MSTA2(X,N,14-i);
            i++;
            if(i==14)
                M=0;
        }
        F0=0.0;
        F1=1.0-100;
        F=1;
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

/**
 * Helper function for bessel_kn
 *
 * Original Fortran code: "Fortran Routines for Computation of Special
 * Functions": jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr)
 */
static void SPHK
(
    int N,
    double X,
    int *NM,
    double *SK,
    double *DK
)
{
    int K;
    double F, F0, F1;

    *NM=N;
    if (X < 1e-20) {
        for (K=0; K<=N; K++) {
            SK[K]=1.0e+300;
            DK[K]=-1.0e+300;
        }
        return;
    }
    SK[0]=0.5*SAF_PId/X*exp(-X);
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

/**
 * Helper function for bessel_in
 *
 * Original Fortran code: "Fortran Routines for Computation of Special
 * Functions": jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr)
 *
 * @note This function has been modified to avoid numerical instability, at the
 *       cost of slightly higher numerical inprecision (but only when such cases
 *       arise)
 */
static void SPHJ
(
    int N,
    double X,
    int *NM,
    double *SJ,
    double *DJ
)
{
    int K, M, i;
    double CS, F, F0, F1, SA, SB;

    *NM=N;
    if (fabs(X) < 1e-80) {
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
        /* I had to add this while loop to avoid NaNs and sacrifice some precision, but only when needed */
        i=0;
        while (M < 0) {
            M=MSTA2(X,N,14-i);
            i++;
            if(i==14)
                M=0;
        }
        F0=0.0;
        F1=1.0-100;
        F=1;
        CS=1;
        for (K=M; K>-1; K--) {
            F=(2.0*K+3.0)*F1/X-F0;
            if (K <= *NM)  SJ[K]=F;
            F0=F1;
            F1=F;
        }
        if (fabs(SA) > fabs(SB))  CS=SA/F;
        if (fabs(SA) <= fabs(SB)) CS=SB/F0;
        for (K=0; K<=*NM; K++) SJ[K] *= CS;
    }
    DJ[0]=(cos(X)-sin(X)/X)/X;
    for (K=1; K<=*NM; K++)
        DJ[K]=SJ[K-1]-(K+1.0)*SJ[K]/X;
}

/**
 * Helper function for bessel_yn
 *
 * Original Fortran code: "Fortran Routines for Computation of Special
 * Functions": jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
static void SPHY
(
    int N,
    double X,
    int *NM,
    double *SY,
    double *DY
)
{
    int K;
    double F, F0, F1;

    *NM=N;
    if (X < 1e-20) {
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

/**
 * Wrapper for the cylindrical Bessel function of the first kind
 */
static double Jn(int n, double z)
{
#ifndef _MSC_VER
    return jn(n,z);
#else
    return _jn(n,z); /* Why Microsoft?! Why?! */
#endif
}

/**
 * Wrapper for the cylindrical Bessel function of the second kind
 */
static double Yn(int n, double z)
{
#ifndef _MSC_VER
    return yn(n,z);
#else
    return _yn(n,z); /* ... */
#endif
}


/* ========================================================================== */
/*                        Cylindrical Bessel Functions                        */
/* ========================================================================== */

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
            for(n=0; n<N+1; n++){
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
            for(n=0; n<N+1; n++){
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
            for(n=0; n<N+1; n++){
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
            for(n=0; n<N+1; n++){
                if(H_n2!=NULL)
                    H_n2[i*(N+1)+n] = cmplx(Jn(n, z[i]), -Yn(n, z[i]));
                if(n==0 && dH_n2!=NULL)
                    dH_n2[i*(N+1)+n] = crmul(ccsub(ccmul(cmplx(Jn(1, z[i]), Yn(1, z[i])), cexp(cmplx(0.0, -SAF_PId))), cmplx(Jn(1, z[i]), -Yn(1, z[i]))), 0.5);
                else if(dH_n2!=NULL)
                    dH_n2[i*(N+1)+n] = crmul(ccsub(cmplx(Jn(n-1, z[i]), -Yn(n-1, z[i])), cmplx(Jn(n+1, z[i]), -Yn(n+1, z[i]))), 0.5);
            }
        }
    }
}


/* ========================================================================== */
/*                         Spherical Bessel Functions                         */
/* ========================================================================== */

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

    j_n_tmp = malloc1d((N+1)*sizeof(double));
    dj_n_tmp = malloc1d((N+1)*sizeof(double));
    *maxN = 1000000000;
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
            for(n=0; n<NM+1; n++){
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
#ifndef NDEBUG
    if(*maxN<N)
        saf_error_print(SAF_WARNING__UNABLE_TO_COMPUTE_BESSEL_FUNCTION_AT_SPECIFIED_ORDER);
#endif

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

    i_n_tmp = malloc1d((N+1)*sizeof(double));
    di_n_tmp = malloc1d((N+1)*sizeof(double));
    *maxN = 1000000000;
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
            for(n=0; n<NM+1; n++){
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
#ifndef NDEBUG
    if(*maxN<N)
        saf_error_print(SAF_WARNING__UNABLE_TO_COMPUTE_BESSEL_FUNCTION_AT_SPECIFIED_ORDER);
#endif

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

    y_n_tmp = malloc1d((N+1)*sizeof(double));
    dy_n_tmp = malloc1d((N+1)*sizeof(double));
    *maxN = 1000000000;
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
            for(n=0; n<NM+1; n++){
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
#ifndef NDEBUG
    if(*maxN<N)
        saf_error_print(SAF_WARNING__UNABLE_TO_COMPUTE_BESSEL_FUNCTION_AT_SPECIFIED_ORDER);
#endif

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

    k_n_tmp = malloc1d((N+1)*sizeof(double));
    dk_n_tmp = malloc1d((N+1)*sizeof(double));
    *maxN = 1000000000;
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
            for(n=0; n<NM+1; n++){
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
#ifndef NDEBUG
    if(*maxN<N)
        saf_error_print(SAF_WARNING__UNABLE_TO_COMPUTE_BESSEL_FUNCTION_AT_SPECIFIED_ORDER);
#endif

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

    j_n_tmp = calloc1d((N+1),sizeof(double));
    dj_n_tmp = calloc1d((N+1),sizeof(double));
    y_n_tmp = calloc1d((N+1),sizeof(double));
    dy_n_tmp = calloc1d((N+1),sizeof(double));
    *maxN = 1000000000;
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
            for(n=0; n<MIN(NM1,NM2)+1; n++){
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
#ifndef NDEBUG
    if(*maxN<N)
        saf_error_print(SAF_WARNING__UNABLE_TO_COMPUTE_BESSEL_FUNCTION_AT_SPECIFIED_ORDER);
#endif

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

    j_n_tmp = calloc1d((N+1),sizeof(double));
    dj_n_tmp = calloc1d((N+1),sizeof(double));
    y_n_tmp = calloc1d((N+1),sizeof(double));
    dy_n_tmp = calloc1d((N+1),sizeof(double));
    *maxN = 1000000000;
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
            for(n=0; n<MIN(NM1,NM2)+1; n++){
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
#ifndef NDEBUG
    if(*maxN<N)
        saf_error_print(SAF_WARNING__UNABLE_TO_COMPUTE_BESSEL_FUNCTION_AT_SPECIFIED_ORDER);
#endif

    free(j_n_tmp);
    free(dj_n_tmp);
    free(y_n_tmp);
    free(dy_n_tmp);
}
