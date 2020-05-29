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
 * @file saf_utility_misc.c
 * @brief Utility: A collection of miscellaneous functions
 *
 * @author Leo McCormack
 * @date 29.01.2020
 */

#include "saf_utility_misc.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * Precomputed factorials for up to !15 (i.e. the "getSH" functions will employ
 * these up to 7th order)
 */
static const long double factorials_15[15] =
{1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0, 479001600.0, 6.2270208e9, 8.71782891e10};

long double factorial(int n)
{
    int i;
    long double ff;
    if(n<15)
        return factorials_15[n];
    else{
        ff = 1.0;
        for(i = 1; i<=n; i++)
            ff *= (long double)i;
        return ff;
    }
}

float matlab_fmodf(float x, float y) {
    float tmp = fmodf(x, y);
    return tmp >= 0 ? tmp : tmp + y;
}

void cxcorr
(
    float* a,
    float* b,
    float* x_ab,
    size_t la,
    size_t lb
)
{
    int m, n, negFLAG, arg;
    size_t len, lim;
    
    len = la + lb - 1;
    memset(x_ab, 0, len*sizeof(float));
    for(m=1; m<=len; m++){
        arg = m-(int)la;
        if(arg<0){
            negFLAG = 1;
            lim = la + arg;
        }
        else{
            negFLAG = 0;
            lim = la - arg;
        }
        for(n=1; n<=lim; n++){
            if(negFLAG == 0)
                x_ab[m-1] += (a[arg+n-1] * b[n-1]);
            else
                x_ab[m-1] += (a[n-1] * b[n-arg-1]);
        }
    }
}

void rand_m1_1
(
    float* vector,
    int length
)
{
    int i;
    for(i=0; i<length; i++)
        vector[i] = (2.0f*rand()/(float)RAND_MAX)-1.0f;
}

void rand_0_1
(
    float* vector,
    int length
)
{
    int i;
    for(i=0; i<length; i++)
        vector[i] = rand()/(float)RAND_MAX;
}

void convd
(
    double* x,
    double* h,
    int len_x,
    int len_h,
    double* y
)
{
    int i, j, h_start, x_start, x_end, len_y;

    len_y = len_h+len_x-1;
    memset(y, 0, len_y*sizeof(double));
    for (i=0; i<len_y; i++) {
        x_start = MAX(0,i-len_h+1);
        x_end   = MIN(i+1,len_x);
        h_start = MIN(i,len_h-1);
        for(j=x_start; j<x_end; j++)
            y[i] += h[h_start--]*x[j];
    }
}

void convz
(
    double_complex* x,
    double_complex* h,
    int len_x,
    int len_h,
    double_complex* y
)
{
    int i, j, h_start, x_start, x_end, len_y;

    len_y = len_h+len_x-1;
    memset(y, 0, len_y*sizeof(double_complex));
    for (i=0; i<len_y; i++) {
        x_start = MAX(0,i-len_h+1);
        x_end   = MIN(i+1,len_x);
        h_start = MIN(i,len_h-1);
        for(j=x_start; j<x_end; j++)
            y[i] = ccadd(y[i], ccmul(h[h_start--], x[j]));
    }
}

void polyd_v
(
    double* x,
    double* poly,
    int len_x
)
{
    int j,i;

    memset(poly, 0, (len_x+1)*sizeof(double));
    poly[0] = 1.0;
    for (j=0; j<len_x; j++){
        for(i=j+1; i>0; i--){
            poly[i] = poly[i] - x[j] * (poly[i-1]);
        }
    } 
}

void polyz_v
(
    double_complex* x,
    double_complex* poly,
    int len_x
)
{
    int j,i;

    memset(poly, 0, (len_x+1)*sizeof(double_complex));
    poly[0] = cmplx(1.0, 0.0);
    for (j=0; j<len_x; j++){
        for(i=j+1; i>0; i--){
            poly[i] = ccsub(poly[i], ccmul(x[j], poly[i-1]));
        }
    }
}

void polyd_m
(
    double* X,
    double_complex* poly,
    int size_x
)
{
    int j,i;
    double_complex* Xcmplx, *e;

    /* Characteristic polynomial */
    Xcmplx = malloc1d(size_x*size_x*sizeof(double_complex));
    e = malloc1d(size_x*sizeof(double_complex));
    for(i=0; i<size_x*size_x; i++)
        Xcmplx[i] = cmplx(X[i], 0.0);
    utility_zeig(Xcmplx, size_x, NULL, NULL, NULL, e);

    /* recursion formula */
    memset(poly, 0, (size_x+1)*sizeof(double_complex));
    poly[0] = cmplx(1.0, 0.0);
    for (j=0; j<size_x; j++){
        for(i=j+1; i>0; i--){
            poly[i] = ccsub(poly[i], ccmul(e[j], poly[i-1]));
        }
    }

    /* clean-up */
    free(Xcmplx);
    free(e);
}
