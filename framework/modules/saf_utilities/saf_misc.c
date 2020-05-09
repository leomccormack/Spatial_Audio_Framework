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
 * @file saf_misc.c
 * @brief Miscellaneous functions
 *
 * @author Leo McCormack
 * @date 29.01.2020
 */

#include "saf_misc.h"
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
