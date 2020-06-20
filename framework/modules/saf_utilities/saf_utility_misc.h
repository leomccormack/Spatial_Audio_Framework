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
 * @file saf_utility_misc.h
 * @brief A collection of miscellaneous functions
 * @author Leo McCormack
 * @date 29.01.2020 
 */

#ifndef SAF_MISC_H_INCLUDED
#define SAF_MISC_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "saf_utilities.h"

/* Cross-platform sleep macro (slightly modified) from:
 https://cboard.cprogramming.com/c-programming/170381-cross-platform-wait-sleep.html */
#ifdef _WIN32
    /* For Windows (32- and 64-bit) */
# include <windows.h>
# define SAF_SLEEP(msecs) Sleep(msecs)
#elif defined(__unix) || defined(__APPLE__)
    /* For linux, OSX, and other unixes */
# ifndef _POSIX_C_SOURCE
#  define _POSIX_C_SOURCE 199309L /* or greater */
# endif
# include <time.h>
# define SAF_SLEEP(msecs) do {     \
   struct timespec ts;             \
   ts.tv_sec = msecs/1000;         \
   ts.tv_nsec = msecs%1000*1000;   \
   nanosleep(&ts, NULL);           \
} while (0)
#else
# error "Unknown system"
#endif
#include <stdio.h>
#include <string.h>
#include <math.h>

/**
 * Factorial, accurate up to n<=25
 *
 * @note The magnitude will still be correct >25, but the precision will be
 *       truncated. The function also returns pre-computed values up to n=15
 *       to make it faster.
 *
 * @param[in] n Order
 * @returns factorial(n)
 */
long double factorial(int n);

/**
 * C fmodf function, which behaves like 'mod' in Matlab (i.e. with the wrap
 * around for negative 'x')
 *
 * @param[in] x Value 'x'
 * @param[in] y Value 'y'
 * @returns fmodf(n)
 */
float matlab_fmodf(float x,  float y);

/**
 * Calculates the cross correlation between two vectors
 *
 * @param[in]  a    Vector a; la x 1
 * @param[in]  b    Vector b; lb x 1
 * @param[in]  la   Length of vector a
 * @param[in]  lb   Length of vector b
 * @param[out] x_ab Cross-correlation between a and b; (la + lb - 1) x 1
 */
void cxcorr(float* a,
            float* b,
            float* x_ab,
            size_t la,
            size_t lb);

/**
 * Generates random numbers between -1 and 1 and stores them in the input vector
 *
 * @param[in,out] vector Vector to populate with random numbers; length x 1
 * @param[in]     length Length of the vector
 */
void rand_m1_1(float* vector,
               int length);

/**
 * Generates random numbers between -1 and 1 and stores them in the input vector
 * for both the real and imaginary parts
 *
 * @param[in,out] vector Vector to populate with random numbers; length x 1
 * @param[in]     length Length of the vector
 */
void rand_cmplx_m1_1(float_complex* vector,
                     int length);

/**
 * Generates random numbers between 0 and 1 and stores them in the input vector
 *
 * @param[in,out] vector Vector to populate with random numbers; length x 1
 * @param[in]     length Length of the vector
 */
void rand_0_1(float* vector,
              int length);

/**
 * Basic 1-D direct convolution in the time-domain (real double precision)
 *
 * @param[in]  x     Input sequence; len_x x 1
 * @param[in]  h     Filter sequence; len_h x 1
 * @param[in]  len_x Length of 'x'
 * @param[in]  len_h Length of 'h'
 * @param[out] y     Output sequence; (len_x+len_h-1) x 1
 */
void convd(double* x,
           double* h,
           int len_x,
           int len_h,
           double* y);

/**
 * Basic 1-D direct convolution in the time-domain (complex double precision)
 *
 * @param[in]  x     Input sequence; len_x x 1
 * @param[in]  h     Filter sequence; len_h x 1
 * @param[in]  len_x Length of 'x'
 * @param[in]  len_h Length of 'h'
 * @param[out] y     Output sequence; (len_x+len_h-1) x 1
 */
void convz(double_complex* x,
           double_complex* h,
           int len_x,
           int len_h,
           double_complex* y);

/**
 * Convert roots of a vector to polynomial (real double precision)
 *
 * @param[in]  x     Input vector; len_x x 1
 * @param[out] poly  Polynomials; (len_x+1) x 1
 * @param[in]  len_x Length of vector 'x'
 */
void polyd_v(double* x,
             double* poly,
             int len_x);

/**
 * Convert roots of a vector to polynomial (complex double precision)
 *
 * @param[in]  x     Input vector; len_x x 1
 * @param[out] poly  Polynomials; (len_x+1) x 1
 * @param[in]  len_x Length of vector 'x'
 */
void polyz_v(double_complex* x,
             double_complex* poly,
             int len_x);

/**
 * Convert roots of a matrix to polynomial (real double precision)
 *
 * @param[in]  X      Square input matrix; size_x x size_x
 * @param[out] poly   Polynomials; (size_x+1) x 1
 * @param[in]  size_x Dimensions of square matrix 'X'
 */
void polyd_m(double* X,
             double_complex* poly,
             int size_x);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_MISC_H_INCLUDED */
