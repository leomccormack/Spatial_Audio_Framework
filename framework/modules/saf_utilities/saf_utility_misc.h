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
 *@addtogroup Utilities
 *@{
 * @file saf_utility_misc.h
 * @brief A collection of miscellaneous functions
 *
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
 * https://cboard.cprogramming.com/c-programming/170381-cross-platform-wait-sleep.html */
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
 * Computes Lagrange interpolation weights of order 'N' for value 'x'
 *
 * @param[in]  N        Order
 * @param[in]  x        Values; len_x x 1
 * @param[in]  len_x    Number of values
 * @param[out] weights  Weights; (ORDER+1) x len_x
 */
void lagrangeWeights(int N,
                     float* x,
                     int len_x,
                     float* weights);

/**
 * This function takes a frequency vector and groups its frequencies into
 * critical bands [Equivalent-Rectangular Bandwidth (ERB)].
 *
 * e.g.
 *   - centerFreq[erb_idx[0]] -> centerFreq[erb_idx[1]] is ERB band 1
 *   - centerFreq[erb_idx[1]] -> centerFreq[erb_idx[2]] is ERB band 2
 *
 * @param[in]  centerFreq Frequency vector; nBands x 1
 * @param[in]  nBands     Number of bins/bands in frequency vector
 * @param[in]  maxFreqLim Past this frequency the bands are grouped into 1 band
 * @param[out] erb_idx    (&) ERB indices; nERBBands x 1
 * @param[out] erb_freqs  (&) ERB frequencies; nERBBands x 1
 * @param[out] nERBBands  (&) Number of ERB bands; 1 x 1
 */
void findERBpartitions(/* Input Arguments */
                       float* centerFreq,
                       int nBands,
                       float maxFreqLim,
                       /* Output Arguments */
                       int** erb_idx,
                       float** erb_freqs,
                       int* nERBBands);

/**
 * Returns the indices required to randomly permute a vector of length 'len'
 */
void randperm(/* Input Arguments */
              int len,
              /* Output Arguments */
              int* randperm_inds);

/**
 * Factorial, accurate up to n<=25
 *
 * @note The magnitude will still be correct above 25, but the precision will be
 *       truncated. The function also returns pre-computed values up to n=15
 *       to make it faster (e.g. for up to 7th order SH computations...).
 *
 * @param[in] n Order
 * @returns factorial(n)
 */
long double factorial(int n);

/**
 * C fmodf function, except it behaves like 'mod' in Matlab (i.e. with the wrap
 * around for negative values of 'x')
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

/**
 * Returns the sum of all values
 */
float sumf(float* values,
           int nValues);

/**
 * Finds the unique values (and their indices) of the input vector
 *
 * @note this is equivalent to using "unique(vals, 'last')" in Matlab
 * @test test__unique_i()
 *
 * @param[in]  input      Input vector; nInputs x 1
 * @param[in]  nInputs    Number of elements in the input vector
 * @param[out] uniqueVals (&) Unique values (set to NULL if not wanted);
 *                        nUnique x 1
 * @param[out] uniqueInds (&) Indices corresponding to Unique values (set to
 *                        NULL if not wanted); nUnique x 1
 * @param[out] nUnique    (&) Number of Unique values; 1 x 1
 *
 */
void unique_i(int* input,
              int nInputs,
              int** uniqueVals,
              int** uniqueInds,
              int* nUnique);


/**
 * Numerically solves first-order, linear, homogeneous differential equation
 * systems, with non-constant coefficients, by generalization of the Pade-
 * approximant method for exponential matrices.
 *
 * The equations are described in matrix form as
 *     Y'(t) = D(t)*Y(t)
 * where D and Y are square-matrix functions of scalar t. The initial condition
 * is Y(0) = I (the identity matrix), and the result is Y(1). For the special
 * case of a constant coefficient matrix D, gexpm is equivalent to the standard
 * matrix exponential (expm).
 *
 * m1: true or false, optional, default should be false.
 *     "minus 1" flag: if m1 = false the generalized exponential is Y; if
 *     true it is Y+I. gexpm is analogous to the expm1 function
 *     ("exponential-minus-1") when m1 is true.
 *
 * @note For both cases of constant and non-constant D, the solution is
 *       determined from a Pade approximation of order 6, using scale-and-square
 *       for constant D and multi-step integration for non-constant D. The form
 *       of the Pade approximant is outlined in the associated document
 *       KJohnson_2015_04_01.pdf. Notes on error control are in the code
 *       comments.
 *
 * This function is based on the Matlab script found (BSD-3-clause license):
 * mathworks.com/matlabcentral/fileexchange/50413-generalized-matrix-exponential
 *
 * Copyright (c) 2015, Kenneth Johnson
 */
void gexpm(float* D,
           int sizeD,
           int m1,
           float* Y);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_MISC_H_INCLUDED */

/**@} */ /* doxygen addtogroup Utilities */
