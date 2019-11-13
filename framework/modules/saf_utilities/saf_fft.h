/*
 * Copyright 2019 Leo McCormack
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

/*
 * Filename: saf_fft.h
 * -------------------
 * Wrappers for optimised fast Fourier transform (FFT) routines. If none are
 * linked, then it employs the highly respectable KissFFT from here
 * (BSD 3-Clause License): https://github.com/mborgerding/kissfft
 * If linking Apple Accelerate: KissFFT is also used in cases where the FFT size
 * is not a power of 2.
 *
 * Dependencies:
 *     Intel MKL, Apple Accelerate, or KissFFT (included in framework)
 * Author, date created:
 *     Leo McCormack, 06.04.2019
 */

#ifndef SAF_FFT_H_INCLUDED
#define SAF_FFT_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "saf_complex.h"
    
/*
 * Example Usage:
 * const int N = 256;                  // FFT size
 * float x_in[N];                      // input buffer (time-domain)
 * x_in[0] = ... x_in[N-1] =           // fill with data
 * float_complex x_out[(N/2+1)];       // output buffer (frequency-domain)
 * float test[N];                      // test (time-domain)
 * void *hFFT;                         // safFFT handle
 *
 * safFFT_create(&hFFT, N);            // creates an instance of safFFT, with size 'N'
 * safFFT_forward(hFFT, x_in, x_out);  // perform forward transform
 * safFFT_backward(hFFT, x_out, test); // perform backwards transform, (here x_in (should) = test)
 * safFFT_destroy(&hFFT);              // destroys an instance of safFFT
 */
    
/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

/*
 * Function: getUniformFreqVector
 * ------------------------------
 * Calcuates the frequencies (in Hz) of uniformly spaced bins, for a given
 * FFT size and sampling rate.
 *
 * Input Arguments:
 *     fftSize    - FFT size
 *     fs         - sampling rate
 * Output Arguments:
 *     freqVector - 0:fs/(fftSize/2):fs/2; (fftSize/2+1) x 1
 */
void getUniformFreqVector(int fftSize,
                          float fs,
                          float* freqVector);
    
 
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */
    
/*
 * Function: safFFT_create
 * -----------------------
 * Creates an instance of safFFT.
 * Note: Only Even FFT sizes are supported.
 *
 * Input Arguments:
 *     phFFT - & address of safFFT handle
 *     N     - FFT size
 */
void safFFT_create(void ** const phFFT,
                   int N);

/*
 * Function: safFFT_destroy
 * ------------------------
 * Destroys an instance of safFFT.
 *
 * Input Arguments:
 *     phFFT - & address of safFFT handle
 */
void safFFT_destroy(void ** const phFFT);

/*
 * Function: safFFT_forward
 * ------------------------
 * Performs the forward-FFT operation.
 *
 * Input Arguments:
 *     hFFT     - safFFT handle
 *     inputTD  - time-domain input; N x 1
 * Output Arguments:
 *     outputFD - frequency-domain output; (N/2 + 1) x 1
 */
void safFFT_forward(void * const hFFT,
                    float* inputTD,
                    float_complex* outputFD);

/*
 * Function: safFFT_backward
 * -------------------------
 * Performs the backward-FFT operation.
 *
 * Input Arguments:
 *     hFFT     - safFFT handle
 *     inputFD  - frequency-domain input; (N/2 + 1) x 1
 * Output Arguments:
 *     outputTD - time-domain output;  N x 1
 */
void safFFT_backward(void * const hFFT,
                     float_complex* inputFD,
                     float* outputTD);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_FFT_H_INCLUDED */
