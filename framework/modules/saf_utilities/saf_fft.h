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

/**
 * @file saf_fft.h
 * @brief Wrappers for optimised fast Fourier transform (FFT) routines
 *
 * @note If none of the supported optimised FFT implementations are linked, then
 *       saf_fft employs the highly respectable KissFFT from here (BSD 3-Clause
 *       License): https://github.com/mborgerding/kissfft
 *
 * @note If linking Apple Accelerate: KissFFT is also used in cases where the
 *       FFT size is not 2^x.
 *
 * ## Dependencies
 *   Intel MKL, Apple Accelerate, or KissFFT (included in framework)
 *
 * @author Leo McCormack
 * @date 06.04.2019
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
#include <limits.h>
#include "saf_complex.h"

/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

/**
 * Calcuates the frequencies (in Hz) of uniformly spaced bins, for a given
 * FFT size and sampling rate.
 *
 * @param[in]  fftSize    FFT size
 * @param[in]  fs         Sampling rate
 * @param[out] freqVector 0:fs/(fftSize/2):fs/2; (fftSize/2+1) x 1
 */
void getUniformFreqVector(int fftSize,
                          float fs,
                          float* freqVector);

/**
 * FFT-based convolution of signal 'x' with filter 'h'.
 *
 * Input channels and filters are zero padded to avoid circular convolution
 * artefacts.
 *
 * @note The output must be of size: nCH x (x_len+h_len-1)
 *
 * @param[in]  x     Input(s); FLAT: nCH x x_len
 * @param[in]  h     Filter(s); FLAT: nCH x h_len
 * @param[in]  x_len Length of input signal, in samples
 * @param[in]  h_len Length of filter, in samples
 * @param[in]  nCH   Number of channels
 * @param[out] y     Output signal(s); FLAT: nCH x (x_len+h_len-1)
 */
void fftconv(float* x,
             float* h,
             int x_len,
             int h_len,
             int nCH,
             float* y);

/**
 * FFT-based convolution for FIR filters.
 *
 * Similar to fftconv, other than only the first x_len samples of y are
 * returned. It has parity with the 'fftfilt' function in Matlab, except it just
 * uses one big FFT (i.e. no overlap-add).
 *
 * @param[in]  x     Input(s); FLAT: nCH x x_len
 * @param[in]  h     Filter(s); FLAT: nCH x h_len
 * @param[in]  x_len Length of input signal, in samples
 * @param[in]  h_len Length of filter, in samples
 * @param[in]  nCH   Number of channels
 * @param[out] y     Output signal(s); FLAT: nCH x x_len
 */
void fftfilt(float* x,
             float* h,
             int x_len,
             int h_len,
             int nCH,
             float* y);

/**
 * Computes the discrete-time analytic signal via the Hilbert transform.
 *
 * The magnitude of the output is the envelope, and imaginary part is the
 * actual Hilbert transform. (Functionally identical to Matlab's 'hilbert'
 * function)
 *
 * @param[in]  x     Input; x_len x 1
 * @param[in]  x_len Length of input signal, in samples
 * @param[out] y     Output analytic signal; x_len x 1
 */
void hilbert(float_complex* x,
             int x_len,
             float_complex* y);
 

/* ========================================================================== */
/*                Real<->Half-Complex (Conjugate-Symmetric) FFT               */
/* ========================================================================== */

/**
 * Creates an instance of saf_rfft; real<->half-complex (conjugate-symmetric)
 * FFT
 *
 * @note Only Even FFT sizes are supported.
 *
 * ## Example Usage
 * \code{.c}
 *   const int N = 256;                    // FFT size
 *   float x_in[N];                        // input buffer (time-domain)
 *   x_in[0] = ... x_in[N-1] =             // fill with data
 *   float_complex x_out[(N/2+1)];         // output buffer (frequency-domain)
 *   float test[N];                        // test (time-domain)
 *   void *hFFT;                           // safFFT handle
 *
 *   saf_rfft_create(&hFFT, N);            // creates instance of safFFT
 *   saf_rfft_forward(hFFT, x_in, x_out);  // perform forward transform
 *   saf_rfft_backward(hFFT, x_out, test); // perform backwards transform
 *   // 'x_in' should equal 'test' (given some numerical error)
 *   saf_rfft_destroy(&hFFT);              // destroys instance of safFFT
 * \endcode
 *
 * @param[in] phFFT (&) address of saf_rfft handle
 * @param[in] N     FFT size
 */
void saf_rfft_create(void ** const phFFT,
                     int N);

/**
 * Destroys an instance of saf_rfft
 *
 * @param[in] phFFT (&) address of saf_rfft handle
 */
void saf_rfft_destroy(void ** const phFFT);

/**
 * Performs the forward-FFT operation; use for real to complex (conjugate
 * symmetric) transformations
 *
 * @note Only the first N/2 + 1 bins are returned in outputFD.
 *
 * @param[in]  hFFT     saf_rfft handle
 * @param[in]  inputTD  Time-domain input; N x 1
 * @param[out] outputFD Frequency-domain output; (N/2 + 1) x 1
 */
void saf_rfft_forward(void * const hFFT,
                      float* inputTD,
                      float_complex* outputFD);

/**
 * Performs the backward-FFT operation; use for complex (conjugate symmetric)
 * to real transformations
 *
 * @note Only the first N/2 + 1 bins are needed to be passed in inputFD.
 *
 * @param[in]  hFFT     saf_rfft handle
 * @param[in]  inputFD  Frequency-domain input; (N/2 + 1) x 1
 * @param[out] outputTD Time-domain output;  N x 1
 */
void saf_rfft_backward(void * const hFFT,
                       float_complex* inputFD,
                       float* outputTD);


/* ========================================================================== */
/*                            Complex<->Complex FFT                           */
/* ========================================================================== */

/**
 * Creates an instance of saf_fft; complex<->complex FFT
 *
 * @note Only Even FFT sizes are supported.
 *
 * @param[in] phFFT (&) address of saf_fft handle
 * @param[in] N     FFT size
 */
void saf_fft_create(void ** const phFFT,
                    int N);

/**
 * Destroys an instance of saf_fft
 *
 * @param[in] phFFT (&) address of saf_fft handle
 */
void saf_fft_destroy(void ** const phFFT);

/**
 * Performs the forward-FFT operation; use for complex to complex
 * transformations
 *
 * @param[in]  hFFT     saf_fft handle
 * @param[in]  inputTD  Time-domain input; N x 1
 * @param[out] outputFD Frequency-domain output; N x 1
 */
void saf_fft_forward(void * const hFFT,
                     float_complex* inputTD,
                     float_complex* outputFD);

/**
 * Performs the backward-FFT operation; use for complex to complex
 * transformations
 *
 * @param[in]  hFFT     saf_fft handle
 * @param[in]  inputFD  Frequency-domain input; N x 1
 * @param[out] outputTD Time-domain output;  N x 1
 */
void saf_fft_backward(void * const hFFT,
                      float_complex* inputFD,
                      float_complex* outputTD);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_FFT_H_INCLUDED */
