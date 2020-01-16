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
#include <limits.h>
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
    
/*
 * Function: fftconv
 * -----------------
 * FFT-based convolution. Input channels and filters are zero padded to avoid
 * circular convolution artefacts.
 * Note: the output must be of size: nCH x (x_len+h_len-1)
 *
 * Input Arguments:
 *     x     - input(s); FLAT: nCH x x_len
 *     h     - filter(s); FLAT: nCH x h_len
 *     x_len - length of input signal, in samples
 *     h_len - length of filter, in samples
 *     nCH   - number of channels
 * Output Arguments:
 *     y     - output signal(s); FLAT: nCH x (x_len+h_len-1)
 */
void fftconv(float* x,
             float* h,
             int x_len,
             int h_len,
             int nCH,
             float* y);

/*
 * Function: fftfilt
 * -----------------
 * FFT-based convolution for FIR filters. Similar to fftconv, other than only
 * the first x_len samples of y are returned. It has parity with the 'fftfilt'
 * function in Matlab, except it just uses one big FFT (i.e. no overlap-add).
 *
 * Input Arguments:
 *     x     - input(s); FLAT: nCH x x_len
 *     h     - filter(s); FLAT: nCH x h_len
 *     x_len - length of input signal, in samples
 *     h_len - length of filter, in samples
 *     nCH   - number of channels
 * Output Arguments:
 *     y     - output signal(s); FLAT: nCH x x_len
 */
void fftfilt(float* x,
             float* h,
             int x_len,
             int h_len,
             int nCH,
             float* y);

/*
 * Function: hilbert
 * -----------------
 * Computes the discrete-time analytic signal via the Hilbert transform.
 * The magnitude of the output is the envelope, and imaginary part is the
 * actual Hilbert transform.
 * (Functionally identical to Matlab's 'hilbert' function)
 *
 * Input Arguments:
 *     x     - input; x_len x 1
 *     x_len - length of input signal, in samples
 * Output Arguments:
 *     y     - output analytic signal; x_len x 1
 */
void hilbert(float_complex* x,
             int x_len,
             float_complex* y);


/* ========================================================================== */
/*                Real<->Half-Complex (Conjugate-Symmetric) FFT               */
/* ========================================================================== */

/*
 * Function: saf_rfft_create
 * -------------------------
 * Creates an instance of saf_rfft.
 * Note: Only Even FFT sizes are supported.
 *
 * Input Arguments:
 *     phFFT - & address of saf_rfft handle
 *     N     - FFT size
 */
void saf_rfft_create(void ** const phFFT,
                     int N);

/*
 * Function: saf_rfft_destroy
 * --------------------------
 * Destroys an instance of saf_rfft.
 *
 * Input Arguments:
 *     phFFT - & address of saf_rfft handle
 */
void saf_rfft_destroy(void ** const phFFT);

/*
 * Function: saf_rfft_forward
 * --------------------------
 * Performs the forward-FFT operation. Use for real to complex (conjugate
 * symmetric) transformations.
 * Note: only the first N/2 + 1 bins are returned in outputFD.
 *
 * Input Arguments:
 *     hFFT     - saf_rfft handle
 *     inputTD  - time-domain input; N x 1
 * Output Arguments:
 *     outputFD - frequency-domain output; (N/2 + 1) x 1
 */
void saf_rfft_forward(void * const hFFT,
                      float* inputTD,
                      float_complex* outputFD);

/*
 * Function: saf_rfft_backward
 * ---------------------------
 * Performs the backward-FFT operation. Use for complex (conjugate symmetric)
 * to real transformations.
 * Note: only the first N/2 + 1 bins are needed to be passed in inputFD.
 *
 * Input Arguments:
 *     hFFT     - saf_rfft handle
 *     inputFD  - frequency-domain input; (N/2 + 1) x 1
 * Output Arguments:
 *     outputTD - time-domain output;  N x 1
 */
void saf_rfft_backward(void * const hFFT,
                       float_complex* inputFD,
                       float* outputTD);


/* ========================================================================== */
/*                            Complex<->Complex FFT                           */
/* ========================================================================== */

/*
 * Function: saf_fft_create
 * ------------------------
 * Creates an instance of saf_fft.
 * Note: Only Even FFT sizes are supported.
 *
 * Input Arguments:
 *     phFFT - & address of saf_fft handle
 *     N     - FFT size
 */
void saf_fft_create(void ** const phFFT,
                    int N);

/*
 * Function: saf_fft_destroy
 * -------------------------
 * Destroys an instance of saf_fft.
 *
 * Input Arguments:
 *     phFFT - & address of saf_fft handle
 */
void saf_fft_destroy(void ** const phFFT);

/*
 * Function: saf_fft_forward
 * -------------------------
 * Performs the forward-FFT operation. Use for complex to complex
 * transformations.
 *
 * Input Arguments:
 *     hFFT     - saf_rfft handle
 *     inputTD  - time-domain input; N x 1
 * Output Arguments:
 *     outputFD - frequency-domain output; N x 1
 */
void saf_fft_forward(void * const hFFT,
                     float_complex* inputTD,
                     float_complex* outputFD);

/*
 * Function: saf_fft_backward
 * ---------------------------
 * Performs the backward-FFT operation. Use for complex to complex
 * transformations.
 *
 * Input Arguments:
 *     hFFT     - saf_rfft handle
 *     inputFD  - frequency-domain input; N x 1
 * Output Arguments:
 *     outputTD - time-domain output;  N x 1
 */
void saf_fft_backward(void * const hFFT,
                      float_complex* inputFD,
                      float_complex* outputTD);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_FFT_H_INCLUDED */
