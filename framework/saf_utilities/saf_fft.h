/*
 Copyright 2019 Leo McCormack
 
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
 *     saf_fft.h
 * Description:
 *     Wrapper for optimised fast Fourier transforms (FFT).
 * Dependencies:
 *     Windows users only: custom Intel MKL '.lib/.dll' files are required.
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library by default.
 *     However, Mac users may elect to use a custom Intel MKL '.dylib' instead.
 *     Further instructions for both Windows/Mac users can be found here:
 *     https://github.com/leomccormack/Spatial_Audio_Framework
 * Author, date created:
 *     Leo McCormack, 06.04.2019
 */

#ifndef SAF_FFT_H_INCLUDED
#define SAF_FFT_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif
    
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "saf_complex.h"
    
/* NOTE: vDSP_fft hasn't been extensively tested, and doesn't seem to return the Nyquist value?! */
/* INTEL MKL Allows for any even FFT size, vDSP must be 2^(int value) */
    
/*
 * Example Usage:
 * const int N = 256;                  // FFT size
 * float x_in[N];                      // input  (time-domain)
 * float_complex x_out[(N/2+1)];       // output (frequency-domain)
 * float test[N];                      // test (time-domain)
 * void *hFFT;                         // safFFT handle
 *
 * safFFT_create(&hFFT, N);            // creates an instance of safFFT, with size 'N'
 * safFFT_forward(hFFT, x_in, x_out);  // perform forward transform
 * safFFT_backward(hFFT, x_out, test); // perform backwards transform, (here x_in = test)
 * safFFT_destroy(&hFFT);              // destroys an instance of safFFT
 */
 
/* creates an instance of safFFT. Only 2^(int value) FFT sizes are supported for vDSP. */
void safFFT_create(void ** const phFFT,         /* address of safFFT handle */
                   int N);                      /* FFT size */

/* destroys an instance ofsafFFT */
void safFFT_destroy(void ** const phFFT);       /* address of safFFT handle */
    
/* performs the forward-FFT operation.  */
void safFFT_forward(void * const hFFT,          /* safFFT handle */
                    float* inputTD,             /* time-domain input; N x 1 */
                    float_complex* outputFD);   /* frequency-domain output; (N/2 + 1) x 1 */
    
/* performs the backward-FFT operation.  */
void safFFT_backward(void * const hFFT,         /* safFFT handle */
                     float_complex* inputFD,    /* frequency-domain input; (N/2 + 1) x 1 */
                     float* outputTD);          /* time-domain output; N x 1 */
    

#ifdef __cplusplus
}/* extern "C" */
#endif

#endif /* SAF_FFT_H_INCLUDED */






