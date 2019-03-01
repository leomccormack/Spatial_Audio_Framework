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
 *     saf_filters.h
 * Description:
 *     Contains a collection of filter design equations.
 * Dependencies:
 *     Windows users only: Intel's MKL must be installed, which can be freely aquired via:
 *     https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library.
 * Author, date created:
 *     Leo McCormack, 01.03.2019
 */

#ifndef SAF_FILTERS_H_INCLUDED
#define SAF_FILTERS_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif
    
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
    
#ifndef M_PI
#define M_PI ( 3.14159265358979323846264338327950288f )
#endif
    
/****************/
/* Enum options */
/****************/

typedef enum _BIQUAD_FILTER_TYPES {
    BIQUAD_FILTER_LPF,         /* low-pass filter */
    BIQUAD_FILTER_HPF,         /* high-pass filter */
    BIQUAD_FILTER_PEAK,        /* peaking filter */
    BIQUAD_FILTER_LOW_SHELF,   /* low-shelving filter */
    BIQUAD_FILTER_HI_SHELF     /* high-shelving filter */
    
}BIQUAD_FILTER_TYPES;


/******************/
/* Main Functions */
/******************/
    
/* calculates 2nd order IIR filter coefficients */
void biQuadCoeffs(/* input arguments */
                  BIQUAD_FILTER_TYPES filterType,     /* see 'BIQUAD_FILTER_TYPES' enum */
                  float fc,                           /* centre frequency, Hz */
                  float fs,                           /* sampling frequency, Hz */
                  float Q,                            /* Q-factor  */
                  float gain_dB,                      /* gain, dB */
                  /* output arguments */
                  float b[3],                         /* b filter coefficients */
                  float a[3]);                        /* a filter coefficients */
    
/* applies biQuad filter to an input signal using the direct form II difference equation:
 * https://en.wikipedia.org/wiki/Digital_biquad_filter
 * input 'signal' is filtered in place (i.e. it becomes the output signal) */
void applyBiQuadFilter(float b[3],                    /* b filter coefficients */
                       float a[3],                    /* a filter coefficients */
                       float w_z_12[2],               /* previous 2 wn samples (init as 0s) */ 
                       float* signal,                 /* signal to be filtered/filtered signal; nSamples x 1 */
                       int nSamples);                 /* number of samples in the signal */
    
/* evaluates the 2nd order IIR transfer function at one or more frequencies, returning the
 * magnitude and/or phase response. */
void evalBiQuadTransferFunction(/* input arguments */
                                float b[3],           /* b filter coefficients */
                                float a[3],           /* a filter coefficients */
                                float* freqs,         /* frequencies, Hz; nFreqs x 1 */
                                int nFreqs,           /* number of frequencies at which to avaluate */
                                float fs,             /* sampling frequency, Hz */
                                /* output arguments */
                                float* magnitude_dB,  /* magnitude, dB, at each frequency (set to NULL of not wanted); nFreqs x 1 */
                                float* phase_rad);    /* phase, radians, at each frequency (set to NULL of not wanted); nFreqs x 1 */


#ifdef __cplusplus
}/* extern "C" */
#endif

#endif /* SAF_FILTERS_H_INCLUDED */






