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
 * Filename: saf_filters.h
 * -----------------------
 * Contains a collection of filter design equations.
 *
 * Dependencies:
 *     none
 * Author, date created:
 *     Leo McCormack, 01.03.2019
 */

#ifndef SAF_FILTERS_H_INCLUDED
#define SAF_FILTERS_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
    
/* ========================================================================== */
/*                                    Enums                                   */
/* ========================================================================== */

/*
 * Enum: BIQUAD_FILTER_TYPES
 * -------------------------
 * Bi-quadratic (second-order) filter design options
 *
 * Options:
 *     BIQUAD_FILTER_LPF       - low-pass filter
 *     BIQUAD_FILTER_HPF       - high-pass filter
 *     BIQUAD_FILTER_PEAK      - peaking filter
 *     BIQUAD_FILTER_LOW_SHELF - low-shelving filter
 *     BIQUAD_FILTER_HI_SHELF  - high-shelving filter
 */
typedef enum _BIQUAD_FILTER_TYPES {
    BIQUAD_FILTER_LPF,
    BIQUAD_FILTER_HPF,
    BIQUAD_FILTER_PEAK,
    BIQUAD_FILTER_LOW_SHELF,
    BIQUAD_FILTER_HI_SHELF
    
}BIQUAD_FILTER_TYPES;


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: biQuadCoeffs
 * ----------------------
 * Calculates 2nd order IIR filter coefficients [1]
 *
 * Input Arguments:
 *     filterType - see 'BIQUAD_FILTER_TYPES' enum
 *     fc         - centre frequency, Hz
 *     fs         - sampling frequency, Hz
 *     Q          - Q-factor
 *     gain_dB    - gain, dB
 * Output Arguments:
 *     b          - b filter coefficients; 3 x 1
 *     a          - a filter coefficients; 3 x 1
 *
 * [1] Zölzer, U. (Ed.). (2011). DAFX: digital audio effects. John Wiley & Sons.
 */
void biQuadCoeffs(/* input arguments */
                  BIQUAD_FILTER_TYPES filterType,
                  float fc,
                  float fs,
                  float Q,
                  float gain_dB,
                  /* output arguments */
                  float b[3],
                  float a[3]);
    
/*
 * Function: applyBiQuadFilter
 * ---------------------------
 * Applies biQuad filter to an input signal using the direct form II difference
 * equation: https://en.wikipedia.org/wiki/Digital_biquad_filter
 * Note: input 'signal' is filtered in place (i.e. it becomes the output signal)
 *
 * Input Arguments:
 *     b        - b filter coefficients; 3 x 1
 *     a        - a filter coefficients; 3 x 1
 *     w_z_12   - previous 2 wn samples (init as 0s); 2 x 1
 *     signal   - signal to be filtered/filtered signal; nSamples x 1
 *     gain_dB  - gain, dB
 *     nSamples - number of samples in the signal
 */
void applyBiQuadFilter(/* Input arguments */
                       float b[3],
                       float a[3],
                       float w_z_12[2],
                       float* signal,
                       int nSamples);
    
/*
 * Function: evalBiQuadTransferFunction
 * ------------------------------------
 * Evaluates the 2nd order IIR transfer function at one or more frequencies,
 * returning its magnitude and/or phase response
 *
 * Input Arguments:
 *     b            - b filter coefficients; 3 x 1
 *     a            - a filter coefficients; 3 x 1
 *     freqs        - frequencies at which to evaluate, Hz; nFreqs x 1
 *     nFreqs       - number of frequencies at which to avaluate
 *     fs           - sampling frequency, Hz
 * Output arguments
 *     magnitude_dB - magnitude, dB, at each frequency (set to NULL of not
 *                    wanted); nFreqs x 1
 *     phase_rad    - phase, radians, at each frequency (set to NULL of not
 *                    wanted); nFreqs x 1
 */
void evalBiQuadTransferFunction(/* Input arguments */
                                float b[3],
                                float a[3],
                                float* freqs,
                                int nFreqs,
                                float fs,
                                /* Output arguments */
                                float* magnitude_dB,
                                float* phase_rad);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_FILTERS_H_INCLUDED */
