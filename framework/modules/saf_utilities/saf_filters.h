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
 * Bi-quadratic (second-order) IIR filter design options.
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

/*
 * Enum: FIR_FILTER_TYPES
 * ----------------------
 * Finite Impulse Response (FIR) filter design options.
 *
 * Options:
 *     FIR_FILTER_LPF - low-pass filter
 *     FIR_FILTER_HPF - high-pass filter
 *     FIR_FILTER_BPF - band-pass filter
 *     FIR_FILTER_BSF - band-stop filter
 */
typedef enum _FIR_FILTER_TYPES {
    FIR_FILTER_LPF,
    FIR_FILTER_HPF,
    FIR_FILTER_BPF,
    FIR_FILTER_BSF
    
}FIR_FILTER_TYPES;

/*
 * Enum: WINDOWING_FUNCTION_TYPES
 * ------------------------------
 * Windowing function types. Symmetric if winlength is odd, and asymmetric if
 * winlength is even. Windows are evaluated: 0 <= n < winlength.
 * Largely taken from: https://en.wikipedia.org/wiki/Window_function
 *
 * Options:
 *     WINDOWING_FUNCTION_RECTANGULAR
 *     WINDOWING_FUNCTION_HAMMING
 *     WINDOWING_FUNCTION_HANN
 *     WINDOWING_FUNCTION_BLACKMAN
 *     WINDOWING_FUNCTION_NUTTALL
 *     WINDOWING_FUNCTION_BLACKMAN_NUTTALL
 *     WINDOWING_FUNCTION_BLACKMAN_HARRIS  
 */
typedef enum _WINDOWING_FUNCTION_TYPES {
    WINDOWING_FUNCTION_RECTANGULAR,
    WINDOWING_FUNCTION_HAMMING,
    WINDOWING_FUNCTION_HANN,
    WINDOWING_FUNCTION_BARTLETT,
    WINDOWING_FUNCTION_BLACKMAN,
    WINDOWING_FUNCTION_NUTTALL,
    WINDOWING_FUNCTION_BLACKMAN_NUTTALL,
    WINDOWING_FUNCTION_BLACKMAN_HARRIS
    
}WINDOWING_FUNCTION_TYPES;


/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

/*
 * Function: getWindowingFunction
 * ------------------------------
 * Computes the weights of a specific windowing function.
 * Weights symmetric if winlength is odd, and asymmetric if winlength is even
 * i.e. if winlength is even:
 *  - index "winlength/2" = 1, and first value!=last value
 * if odd:
 *  - index "(winlength-1)/2" = 1, and first value==last value
 *
 * Input Arguments:
 *     type      - see 'WINDOWING_FUNCTION_TYPES' enum
 *     winlength - window length in samples
 * Output Arguments:
 *     win       - windowing function; winlength x 1
 */
void getWindowingFunction(WINDOWING_FUNCTION_TYPES type,
                          int winlength,
                          float* win);

/*
 * Function: getOctaveBandCutoffFreqs
 * ----------------------------------
 * Converts octave band CENTRE frequencies into CUTOFF frequencies.
 * Note: the lower and upper CENTRE frequencies only have their upper and lower
 * CUTOFF frequencies computed, respectively. e.g.:
 *   centreFreqs[6] = { 125, 250, 500, 1000, 2000, 4000 }, becomes:
 *   cutoffFreqs[5] = { 176, 354, 707, 1410, 2830 }
 * Passing cutoffFreqs[5] to "FIRFilterbank", will give filter coefficients for
 * the following:
 *    Band1: LPF @ 176Hz
 *    Band2: BFP @ 176-354Hz
 *    Band3: BFP @ 354-707Hz
 *    Band4: BFP @ 707-1410Hz
 *    Band5: BFP @ 1410-2830Hz
 *    Band6: HPF @ 2830Hz
 * (Basically, band 125Hz also encapsulates everything down to DC, and band 4kHz
 * also encapsulates everything up to Nyquist)
 * Also note: cutoffFreqs vector is shorter than centreFreqs by 1 element.
 *
 * Input Arguments:
 *     centreFreqs  - centre frequencies (octave bands); nCentreFreqs x 1
 *     nCentreFreqs - number of centre frequencies
 * Output Arguments:
 *     cutoffFreqs  - cutoff frequencies, which encapsulate the specified centre
 *                    frequencies by 1 octave; (nCentreFreqs-1) x 1
 */
void getOctaveBandCutoffFreqs(float* centreFreqs,
                              int nCentreFreqs,
                              float* cutoffFreqs);

/*
 * Function: flattenMinphase
 * -------------------------
 * Equalises input sequence by its minimum phase form, in order to bring its
 * magnitude response to unity.
 *
 * Input/output Arguments:
 *     x   - input; len x 1
 * Input Arguments:
 *     len - length of input.
 */
void flattenMinphase(float* x,
                     int len);


/* ========================================================================== */
/*                              Bi-Quad Functions                             */
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


/* ========================================================================== */
/*                            FIR Filter Functions                            */
/* ========================================================================== */

/*
 * Function: FIRCoeffs
 * -------------------
 * FIR filters by windowing. When using the Hamming window, and scalingFLAG=1,
 * the function is numerically identical to the default 'fir1' function in
 * Matlab (when using it in single precision mode) [1].
 * Note: input argument 'order' cannot be odd valued.
 *
 * Some guidelines regarding the approx order (N) for certain filters. i.e.
 * the orders where you actually get the expected -6dB attenuation at the cutoff
 * frequency specified (fs=48kHz, Hamming window, scalingFLAG=1).
 * (Use these figures only just to get a rough idea)
 *  - LPF @ 100Hz - N~1400
 *  - LPF @ 250Hz - N~550
 *  - LPF @ 1kHz  - N~150
 *  - LPF @ 4kHz  - N~40
 *  - BPF @ 88-176Hz   - N~2500
 *  - BPF @ 176-354Hz  - N~1600
 *  - BPF @ 707-1410Hz - N~400
 *  - HPF @ 200Hz - N~450
 *  - HPF @ 4kHz  - N~60
 *
 * Input Arguments:
 *     filterType  - see 'FIR_FILTER_TYPES' enum
 *     order       - filter order (N). Must be even.
 *     cutoff1     - filter1 cutoff in Hz, for LPF/HPF, and lower cutoff for
 *                   BPF/BSF
 *     cutoff2     - filter2 cutoff in Hz, not needed for LPF/HPF, this is the
 *                   upper cutoff for BPF/BSF
 *     samplerate  - sampling rate in Hz
 *     windowType  - see 'WINDOWING_FUNCTION_TYPES' enum
 *     scalingFLAG - 0: none, 1: scaling applied to ensure passband is at 0dB
 * Output Arguments:
 *     filter      - filter coefficients/weights/taps; (order+1) x 1
 *
 * [1] "Programs for Digital Signal Processing", IEEE Press John Wiley & Sons,
 *     1979, pg. 5.2-1.
 */
void FIRCoeffs(FIR_FILTER_TYPES filterType,
               int order,
               float cutoff1,
               float cutoff2,
               float sampleRate,
               WINDOWING_FUNCTION_TYPES windowType,
               int scalingFLAG,
               float* filter);

/*
 * Function: FIRFilterbank
 * -----------------------
 * Returns a bank of FIR filter coefficients required to divide a signal into
 * frequency bands. Provided the order is sufficient, the sum of the bands
 * should recontruct the original (although, shifted in time due to group delay)
 * e.g fc[1] = { 1000 };
 *   Band1, &filter[0*(order+1)] : LPF @ 1kHz
 *   Band2, &filter[1*(order+1)] : HPF @ 1kHz
 * e.g fc[3] = { 1000, 2000, 4000 };
 *   Band1, &filter[0*(order+1)] : LPF @ 1kHz
 *   Band2, &filter[1*(order+1)] : BPF @ 1-2kHz
 *   Band3, &filter[2*(order+1)] : BPF @ 2-4kHz
 *   Band4, &filter[3*(order+1)] : HPF @ 4kHz
 *
 * Input Arguments:
 *     order        - filter order. Must be even.
 *     fc           - vector of cutoff frequencies; nCutoffFreqs x 1
 *     nCutoffFreqs - number of cutoff frequencies in vector 'fc'.
 *     samplerate   - sampling rate in Hz
 *     windowType   - see 'WINDOWING_FUNCTION_TYPES' enum
 *     scalingFLAG  - 0: none, 1: scaling applied to ensure passbands are at 0dB
 * Output Arguments:
 *     filter       - filter coefficients/weights/taps;
 *                    FLAT: (nCutoffFreqs+1) x (order+1)
 */
void FIRFilterbank(int order,
                   float* fc,
                   int nCutoffFreqs,
                   float sampleRate,
                   WINDOWING_FUNCTION_TYPES windowType,
                   int scalingFLAG,
                   float* filterbank);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_FILTERS_H_INCLUDED */
