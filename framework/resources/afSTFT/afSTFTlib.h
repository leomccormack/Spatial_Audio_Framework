/*
 Copyright (c) 2015 Juha Vilkamo
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

/**
 * @file afSTFTlib.h
 * @brief Slightly modified version of afSTFTlib
 *
 * The original afSTFT code, written by Juha Vilkamo, can be found here:
 * https://github.com/jvilkamo/afSTFT
 * This version is slightly modified. It adds a function to change the number of
 * channels on the fly and includes vectors for the hybrid mode centre
 * frequencies @44.1kHz/48kHz with 128 hop size for convenience.
 * It also supports the use of SAF utilities (for the vectorisation and FFT).
 *
 * Note that the design is also detailed in chapter 1 of [1]
 *
 * @see [1] Pulkki, V., Delikaris-Manias, S. and Politis, A. 2018. Parametric
 *          time--frequency domain spatial audio. John Wiley & Sons,
 *          Incorporated.
 *
 * @author Juha Vilkamo
 * @date 08.04.2015
 */

#ifndef __afSTFTlib_INCLUDED__
#define __afSTFTlib_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
   
/**
 * Remove the "AFSTFT_USE_SAF_UTILITIES" definition, and add the vecTools.h/c
 * and fft4g.h/c to your project, if you want to use the original afSTFT vector
 * code. Note the vecTools.h/c and fft4g.h/c files, may be found here:
 *     https://github.com/jvilkamo/afSTFT
 */
#define AFSTFT_USE_SAF_UTILITIES
#ifdef AFSTFT_USE_SAF_UTILITIES
# include "../../modules/saf_utilities/saf_utilities.h"
#endif
    
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Options for how the frequency domain data is permuted when using afSTFT */
typedef enum _AFSTFT_FDDATA_FORMAT{
    AFSTFT_BANDS_CH_TIME, /**< nBands x nChannels x nTimeHops */
    AFSTFT_TIME_CH_BANDS  /**< nTimeHops x nChannels x nBands */

}AFSTFT_FDDATA_FORMAT;
    
/**
 * afSTFT centre frequencies for 128 hop size and hybrid-mode enabled (48kHz) */
extern const double __afCenterFreq48e3[133];
/**
 * afSTFT centre frequencies for 128 hop size and hybrid-mode enabled (44.1kHz)
 */
extern const double __afCenterFreq44100[133];

/**
 * Creates an instance of the qmf filterbank
 *
 * @test test__afSTFT()
 *
 * @param[in] phQMF      (&) address of qmf handle
 * @param[in] nCHin      Number of input channels
 * @param[in] nCHout     Number of output channels
 * @param[in] hopsize    Hop size, in samples
 * @param[in] hybridmode 0: disabled, 1: hybrid-filtering enabled
 * @param[in] format     frequency-domain frame format, see #_QMF_FDDATA_FORMAT
 *                       enum
 */
// * @test test__afSTFT()
void afSTFT_create(void ** const phSTFT,
                   int nCHin,
                   int nCHout,
                   int hopsize,
                   int lowDelayMode,
                   int hybridmode,
                   AFSTFT_FDDATA_FORMAT format);

void afSTFT_destroy(void ** const phSTFT);

void afSTFT_forward(void * const hSTFT,
                    float** dataTD,
                    int framesize,
                    float_complex*** dataFD);

void afSTFT_backward(void * const hSTFT,
                     float_complex*** dataFD,
                     int framesize,
                     float** dataTD);


/**
 * Re-allocates memory to support a change in the number of input/output
 * channels */
void afSTFT_channelChange(void * const hSTFT, int new_nCHin, int new_nCHout);

/**
 * Flushes time-domain buffers with zeros. */
void afSTFT_clearBuffers(void * const hSTFT);

int afSTFT_getNBands(void * const hSTFT);

int afSTFT_getProcDelay(void * const hSTFT);


void afSTFT_getCentreFreqs(void * const hSTFT,
                           float fs,
                           int nBands,
                           float* freqVector);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* __afSTFTlib_INCLUDED__ */
