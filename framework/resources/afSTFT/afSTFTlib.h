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
 * The original afSTFT code written by Juha Vilkamo can be found here:
 * https://github.com/jvilkamo/afSTFT
 * This version is slightly modified to be more in-line with how the rest of SAF
 * is structured.
 * The files afSTFTlib.h/.c act as the interface to afSTFT, which is then
 * implemented in afSTFT_internal.h/.c.
 *
 * This version also adds functionality to change the number of channels on the
 * fly, flush the run-time buffers with zeros, return the current frequency
 * vector and the current processing delay.
 * It also incorporates SAF utilities (for the vectorisation and FFT).
 *
 * Note that the afSTFT design is layed out in detail in chapter 1 of [1]
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
 * Options for how the frequency domain data is permuted when using afSTFT
 */
typedef enum _AFSTFT_FDDATA_FORMAT{
    AFSTFT_BANDS_CH_TIME, /**< nBands x nChannels x nTimeHops */
    AFSTFT_TIME_CH_BANDS  /**< nTimeHops x nChannels x nBands */

}AFSTFT_FDDATA_FORMAT;
    
/**
 * Creates an instance of afSTFT
 *
 * @test test__afSTFT()
 *
 * @param[in] phSTFT       (&) address of afSTFT handle
 * @param[in] nCHin        Number of input channels
 * @param[in] nCHout       Number of output channels
 * @param[in] hopsize      Hop size, in samples
 * @param[in] lowDelayMode 0: disabled, 1: low-delay mode enabled
 * @param[in] hybridmode   0: disabled, 1: hybrid-filtering enabled
 * @param[in] format       Frequency-domain frame format, see
 *                         #_AFSTFT_FDDATA_FORMAT enum
 */
void afSTFT_create(void ** const phSTFT,
                   int nCHin,
                   int nCHout,
                   int hopsize,
                   int lowDelayMode,
                   int hybridmode,
                   AFSTFT_FDDATA_FORMAT format);

/**
 * Destroys an instance of afSTFT
 *
 * @param[in] phSTFT  (&) address of afSTFT handle
 */
void afSTFT_destroy(void ** const phSTFT);

/**
 * Performs forward afSTFT transform
 *
 * @param[in]  hSTFT     afSTFT handle
 * @param[in]  dataTD    Time-domain input; nCHin x framesize
 * @param[in]  framesize Frame size of time-domain data
 * @param[out] dataFD    Frequency-domain output; see #_AFSTFT_FDDATA_FORMAT
 *                       enum
 */
void afSTFT_forward(void * const hSTFT,
                    float** dataTD,
                    int framesize,
                    float_complex*** dataFD);

/**
 * Performs backward afSTFT transform
 *
 * @param[in]  hSTFT     afSTFT handle
 * @param[in]  dataFD    Frequency-domain input; see #_AFSTFT_FDDATA_FORMAT enum
 * @param[in]  framesize Frame size of time-domain data
 * @param[out] dataTD    Time-domain output;  nCHout x framesize
 */
void afSTFT_backward(void * const hSTFT,
                     float_complex*** dataFD,
                     int framesize,
                     float** dataTD);

/**
 * Re-allocates memory to support a change in the number of input/output
 * channels
 *
 * @param[in] hSTFT      afSTFT handle
 * @param[in] new_nCHin  New number of input channels
 * @param[in] new_nCHout New number of output channels
 */
void afSTFT_channelChange(void * const hSTFT,
                          int new_nCHin,
                          int new_nCHout);

/**
 * Flushes time-domain buffers with zeros.
 */
void afSTFT_clearBuffers(void * const hSTFT);

/**
 * Returns number of frequency bands
 */
int afSTFT_getNBands(void * const hSTFT);

/**
 * Returns current processing delay, in samples
 *
 * @note The base delay is 9*hopsize, which is increased to 12*hopsize when the
 *       hybrid filtering mode is enabled.
 * @warning Currently only correct when low delay mode is disabled!
 */
int afSTFT_getProcDelay(void * const hSTFT);

/**
 * Returns current frequency vector
 */
void afSTFT_getCentreFreqs(void * const hSTFT,
                           float fs,
                           int nBands,
                           float* freqVector);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* __afSTFTlib_INCLUDED__ */
