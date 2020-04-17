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
 */

#ifndef __afSTFTlib_tester__afSTFTlib__
#define __afSTFTlib_tester__afSTFTlib__

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
    
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
    
/**
 * afSTFT centre frequencies for 128 hop size and hybrid-mode enabled (48kHz) */
extern const double __afCenterFreq48e3[133];
/**
 * afSTFT centre frequencies for 128 hop size and hybrid-mode enabled (44.1kHz)
 */
extern const double __afCenterFreq44100[133];
/**
 * Complex data type used by afSTFTlib
 */
typedef struct {
    float *re;
    float *im;
} complexVector;

/**
 * Initialises an instance of afSTFTlib [1]
 *
 * @param[in] handle      (&) afSTFTlib handle
 * @param[in] hopSize     Hop size, in samples
 * @param[in] inChannels  Number of input channels
 * @param[in] outChannels Number of output channels
 * @param[in] LDmode      '0' disable low-delay mode, '1' enable
 * @param[in] hybridMode  '0' disable hybrid-mode, '1' enable
 *
 * @see [1] Vilkamo, J., & Backstrom, T. (2018). Time--Frequency Processing:
 *          Methods and Tools. In Parametric Time-Frequency Domain Spatial
 *          Audio. John Wiley & Sons.
 */
void afSTFTinit(void** handle,
                int hopSize,
                int inChannels,
                int outChannels,
                int LDmode,
                int hybridMode);

/**
 * Re-allocates memory to support a change in the number of input/output
 * channels
 *
 * @note Not thread safe. So do not call in the middle of a real-time loop.
 *
 * @param[in] handle          afSTFTlib handle
 * @param[in] new_inChannels  New number of input channels
 * @param[in] new_outChannels New number of output channels
 */
void afSTFTchannelChange(void* handle, int new_inChannels, int new_outChannels);

/**
 * Flushes time-domain buffers with zeros.
 *
 * @param[in] handle afSTFTlib handle
 */
void afSTFTclearBuffers(void* handle);

/**
 * Applies the forward afSTFT transform.
 *
 * @param[in] handle afSTFTlib handle
 * @param[in] inTD   input time-domain signals; inChannels x hopSize
 * @param[in] outFD  input time-frequency domain signals; inChannels x nBands
 */
void afSTFTforward(void* handle, float** inTD, complexVector* outFD);

/**
 * Applies the backward afSTFT transform.
 *
 * @param[in] handle afSTFTlib handle
 * @param[in] inFD   output time-domain signals; outChannels x hopSize
 * @param[in] outTD  output time-frequency domain signals; outChannels x nBands
 */
void afSTFTinverse(void* handle, complexVector* inFD, float** outTD);

/**
 * Destroys an instance of afSTFTlib
 *
 * @param[in] handle (&) afSTFTlib handle
 */
void afSTFTfree(void* handle);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* __afSTFTlib_tester__afSTFTlib__ */ 
