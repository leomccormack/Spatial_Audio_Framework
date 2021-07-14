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
 * @file afSTFT_internal.h
 * @brief A modified version of afSTFTlib
 *
 * The original afSTFT code (by Juha Vilkamo) can be found here:
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
 * The afSTFT design is also described in more detail in [1]
 *
 * @see [1] Vilkamo, J., & Ba"ckstro"m, T. (2018). Time--Frequency Processing:
 *          Methods and Tools. In Parametric Time--Frequency Domain Spatial
 *          Audio. John Wiley & Sons.
 *
 * @author Juha Vilkamo
 * @date 08.04.2015
 * @license MIT
 */

#ifndef __afSTFTlib_INTERNAL_INCLUDED__
#define __afSTFTlib_INTERNAL_INCLUDED__

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
#include "saf_externals.h" 
#else
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
#endif
 
/* Coefficients for a half-band filter, i.e., the "hybrid filter" applied optionally at the bands 1--4. */
#define COEFF1 0.031273141818515176604f /**< Filter coefficient 0 for hybrid filtering */
#define COEFF2 0.28127313041521179171f  /**< Filter coefficient 1 for hybrid filtering */
#define COEFF3 0.5f /**< Filter coefficient 3 for hybrid filtering */


/* ========================================================================== */
/*                            Internal structures                             */
/* ========================================================================== */

/**
 * Complex data type used by afSTFTlib
 */
typedef struct {
    float *re;
    float *im;
} complexVector;

/**
 * Main data structure for afSTFTlib
 */
typedef struct{
    int inChannels;
    int outChannels;
    int hopSize;
    int hLen;
    int LDmode;
    int hopIndexIn;
    int hopIndexOut;
    int totalHops;
    float *protoFilter;
    float *protoFilterI;
    float **inBuffer;
    float *fftProcessFrameTD;
    float **outBuffer;
#ifdef AFSTFT_USE_SAF_UTILITIES
    void* hSafFFT;
    float_complex *fftProcessFrameFD;
    float* tempHopBuffer;
#else
    int pr;
    int log2n;
    void *vtFFT;
    float *fftProcessFrameFD;
#endif
    void *h_afHybrid;
    int hybridMode;
    
} afSTFTlib_internal_data;

/**
 * Data structure for the hybrid filtering employed by afSTFTlib.
 *
 * The purpose of this filtering is to further divide the 4 lowest FFT-bins,
 * to improve the frequency resolution at low-frequencies. For example, 129
 * bins would become 133 hybrid-bins.
 */
typedef struct{
    int inChannels;
    int outChannels;
    int hopSize;
    float hybridCoeffs[3];
    complexVector **analysisBuffer;
    int loopPointer;
} afHybrid;


/* ========================================================================== */
/*                            Internal functions                              */
/* ========================================================================== */

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
void afSTFTlib_init(void** handle,
                    int hopSize,
                    int inChannels,
                    int outChannels,
                    int LDmode,
                    int hybridMode);

/**
 * Re-allocates memory to support a change in the number of input/output
 * channels */
void afSTFTlib_channelChange(void* handle,
                             int new_inChannels,
                             int new_outChannels);

/** Flushes time-domain buffers with zeros */
void afSTFTlib_clearBuffers(void* handle);

/** Applies the forward afSTFT transform */
void afSTFTlib_forward(void* handle,
                       float** inTD,
                       complexVector* outFD);

/** Applies the backward afSTFT transform. */
void afSTFTlib_inverse(void* handle,
                       complexVector* inFD,
                       float** outTD);

/** Destroys an instance of afSTFTlib */
void afSTFTlib_free(void* handle);

/** Creates and initialises an instance of the afHybrid filtering structure */
void afHybridInit(void** handle,
                  int hopSize,
                  int inChannels,
                  int outChannels);

/** Forward hybrid-filtering transform */
void afHybridForward(void* handle,
                     complexVector* FD);

/** Inverse hybrid-filtering transform */
void afHybridInverse(void* handle,
                     complexVector* FD);

/** Frees an instnce of the afHybrid filtering structure */
void afHybridFree(void* handle);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* __afSTFTlib_INTERNAL_INCLUDED__ */ 
