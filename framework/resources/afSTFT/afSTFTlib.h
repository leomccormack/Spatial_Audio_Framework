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

#if defined(SAF_ENABLE_AFSTFT) && !defined(__afSTFTlib_tester__afSTFTlib__)
#define __afSTFTlib_tester__afSTFTlib__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
   
/* The original afSTFT code, written by Juha Vilkamo, can be found here:
 * https://github.com/jvilkamo/afSTFT
 * This version is slightly modified. It adds a function to change the number of
 * channels on the fly and includes vectors for the hybrid mode centre
 * frequencies @44.1kHz/48kHz with 128 hop size for convenience.
 * It also supports the use of SAF utilities (for the vectorisation and FFT).
 *
 * Remove the "AFSTFT_USE_SAF_UTILITIES" definition, and add the vecTools.h/c
 * and fft4g.h/c to your project, if you want to use the original afSTFT code.
 * Note the vecTools.h/c and fft4g.h/c files, may be found here:
 *     https://github.com/jvilkamo/afSTFT */
#define AFSTFT_USE_SAF_UTILITIES
    
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
    
extern const double __afCenterFreq48e3[133];
extern const double __afCenterFreq44100[133];

typedef struct {
    float *re;
    float *im;
} complexVector;

/* Main Functions: */

void afSTFTinit(void** handle, int hopSize, int inChannels, int outChannels, int LDmode, int hybridMode);

/* Re-allocates memory to support a channel change. Do not call in the real-time loop. */
void afSTFTchannelChange(void* handle, int new_inChannels, int new_outChannels);

/* Flushes time-domain buffers with zeros. Do not call in the real-time loop. */
void afSTFTclearBuffers(void* handle);

void afSTFTforward(void* handle, float** inTD, complexVector* outFD);

void afSTFTinverse(void* handle, complexVector* inFD, float** outTD);

void afSTFTfree(void* handle);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* defined(SAF_ENABLE_AFSTFT) && !defined(__afSTFTlib_tester__afSTFTlib__) */ 
