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

#ifndef __afSTFTlib_tester__afSTFTlib__
#define __afSTFTlib_tester__afSTFTlib__

/* The original afSTFT code, written by Juha Vilkamo, can be found here:
 * https://github.com/jvilkamo/afSTFT
 * This version is slightly modified. It adds a function to change the number of
 * channels on the fly and includes vectors for the hybrid mode centre
 * frequencies @44.1kHz/48kHz with 128 hop size.
 * It also supports Intel's MKL FFT if "SAF_USE_INTEL_MKL" is defined. */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vecTools.h"

extern const double __afCenterFreq48e3[133];
extern const double __afCenterFreq44100[133];

typedef struct {
    float *re;
    float *im;
} complexVector;

typedef struct{
    int inChannels;
    int outChannels;
    int maxChannels;
    int hopSize;
    int hLen;
    int pr;
    int LDmode;
    int hopIndexIn;
    int hopIndexOut;
    int totalHops;
    float *protoFilter;
    float *protoFilterI;
    float **inBuffer;
    float *fftProcessFrameTD;
    float *fftProcessFrameFD;
    float **outBuffer;
    int log2n;
    void *vtFFT;
    void *h_afHybrid;
    int hybridMode;
} afSTFT;

typedef struct{
    int inChannels;
    int outChannels;
    int hopSize;
    float hybridCoeffs[3];
    complexVector **analysisBuffer;
    int loopPointer;

} afHybrid;


/* Call these */

void afSTFTinit(void** handle, int hopSize, int inChannels, int outChannels, int LDmode, int hybridMode);

/* Re-allocates memory to support a channel change. Do not call in the real-time loop. */
void afSTFTchannelChange(void* handle, int new_inChannels, int new_outChannels);

/* Flushes time-domain buffers with zeros. Do not call in the real-time loop. */
void afSTFTclearBuffers(void* handle);

void afSTFTforward(void* handle, float** inTD, complexVector* outFD);

void afSTFTinverse(void* handle, complexVector* inFD, float** outTD);

void afSTFTfree(void* handle);


/* Internal functions */

void afHybridInit(void** handle, int hopSize, int inChannels, int outChannels);

void afHybridForward(void* handle, complexVector* FD);

void afHybridInverse(void* handle, complexVector* FD);

void afHybridFree(void* handle);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* defined(__afSTFTlib_tester__afSTFTlib__) */
