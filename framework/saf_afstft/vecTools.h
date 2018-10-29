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

#ifndef __afSTFT_MEXfile__vecTools__
#define __afSTFT_MEXfile__vecTools__

#include "saf_utilities.h"
#ifdef __APPLE__
  #ifdef NDEBUG
    #define VDSP 1
  #endif
#elif defined(_MSC_VER) && defined(INTEL_MKL_VERSION)
  /* Real to Complex FFT is not currently supported by intel! WTF. */
  //#define MKL_FFT 1
#endif

#ifdef VDSP
#include <Accelerate/Accelerate.h>
#elif MKL_FFT
#include "mkl_dfti.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fft4g.h"


typedef struct {
    float *timeData;
    float *frequencyData;
    int N;
    int log2n;
#ifdef VDSP
    FFTSetup FFT; 
    DSPSplitComplex VDSP_split;
#elif MKL_FFT
	float_complex *timeData_mkl;
	float_complex *frequencyData_mkl;
	DFTI_DESCRIPTOR_HANDLE my_desc1_handle;
	DFTI_DESCRIPTOR_HANDLE my_desc2_handle;
	MKL_LONG status;
#else
    float *a,*w;
    int *ip;
#endif
} vtFFT;

void vtClr(float* vec, int N);

void vtVma(float* vec1, float* vec2, float* vec3, int N);


#endif /* defined(__afSTFT_MEXfile__vecTools__) */
