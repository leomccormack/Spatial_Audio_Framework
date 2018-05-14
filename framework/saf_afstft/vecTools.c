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

#include "vecTools.h"


/* VECTOR FLUSH */
void vtClr(float* vec, int N)
{
#ifdef VDSP
    vDSP_vclr(vec,1, N);
#else
    float *p1 = vec;
    for (int k=0;k<N;k++)
    {
        *p1=0.0f;
        p1++;
    }
#endif
}

/* VECTOR MUL ADD */
void vtVma(float* vec1, float* vec2, float* vec3, int N)
{
#ifdef VDSP
    vDSP_vma(vec1,1, vec2,1, vec3,1, vec3,1,N);
#else
    for (int k=0;k<N;k++)
    {
        vec3[k] += vec1[k]*vec2[k];
    }
#endif
}


/* FFT INITIALIZATION */
void vtInitFFT(void** planPr, float* timeData, float* frequencyData, int log2n)
{
    *planPr = (void*)malloc(sizeof(vtFFT));
    vtFFT *h = (vtFFT*)(*planPr);
    h->timeData = timeData;
    h->frequencyData = frequencyData;
    h->N = (int)pow(2,log2n);
    h->log2n = log2n;
#ifdef VDSP
    h->FFT = (void*)vDSP_create_fftsetup( h->log2n, FFT_RADIX2);
    h->VDSP_split.realp = frequencyData;
    h->VDSP_split.imagp = &(frequencyData[(h->N)/2]);
#else
    /* Ooura */
    h->w = (float*)malloc(sizeof(float)*(h->N)/2);
    h->ip = (int*)malloc(sizeof(int)*(2+h->N));
    h->a = (float*)malloc(sizeof(float)*(h->N));
    h->ip[0]=0;
    rdft(h->N,1,h->a, h->ip, h->w);
#endif
}

/* FFT FREE */
void vtFreeFFT(void* planPr)
{
    vtFFT *h = (vtFFT*)(planPr);
#ifdef VDSP
    vDSP_destroy_fftsetup(h->FFT);
#else
    free(h->w);
    free(h->ip);
    free(h->a);
#endif
    free(planPr);
}



/* FFT RUN */
void vtRunFFT(void* planPr, int positiveForForwardTransform)
{
    vtFFT *h = (vtFFT*)planPr;
#ifdef VDSP
    if (positiveForForwardTransform > 0) /* FORWARD FFT */
    {
        vDSP_ctoz((DSPComplex*)(h->timeData), 2, &(h->VDSP_split), 1, (h->N)/2);
        vDSP_fft_zrip((FFTSetup)(h->FFT),&(h->VDSP_split),1, h->log2n, FFT_FORWARD);
    }
    else /* INVERSE FFT */
    {
        vDSP_fft_zrip(h->FFT,&(h->VDSP_split),1, h->log2n, FFT_INVERSE);
        vDSP_ztoc(&(h->VDSP_split),1, (DSPComplex*)h->timeData, 2, (h->N)/2);
    }
#else
    /* Note (A): The phase is conjugated below for Ooura's FFT to produce the same output than that of the vDSP FFT. */
    int k;
    if (positiveForForwardTransform > 0) //
    {
        memcpy(h->a,h->timeData,sizeof(float)*h->N);
        rdft(h->N, 1, h->a, h->ip, h->w);
        for (k=0;k<(h->N)/2;k++)
        {
            h->frequencyData[k] = h->a[2*k];
            h->frequencyData[k+(h->N)/2] = - h->a[2*k+1]; /* Check note (A) above */
        }
        h->frequencyData[(h->N)/2] *= -1.0f; /* Check note (A) above */
    }
    else //
    {
        for (k=0;k<(h->N)/2;k++)
        {
            h->a[2*k] = 4.0f*h->frequencyData[k];
            h->a[2*k+1] = -4.0f*h->frequencyData[k+(h->N)/2]; /* Check note (A) above */
        }
        h->a[1] *= -1.0f; /* Check note (A) above */
        rdft(h->N, -1, h->a, h->ip, h->w);
        memcpy(h->timeData,h->a,sizeof(float)*h->N);
    }
#endif
}
