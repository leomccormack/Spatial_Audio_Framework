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
    memset(vec, 0, N*sizeof(float));
#endif
}

/* VECTOR MUL ADD */
void vtVma(float* vec1, float* vec2, float* vec3, int N)
{
#ifdef VDSP
    vDSP_vma(vec1,1, vec2,1, vec3,1, vec3,1,N);
#else
    int k;
    for (k=0;k<N;k++) {
        vec3[k] += vec1[k]*vec2[k];
    }
#endif
}

/* FFT INITIALISATION */
void vtInitFFT(void** planPr, float* timeData, float* frequencyData, int log2n)
{
    *planPr = (void*)malloc(sizeof(vtFFT));
    vtFFT *h = (vtFFT*)(*planPr);
    h->timeData = timeData;
    h->frequencyData = frequencyData;
    h->N = powf(2,log2n);
    h->log2n = log2n;
#if defined(VDSP)
    h->FFT = (void*)vDSP_create_fftsetup( h->log2n, FFT_RADIX2);
    h->VDSP_split.realp = frequencyData;
    h->VDSP_split.imagp = &(frequencyData[(h->N)/2]);
#elif defined(MKL_FFT)
    float  Scale;
    const int number_of_channels = 1; /* hard coded here for 1 channel */
    MKL_LONG input_strides[2], output_strides[2], Status;
    h->mkl_fft_out = calloc(((h->N)/2+1)*2, sizeof(float)); 
    /* create handle */
    Status = DftiCreateDescriptor(&(h->MKL_FFT_Handle), DFTI_SINGLE, DFTI_REAL, 1, h->N); /* 1-D, single precision, real_input->fft->half_complex->ifft->real_output */
    /* Configure handle */
    Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE); /* Not inplace, i.e. output has its own dedicated memory */
    /* specify output format as complex conjugate-symmetric data. This is the same as MatLab, except only the
     * first N/2+1 elements are returned. The inverse transform will automatically symmetrically+conjugate
     * replicate these elements, in order to get the required N elements internally. */
    Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if(number_of_channels > 1) /* only required for multiple channels */
        Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_NUMBER_OF_TRANSFORMS, number_of_channels);
    Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_INPUT_DISTANCE, 1);  /* strides between samples (default=1) */
    Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_OUTPUT_DISTANCE, 1); /* strides between samples (default=1) */
    input_strides[0]  = 0; input_strides[1]  = 1;
    output_strides[0] = 0; output_strides[1] = 1;
    Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_INPUT_STRIDES, input_strides);   /* strides between channels (default=[0,1]) */
    Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_OUTPUT_STRIDES, output_strides); /* strides between channels (default=[0,1]) */
    /* Configuration parameters for backward-FFT */
    /* Since vDSP uses an energy preserving scaling for their FFT, I've chosen here to scale the output by 2 to keep them
     * consistent, but it doesn't really matter. The 1/N scaling is done internally.
     *
     * The Logic for this: "for real signals, half the energy of any spectrum are in thier complex conjugate bins of the FFT result.
     *     You can add the two equal magnitude bins, or more commonly, just scale one of them by 2X and ignore the other half of
     *     the FFT result array"*/
    Scale = 2.0f;
    Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_BACKWARD_SCALE, Scale); /* scalar applied after ifft */
    /* commit these chosen parameters */
    Status = DftiCommitDescriptor(h->MKL_FFT_Handle);
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
#if defined(VDSP)
    vDSP_destroy_fftsetup(h->FFT);
#elif defined(MKL_FFT)
    MKL_LONG Status;
    Status = DftiFreeDescriptor(&(h->MKL_FFT_Handle));
    free(h->mkl_fft_out); 
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
#if defined(VDSP)
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
#elif defined(MKL_FFT)
    int i;
    MKL_LONG Status;
    if (positiveForForwardTransform > 0) /* FORWARD FFT */
    {
        Status = DftiComputeForward(h->MKL_FFT_Handle, h->timeData, h->mkl_fft_out);
        for(i=0; i<(h->N)/2; i++){
            h->frequencyData[i] = h->mkl_fft_out[2*i];
            h->frequencyData[i+(h->N)/2] = h->mkl_fft_out[2*i+1];
        }
    }
    else /* INVERSE FFT */
    {
        for(i=0; i<(h->N)/2; i++){
            h->mkl_fft_out[2*i] = h->frequencyData[i];
            h->mkl_fft_out[2*i+1] = h->frequencyData[i+(h->N)/2];
        }
        /* not 100% why, but if I don't zero this value here, then small numerical
         * distortions creep into the output signals of afSTFT.
         * This also makes it consistent with the other 2 FFT options. */
        h->mkl_fft_out[(h->N)] = 0.0f;
        h->mkl_fft_out[(h->N)+1] = 0.0f;
        Status = DftiComputeBackward(h->MKL_FFT_Handle, h->mkl_fft_out, h->timeData);
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
