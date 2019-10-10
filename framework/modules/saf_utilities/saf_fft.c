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
 * Filename: saf_fft.c
 * -------------------
 * Wrapper for optimised fast Fourier transform (FFT) routines. If none are
 * linked, then it employs the highly respectable KissFFT from here
 * (BSD 3-Clause License): https://github.com/mborgerding/kissfft
 * If linking Apple Accelerate: KissFFT is also used in cases where the FFT size
 * is not a power of 2.
 *
 * Dependencies:
 *     Intel MKL, Apple Accelerate, or KissFFT (included in framework)
 * Author, date created:
 *     Leo McCormack, 06.04.2019
 */

#include "saf_utilities.h"
#include "saf_fft.h"

typedef struct _safFFT_data {
    int N;
    float  Scale;
#if defined(__ACCELERATE__)
    int log2n;
    FFTSetup FFT;
    DSPSplitComplex VDSP_split;
#elif defined(INTEL_MKL_VERSION)
    DFTI_DESCRIPTOR_HANDLE MKL_FFT_Handle;
    MKL_LONG input_strides[2], output_strides[2], Status;
#endif
    int useKissFFT_flag;
    kiss_fftr_cfg kissFFThandle_fwd;
    kiss_fftr_cfg kissFFThandle_bkw;
    
}safFFT_data;

void safFFT_create
(
    void ** const phFFT,
    int N
)
{
    *phFFT = malloc1d(sizeof(safFFT_data));
    safFFT_data *h = (safFFT_data*)(*phFFT);
    
    h->N = N;
    h->Scale = 1.0f/(float)N; /* output scaling after ifft */
    assert(N>=2); /* only even (non zero) FFT sizes allowed */
#if defined(__ACCELERATE__)
    if(ceilf(log2f(N)) == floorf(log2f(N))) /* true if N is 2 to the power of some integer number */
        h->useKissFFT_flag = 0;
    else
        h->useKissFFT_flag = 1;
    /* Apple Accelerate only supports 2^x FFT sizes */
    if(!h->useKissFFT_flag){
        h->log2n = (int)(log2f((float)N)+0.1f);
        h->FFT = (void*)vDSP_create_fftsetup(h->log2n, FFT_RADIX2);
        h->VDSP_split.realp = malloc1d((h->N/2)*sizeof(float));
        h->VDSP_split.imagp = malloc1d((h->N/2)*sizeof(float));
    }
#elif defined(INTEL_MKL_VERSION)
    h->useKissFFT_flag = 0;
    h->MKL_FFT_Handle = 0;
    h->Status = DftiCreateDescriptor(&(h->MKL_FFT_Handle), DFTI_SINGLE, DFTI_REAL, 1, h->N); /* 1-D, single precision, real_input->fft->half_complex->ifft->real_output */
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE); /* Not inplace, i.e. output has its own dedicated memory */
    /* specify output format as complex conjugate-symmetric data. This is the same as MatLab, except only the
     * first N/2+1 elements are returned. The inverse transform will automatically symmetrically+conjugate
     * replicate these elements, in order to get the required N elements internally. */
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    const int number_of_channels = 1; /* hard coded here for 1 channel */
    if(number_of_channels > 1) /* only required for multiple channels */
        h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_NUMBER_OF_TRANSFORMS, number_of_channels);
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_INPUT_DISTANCE, 1);  /* strides between samples (default=1) */
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_OUTPUT_DISTANCE, 1); /* strides between samples (default=1) */
    h->input_strides[0]  = 0; h->input_strides[1]  = 1; /* hard coded here for 1 channel */
    h->output_strides[0] = 0; h->output_strides[1] = 1; /* hard coded here for 1 channel */
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_INPUT_STRIDES, h->input_strides);   /* strides between channels (default=[0,1]) */
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_OUTPUT_STRIDES, h->output_strides); /* strides between channels (default=[0,1]) */
    /* Configuration parameters for backward-FFT */
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_BACKWARD_SCALE, h->Scale);      /* scalar applied after ifft */
    /* commit these chosen parameters */
    h->Status = DftiCommitDescriptor(h->MKL_FFT_Handle);
#else
    h->useKissFFT_flag = 1;
#endif
    if(h->useKissFFT_flag){
       h->kissFFThandle_fwd = kiss_fftr_alloc(h->N, 0, NULL, NULL);
       h->kissFFThandle_bkw = kiss_fftr_alloc(h->N, 1, NULL, NULL);
    }
}

void safFFT_destroy
(
    void ** const phFFT
)
{
    safFFT_data *h = (safFFT_data*)(*phFFT);
    if(h!=NULL){
#if defined(__ACCELERATE__)
        if(!h->useKissFFT_flag){
            vDSP_destroy_fftsetup(h->FFT);
            free(h->VDSP_split.realp);
            free(h->VDSP_split.imagp);
        }
#elif defined(INTEL_MKL_VERSION)
        h->Status = DftiFreeDescriptor(&(h->MKL_FFT_Handle));
#endif
        if(h->useKissFFT_flag){
            kiss_fftr_free(h->kissFFThandle_fwd);
            kiss_fftr_free(h->kissFFThandle_bkw);
        }
        free(h);
        h=NULL;
    }
}

void safFFT_forward
(
    void * const hFFT,
    float* inputTD,
    float_complex* outputFD
)
{
    safFFT_data *h = (safFFT_data*)(hFFT);
#if defined(__ACCELERATE__)
    int i;
    if(!h->useKissFFT_flag){
        vDSP_ctoz((DSPComplex*)inputTD, 2, &(h->VDSP_split), 1, (h->N)/2);
        vDSP_fft_zrip((FFTSetup)(h->FFT),&(h->VDSP_split), 1, h->log2n, FFT_FORWARD);
        /* DC */
        outputFD[0] = cmplxf(h->VDSP_split.realp[0]/2.0f, 0.0f);
        /* Note: the output is scaled by 2, because vDSP_fft automatically compensates for the loss of energy
         * when removing the symmetric/conjugate (N/2+2:N) bins. However, this is dumb, so the 2x scaling
         * is removed here, so it has parity with the other FFT implementations supported by SAF. */
        for(i=1; i<h->N/2; i++)
            outputFD[i] = cmplxf(h->VDSP_split.realp[i]/2.0f, h->VDSP_split.imagp[i]/2.0f);
        /* the real part of the Nyquist value is the imaginary part of DC. (At least for this fucked up library it is). */
        outputFD[h->N/2] = cmplxf(h->VDSP_split.imagp[0]/2.0f, 0.0f);
        /* https://stackoverflow.com/questions/43289265/implementing-an-fft-using-vdsp */
    }
#elif defined(INTEL_MKL_VERSION)
    h->Status = DftiComputeForward(h->MKL_FFT_Handle, inputTD, outputFD);
#endif
    if(h->useKissFFT_flag)
        kiss_fftr(h->kissFFThandle_fwd, inputTD, (kiss_fft_cpx*)outputFD);
}

void safFFT_backward
(
    void * const hFFT,
    float_complex* inputFD,
    float* outputTD
)
{
    safFFT_data *h = (safFFT_data*)(hFFT);
    int i;
#if defined(__ACCELERATE__)
    if(!h->useKissFFT_flag){
        h->VDSP_split.realp[0] = crealf(inputFD[0]);
        h->VDSP_split.imagp[0] = crealf(inputFD[h->N/2]);
        for(i=1; i<h->N/2; i++){
            h->VDSP_split.realp[i] = crealf(inputFD[i]);
            h->VDSP_split.imagp[i] = cimagf(inputFD[i]);
        }
        vDSP_fft_zrip(h->FFT, &(h->VDSP_split), 1, h->log2n, FFT_INVERSE);
        vDSP_ztoc(&(h->VDSP_split), 1, (DSPComplex*)outputTD, 2, (h->N)/2);
        vDSP_vsmul(outputTD, 1, &(h->Scale), outputTD, 1, h->N);
    }
#elif defined(INTEL_MKL_VERSION)
    h->Status = DftiComputeBackward(h->MKL_FFT_Handle, inputFD, outputTD);
#endif
    if(h->useKissFFT_flag){
        kiss_fftri(h->kissFFThandle_bkw, (kiss_fft_cpx*)inputFD, outputTD);
        for(i=0; i<h->N; i++)
            outputTD[i] /= (float)(h->N);
    }
}






 
