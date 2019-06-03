/*
 Copyright 2019 Leo McCormack
 
 Permission to use, copy, modify, and/or distribute this software for any purpose with or
 without fee is hereby granted, provided that the above copyright notice and this permission
 notice appear in all copies.
 
 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 OR PERFORMANCE OF THIS SOFTWARE.
*/
/*
 * Filename:
 *     saf_fft.c
 * Description:
 *     Wrapper for optimised fast Fourier transforms (FFT).
 * Dependencies:
 *     Windows users only: custom Intel MKL '.lib/.dll' files are required.
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library by default.
 *     However, Mac users may elect to use a custom Intel MKL '.dylib' instead.
 *     Further instructions for both Windows/Mac users can be found here:
 *     https://github.com/leomccormack/Spatial_Audio_Framework
 * Author, date created:
 *     Leo McCormack, 06.04.2019
 */

#include "saf_utilities.h"
#include "saf_fft.h"

/* NOTE: vDSP_fft hasn't been extensively tested, and doesn't seem to return the Nyquist value?! */

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
    
}safFFT_data;


void safFFT_create
(
    void ** const phFFT,
    int N
)
{
    *phFFT = malloc(sizeof(safFFT_data));
    safFFT_data *h = (safFFT_data*)(*phFFT);
    
    h->N = N;
    h->Scale = 1.0f/(float)N; /* output scaling after ifft */
#if defined(__ACCELERATE__)
    h->log2n = (int)(log2f((float)N)+0.1f);
    h->FFT = (void*)vDSP_create_fftsetup(h->log2n, FFT_RADIX2);
    h->VDSP_split.realp = malloc((h->N/2+1)*sizeof(float));
    h->VDSP_split.imagp = malloc((h->N/2+1)*sizeof(float));
#elif defined(INTEL_MKL_VERSION)
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
#endif
}

void safFFT_destroy
(
    void ** const phFFT
)
{
    safFFT_data *h = (safFFT_data*)(*phFFT);
#if defined(__ACCELERATE__)
    vDSP_destroy_fftsetup(h->FFT);
    free(h->VDSP_split.realp);
    free(h->VDSP_split.imagp);
#elif defined(INTEL_MKL_VERSION)
    h->Status = DftiFreeDescriptor(&(h->MKL_FFT_Handle));
#endif
    free(h);
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
    vDSP_ctoz((DSPComplex*)inputTD, 2, &(h->VDSP_split), 1, (h->N)/2);
    vDSP_fft_zrip((FFTSetup)(h->FFT),&(h->VDSP_split),1, h->log2n, FFT_FORWARD);
    /* the output is scaled by 2, because vDSP_fft automatically compensates for the loss of energy
     * when removing the symmetric/conjugate (N/2+2:N) bins. However, this is dumb, so the 2x scaling
     * is removed here. */
    for(i=0; i<h->N/2+1; i++)
        outputFD[i] = cmplxf(h->VDSP_split.realp[i]/2.0f, h->VDSP_split.imagp[i]/2.0f);
#elif defined(INTEL_MKL_VERSION)
    h->Status = DftiComputeForward(h->MKL_FFT_Handle, inputTD, outputFD);
#endif
}

void safFFT_backward
(
    void * const hFFT,
    float_complex* inputFD,
    float* outputTD
)
{
    safFFT_data *h = (safFFT_data*)(hFFT);
#if defined(__ACCELERATE__)
    int i;
    for(i=0; i<h->N/2+1; i++){
        h->VDSP_split.realp[i] = crealf(inputFD[i]);
        h->VDSP_split.imagp[i] = cimagf(inputFD[i]);
    }
    vDSP_fft_zrip(h->FFT,&(h->VDSP_split),1, h->log2n, FFT_INVERSE);
    vDSP_ztoc(&(h->VDSP_split),1, (DSPComplex*)outputTD, 2, (h->N)/2);
    utility_svsmul(outputTD, &(h->Scale), h->N, NULL);
#elif defined(INTEL_MKL_VERSION)
    h->Status = DftiComputeBackward(h->MKL_FFT_Handle, inputFD, outputTD);
#endif
}






 
