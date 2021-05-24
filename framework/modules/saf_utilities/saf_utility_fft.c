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

/**
 * @file saf_utility_fft.c
 * @ingroup Utilities
 * @brief Wrappers for optimised fast Fourier transform (FFT) routines
 *
 * @note If none of the supported optimised FFT implementations are linked, then
 *       saf_fft employs the highly respectable KissFFT from here (BSD 3-Clause
 *       License): https://github.com/mborgerding/kissfft
 * @note If linking Apple Accelerate: KissFFT is also used in cases where the
 *       FFT size is not 2^x.
 * @note The developers are aware of the existance of FFTW. However, due to how
 *       cumbersome it is to install for multiple platforms, and since it is not
 *       as fast as Intel IPP/MKL in our experience, we haven't yet included it
 *       in this wrapper. However, if you wish to add support for it, then we
 *       are happy to accept pull requests.
 *
 * ## Dependencies
 *   Intel MKL, Apple Accelerate, or KissFFT (included in the framework)
 *
 * @author Leo McCormack
 * @date 06.04.2019
 */

#include "saf_utilities.h"
#include "saf_externals.h"

/**
 * Data structure for short-time Fourier transform
 */
typedef struct _saf_stft_data {
    int winsize, hopsize, fftsize, nCHin, nCHout, nBands;
    void* hFFT;
    int numOvrlpAddBlocks, bufferlength, nPrevHops;
    float* window, *insig_win, *outsig_win;
    float** overlapAddBuffer;
    float*** prev_inhops;
    float_complex* tmp_fft;
    SAF_STFT_FDDATA_FORMAT FDformat;

}saf_stft_data;

/**
 * Data structure for real-(half)complex FFT transforms
 */
typedef struct _saf_rfft_data {
    int N;
    float  Scale;
#if defined(SAF_USE_INTEL_IPP)
    int useIPPfft_FLAG;
    int specSize, specBufferSize, bufferSize, log2n;
    IppsDFTSpec_R_32f* hDFTspec;
    IppsFFTSpec_R_32f* hFFTspec;
    Ipp8u* memSpec;
    Ipp8u* buffer;
    Ipp8u* memInit;
#elif defined(__ACCELERATE__)
    vDSP_DFT_Setup DFT_fwd;
    vDSP_DFT_Setup DFT_bwd;
    DSPSplitComplex VDSP_split;
    DSPSplitComplex VDSP_split_tmp;
#elif defined(INTEL_MKL_VERSION)
    DFTI_DESCRIPTOR_HANDLE MKL_FFT_Handle;
    MKL_LONG input_strides[2], output_strides[2], Status;
#else /* DEFAULT: */
    kiss_fftr_cfg kissFFThandle_fwd;
    kiss_fftr_cfg kissFFThandle_bkw;
#endif
    
}saf_rfft_data;

/**
 * Data structure for complex-complex FFT transforms
 */
typedef struct _saf_fft_data {
    int N;
    float  Scale;
#if defined(SAF_USE_INTEL_IPP)
    int useIPPfft_FLAG;
    int specSize, specBufferSize, bufferSize, log2n;
    IppsDFTSpec_C_32fc* hDFTspec;
    IppsFFTSpec_C_32fc* hFFTspec;
    Ipp8u* memSpec;
    Ipp8u* buffer;
    Ipp8u* memInit;
#elif defined(__ACCELERATE__)
    vDSP_DFT_Setup DFT_fwd;
    vDSP_DFT_Setup DFT_bwd;
    DSPSplitComplex VDSP_split;
    DSPSplitComplex VDSP_split_tmp;
#elif defined(INTEL_MKL_VERSION)
    DFTI_DESCRIPTOR_HANDLE MKL_FFT_Handle;
    MKL_LONG Status;
#else
    kiss_fft_cfg kissFFThandle_fwd;
    kiss_fft_cfg kissFFThandle_bkw;
#endif

}saf_fft_data;


/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

void getUniformFreqVector
(
    int fftSize,
    float fs,
    float* freqVector
)
{
    int k;
    for(k=0; k<fftSize/2+1; k++)
        freqVector[k] = (float)k * fs/(float)fftSize;
}

void fftconv
(
    float* x,
    float* h,
    int x_len,
    int h_len,
    int nCH,
    float* y
)
{
    int i, y_len, fftSize, nBins;
    float* h0, *x0, *y0;
    float_complex* H, *X, *Y;
    void* hfft;
    
    /* prep */
    y_len = x_len + h_len - 1;
    fftSize =  (int)((float)nextpow2(y_len)+0.5f);
    nBins = fftSize/2+1;
    h0 = calloc1d(fftSize, sizeof(float));
    x0 = calloc1d(fftSize, sizeof(float));
    y0 = malloc1d(fftSize * sizeof(float));
    H = malloc1d(nBins*sizeof(float_complex));
    X = malloc1d(nBins*sizeof(float_complex));
    Y = malloc1d(nBins*sizeof(float_complex));
    saf_rfft_create(&hfft, fftSize);
    
    /* apply convolution per channel */
    for(i=0; i<nCH; i++){
        /* zero pad to avoid circular convolution artefacts, prior to fft */
        memcpy(h0, &h[i*h_len], h_len*sizeof(float));
        memcpy(x0, &x[i*x_len], x_len*sizeof(float));
        saf_rfft_forward(hfft, x0, X);
        saf_rfft_forward(hfft, h0, H);
        
        /* multiply the two spectra */
        utility_cvvmul(X, H, nBins, Y);
        
        /* ifft, truncate and store to output */
        saf_rfft_backward(hfft, Y, y0);
        memcpy(&y[i*y_len], y0, y_len*sizeof(float));
    }
    
    /* tidy up */
    saf_rfft_destroy(&hfft);
    free(h0);
    free(x0);
    free(y0);
    free(H);
    free(X);
    free(Y);
}

void fftfilt
(
    float* x,
    float* h,
    int x_len,
    int h_len,
    int nCH,
    float* y
)
{
    int i;
    float* y_tmp;
    
    y_tmp = malloc1d(nCH*(x_len+h_len-1)*sizeof(float));
    fftconv(x, h, x_len, h_len, nCH, y_tmp);
    for(i=0; i<nCH; i++)
        memcpy(&y[i*x_len], &y_tmp[i*(x_len+h_len-1)], x_len*sizeof(float));
    free(y_tmp);
}

void hilbert
(
    float_complex* x,
    int x_len,
    float_complex* y
)
{
    int i;
    float_complex *xfft, *h, *xhfft; 
    void* hfft;
    
    saf_fft_create(&hfft, x_len);
    xfft = malloc1d(x_len*sizeof(float_complex));
    h = malloc1d(x_len*sizeof(float_complex));
    xhfft = malloc1d(x_len*sizeof(float_complex));
    
    /* Forward fft */
    saf_fft_forward(hfft, x, xfft);
    
    /* define vector h */
    memset(h, 0, sizeof(float_complex)*x_len);
    if(x_len % 2 == 0){
        /* even */
        h[0] = cmplxf(1.0f, 0.0f);
        h[x_len/2] = cmplxf(1.0f, 0.0f);
        for(i=1;i<x_len/2;i++)
            h[i] = cmplxf(2.0f, 0.0f);
    }
    else{
        saf_assert(x_len % 2 == 0, "Uneven lengths are not supported by saf_fft"); 
        /* odd */
        h[0] = cmplxf(1.0f, 0.0f);
        for(i=1;i<(x_len+1)/2;i++)
            h[i] = cmplxf(2.0f, 0.0f);
    }
    
    /* apply h, and ifft */
    utility_cvvmul(xfft, h, x_len, xhfft);
    saf_fft_backward(hfft, xhfft, y);
    
    /* tidy up */
    saf_fft_destroy(&hfft); 
    free(xfft);
    free(h);
    free(xhfft);
}

/* ========================================================================== */
/*                     Short-time Fourier Transform (STFT)                    */
/* ========================================================================== */

void saf_stft_create
(
    void ** const phSTFT,
    int winsize,
    int hopsize,
    int nCHin,
    int nCHout,
    SAF_STFT_FDDATA_FORMAT FDformat
)
{
    *phSTFT = malloc1d(sizeof(saf_stft_data));
    saf_stft_data *h = (saf_stft_data*)(*phSTFT);

    h->nCHin = nCHin;
    h->nCHout = nCHout;
    h->winsize = winsize;
    h->hopsize = hopsize;
    h->nBands = winsize+1;
    h->FDformat = FDformat;

    /* set-up FFT */
    h->fftsize = 2*winsize;
    saf_rfft_create(&(h->hFFT), h->fftsize);
    h->insig_win = calloc1d(h->fftsize, sizeof(float));

    /* Intermediate buffers */
    h->tmp_fft = malloc1d(h->nBands * sizeof(float_complex));
    h->outsig_win = malloc1d(h->fftsize*sizeof(float));
    h->nPrevHops = winsize/hopsize-1;
    if (h->nPrevHops>0)
        h->prev_inhops = (float***)calloc3d(h->nPrevHops, nCHin, hopsize, sizeof(float));
    else
        h->prev_inhops = NULL;

    /* Windowing function */
    if(winsize==hopsize)
        h->window = NULL;
    else{
        h->window = malloc1d(winsize*sizeof(float));
        getWindowingFunction(WINDOWING_FUNCTION_HANN, winsize, h->window);
    }

    /* Overlap-add buffer */
    h->numOvrlpAddBlocks = winsize/hopsize;
    h->bufferlength = h->numOvrlpAddBlocks * (h->fftsize);
    h->overlapAddBuffer = (float**)calloc2d(nCHout, h->bufferlength, sizeof(float));
}

void saf_stft_destroy
(
    void ** const phSTFT
)
{
    saf_stft_data *h = (saf_stft_data*)(*phSTFT);
    if(h!=NULL){
        saf_rfft_destroy(&(h->hFFT));
        free(h->window);
        free(h->overlapAddBuffer);
        free(h->insig_win);
        free(h->tmp_fft);
        free(h->prev_inhops);
        free(h);
        h=NULL;
        *phSTFT = NULL;
    }
}

void saf_stft_forward
(
    void * const hSTFT,
    float** dataTD,
    int framesize,
    float_complex*** dataFD

)
{
    saf_stft_data *h = (saf_stft_data*)(hSTFT);
    int ch, j, nHops, t, band, idx, hIdx;

    saf_assert(framesize % h->hopsize == 0, "framesize must be multiple of hopsize");  
    nHops = framesize/h->hopsize;

    /* For linear time-invariant (LTI) operation (i.e. no previous hops are
     * required) */
    if(h->winsize==h->hopsize){
        for (t = 0; t<nHops; t++){
            for(ch=0; ch < h->nCHin; ch++){
                /* Window input signal (Rectangular) */
                memcpy(h->insig_win, &dataTD[ch][t*(h->hopsize)], h->winsize*sizeof(float));

                /* Apply FFT and copy data to output dataFD buffer */
                switch(h->FDformat){
                    case SAF_STFT_TIME_CH_BANDS:
                        saf_rfft_forward(h->hFFT, h->insig_win, dataFD[t][ch]);
                        break;

                    case SAF_STFT_BANDS_CH_TIME:
                        saf_rfft_forward(h->hFFT, h->insig_win, h->tmp_fft);
                        for(band=0; band<h->nBands; band++)
                            dataFD[band][ch][t] = h->tmp_fft[band];
                        break;
                }
            }
        }
    }
    /* For oversampled TF transforms */
    else{
        idx = 0;
        for (t = 0; t<nHops; t++){
            for(ch=0; ch < h->nCHin; ch++){
                hIdx = 0;
                /* Window input signal */
                while (hIdx < h->winsize){
                    memcpy(&(h->insig_win[hIdx]), h->prev_inhops[0][ch], h->hopsize*sizeof(float));
                    for(j=0; j< h->nPrevHops-1; j++)
                        memcpy(h->prev_inhops[j][ch], h->prev_inhops[j+1][ch], h->hopsize*sizeof(float));
                    memcpy(h->prev_inhops[h->nPrevHops-1][ch], &dataTD[ch][idx], h->hopsize*sizeof(float));
                    hIdx += h->hopsize;
                }
                utility_svvmul(h->insig_win, h->window, h->winsize, h->insig_win);

                /* Apply FFT and copy data to output dataFD buffer */
                switch(h->FDformat){
                    case SAF_STFT_TIME_CH_BANDS:
                        saf_rfft_forward(h->hFFT, h->insig_win, dataFD[t][ch]);
                        break;

                    case SAF_STFT_BANDS_CH_TIME:
                        saf_rfft_forward(h->hFFT, h->insig_win, h->tmp_fft);
                        for(band=0; band<h->nBands; band++)
                            dataFD[band][ch][t] = h->tmp_fft[band];
                        break;
                }
            }
            idx += h->hopsize;
        } 
    }
}

void saf_stft_backward
(
    void * const hSTFT,
    float_complex*** dataFD,
    int framesize,
    float** dataTD
)
{
    saf_stft_data *h = (saf_stft_data*)(hSTFT);
    int t, ch, nHops, band;

    saf_assert(framesize % h->hopsize == 0, "framesize must be multiple of hopsize");
    nHops = framesize/h->hopsize;

    for (t = 0; t<nHops; t++){
        for(ch=0; ch < h->nCHout; ch++){
            /* Shift data down */
            memcpy(h->overlapAddBuffer[ch], &h->overlapAddBuffer[ch][h->hopsize], (h->numOvrlpAddBlocks-1)*(h->hopsize)*sizeof(float));

            /* Append with zeros */
            memset(&h->overlapAddBuffer[ch][(h->numOvrlpAddBlocks-1)*(h->hopsize)], 0, h->hopsize*sizeof(float));

            /* Apply inverse FFT */
            switch(h->FDformat){
                case SAF_STFT_TIME_CH_BANDS:
                    saf_rfft_backward(h->hFFT, dataFD[t][ch], h->outsig_win);
                    break;
                case SAF_STFT_BANDS_CH_TIME:
                    for(band=0; band<h->nBands; band++)
                        h->tmp_fft[band] = dataFD[band][ch][t];
                    saf_rfft_backward(h->hFFT, h->tmp_fft, h->outsig_win);
                    break;
            }

            /* Overlap-Add and copy 1:hopsize to output buffer */
            utility_svvadd(h->overlapAddBuffer[ch], h->outsig_win, h->fftsize, h->overlapAddBuffer[ch]);
            memcpy(&dataTD[ch][t*(h->hopsize)], h->overlapAddBuffer[ch], h->hopsize*sizeof(float)); 
        }
    }
}

void saf_stft_flushBuffers
(
    void * const hSTFT
)
{
    saf_stft_data *h = (saf_stft_data*)(hSTFT);
    if(h->nPrevHops > 0)
        memset(FLATTEN3D(h->prev_inhops), 0, h->nPrevHops * (h->nCHin) * (h->hopsize) * sizeof(float));
    memset(FLATTEN2D(h->overlapAddBuffer), 0, h->nCHout * (h->bufferlength) * sizeof(float));
}

void saf_stft_channelChange
(
    void * const hSTFT,
    int new_nCHin,
    int new_nCHout
)
{
    saf_stft_data *h = (saf_stft_data*)(hSTFT);
    int i, ch;

    if(new_nCHin != h->nCHin && h->nPrevHops > 0){
        /* Reallocate memory while retaining previous values (which will be
         * truncated if new_nCHin < nCHin) */
        h->prev_inhops = (float***)realloc3d_r((void***)h->prev_inhops, h->nPrevHops, new_nCHin, h->hopsize,
                                               h->nPrevHops, h->nCHin, h->hopsize, sizeof(float));
        /* Zero new channels if new_nCHin > nCHin */
        for(i=0; i< h->nPrevHops; i++)
            for(ch=h->nCHin; ch<new_nCHin; ch++)
                memset(h->prev_inhops[i][ch], 0, h->hopsize * sizeof(float));

        h->nCHin = new_nCHin;
    }

    if(new_nCHout != h->nCHout){
        /* Reallocate memory while retaining previous values (which will be
         * truncated if new_nCHout < nCHout) */
        h->overlapAddBuffer = (float**)realloc2d_r((void**)h->overlapAddBuffer, new_nCHout, h->bufferlength,
                                                   h->nCHout, h->bufferlength, sizeof(float));
        /* Zero new channels if new_nCHout > nCHout */
        for(ch=h->nCHout; ch<new_nCHout; ch++)
            memset(h->overlapAddBuffer[ch], 0, h->bufferlength * sizeof(float));

        h->nCHout = new_nCHout;
    }
}


/* ========================================================================== */
/*                Real<->Half-Complex (Conjugate-Symmetric) FFT               */
/* ========================================================================== */

void saf_rfft_create
(
    void ** const phFFT,
    int N
)
{
    *phFFT = malloc1d(sizeof(saf_rfft_data));
    saf_rfft_data *h = (saf_rfft_data*)(*phFFT);
    
    h->N = N;
    h->Scale = 1.0f/(float)N; /* output scaling after ifft */
    saf_assert(N>=2 && ISEVEN(N), "Only even (non zero) FFT sizes are supported");
#if defined(SAF_USE_INTEL_IPP) 
    /* Use ippsFFT if N is 2^x, otherwise, use ippsDFT */
    if(ceilf(log2f(N)) == floorf(log2f(N))){
        h->useIPPfft_FLAG = 1;
        h->log2n = (int)(log2f((float)N)+0.1f);
        ippsFFTGetSize_R_32f(h->log2n, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &(h->specSize), &(h->specBufferSize), &(h->bufferSize));
        h->hFFTspec = NULL;
        h->memSpec = (Ipp8u*) ippMalloc(h->specSize);
        h->buffer  = (Ipp8u*) ippMalloc(h->bufferSize);
        h->memInit = (Ipp8u*) ippMalloc(h->specBufferSize);
        ippsFFTInit_R_32f(&(h->hFFTspec), h->log2n, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, h->memSpec, h->memInit);
    }
    else{
        h->useIPPfft_FLAG = 0;
        ippsDFTGetSize_R_32f(N, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &(h->specSize), &(h->specBufferSize), &(h->bufferSize));
        h->hDFTspec = (IppsDFTSpec_R_32f*) ippMalloc(h->specSize);
        h->buffer  = (Ipp8u*) ippMalloc(h->bufferSize);
        h->memInit = (Ipp8u*) ippMalloc(h->specBufferSize);
        ippsDFTInit_R_32f(N, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, h->hDFTspec, h->memInit);
    }
    if (h->memInit)
        ippFree(h->memInit);
#elif defined(__ACCELERATE__)
    h->DFT_fwd = vDSP_DFT_zrop_CreateSetup(0, N, vDSP_DFT_FORWARD);
    h->DFT_bwd = vDSP_DFT_zrop_CreateSetup(0, N, vDSP_DFT_INVERSE);
    /* Note that DFT lengths must satisfy: f * 2.^g, where f is 1, 3, 5, or 15, and g >=4 */
    saf_assert(h->DFT_fwd!=0 && h->DFT_bwd!=0, "Failed to create vDSP DFT");
    h->VDSP_split_tmp.realp = malloc1d((h->N/2)*sizeof(float));
    h->VDSP_split_tmp.imagp = malloc1d((h->N/2)*sizeof(float));
    h->VDSP_split.realp = malloc1d((h->N/2)*sizeof(float));
    h->VDSP_split.imagp = malloc1d((h->N/2)*sizeof(float));
#elif defined(INTEL_MKL_VERSION)
    h->MKL_FFT_Handle = 0;
    h->Status = DftiCreateDescriptor(&(h->MKL_FFT_Handle), DFTI_SINGLE, DFTI_REAL, 1, h->N); /* 1-D, single precision, real_input->fft->half_complex->ifft->real_output */
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE); /* Not inplace, i.e. output has its own dedicated memory */
    /* specify output format as complex conjugate-symmetric data. This is the same as MatLab, except only the
     * first N/2+1 elements are returned. The inverse transform will automatically symmetrically+conjugate
     * replicate these elements, in order to get the required N elements internally. */
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    /* Configuration parameters for backward-FFT */
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_BACKWARD_SCALE, h->Scale);      /* scalar applied after ifft */
    /* commit these chosen parameters */
    h->Status = DftiCommitDescriptor(h->MKL_FFT_Handle);
#else
    h->kissFFThandle_fwd = kiss_fftr_alloc(h->N, 0, NULL, NULL);
    h->kissFFThandle_bkw = kiss_fftr_alloc(h->N, 1, NULL, NULL);
#endif
}

void saf_rfft_destroy
(
    void ** const phFFT
)
{
    saf_rfft_data *h = (saf_rfft_data*)(*phFFT);
    if(h!=NULL){
#if defined(SAF_USE_INTEL_IPP)
        if(h->useIPPfft_FLAG){
            if(h->memSpec)
                ippFree(h->memSpec);
        }
        else {
            if(h->hDFTspec)
                ippFree(h->hDFTspec);
        }
        if(h->buffer)
            ippFree(h->buffer);
#elif defined(__ACCELERATE__)
        vDSP_DFT_DestroySetup(h->DFT_fwd);
        vDSP_DFT_DestroySetup(h->DFT_bwd);
        free(h->VDSP_split_tmp.realp);
        free(h->VDSP_split_tmp.imagp);
        free(h->VDSP_split.realp);
        free(h->VDSP_split.imagp);
#elif defined(INTEL_MKL_VERSION)
        h->Status = DftiFreeDescriptor(&(h->MKL_FFT_Handle));
#else
        kiss_fftr_free(h->kissFFThandle_fwd);
        kiss_fftr_free(h->kissFFThandle_bkw);
#endif
        free(h);
        h=NULL;
    }
}

void saf_rfft_forward
(
    void * const hFFT,
    float* inputTD,
    float_complex* outputFD
)
{
    saf_rfft_data *h = (saf_rfft_data*)(hFFT);

#if defined(SAF_USE_INTEL_IPP)
    if(h->useIPPfft_FLAG)
        ippsFFTFwd_RToCCS_32f((Ipp32f*)inputTD, (Ipp32f*)outputFD, h->hFFTspec, h->buffer);
    else
        ippsDFTFwd_RToCCS_32f((Ipp32f*)inputTD, (Ipp32f*)outputFD, h->hDFTspec, h->buffer);
#elif defined(__ACCELERATE__)
    vDSP_ctoz((DSPComplex*)inputTD, 2, &(h->VDSP_split_tmp), 1, (h->N)/2);
    vDSP_DFT_Execute(h->DFT_fwd, h->VDSP_split_tmp.realp, h->VDSP_split_tmp.imagp, h->VDSP_split.realp, h->VDSP_split.imagp);
    /* DC */
    outputFD[0] = cmplxf(h->VDSP_split.realp[0], 0.0f);
    cblas_scopy(h->N/2-1, &h->VDSP_split.realp[1], 1, &((float*)(outputFD))[2], 2);
    cblas_scopy(h->N/2-1, &h->VDSP_split.imagp[1], 1, &((float*)(outputFD))[3], 2);
    /* the real part of the Nyquist value is the imaginary part of DC. */
    /* https://stackoverflow.com/questions/43289265/implementing-an-fft-using-vdsp */
    outputFD[h->N/2] = cmplxf(h->VDSP_split.imagp[0], 0.0f);
    /* Note: the output is scaled by 2, because vDSP_fft automatically compensates for the loss of energy
     * when removing the symmetric/conjugate (N/2+2:N) bins. However, this is dumb... so the 2x scaling
     * is removed here; so it has parity with the other FFT implementations supported by SAF. */
    cblas_sscal(2*(h->N/2+1), 0.5f, (float*)outputFD, 1);
#elif defined(INTEL_MKL_VERSION)
    h->Status = DftiComputeForward(h->MKL_FFT_Handle, inputTD, outputFD);
#else
    kiss_fftr(h->kissFFThandle_fwd, inputTD, (kiss_fft_cpx*)outputFD);
#endif
}

void saf_rfft_backward
(
    void * const hFFT,
    float_complex* inputFD,
    float* outputTD
)
{
    saf_rfft_data *h = (saf_rfft_data*)(hFFT); 
    
#if defined(SAF_USE_INTEL_IPP)
    if(h->useIPPfft_FLAG)
        ippsFFTInv_CCSToR_32f((Ipp32f*)inputFD, (Ipp32f*)outputTD, h->hFFTspec, h->buffer);
    else
        ippsDFTInv_CCSToR_32f((Ipp32f*)inputFD, (Ipp32f*)outputTD, h->hDFTspec, h->buffer);
#elif defined(__ACCELERATE__)
    h->VDSP_split_tmp.realp[0] = crealf(inputFD[0]);
    h->VDSP_split_tmp.imagp[0] = crealf(inputFD[h->N/2]);
    cblas_scopy(h->N/2-1, &((float*)(inputFD))[2], 2, &h->VDSP_split_tmp.realp[1], 1);
    cblas_scopy(h->N/2-1, &((float*)(inputFD))[3], 2, &h->VDSP_split_tmp.imagp[1], 1);
    vDSP_DFT_Execute(h->DFT_bwd, h->VDSP_split_tmp.realp, h->VDSP_split_tmp.imagp, h->VDSP_split.realp, h->VDSP_split.imagp);
    vDSP_ztoc(&(h->VDSP_split), 1, (DSPComplex*)outputTD, 2, (h->N)/2);
    vDSP_vsmul(outputTD, 1, &(h->Scale), outputTD, 1, h->N);
#elif defined(INTEL_MKL_VERSION)
    h->Status = DftiComputeBackward(h->MKL_FFT_Handle, inputFD, outputTD);
#else
    kiss_fftri(h->kissFFThandle_bkw, (kiss_fft_cpx*)inputFD, outputTD);
    cblas_sscal(h->N, 1.0f/(float)(h->N), outputTD, 1);
#endif
}


/* ========================================================================== */
/*                            Complex<->Complex FFT                           */
/* ========================================================================== */

void saf_fft_create
(
    void ** const phFFT,
    int N
)
{
    *phFFT = malloc1d(sizeof(saf_fft_data));
    saf_fft_data *h = (saf_fft_data*)(*phFFT);
    
    h->N = N;
    h->Scale = 1.0f/(float)N; /* output scaling after ifft */
    saf_assert(N>=2, "Only even (non zero) FFT sizes are supported");
#if defined(SAF_USE_INTEL_IPP)
    /* Use ippsFFT if N is 2^x, otherwise, use ippsDFT */
    if(ceilf(log2f(N)) == floorf(log2f(N))){
        h->useIPPfft_FLAG = 1;
        h->log2n = (int)(log2f((float)N)+0.1f);
        ippsFFTGetSize_C_32f(h->log2n, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &(h->specSize), &(h->specBufferSize), &(h->bufferSize));
        h->hFFTspec = NULL;
        h->memSpec = (Ipp8u*) ippsMalloc_8u(h->specSize);
        h->buffer  = (Ipp8u*) ippsMalloc_8u(h->bufferSize);
        h->memInit = (Ipp8u*) ippsMalloc_8u(h->specBufferSize);
        ippsFFTInit_C_32fc(&(h->hFFTspec), h->log2n, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, h->memSpec, h->memInit);
    }
    else{
        h->useIPPfft_FLAG = 0;
        ippsDFTGetSize_C_32f(N, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &(h->specSize), &(h->specBufferSize), &(h->bufferSize));
        h->hDFTspec = (IppsDFTSpec_C_32fc*) ippsMalloc_8u(h->specSize);
        h->buffer  = (Ipp8u*) ippsMalloc_8u(h->bufferSize);
        h->memInit = (Ipp8u*) ippsMalloc_8u(h->specBufferSize);
        ippsDFTInit_C_32fc(N, IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, h->hDFTspec, h->memInit);
    }
    if (h->memInit)
        ippFree(h->memInit);
#elif defined(__ACCELERATE__)
    h->DFT_fwd = vDSP_DFT_zop_CreateSetup(0, N, vDSP_DFT_FORWARD);
    h->DFT_bwd = vDSP_DFT_zop_CreateSetup(0, N, vDSP_DFT_INVERSE);
    /* Note that DFT lengths must satisfy: f * 2.^g, where f is 1, 3, 5, or 15, and g >=3 */
    saf_assert(h->DFT_fwd!=0 && h->DFT_bwd!=0, "Failed to create vDSP DFT");
    h->VDSP_split_tmp.realp = malloc1d((h->N)*sizeof(float));
    h->VDSP_split_tmp.imagp = malloc1d((h->N)*sizeof(float));
    h->VDSP_split.realp = malloc1d((h->N)*sizeof(float));
    h->VDSP_split.imagp = malloc1d((h->N)*sizeof(float));
#elif defined(INTEL_MKL_VERSION)
    h->MKL_FFT_Handle = 0;
    h->Status = DftiCreateDescriptor( &(h->MKL_FFT_Handle), DFTI_SINGLE,
                                  DFTI_COMPLEX, 1, h->N); /* 1-D, single precision, complex_input_td->fft->complex_input_fd->ifft->complex_output_td */
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE); /* Not inplace, i.e. output has its own dedicated memory */
    /* Configuration parameters for backward-FFT */
    h->Status = DftiSetValue(h->MKL_FFT_Handle, DFTI_BACKWARD_SCALE, h->Scale);      /* scalar applied after ifft */
    /* commit these chosen parameters */
    h->Status = DftiCommitDescriptor(h->MKL_FFT_Handle);
#else
    h->kissFFThandle_fwd = kiss_fft_alloc(h->N, 0, NULL, NULL);
    h->kissFFThandle_bkw = kiss_fft_alloc(h->N, 1, NULL, NULL);
#endif
}

void saf_fft_destroy
(
    void ** const phFFT
)
{
    saf_fft_data *h = (saf_fft_data*)(*phFFT);
    
    if(h!=NULL){
#if defined(SAF_USE_INTEL_IPP)
        if(h->useIPPfft_FLAG){
            if(h->memSpec)
                ippFree(h->memSpec);
        }
        else {
            if(h->hDFTspec)
                ippFree(h->hDFTspec);
        }
        if(h->buffer)
            ippFree(h->buffer);
#elif defined(__ACCELERATE__)
        vDSP_DFT_DestroySetup(h->DFT_fwd);
        vDSP_DFT_DestroySetup(h->DFT_bwd);
        free(h->VDSP_split_tmp.realp);
        free(h->VDSP_split_tmp.imagp);
        free(h->VDSP_split.realp);
        free(h->VDSP_split.imagp);
#elif defined(INTEL_MKL_VERSION)
        h->Status = DftiFreeDescriptor(&(h->MKL_FFT_Handle));
#else
        kiss_fft_free(h->kissFFThandle_fwd);
        kiss_fft_free(h->kissFFThandle_bkw);
#endif
        free(h);
        h=NULL;
    }
}

void saf_fft_forward
(
    void * const hFFT,
    float_complex* inputTD,
    float_complex* outputFD
)
{
    saf_fft_data *h = (saf_fft_data*)(hFFT);
    
#if defined(SAF_USE_INTEL_IPP)
    if(h->useIPPfft_FLAG)
        ippsFFTFwd_CToC_32fc((Ipp32fc*)inputTD, (Ipp32fc*)outputFD, h->hFFTspec, h->buffer);
    else
        ippsDFTFwd_CToC_32fc((Ipp32fc*)inputTD, (Ipp32fc*)outputFD, h->hDFTspec, h->buffer);
#elif defined(__ACCELERATE__)
    cblas_scopy(h->N, &((float*)(inputTD))[0], 2, h->VDSP_split_tmp.realp, 1);
    cblas_scopy(h->N, &((float*)(inputTD))[1], 2, h->VDSP_split_tmp.imagp, 1);
    vDSP_DFT_Execute(h->DFT_fwd, h->VDSP_split_tmp.realp, h->VDSP_split_tmp.imagp, h->VDSP_split.realp, h->VDSP_split.imagp);
    cblas_scopy(h->N, h->VDSP_split.realp, 1, &((float*)(outputFD))[0], 2);
    cblas_scopy(h->N, h->VDSP_split.imagp, 1, &((float*)(outputFD))[1], 2);
#elif defined(INTEL_MKL_VERSION)
    h->Status = DftiComputeForward(h->MKL_FFT_Handle, inputTD, outputFD);
#else
    kiss_fft(h->kissFFThandle_fwd, (kiss_fft_cpx*)inputTD, (kiss_fft_cpx*)outputFD);
#endif
}

void saf_fft_backward
(
    void * const hFFT,
    float_complex* inputFD,
    float_complex* outputTD
)
{
    saf_fft_data *h = (saf_fft_data*)(hFFT);
#if defined(SAF_USE_INTEL_IPP)
    if(h->useIPPfft_FLAG)
        ippsFFTInv_CToC_32fc((Ipp32fc*)inputFD, (Ipp32fc*)outputTD, h->hFFTspec, h->buffer);
    else
        ippsDFTInv_CToC_32fc((Ipp32fc*)inputFD, (Ipp32fc*)outputTD, h->hDFTspec, h->buffer);
#elif defined(__ACCELERATE__)
    cblas_scopy(h->N, &((float*)(inputFD))[0], 2, h->VDSP_split_tmp.realp, 1);
    cblas_scopy(h->N, &((float*)(inputFD))[1], 2, h->VDSP_split_tmp.imagp, 1);
    vDSP_DFT_Execute(h->DFT_bwd, h->VDSP_split_tmp.realp, h->VDSP_split_tmp.imagp, h->VDSP_split.realp, h->VDSP_split.imagp);
    cblas_scopy(h->N, h->VDSP_split.realp, 1, &((float*)(outputTD))[0], 2);
    cblas_scopy(h->N, h->VDSP_split.imagp, 1, &((float*)(outputTD))[1], 2);
    cblas_sscal(2*(h->N), 1.0f/(float)(h->N), (float*)outputTD, 1);
#elif defined(INTEL_MKL_VERSION)
    h->Status = DftiComputeBackward(h->MKL_FFT_Handle, inputFD, outputTD);
#else
    kiss_fft(h->kissFFThandle_bkw, (kiss_fft_cpx*)inputFD, (kiss_fft_cpx*)outputTD);
    cblas_sscal(2*(h->N), 1.0f/(float)(h->N), (float*)outputTD, 1);
#endif
}
