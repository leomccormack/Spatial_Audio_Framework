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
 * @file saf_utility_matrixConv.c
 * @brief Utility: Matrix convolver functions
 *
 * @author Leo McCormack
 * @date 06.04.2019
 */

#include "saf_utilities.h"
#include "saf_utility_matrixConv.h"


/* ========================================================================== */
/*                              Matrix Convolver                              */
/* ========================================================================== */

/**
 * Data structure for the matrix convolver.
 */
typedef struct _safMatConv_data {
    int hopSize, fftSize, nBins;
    int length_h, nCHin, nCHout;
    int numFilterBlocks, numOvrlpAddBlocks;
    int usePartFLAG;
    void* hFFT;
    float* x_pad, *y_pad, *hx_n, *z_n, *ovrlpAddBuffer, *y_n_overlap;
    float_complex* H_f, *X_n, *HX_n;
    float_complex** Hpart_f;
    
}safMatConv_data;
 
void  saf_matrixConv_create
(
    void ** const phMC,
    int hopSize,
    float* H,         /* nCHout x nCHin x length_h */
    int length_h,
    int nCHin,
    int nCHout,
    int usePartFLAG
)
{
    *phMC = malloc1d(sizeof(safMatConv_data));
    safMatConv_data *h = (safMatConv_data*)(*phMC);
    int no, ni, nb;
    float* h_pad, *h_pad_2hops;
    
    h->hopSize = hopSize;
    h->length_h = length_h;
    h->nCHin = nCHin;
    h->nCHout = nCHout;
    h->usePartFLAG = usePartFLAG;
    if(hopSize>length_h && h->usePartFLAG)
        h->usePartFLAG = 0; /* no benefit in partitioning in this case */
    
    if(!h->usePartFLAG){
        /* intialise non-partitioned convolution mode */
        h->numOvrlpAddBlocks = (int)(ceilf((float)(hopSize+length_h-1)/(float)hopSize)+0.1f);
        h->fftSize = (h->numOvrlpAddBlocks)*hopSize;
        h->nBins = h->fftSize/2 + 1;
        
        /* Allocate memory for buffers and perform fft on H */
        h->ovrlpAddBuffer = calloc1d(nCHout*(h->fftSize), sizeof(float));
        h->x_pad = calloc1d((h->nCHin)*(h->fftSize), sizeof(float)); // CALLOC
        h->y_pad = malloc1d((h->nCHout)*(h->fftSize)*sizeof(float));
        h->hx_n = malloc1d((h->fftSize) * sizeof(float));
        h->H_f = malloc1d((h->nCHout)*(h->nCHin)*(h->nBins)*sizeof(float_complex));
        h->X_n = malloc1d((h->nCHin)*(h->nBins)*sizeof(float_complex));
        h->HX_n = malloc1d((h->nBins)*sizeof(float_complex));
        h->z_n = malloc1d((h->fftSize) * sizeof(float));
        saf_rfft_create(&(h->hFFT), h->fftSize);
        h_pad = calloc1d(h->fftSize, sizeof(float));
        for(no=0; no<nCHout; no++){
            for(ni=0; ni<nCHin; ni++){
                memcpy(h_pad, &(H[no*nCHin*length_h+ni*length_h]), length_h*sizeof(float));
                saf_rfft_forward(h->hFFT, h_pad, &(h->H_f[no*nCHin*(h->nBins)+ni*(h->nBins)]));
            }
        }
        free(h_pad);
    }
    else{
        /* intialise partitioned convolution mode */
        h->fftSize = 2*(h->hopSize);
        h->nBins = hopSize+1;
        h->numFilterBlocks = (int)ceilf((float)length_h/(float)hopSize); /* number of partitions */
        assert(h->numFilterBlocks>=1);
        
        /* Allocate memory for buffers and perform fft on partitioned H */
        h_pad = calloc1d(h->numFilterBlocks * hopSize, sizeof(float));
        h_pad_2hops = calloc1d(2 * hopSize, sizeof(float));
        h->Hpart_f = malloc1d(nCHout*sizeof(float_complex*));
        h->X_n = calloc1d(h->numFilterBlocks * nCHin * (h->nBins), sizeof(float_complex));
        h->HX_n = malloc1d(h->numFilterBlocks * nCHin * (h->nBins) * sizeof(float_complex));
        h->x_pad = calloc1d(2 * hopSize, sizeof(float));
        h->hx_n = malloc1d(h->numFilterBlocks*nCHin*(h->fftSize)*sizeof(float));
        h->y_n_overlap = calloc1d(nCHout*hopSize, sizeof(float));
        h->z_n = malloc1d((h->fftSize) * sizeof(float));
        saf_rfft_create(&(h->hFFT), h->fftSize);
        for(no=0; no<nCHout; no++){
            h->Hpart_f[no] = malloc1d(h->numFilterBlocks*nCHin*(h->nBins)*sizeof(float_complex));
            for(ni=0; ni<nCHin; ni++){
                memcpy(h_pad, &H[no*nCHin*length_h+ni*length_h], length_h*sizeof(float)); /* zero pad filter, to be multiple of hopsize */
                for (nb=0; nb<h->numFilterBlocks; nb++){
                    memcpy(h_pad_2hops, &(h_pad[nb*hopSize]), hopSize*sizeof(float));
                    saf_rfft_forward(h->hFFT, h_pad_2hops, &(h->Hpart_f[no][nb*nCHin*(h->nBins)+ni*(h->nBins)]));
                }
            }
        }
        
        free(h_pad);
        free(h_pad_2hops);
    }
}

void saf_matrixConv_destroy
(
    void ** const phMC
)
{
    safMatConv_data *h = (safMatConv_data*)(*phMC);
    int no;
    
    if(h!=NULL){
        saf_rfft_destroy(&(h->hFFT));
        free(h->X_n);
        free(h->x_pad);
        free(h->z_n);
        free(h->hx_n);
        free(h->HX_n);
        if(!h->usePartFLAG){
            free(h->ovrlpAddBuffer);
            free(h->y_pad);
            free(h->H_f);
        }
        else{
            free(h->y_n_overlap);
            for(no=0; no<h->nCHout; no++)
                free(h->Hpart_f[no]);
            free(h->Hpart_f);
        }
        free(h);
        h=NULL;
    }
}

void saf_matrixConv_apply
(
    void * const hMC,
    float* inputSig,
    float* outputSig
)
{
    safMatConv_data *h = (safMatConv_data*)(hMC);
    int ni, no, nb;
    
    /* apply non-partitioned convolution */
    if(!h->usePartFLAG){
        /* zero-pad input signals and perform fft */
        for(ni=0; ni<h->nCHin; ni++){
            memcpy(&(h->x_pad[ni*(h->fftSize)]), &inputSig[ni*(h->hopSize)], h->hopSize *sizeof(float));
            saf_rfft_forward(h->hFFT, &(h->x_pad[ni*(h->fftSize)]), &(h->X_n[ni*(h->nBins)]));
        }
        
        for(no=0; no<h->nCHout; no++){
            /* Apply filter and perform ifft, and sum over input channels */
            memset(h->z_n, 0, (h->fftSize) * sizeof(float));
            for(ni=0; ni<h->nCHin; ni++){



                /* This is the bulk of the CPU work (90.3%): */
                utility_cvvmul(&(h->H_f[no*(h->nCHin)*(h->nBins)+ni*(h->nBins)]), &(h->X_n[ni*(h->nBins)]), h->nBins, h->HX_n);
                saf_rfft_backward(h->hFFT, h->HX_n, h->hx_n);
                utility_svvadd(h->z_n,  h->hx_n, h->fftSize, h->z_n); /* 7.8% */
            }
            
            /* over-lap add buffer */
            memcpy(&(h->ovrlpAddBuffer[no*(h->fftSize)]), &(h->ovrlpAddBuffer[no*(h->fftSize)+(h->hopSize)]), (h->numOvrlpAddBlocks-1)*(h->hopSize)*sizeof(float));
            memset(&(h->ovrlpAddBuffer[no*(h->fftSize)+(h->numOvrlpAddBlocks-1)*(h->hopSize)]), 0, (h->hopSize)*sizeof(float));

            /* sum with overlap buffer and copy the result to the output buffer */
            utility_svvadd(&(h->ovrlpAddBuffer[no*(h->fftSize)]),  h->z_n, (h->fftSize), &(h->ovrlpAddBuffer[no*(h->fftSize)]));

            /* truncate buffer and output */
            memcpy(&(outputSig[no*(h->hopSize)]), &(h->ovrlpAddBuffer[no*(h->fftSize)]), h->hopSize*sizeof(float));
        }
    }
    /* apply partitioned convolution */
    else{
        /* zero-pad input signals and perform fft. Store in partition slot 1. */
        memcpy(&(h->X_n[1*(h->nCHin)*(h->nBins)]), h->X_n, (h->numFilterBlocks-1)*(h->nCHin)*(h->nBins)*sizeof(float_complex)); /* shuffle */
        for(ni=0; ni<h->nCHin; ni++){
            memcpy(h->x_pad, &(inputSig[ni*(h->hopSize)]), h->hopSize *sizeof(float));
            saf_rfft_forward(h->hFFT, h->x_pad, &(h->X_n[0*(h->nCHin)*(h->nBins)+ni*(h->nBins)]));
        }
        
        /* apply convolution and inverse fft */
        for(no=0; no<h->nCHout; no++){
            utility_cvvmul(h->Hpart_f[no], h->X_n, h->numFilterBlocks * (h->nCHin) * (h->nBins), h->HX_n); /* This is the bulk of the CPU work */
            for(nb=0; nb<h->numFilterBlocks; nb++)
                for(ni=0; ni<h->nCHin; ni++)
                    saf_rfft_backward(h->hFFT, &(h->HX_n[nb*(h->nCHin)*(h->nBins)+ni*(h->nBins)]), &(h->hx_n[nb*(h->nCHin)*(h->fftSize)+ni*(h->fftSize)]));
            
            /* output frame for this channel is the sum over all partitions and input channels */
            memset(h->z_n, 0, (h->fftSize) * sizeof(float));
            for(nb=0; nb<h->numFilterBlocks*(h->nCHin); nb++)
                utility_svvadd(h->z_n, (const float*)&(h->hx_n[nb*(h->fftSize)]), h->fftSize, h->z_n);

            /* sum with overlap buffer and copy the result to the output buffer */
            utility_svvadd(h->z_n, (const float*)&(h->y_n_overlap[no*(h->hopSize)]), h->hopSize, &(outputSig[no*(h->hopSize)]));

            /* for next iteration: */
            memcpy(&(h->y_n_overlap[no*(h->hopSize)]), &(h->z_n[h->hopSize]), h->hopSize*sizeof(float));
        }
    }
}


/* ========================================================================== */
/*                           Multi-Channel Convolver                          */
/* ========================================================================== */

/**
 * Data structure for the multi-channel convolver.
 */
typedef struct _safMulConv_data {
    int hopSize, fftSize, nBins;
    int length_h, nCH;
    int numOvrlpAddBlocks, numFilterBlocks;
    int usePartFLAG;
    void* hFFT;
    float* x_pad, *z_n, *ovrlpAddBuffer, *hx_n, *y_n_overlap;
    float_complex* X_n, *HX_n, *Z_n, *H_f, *Hpart_f;
    
}safMulConv_data;

void saf_multiConv_create
(
    void ** const phMC,
    int hopSize,
    float* H,         /* nCH x length_h */
    int length_h,
    int nCH,
    int usePartFLAG
)
{
    *phMC = malloc1d(sizeof(safMulConv_data));
    safMulConv_data *h = (safMulConv_data*)(*phMC);
    int nc, nb;
    float* h_pad, *h_pad_2hops;
    
    h->hopSize = hopSize;
    h->length_h = length_h;
    h->nCH = nCH;
    h->usePartFLAG = usePartFLAG;
    if(hopSize>length_h && h->usePartFLAG)
        h->usePartFLAG = 0; /* no benefit in partitioning in this case */
    
    if(!h->usePartFLAG){
        /* intialise non-partitioned convolution mode */
        h->numOvrlpAddBlocks = (int)(ceilf((float)(hopSize+length_h-1)/(float)hopSize)+0.1f);
        h->fftSize = (h->numOvrlpAddBlocks*hopSize);
        h->nBins = h->fftSize/2 + 1;
        
        /* Allocate memory for buffers and perform fft on partitioned H */
        h->ovrlpAddBuffer = calloc1d(nCH*h->fftSize, sizeof(float));
        h_pad = calloc1d(h->fftSize, sizeof(float));
        h->H_f = malloc1d(nCH*(h->nBins)*sizeof(float_complex));
        h->X_n = calloc1d(nCH * (h->nBins), sizeof(float_complex));
        h->Z_n = malloc1d(nCH * (h->nBins) * sizeof(float_complex));
        h->x_pad = calloc1d(h->fftSize, sizeof(float));
        h->z_n = malloc1d(nCH*(h->fftSize)*sizeof(float));
        saf_rfft_create(&(h->hFFT), h->fftSize);
        for(nc=0; nc<nCH; nc++){
            memcpy(h_pad, &H[nc*length_h], length_h*sizeof(float)); /* zero pad filter, to be multiple of hopsize */
            saf_rfft_forward(h->hFFT, h_pad, &(h->H_f[nc*(h->nBins)]));
        }
        
        free(h_pad);
    }
    else{
        /* intialise partitioned convolution mode */
        h->fftSize = 2*(h->hopSize);
        h->nBins = hopSize+1;
        h->numFilterBlocks = (int)ceilf((float)length_h/(float)hopSize); /* number of partitions */
        assert(h->numFilterBlocks>=1);
        
        /* Allocate memory for buffers and perform fft on partitioned H */
        h_pad = calloc1d(h->numFilterBlocks * hopSize, sizeof(float));
        h_pad_2hops = calloc1d(2 * hopSize, sizeof(float));
        h->Hpart_f = malloc1d(h->numFilterBlocks*nCH*(h->nBins)*sizeof(float_complex));
        h->X_n = calloc1d(h->numFilterBlocks * nCH * (h->nBins), sizeof(float_complex));
        h->HX_n = calloc1d(h->numFilterBlocks * nCH * (h->nBins), sizeof(float_complex));
        h->x_pad = calloc1d(2 * hopSize, sizeof(float));
        h->hx_n = malloc1d(h->numFilterBlocks*nCH*(h->fftSize)*sizeof(float));
        h->z_n = calloc1d(h->fftSize, sizeof(float));
        h->y_n_overlap = calloc1d(nCH*hopSize, sizeof(float));
        saf_rfft_create(&(h->hFFT), h->fftSize);
        for(nc=0; nc<nCH; nc++){
            memcpy(h_pad, &H[nc*length_h], length_h*sizeof(float)); /* zero pad filter, to be multiple of hopsize */
            for (nb=0; nb<h->numFilterBlocks; nb++){
                memcpy(h_pad_2hops, &(h_pad[nb*hopSize]), hopSize*sizeof(float));
                saf_rfft_forward(h->hFFT, h_pad_2hops, &(h->Hpart_f[nb*nCH*(h->nBins)+nc*(h->nBins)]));
            }
        }
        
        free(h_pad);
        free(h_pad_2hops);
    }
}

void saf_multiConv_destroy
(
    void ** const phMC
)
{
    safMulConv_data *h = (safMulConv_data*)(*phMC);
    
    if(h!=NULL){
        saf_rfft_destroy(&(h->hFFT));
        free(h->X_n);
        free(h->x_pad);
        free(h->z_n);
        if(!h->usePartFLAG){
            free(h->Z_n);
            free(h->H_f);
        }
        else{
            free(h->HX_n);
            free(h->hx_n);
            free(h->y_n_overlap);
            free(h->Hpart_f);
        }
        free(h);
        h=NULL;
    }
}

void saf_multiConv_apply
(
    void * const hMC,
    float* inputSig,
    float* outputSig
)
{
    safMulConv_data *h = (safMulConv_data*)(hMC);
    int nc, nb;
    
    /* apply non-partitioned convolution */
    if(!h->usePartFLAG){
        /* zero-pad input signals and perform fft. */
        for(nc=0; nc<h->nCH; nc++){
            memcpy(h->x_pad, &(inputSig[nc*(h->hopSize)]), h->hopSize *sizeof(float));
            saf_rfft_forward(h->hFFT, h->x_pad, &(h->X_n[nc*(h->nBins)]));
        }
        
        /* apply convolution and inverse fft */
        utility_cvvmul(h->H_f, h->X_n, (h->nCH) * (h->nBins), h->Z_n); /* This is the bulk of the CPU work */
        for(nc=0; nc<h->nCH; nc++){
            saf_rfft_backward(h->hFFT, &(h->Z_n[nc*(h->nBins)]), &(h->z_n[nc*(h->fftSize)]));
            
            /* sum with overlap buffer and copy the result to the output buffer */
            utility_svvcopy(&(h->ovrlpAddBuffer[nc*(h->fftSize)+(h->hopSize)]), (h->numOvrlpAddBlocks-1)*(h->hopSize), &(h->ovrlpAddBuffer[nc*(h->fftSize)]));
            memset(&(h->ovrlpAddBuffer[nc*(h->fftSize)+(h->numOvrlpAddBlocks-1)*(h->hopSize)]), 0, (h->hopSize)*sizeof(float));
            utility_svvadd(&(h->ovrlpAddBuffer[nc*(h->fftSize)]),  &(h->z_n[nc*(h->fftSize)]), (h->fftSize), &(h->ovrlpAddBuffer[nc*(h->fftSize)]));
            utility_svvcopy(&(h->ovrlpAddBuffer[nc*(h->fftSize)]), h->hopSize, &(outputSig[nc*(h->hopSize)]));
        }
    }
    /* apply partitioned convolution */
    else{
        /* zero-pad input signals and perform fft. Store in partition slot 1. */
        memcpy(&(h->X_n[1*(h->nCH)*(h->nBins)]), h->X_n, (h->numFilterBlocks-1)*(h->nCH)*(h->nBins)*sizeof(float_complex));
        for(nc=0; nc<h->nCH; nc++){
            memcpy(h->x_pad, &(inputSig[nc*(h->hopSize)]), h->hopSize * sizeof(float));
            saf_rfft_forward(h->hFFT, h->x_pad, &(h->X_n[0*(h->nCH)*(h->nBins)+nc*(h->nBins)]));
        }
        
        /* apply convolution and inverse fft */
        utility_cvvmul(h->Hpart_f, h->X_n, h->numFilterBlocks * (h->nCH) * (h->nBins), h->HX_n); /* This is the bulk of the CPU work */
        for(nc=0; nc<h->nCH; nc++){
            for(nb=0; nb<h->numFilterBlocks; nb++)
                saf_rfft_backward(h->hFFT, &(h->HX_n[nb*(h->nCH)*(h->nBins)+nc*(h->nBins)]), &(h->hx_n[nb*(h->nCH)*(h->fftSize)+nc*(h->fftSize)]));
            
            /* output frame for this channel is the sum over all partitions */
            memset(h->z_n, 0, h->fftSize*sizeof(float));
            for(nb=0; nb<h->numFilterBlocks; nb++)
                utility_svvadd(h->z_n, (const float*)&(h->hx_n[nb*(h->nCH)*(h->fftSize)+nc*(h->fftSize)]), h->fftSize, h->z_n);
            
            /* sum with overlap buffer and copy the result to the output buffer */
            utility_svvadd(h->z_n, (const float*)&(h->y_n_overlap[nc*(h->hopSize)]), h->hopSize, &(outputSig[nc* (h->hopSize)]));
            
            /* for next iteration: */
            memcpy(&(h->y_n_overlap[nc*(h->hopSize)]), &(h->z_n[h->hopSize]), h->hopSize*sizeof(float));
        }
    }
}

