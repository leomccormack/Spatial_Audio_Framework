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
 * Filename: saf_matrixConv.c
 * --------------------------
 * Matrix convolver. UNTESTED!
 *
 * Dependencies:
 *     Currently, only Intel MKL
 * Author, date created:
 *     Leo McCormack, 06.04.2019
 */

#include "saf_utilities.h"
#include "saf_matrixConv.h"

#if defined(INTEL_MKL_VERSION)

typedef struct _safMatConv_data {
    int hopSize, fftSize, nBins;
    int length_h, length_h_pad, nCHin, nCHout;
    int numOvrlpAddBlocks, length_h_fft;
    void* hFFT;
    float* x_pad, *y_pad, *hx_n, *z_n;
    float* ovrlpAddBuffer;
    float_complex* H_f, *X_n, *HX_n;
    
}safMatConv_data;
 
void matrixConv_create
(
    void ** const phMC,
    int hopSize,
    float* H,         /* nCHout x nCHin x length_h */
    int length_h,
    int nCHin,
    int nCHout
)
{
    *phMC = malloc1d(sizeof(safMatConv_data));
    safMatConv_data *h = (safMatConv_data*)(*phMC);
    int i, j;
    float* h_pad;
    
    h->hopSize = hopSize;
    h->fftSize = 2*(h->hopSize);
    h->nBins = hopSize+1;
    h->length_h = length_h;
    h->nCHin = nCHin;
    h->nCHout = nCHout;
    h->numOvrlpAddBlocks = (int)(ceilf((float)(hopSize+length_h-1)/(float)hopSize)+0.1f);
    h->length_h_pad = hopSize * (h->numOvrlpAddBlocks);
    h->length_h_fft = h->length_h_pad/2+1;
    
    /* Allocate memory for buffers and perform fft on H */
    h->ovrlpAddBuffer = malloc1d(nCHout*(h->length_h_pad)*sizeof(float));
    h->x_pad = calloc1d((h->nCHin)*(h->length_h_pad), sizeof(float)); // CALLOC
    h->y_pad = malloc1d((h->nCHout)*(h->length_h_pad)*sizeof(float));
    h->hx_n = malloc1d((h->length_h_pad) * sizeof(float));
    h->z_n = malloc1d((h->length_h_pad) * sizeof(float));
    h->H_f = malloc1d((h->nCHout)*(h->nCHin)*(h->length_h_fft)*sizeof(float_complex));
    h->X_n = malloc1d((h->nCHin)*(h->length_h_fft)*sizeof(float_complex));
    h->HX_n = malloc1d((h->length_h_fft)*sizeof(float_complex));
    safFFT_create(&(h->hFFT), h->length_h_pad);
    h_pad = calloc1d(h->length_h_pad, sizeof(float));
    for(i=0; i<nCHout; i++){
        for(j=0; j<nCHin; j++){
            memcpy(h_pad, &(H[i*nCHin*length_h+j*length_h]), length_h*sizeof(float));
            safFFT_forward(h->hFFT, h_pad, &(h->H_f[i*nCHin*(h->length_h_fft)+j*(h->length_h_fft)]));
        }
    }
    free(h_pad);
}

void matrixConv_destroy
(
    void ** const phMC
)
{
    safMatConv_data *h = (safMatConv_data*)(*phMC);
    safFFT_destroy(&(h->hFFT));
    free(h->ovrlpAddBuffer);
    free(h->x_pad);
    free(h->y_pad);
    free(h->hx_n);
    free(h->z_n);
    free(h->H_f);
    free(h->X_n);
    free(h->HX_n);
    free(h);
}

void matrixConv_apply
(
    void * const hMC,
    float* inputSig,
    float* outputSig
)
{
    safMatConv_data *h = (safMatConv_data*)(hMC);
    int i, j, k;
    
    /* zero-pad input signals and perform fft */
    for(i=0; i<h->nCHin; i++){
        memcpy(&(h->x_pad[i*(h->length_h_pad)]), &inputSig[i*(h->hopSize)], h->hopSize *sizeof(float));
        safFFT_forward(h->hFFT, &(h->x_pad[i*(h->length_h_pad)]), &(h->X_n[i*(h->length_h_fft)]));
    }
    
    for(i=0; i<h->nCHout; i++){
        /* Apply filter and perform ifft, and sum over input channels */
        memset(h->z_n, 0, (h->length_h_pad) * sizeof(float));
        for(j=0; j<h->nCHin; j++){
            utility_cvvmul(&(h->H_f[i*(h->nCHin)*(h->length_h_fft)+j*(h->length_h_fft)]), &(h->X_n[j*(h->length_h_fft)]), h->length_h_fft, h->HX_n);
            safFFT_backward(h->hFFT, h->HX_n, h->hx_n);
            for(k=0; k<h->length_h_pad; k++)
                h->z_n[k] += h->hx_n[k];
        }
        
        /* over-lap add buffer */
        for(k=0; k<(h->numOvrlpAddBlocks-1)*(h->hopSize); k++)
            h->ovrlpAddBuffer[i*(h->length_h_pad)+k] = h->ovrlpAddBuffer[i*(h->length_h_pad)+k+(h->hopSize)];
        for(; k<(h->length_h_pad); k++)
            h->ovrlpAddBuffer[i*(h->length_h_pad)+k] = 0.0f;
        for(k=0; k<(h->length_h_pad); k++)
            h->ovrlpAddBuffer[i*(h->length_h_pad)+k] += h->z_n[k];
    }
    
    /* truncate buffer and output */
    for(i=0; i<h->nCHout; i++)
        for(k=0; k<h->hopSize; k++)
            outputSig[i*(h->hopSize)+k] = h->ovrlpAddBuffer[i*(h->length_h_pad)+k];
}

#endif
