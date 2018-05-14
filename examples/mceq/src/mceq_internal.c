/*
 Copyright 2018 Leo McCormack
 
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
 *     mceq_internal.c
 * Description:
 *     Multi-channel equaliser.
 * Dependencies:
 *     saf_utilities, afSTFTlib,
 * Author, date created:
 *     Leo McCormack, 21.03.2018
 */

#include "mceq_internal.h"

void mceq_initTFT
(
    void* const hMEQ
)
{
    mceq_data *pData = (mceq_data*)(hMEQ);
    int t, ch;
    
    if (pData->hSTFT != NULL){
        afSTFTfree(pData->hSTFT);
        pData->hSTFT = NULL;
        for (t = 0; t<TIME_SLOTS; t++) {
            for (ch = 0; ch< pData->nChannels; ch++) {
                free(pData->STFTInputFrameTF[t][ch].re);
                free(pData->STFTInputFrameTF[t][ch].im);
                free(pData->STFTOutputFrameTF[t][ch].re);
                free(pData->STFTOutputFrameTF[t][ch].im);
            }
        }
        free2d((void**)pData->STFTInputFrameTF, TIME_SLOTS);
        free2d((void**)pData->STFTOutputFrameTF, TIME_SLOTS);
        free2d((void**)pData->tempHopFrameTD, pData->nChannels);
    }
    if (pData->hSTFT == NULL){
        afSTFTinit(&(pData->hSTFT), HOP_SIZE, pData->new_nChannels, pData->new_nChannels, 0, 0);
        pData->STFTInputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, pData->new_nChannels, sizeof(complexVector));
        pData->STFTOutputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, pData->new_nChannels, sizeof(complexVector));
        for(t=0; t<TIME_SLOTS; t++) {
            for(ch=0; ch< pData->new_nChannels; ch++) {
                pData->STFTInputFrameTF[t][ch].re = (float*)calloc(NUM_BANDS, sizeof(float));
                pData->STFTInputFrameTF[t][ch].im = (float*)calloc(NUM_BANDS, sizeof(float));
                pData->STFTOutputFrameTF[t][ch].re = (float*)calloc(NUM_BANDS, sizeof(float));
                pData->STFTOutputFrameTF[t][ch].im = (float*)calloc(NUM_BANDS, sizeof(float));
            }
        }
        pData->tempHopFrameTD = (float**)malloc2d( pData->new_nChannels, HOP_SIZE, sizeof(float));
        pData->nChannels = pData->new_nChannels;
    }
}


void mceq_initFilter
(
    filter* f,
    float freqVector_n[NUM_BANDS],
    float disp_freqVector_n[NUM_DISPLAY_FREQS],
    float fs
)
{
    int band;
    float w, K, KK, D, V0;
    float_complex Hw_num, Hw_denum;
    
    /* calculate the IIR filter coefficients */
    switch (f->type){
        case FILTER_LPF:
            /* Filter design equations - DAFX (2nd ed) p50 */
            K = tanf(M_PI * f->fc/fs);
            KK = K * K;
            D = KK * f->Q + K + f->Q;
            f->b[0] = (KK * f->Q)/D;
            f->b[1] = (2.0f * KK * f->Q)/D;
            f->b[2] = f->b[0];
            f->a[1] = (2.0f * f->Q * (KK - 1.0f))/D;
            f->a[2] = (KK * f->Q - K + f->Q)/D;
            break;
            
        case FILTER_HPF:
            /* Filter design equations - DAFX (2nd ed) p50 */
            K = tanf(M_PI * f->fc/fs);
            KK = K * K;
            D = KK * f->Q + K + f->Q;
            f->b[0] = (f->Q)/D;
            f->b[1] = -(2.0f * f->Q)/D;
            f->b[2] = f->b[0];
            f->a[1] = (2.0f * f->Q * (KK - 1.0f))/D;
            f->a[2] = (KK * f->Q - K + f->Q)/D;
            break;
        case FILTER_LO_SHELF:
            /* Filter design equations - DAFX (2nd ed) p64 */
            K = tanf(M_PI * f->fc/fs);
            V0 = powf(10.0f, (f->G/20.0f));
            if (V0 < 1.0f)
                V0 = 1.0f/V0;
            KK = K * K;
            if (f->G > 0.0f){
                D = 1.0f + sqrtf(2.0f) * K + KK;
                f->b[0] = (1.0f + sqrtf(2.0f * V0) * K + V0 * KK)/D;
                f->b[1] = (2.0f*(V0*KK - 1.0f))/D;
                f->b[2] = (1.0f - sqrtf(2.0f * V0) * K + V0 * KK)/D;
                f->a[1] = (2.0f * (KK - 1.0f))/D;
                f->a[2] = (1.0f - sqrtf(2.0f) * K + KK)/D;
            }
            else{
                D = V0 + sqrtf(2.0f*V0)*K + KK;
                f->b[0] = (V0*(1.0f + sqrtf(2.0f)*K + KK))/D;
                f->b[1] = (2.0f*V0*(KK - 1.0f))/D;
                f->b[2] = (V0*(1.0f - sqrtf(2.0f)*K + KK))/D;
                f->a[1] = (2.0f * (KK - V0))/D;
                f->a[2] = (V0 - sqrtf(2.0f*V0)*K + KK)/D;
            }
            break;
        case FILTER_HI_SHELF:
            /* Filter design equations - DAFX (2nd ed) p64 */
            K = tanf(M_PI * f->fc/fs);
            V0 = powf(10.0f, (f->G/20.0f));
            if (V0 < 1.0f)
                V0 = 1.0f/V0;
            KK = K * K;
            if (f->G > 0.0f){
                D = 1.0f + sqrtf(2.0f) * K + KK;
                f->b[0] = (V0 + sqrtf(2.0f * V0) * K + KK)/D;
                f->b[1] = (2.0f*(KK - V0))/D;
                f->b[2] = (V0 - sqrtf(2.0 * V0) * K + KK)/D;
                f->a[1] = (2.0f*(KK - 1.0f))/D;
                f->a[2] = (1.0f - sqrtf(2.0f) * K + KK)/D;
            }
            else{
                D = 1.0f + sqrtf(2.0f*V0) * K + V0*KK;
                f->b[0] = (V0*(1.0f + sqrtf(2.0f)*K + KK))/D;
                f->b[1] = (2.0f*V0*(KK - 1.0f))/D;
                f->b[2] = (V0*(1.0f - sqrtf(2.0f)*K + KK))/D;
                f->a[1] = (2.0f * (V0*KK - 1.0f))/D;
                f->a[2] = (1.0f - sqrtf(2.0f*V0)*K + V0*KK)/D;
            }
            break;
        case FILTER_PEAK:
            /* Filter design equations - DAFX (2nd ed) p66 */
            K = tanf(M_PI * f->fc/fs);
            V0 = powf(10.0f, (f->G/20.0f));
            KK = K * K;
            if (f->G > 0.0f){
                D = 1.0f + (K/f->Q) + KK;
                f->b[0] = (1.0f + (V0/f->Q) * K + KK)/D;
                f->b[1] = (2.0f*(KK - 1.0f))/D;
                f->b[2] = (1.0f - (V0/f->Q) * K + KK)/D;
                f->a[1] = f->b[1];
                f->a[2] = (1.0f - (K/f->Q) + KK)/D;
            }
            else {
                D = 1.0f + (K/(V0*f->Q)) + KK;
                f->b[0] = (1.0f + (K/f->Q) + KK)/D;
                f->b[1] = (2.0f*(KK - 1.0f))/D;
                f->b[2] = (1.0f - (K/f->Q) + KK)/D;
                f->a[1] = f->b[1];
                f->a[2] = (1.0f - (K/(V0*f->Q)) + KK)/D;
            }
            break;
    }
    /* extract only the magnitude response from the IIR filter */
    for(band=0; band<NUM_BANDS; band++){
        w = freqVector_n[band];
        Hw_num =   cmplxf(f->b[0] + f->b[1]*cosf(w) + f->b[2]*cosf(2.0f*w),
                          -f->b[1]*sinf(w) - f->b[2] * sinf(2.0f*w));
        Hw_denum = cmplxf(1.0f + f->a[1]*cosf(w) + f->a[2]*cosf(2.0f*w),
                          -f->a[1]*sinf(w) - f->a[2] * sinf(2.0f*w));
        f->FBmag[band] = cabsf(ccdivf(Hw_num, Hw_denum));
    }
    /* same for plotting */
    for(band=0; band<NUM_DISPLAY_FREQS; band++){
        w = freqVector_n[band];
        Hw_num =   cmplxf(f->b[0] + f->b[1]*cosf(w) + f->b[2]*cosf(2.0f*w),
                          -f->b[1]*sinf(w) - f->b[2] * sinf(2.0f*w));
        Hw_denum = cmplxf(1.0f + f->a[1]*cosf(w) + f->a[2]*cosf(2.0f*w),
                          -f->a[1]*sinf(w) - f->a[2] * sinf(2.0f*w));
        f->disp_mags[band] = cabsf(ccdivf(Hw_num, Hw_denum));
    }
}

 











 
