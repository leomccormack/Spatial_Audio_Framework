    /*
 Copyright 2017-2018 Leo McCormack
 
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
 *     ambi_drc_internal.c
 * Description:
 *     A frequency-dependent spherical harmonic domain dynamic range compressor (DRC). The
 *     implementation can also keep track of the frequency-dependent gain factors for
 *     the omnidirectional component over time, for optional plotting. The design utilises
 *     a similar approach as in:
 *         McCormack, L., & Välimäki, V. (2017). "FFT-Based Dynamic Range Compression". in
 *         Proceedings of the 14th Sound and Music Computing Conference, July 5-8, Espoo,
 *         Finland.
 *     The DRC gain factors are determined based on analysing the omnidirectional component.
 *     These gain factors are then applied to the higher-order components, in a such a manner
 *     as to retain the spatial information within them.
 * Dependencies:
 *     saf_utilities, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 07.01.2017
 */

#include "ambi_drc.h"
#include "ambi_drc_internal.h"

/* Adapted from:
 * D. Giannoulis, M. Massberg, and J. D. Reiss, “Digital dynamic range compressor design: Tutorial and analysis,”
 * Journal of the Audio Engineering Society, vol. 60, no. 6, pp. 399–408, June 2012. */
float ambi_drc_gainComputer
(
    float xG,
    float T,
    float R,
    float W
)
{
    float yG;
    if (2.0f*(xG - T) < -W)
        yG = xG;
    else if (2.0f*(fabsf(xG - T)) <= W)
        yG = xG + (1.0f / R - 1.0f) * powf(xG - T + W / 2.0f, 2.0f) / (2.0f*W);
    else if (2.0f*(xG - T) > W)
        yG = T + (xG - T) / R;
    else
        yG = 0.0f;

    return yG;
}

/* Adapted from:
 * D. Giannoulis, M. Massberg, and J. D. Reiss, “Digital dynamic range compressor design: Tutorial and analysis,”
 * Journal of the Audio Engineering Society, vol. 60, no. 6, pp. 399–408, June 2012. */
float ambi_drc_smoothPeakDetector
(
    float xL,
    float yL_z1,
    float alpha_a,
    float alpha_r
)
{
    float yL;
    if (xL > yL_z1)
        yL = alpha_a*yL_z1 + (1.0f - alpha_a) * xL;
    else
        yL = alpha_r*yL_z1 + (1.0f - alpha_r) * xL;

    return yL;
}

void ambi_drc_initTFT
(
    void* const hAmbi
)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    int t, ch;
    
    /* free afSTFT + buffers, if already allocated */
    if (pData->hSTFT != NULL){
        afSTFTfree(pData->hSTFT);
        pData->hSTFT = NULL;
        for (t = 0; t<TIME_SLOTS; t++) {
            for (ch = 0; ch< pData->nSH; ch++) {
                free(pData->STFTFrameTF[t][ch].re);
                free(pData->STFTFrameTF[t][ch].im);
            }
        }
        free2d((void**)pData->tempHopFrameTD, MAX(MAX_NUM_SH_SIGNALS, pData->nSH));
        free2d((void**)pData->STFTFrameTF, TIME_SLOTS);
    }
    
    /* reallocate afSTFT + buffers */
    if (pData->hSTFT == NULL){
    
        afSTFTinit(&(pData->hSTFT), HOP_SIZE, pData->new_nSH, pData->new_nSH, 0, 1);
        pData->STFTFrameTF = (complexVector**)malloc2d(TIME_SLOTS, pData->new_nSH, sizeof(complexVector));
        for(t=0; t<TIME_SLOTS; t++) {
            for(ch=0; ch< pData->new_nSH; ch++) {
                pData->STFTFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
                pData->STFTFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
            }
        }
        pData->tempHopFrameTD = (float**)malloc2d(pData->new_nSH, HOP_SIZE, sizeof(float));
        pData->nSH = pData->new_nSH;
    }
}


void ambi_drc_setInputOrder(INPUT_ORDER inOrder, int* nSH)
{
    switch(inOrder){
        case INPUT_OMNI:
            (*nSH) = 1;
            break;
        case INPUT_ORDER_1:
            (*nSH) = 4;
            break;
        case INPUT_ORDER_2:
            (*nSH) = 9;
            break;
        case INPUT_ORDER_3:
            (*nSH) = 16;
            break;
        case INPUT_ORDER_4:
            (*nSH) = 25;
            break;
        case INPUT_ORDER_5:
            (*nSH) = 36;
            break;
        case INPUT_ORDER_6:
            (*nSH) = 49;
            break;
        case INPUT_ORDER_7:
            (*nSH) = 64;
            break;
    }
}
