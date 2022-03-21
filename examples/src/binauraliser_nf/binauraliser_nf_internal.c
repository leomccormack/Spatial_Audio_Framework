/*
 * Copyright 2022 Michael McCrea, Leo McCormack
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
 * @file: binauraliser_nf_internal.c
 * @brief Convolves input audio (up to 64 channels) with interpolated HRTFs in
 *        the time-frequency domain, and applies optional near-field binaural
 *        filtering, as described in [1].
 *
 * The HRTFs are interpolated by applying amplitude-preserving VBAP gains to the
 * HRTF magnitude responses and inter-aural time differences (ITDs)
 * individually, before being re-combined. The example also allows the user to
 * specify an external SOFA file for the convolution, and rotations of the
 * source directions to accomodate head-tracking.
 *
 * @see [1] S. Spagnol, E. Tavazzi, and F. Avanzini, “Distance rendering and
 *          perception of nearby virtual sound sources with a near-field filter
 *          model,” Applied Acoustics, vol. 115, pp. 61–73, Jan. 2017,
 *          doi: 10.1016/j.apacoust.2016.08.015.
 *
 * @author Michael McCrea, Leo McCormack
 * @date 22.02.2022
 * @license ISC
 */

#include "binauraliser_nf_internal.h"

void binauraliserNF_initTFT
(
    void* const hBin
)
{
    binauraliserNF_data *pData = (binauraliserNF_data*)(hBin);
 
    if(pData->hSTFT==NULL)
        afSTFT_create(&(pData->hSTFT), pData->new_nSources, pData->new_nSources * NUM_EARS, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    else if(pData->new_nSources!=pData->nSources){
        afSTFT_channelChange(pData->hSTFT, pData->new_nSources, pData->new_nSources * NUM_EARS);
        afSTFT_clearBuffers(pData->hSTFT);
    }
    pData->nSources = pData->new_nSources;
}

void binauraliserNF_resetSourceDistances(void* const hBin)
{
    binauraliserNF_data *pData = (binauraliserNF_data*)(hBin);
    
    for(int i=0; i<MAX_NUM_INPUTS; i++)
        pData->src_dists_m[i] = pData->farfield_thresh_m * pData->farfield_headroom;
}







 
