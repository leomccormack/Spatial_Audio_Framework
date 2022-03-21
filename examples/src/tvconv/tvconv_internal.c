/*
 * Copyright 2020 Leo McCormack
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
 * @file tvconv_internal.c
 * @brief A time-varying multi-channel convolver
 * @author Rapolas Daugintis
 * @date 13.07.2021
 */

#include "tvconv_internal.h"

void tvconv_setCodecStatus(void* const hTVCnv, CODEC_STATUS newStatus)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    if(newStatus==CODEC_STATUS_NOT_INITIALISED){
        /* Pause until current initialisation is complete */
        while(pData->codecStatus == CODEC_STATUS_INITIALISING)
            SAF_SLEEP(10);
    }
    pData->codecStatus = newStatus;
}

void tvconv_findNearestNeigbour(void* const hTVCnv)
{
    float dist = 0, minDist = 0;
    int i, d, min_idx = 0;
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    if (pData->nListenerPositions > 0 && pData->listenerPositions != NULL) {
        for(i = 0; i < pData->nListenerPositions; i++){
            for(d = 0; d < NUM_DIMENSIONS; d++)
                dist += (pData->targetPosition[d] - pData->listenerPositions[i][d]) *
                        (pData->targetPosition[d] - pData->listenerPositions[i][d]);
            
            if(dist < minDist || i == 0){
                minDist = dist;
                min_idx = i;
            }
            dist = 0;
        }
    }

    pData->position_idx = min_idx;
}

void tvconv_setMinMaxDimensions(void* const hTVCnv)
{
    int i, d;
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    if(pData->listenerPositions != NULL){
        for(d = 0; d < NUM_DIMENSIONS; d++){
            pData->minDimensions[d] = pData->listenerPositions[0][d];
            pData->maxDimensions[d] = pData->listenerPositions[0][d];
            for(i = 1; i<pData->nListenerPositions; i++){
                if(pData->listenerPositions[i][d] < pData->minDimensions[d])
                    pData->minDimensions[d] = pData->listenerPositions[i][d];
                else if (pData->listenerPositions[i][d] > pData->maxDimensions[d])
                    pData->maxDimensions[d] = pData->listenerPositions[i][d];
            }
        }
        /* resetting current position to the start */
        for(d = 0; d < NUM_DIMENSIONS; d++)
            pData->targetPosition[d] = pData->minDimensions[d];
    }
}
