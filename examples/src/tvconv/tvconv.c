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
 * @file tvconv.c
 * @brief A time-varying multi-channel convolver
 * @author Rapolas Daugintis
 * @date 18.11.2020
 */

#include "tvconv.h"
#include "tvconv_internal.h"

void tvconv_create
(
    void** const phTVCnv
)
{
    tvconv_data* pData = (tvconv_data*)malloc1d(sizeof(tvconv_data));
    *phTVCnv = (void*)pData;

    printf(SAF_VERSION_LICENSE_STRING);

    /* Default user parameters */
    pData->nChannels = 1;
    pData->enablePartitionedConv = 1;
    
    /* internal values */
}

void tvconv_test(void* const hTVCnv)
{
    tvconv_data* pData = (tvconv_data*) hTVCnv;
//    pData->position = malloc1d(sizeof(vectorND));
    pData->npositions = 3;
    pData->positions = malloc1d(sizeof(vectorND)*pData->npositions);
    for (unsigned int i = 0; i < pData->npositions; i++){
        pData->positions[i][0] = 0.1*i;
        pData->positions[i][1] = 0.1*i;
        pData->positions[i][2] = 0.1*i;
        
        pData->position[i] = 0.04;
        }
    
    tvconv_findNearestNeigbour(hTVCnv);
    printf("Nearest neighbour index: %i\n", pData->position_idx);
}
