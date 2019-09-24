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
 * Filename: multiconv.c
 * ---------------------
 * A multi-channel convolver
 *
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 23.09.2019
 */

#include "multiconv.h"
#include "multiconv_internal.h"

void multiconv_create
(
    void ** const phMCnv
)
{
    multiconv_data* pData = (multiconv_data*)malloc1d(sizeof(multiconv_data));
    *phMCnv = (void*)pData;
    
    /* internal values */
    pData->hostBlockSize = -1; /* force initialisation */
    pData->inputFrameTD = NULL;
    pData->outputFrameTD = NULL;
    pData->hMultiConv = NULL;
    
    /* Default user parameters */
    pData->numChannels = 64;
}

void multiconv_destroy
(
    void ** const phMCnv
)
{
    multiconv_data *pData = (multiconv_data*)(*phMCnv);
    
    if (pData != NULL) {
        free(pData->inputFrameTD);
        free(pData->outputFrameTD);
        multiConvPart_destroy(&(pData->hMultiConv));
        free(pData);
        pData = NULL;
    }
}

void multiconv_init
(
    void * const hMCnv,
    int          sampleRate,
    int          hostBlockSize
)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    const int length_h = 9000;
    float* H; //[64][length_h];
    
    
    int i,j;
    if(pData->hostBlockSize != hostBlockSize){
        H = calloc1d(64*length_h,sizeof(float));
        pData->hostBlockSize = hostBlockSize;
        pData->inputFrameTD  = (float**)realloc2d((void**)pData->inputFrameTD, MAX_NUM_CHANNELS, hostBlockSize, sizeof(float));
        pData->outputFrameTD = (float**)realloc2d((void**)pData->outputFrameTD, MAX_NUM_CHANNELS, hostBlockSize, sizeof(float));
        memset(ADR2D(pData->inputFrameTD), 0, MAX_NUM_CHANNELS*hostBlockSize*sizeof(float));
   
//        for(i=0; i<64; i++)
//            for(j=0; j<length_h; j++)
//                H[i][j] =
        multiConvPart_destroy(&(pData->hMultiConv));
        multiConvPart_create(&(pData->hMultiConv), hostBlockSize, (float*)H, length_h, pData->numChannels);
        free(H);
    }
    
    /* starting values */
    
} 

void multiconv_process
(
    void  *  const hMCnv,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples,
    int            isPlaying
)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    int i, j, n;
    int numChannels;
 
    if (nSamples == pData->hostBlockSize) {
        /* prep */
        numChannels = pData->numChannels;
        
//        /* Load time-domain data */
//        for(i=0; i < MIN(numChannels, nInputs); i++)
//            utility_svvcopy(inputs[i], pData->hostBlockSize, pData->inputFrameTD[i]);
//        for(; i<numChannels; i++)
//            memset(pData->inputFrameTD[i], 0, pData->hostBlockSize * sizeof(float)); /* fill remaining channels with zeros */

        multiConvPart_apply(pData->hMultiConv, ADR2D(pData->inputFrameTD), ADR2D(pData->outputFrameTD));
        
        n=0;
        
        
//        if (isPlaying){
//
//
//        }
//        else
//            memset(pData->outputFrameTD, 0, MAX_NUM_CHANNELS*FRAME_SIZE*sizeof(float));
//
//        /* copy signals to output buffer */
//        for (i = 0; i < MIN(numChannels, nOutputs); i++)
//            utility_svvcopy(pData->outputFrameTD[i], FRAME_SIZE, outputs[i]);
//        for (; i < nOutputs; i++)
//            memset(outputs[i], 0, FRAME_SIZE*sizeof(float));
    }
    else{
//        for (i = 0; i < nOutputs; i++)
//            memset(outputs[i], 0, FRAME_SIZE*sizeof(float));
    }
}


/*gets*/
