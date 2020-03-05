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
 * @file multiconv.c
 * @brief A multi-channel convolver
 * @author Leo McCormack
 * @date 23.09.2019
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
    pData->filters = NULL;
    pData->reInitFilters = 1;
    pData->nfilters = 0;
    pData->filter_length = 0;
    pData->filter_fs = 0;
    
    /* Default user parameters */
    pData->nChannels = 1;
    pData->enablePartitionedConv = 0;
}

void multiconv_destroy
(
    void ** const phMCnv
)
{
    multiconv_data *pData = (multiconv_data*)(*phMCnv);
    
    if (pData != NULL) {
        free2d((void***)&(pData->inputFrameTD));
        free2d((void***)&(pData->outputFrameTD));
        free1d((void**)&(pData->filters));
        saf_multiConv_destroy(&(pData->hMultiConv));
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
    
    pData->host_fs = sampleRate;
    if(pData->hostBlockSize != hostBlockSize){
        pData->hostBlockSize = hostBlockSize;
        pData->inputFrameTD  = (float**)realloc2d((void**)pData->inputFrameTD, MAX_NUM_CHANNELS, hostBlockSize, sizeof(float));
        pData->outputFrameTD = (float**)realloc2d((void**)pData->outputFrameTD, MAX_NUM_CHANNELS, hostBlockSize, sizeof(float));
        memset(ADR2D(pData->inputFrameTD), 0, MAX_NUM_CHANNELS*hostBlockSize*sizeof(float));
        pData->reInitFilters = 1;
    }
    
    multiconv_checkReInit(hMCnv);
} 


void multiconv_process
(
    void  *  const hMCnv,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples
)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    int i;
    int numChannels, nFilters;
 
    multiconv_checkReInit(hMCnv);
    
    if (nSamples == pData->hostBlockSize && pData->reInitFilters == 0) {
        /* prep */
        nFilters = pData->nfilters;
        numChannels = pData->nChannels;
        
        /* Load time-domain data */
        for(i=0; i < MIN(MIN(nFilters,numChannels), nInputs); i++)
            utility_svvcopy(inputs[i], pData->hostBlockSize, pData->inputFrameTD[i]);
        for(; i<MAX(nFilters,numChannels); i++)
            memset(pData->inputFrameTD[i], 0, pData->hostBlockSize * sizeof(float)); /* fill remaining channels with zeros */
 
        /* Apply convolution */
        if(pData->hMultiConv != NULL)
            saf_multiConv_apply(pData->hMultiConv, ADR2D(pData->inputFrameTD), ADR2D(pData->outputFrameTD));
        else
            memcpy(ADR2D(pData->outputFrameTD), ADR2D(pData->inputFrameTD), MAX(nFilters,numChannels) * (pData->hostBlockSize)*sizeof(float));
        
        /* copy signals to output buffer */
        for (i = 0; i < MIN(MAX(nFilters,numChannels), nOutputs); i++)
            utility_svvcopy(pData->outputFrameTD[i], pData->hostBlockSize, outputs[i]);
        for (; i < nOutputs; i++)
            memset(outputs[i], 0, pData->hostBlockSize*sizeof(float));
    }
    else{
        for (i = 0; i < nOutputs; i++)
            memset(outputs[i], 0, nSamples*sizeof(float));
    }
}


/*sets*/

void multiconv_refreshParams(void* const hMCnv)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    pData->reInitFilters = 1;
}

void multiconv_checkReInit(void* const hMCnv)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    
    /* reinitialise if needed */
    if ((pData->reInitFilters == 1) && (pData->filters !=NULL)) {
        pData->reInitFilters = 2;
        saf_multiConv_destroy(&(pData->hMultiConv));
        saf_multiConv_create(&(pData->hMultiConv), pData->hostBlockSize, pData->filters, pData->filter_length, pData->nfilters, pData->enablePartitionedConv);
        pData->reInitFilters = 0;
    }
}

void multiconv_setFilters
(
    void* const hMCnv,
    const float** H,
    int numChannels,
    int numSamples,
    int sampleRate
)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    int i;
    
    pData->filters = realloc1d(pData->filters, numChannels*numSamples*sizeof(float));
    pData->nfilters = numChannels;
    pData->filter_length = numSamples;
    for(i=0; i<numChannels; i++)
        memcpy(&(pData->filters[i*numSamples]), H[i], numSamples*sizeof(float));
    pData->filter_fs = sampleRate;
    pData->reInitFilters = 1;
}

void multiconv_setEnablePart(void* const hMCnv, int newState)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    if(pData->enablePartitionedConv!=newState){
        pData->enablePartitionedConv = newState;
        pData->reInitFilters = 1;
    }
}

void multiconv_setNumChannels(void* const hMCnv, int newValue)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    pData->nChannels = CLAMP(newValue, 1, MULTICONV_MAX_NUM_CHANNELS);
}

/*gets*/

int multiconv_getEnablePart(void* const hMCnv)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    return pData->enablePartitionedConv;
}

int multiconv_getNumChannels(void* const hMCnv)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    return pData->nChannels;
}

int multiconv_getHostBlockSize(void* const hMCnv)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    return pData->hostBlockSize;
}

int multiconv_getNfilters(void* const hMCnv)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    return pData->nfilters;
}

int multiconv_getFilterLength(void* const hMCnv)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    return pData->filter_length;
}

int multiconv_getFilterFs(void* const hMCnv)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    return pData->filter_fs;
}

int multiconv_getHostFs(void* const hMCnv)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    return pData->host_fs;
}
