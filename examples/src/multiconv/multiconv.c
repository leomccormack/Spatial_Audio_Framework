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
 * @license ISC
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

    /* Default user parameters */
    pData->nChannels = 1;
    pData->enablePartitionedConv = 0;
    
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

    /* set FIFO buffers */
    pData->FIFO_idx = 0;
    memset(pData->inFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
    memset(pData->outFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
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
        free(pData->filters);
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
        pData->hostBlockSize_clamped = SAF_CLAMP(pData->hostBlockSize, MIN_FRAME_SIZE, MAX_FRAME_SIZE);
        pData->reInitFilters = 1;
    }
    
    multiconv_checkReInit(hMCnv);
} 


void multiconv_process
(
    void        *  const hMCnv,
    const float *const * inputs,
    float       ** const outputs,
    int                  nInputs,
    int                  nOutputs,
    int                  nSamples
)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    int s, ch, i;
    int numChannels;
 
    multiconv_checkReInit(hMCnv);

    /* prep */
    numChannels = pData->nChannels;

    for(s=0; s<nSamples; s++){
        /* Load input signals into inFIFO buffer */
        for(ch=0; ch<SAF_MIN(SAF_MIN(nInputs,numChannels),MAX_NUM_CHANNELS); ch++)
            pData->inFIFO[ch][pData->FIFO_idx] = inputs[ch][s];
        for(; ch<numChannels; ch++) /* Zero any channels that were not given */
            pData->inFIFO[ch][pData->FIFO_idx] = 0.0f;

        /* Pull output signals from outFIFO buffer */
        for(ch=0; ch<SAF_MIN(SAF_MIN(nOutputs, numChannels),MAX_NUM_CHANNELS); ch++)
            outputs[ch][s] = pData->outFIFO[ch][pData->FIFO_idx];
        for(; ch<nOutputs; ch++) /* Zero any extra channels */
            outputs[ch][s] = 0.0f;

        /* Increment buffer index */
        pData->FIFO_idx++;

        /* Process frame if inFIFO is full and filters are loaded and saf_matrixConv_apply is ready for it */
        if (pData->FIFO_idx >= pData->hostBlockSize_clamped && pData->reInitFilters == 0 ) {
            pData->FIFO_idx = 0;

            /* Load time-domain data */
            for(i=0; i < numChannels; i++)
                utility_svvcopy(pData->inFIFO[i], pData->hostBlockSize_clamped, pData->inputFrameTD[i]);

            /* Apply convolution */
            if(pData->hMultiConv != NULL)
                saf_multiConv_apply(pData->hMultiConv, FLATTEN2D(pData->inputFrameTD), FLATTEN2D(pData->outputFrameTD));
            else
                memset(FLATTEN2D(pData->outputFrameTD), 0, MAX_NUM_CHANNELS * (pData->hostBlockSize_clamped)*sizeof(float));

            /* copy signals to output buffer */
            for (i = 0; i < SAF_MIN(numChannels, MAX_NUM_CHANNELS); i++)
                utility_svvcopy(pData->outputFrameTD[i], pData->hostBlockSize_clamped, pData->outFIFO[i]);
        }
        else if(pData->FIFO_idx >= pData->hostBlockSize_clamped){
            /* clear outFIFO if codec was not ready */
            pData->FIFO_idx = 0;
            memset(pData->outFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
        }
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
        pData->hostBlockSize_clamped = SAF_CLAMP(pData->hostBlockSize, MIN_FRAME_SIZE, MAX_FRAME_SIZE);
        saf_multiConv_create(&(pData->hMultiConv),
                             pData->hostBlockSize_clamped,
                             pData->filters,
                             pData->filter_length,
                             pData->nfilters,
                             pData->enablePartitionedConv);

        /* Resize buffers */
        pData->inputFrameTD  = (float**)realloc2d((void**)pData->inputFrameTD, MAX_NUM_CHANNELS, pData->hostBlockSize_clamped, sizeof(float));
        pData->outputFrameTD = (float**)realloc2d((void**)pData->outputFrameTD, MAX_NUM_CHANNELS, pData->hostBlockSize_clamped, sizeof(float));
        memset(FLATTEN2D(pData->inputFrameTD), 0, MAX_NUM_CHANNELS*(pData->hostBlockSize_clamped)*sizeof(float));

        /* reset FIFO buffers */
        pData->FIFO_idx = 0;
        memset(pData->inFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
        memset(pData->outFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));

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
    pData->nChannels = SAF_CLAMP(newValue, 1, MAX_NUM_CHANNELS);
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

int multiconv_getProcessingDelay(void* const hMCnv)
{
    multiconv_data *pData = (multiconv_data*)(hMCnv);
    return pData->hostBlockSize_clamped;
}
