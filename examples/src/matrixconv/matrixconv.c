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
 * @file matrixconv.c
 * @brief A standard matrix convolver
 * @author Leo McCormack
 * @date 30.09.2019
 * @license ISC
 */

#include "matrixconv.h"
#include "matrixconv_internal.h"

void matrixconv_create
(
    void ** const phMCnv
)
{
    matrixconv_data* pData = (matrixconv_data*)malloc1d(sizeof(matrixconv_data));
    *phMCnv = (void*)pData;

    /* Default user parameters */
    pData->nInputChannels = 1;
    pData->enablePartitionedConv = 0;

    /* internal values */
    pData->hostBlockSize = -1; /* force initialisation */
    pData->inputFrameTD = NULL;
    pData->outputFrameTD = NULL;
    pData->hMatrixConv = NULL;
    pData->filters = NULL;
    pData->reInitFilters = 1;
    pData->nfilters = 0;
    pData->filter_length = 0;
    pData->filter_fs = 0;
    pData->input_wav_length = 0;
    pData->nOutputChannels = 0;

    /* set FIFO buffers */
    pData->FIFO_idx = 0;
    memset(pData->inFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
    memset(pData->outFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
}

void matrixconv_destroy
(
    void ** const phMCnv
)
{
    matrixconv_data *pData = (matrixconv_data*)(*phMCnv);
    
    if (pData != NULL) {
        free(pData->inputFrameTD);
        free(pData->outputFrameTD);
        free(pData->filters);
        saf_matrixConv_destroy(&(pData->hMatrixConv));
        free(pData);
        pData = NULL;
    }
}

void matrixconv_init
(
    void * const hMCnv,
    int          sampleRate,
    int          hostBlockSize
)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);

    pData->host_fs = sampleRate;

    if(pData->hostBlockSize != hostBlockSize){
        pData->hostBlockSize = hostBlockSize;
        pData->hostBlockSize_clamped = SAF_CLAMP(pData->hostBlockSize, MIN_FRAME_SIZE, MAX_FRAME_SIZE);
        pData->reInitFilters = 1;
    }
    
    matrixconv_checkReInit(hMCnv);
} 

void matrixconv_process
(
    void        *  const hMCnv,
    const float *const * inputs,
    float       ** const outputs,
    int                  nInputs,
    int                  nOutputs,
    int                  nSamples
)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    int s, ch, i;
    int numInputChannels, numOutputChannels;
 
    matrixconv_checkReInit(hMCnv);

    /* prep */
    numInputChannels = pData->nInputChannels;
    numOutputChannels = pData->nOutputChannels;

    for(s=0; s<nSamples; s++){
        /* Load input signals into inFIFO buffer */
        for(ch=0; ch<SAF_MIN(SAF_MIN(nInputs,numInputChannels),MAX_NUM_CHANNELS); ch++)
            pData->inFIFO[ch][pData->FIFO_idx] = inputs[ch][s];
        for(; ch<numInputChannels; ch++) /* Zero any channels that were not given */
            pData->inFIFO[ch][pData->FIFO_idx] = 0.0f;

        /* Pull output signals from outFIFO buffer */
        for(ch=0; ch<SAF_MIN(SAF_MIN(nOutputs, numOutputChannels),MAX_NUM_CHANNELS); ch++)
            outputs[ch][s] = pData->outFIFO[ch][pData->FIFO_idx];
        for(; ch<nOutputs; ch++) /* Zero any extra channels */
            outputs[ch][s] = 0.0f;

        /* Increment buffer index */
        pData->FIFO_idx++;

        /* Process frame if inFIFO is full and filters are loaded and saf_matrixConv_apply is ready for it */
        if (pData->FIFO_idx >= pData->hostBlockSize_clamped && pData->reInitFilters == 0 ) {
            pData->FIFO_idx = 0;

            /* Load time-domain data */
            for(i=0; i < numInputChannels; i++)
                utility_svvcopy(pData->inFIFO[i], pData->hostBlockSize_clamped, pData->inputFrameTD[i]);

            /* Apply matrix convolution */
            if(pData->hMatrixConv != NULL && pData->filter_length>0)
                saf_matrixConv_apply(pData->hMatrixConv, FLATTEN2D(pData->inputFrameTD), FLATTEN2D(pData->outputFrameTD));
            /* if the matrix convolver handle has not been initialised yet (i.e. no filters have been loaded) then zero the output */
            else
                memset(FLATTEN2D(pData->outputFrameTD), 0, MAX_NUM_CHANNELS * (pData->hostBlockSize_clamped)*sizeof(float));

            /* copy signals to output buffer */
            for (i = 0; i < SAF_MIN(numOutputChannels, MAX_NUM_CHANNELS); i++)
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

void matrixconv_refreshParams(void* const hMCnv)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    pData->reInitFilters = 1;
}

void matrixconv_checkReInit(void* const hMCnv)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    
    /* reinitialise if needed */
    if ((pData->reInitFilters == 1) && (pData->filters != NULL)) {
        pData->reInitFilters = 2;
        saf_matrixConv_destroy(&(pData->hMatrixConv));
        pData->hMatrixConv = NULL;
        
        /* if length of the loaded wav file was not divisable by the specified number of inputs, then the handle remains NULL,
         * and no convolution is applied */
        pData->hostBlockSize_clamped = SAF_CLAMP(pData->hostBlockSize, MIN_FRAME_SIZE, MAX_FRAME_SIZE);
        if(pData->filter_length>0){
            saf_matrixConv_create(&(pData->hMatrixConv),
                                  pData->hostBlockSize_clamped, /*pData->hostBlockSize,*/
                                  pData->filters,
                                  pData->filter_length,
                                  pData->nInputChannels,
                                  pData->nOutputChannels,
                                  pData->enablePartitionedConv);
        }

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

void matrixconv_setFilters
(
    void* const hMCnv,
    const float** H,
    int numChannels,
    int numSamples,
    int sampleRate
)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    int i;
    saf_assert(numChannels<=MAX_NUM_CHANNELS_FOR_WAV && numChannels > 0 && numSamples > 0, "WAV is limited to 1024 channels");
    
    pData->nOutputChannels = SAF_MIN(numChannels, MAX_NUM_CHANNELS);
    pData->input_wav_length = numSamples;
    pData->nfilters = (pData->nOutputChannels) * (pData->nInputChannels);
    
    /* store the loaded filters */
    pData->filters = realloc1d(pData->filters, numChannels * numSamples * sizeof(float));
    for(i=0; i<numChannels; i++)
        memcpy(&(pData->filters[i*numSamples]), H[i], numSamples * sizeof(float));
    pData->filter_fs = sampleRate;
    
    /* if the number of samples in loaded data is not divisable by the currently specified number of
     * inputs, then the filter length is set to 0 and no further processing is conducted. */
    if(pData->input_wav_length % pData->nInputChannels == 0)
        pData->filter_length = (pData->input_wav_length) / (pData->nInputChannels);
    else
        pData->filter_length = 0;

    pData->reInitFilters = 1;
}

void matrixconv_setEnablePart(void* const hMCnv, int newState)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    if(pData->enablePartitionedConv!=newState){
        pData->enablePartitionedConv = newState;
        pData->reInitFilters = 1;
    }
}

void matrixconv_setNumInputChannels(void* const hMCnv, int newValue)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    pData->nInputChannels = SAF_CLAMP(newValue, 1, MAX_NUM_CHANNELS);
    pData->nfilters = (pData->nOutputChannels) * (pData->nInputChannels);
    if((pData->nOutputChannels > 0) && (pData->input_wav_length % pData->nInputChannels == 0))
        pData->filter_length = (pData->input_wav_length) / (pData->nInputChannels);
    else
        pData->filter_length = 0;
    pData->reInitFilters = 1;
}


/*gets*/

int matrixconv_getEnablePart(void* const hMCnv)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    return pData->enablePartitionedConv;
}

int matrixconv_getNumInputChannels(void* const hMCnv)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    return pData->nInputChannels;
}

int matrixconv_getNumOutputChannels(void* const hMCnv)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    return pData->nOutputChannels;
}

int matrixconv_getHostBlockSize(void* const hMCnv)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    return pData->hostBlockSize;
}

int matrixconv_getNfilters(void* const hMCnv)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    return pData->nfilters;
}

int matrixconv_getFilterLength(void* const hMCnv)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    return pData->filter_length;
}

int matrixconv_getFilterFs(void* const hMCnv)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    return pData->filter_fs;
}

int matrixconv_getHostFs(void* const hMCnv)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    return pData->host_fs;
}

int matrixconv_getProcessingDelay(void* const hMCnv)
{
    matrixconv_data *pData = (matrixconv_data*)(hMCnv);
    return pData->hostBlockSize_clamped;
}

