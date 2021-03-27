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
    pData->nInputChannels = 1;
    pData->enablePartitionedConv = 0;

    /* internal values */
    pData->hostBlockSize = -1; /* force initialisation */
    pData->inputFrameTD = NULL;
    pData->outputFrameTD = NULL;
    pData->hMatrixConv = NULL;
    pData->irs = NULL;
    pData->reInitFilters = 1;
    pData->nIrChannels = 0;
    pData->ir_length = 0;
    pData->ir_fs = 0;
    //pData->input_wav_length = 0;
    pData->nOutputChannels = 0;

    /* set FIFO buffers */
    pData->FIFO_idx = 0;
    memset(pData->inFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
    memset(pData->outFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
    
    /* positions */
    pData->positions = NULL;
    pData->positions_Last = NULL;
    pData->position_idx = 0;
    pData->position_idx_Last = 0;
    pData->position_idx_Last2 = 0;
}

void tvconv_destroy
(
    void** const phTVCnv
)
{
    tvconv_data* pData = (tvconv_data*)(*phTVCnv);
    
    if (pData != NULL){
        free(pData->inputFrameTD);
        free(pData->outputFrameTD);
        free(pData->irs);
        free(pData->positions);
        free(pData->positions_Last);
        saf_matrixConv_destroy(&(pData->hMatrixConv));
        free(pData);
        pData = NULL;
    }
}

void tvconv_init
(
    void* const hTVCnv,
    int sampleRate,
    int hostBlockSize
)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);

    pData->host_fs = sampleRate;

    if(pData->hostBlockSize != hostBlockSize){
        pData->hostBlockSize = hostBlockSize;
        pData->hostBlockSize_clamped = CLAMP(pData->hostBlockSize, MIN_FRAME_SIZE, MAX_FRAME_SIZE);
        pData->reInitFilters = 1;
    }
    
    tvconv_checkReInit(hTVCnv);
}

void tvconv_process
(
    void  *  const hTVCnv,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples
)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    int s, ch, i;
    int numInputChannels, numOutputChannels;
 
    tvconv_checkReInit(hTVCnv);

    /* prep */
    numInputChannels = pData->nInputChannels;
    numOutputChannels = pData->nOutputChannels;

    for(s=0; s<nSamples; s++){
        /* Load input signals into inFIFO buffer */
        for(ch=0; ch<MIN(MIN(nInputs,numInputChannels),MAX_NUM_CHANNELS); ch++)
            pData->inFIFO[ch][pData->FIFO_idx] = inputs[ch][s];
        for(; ch<numInputChannels; ch++) /* Zero any channels that were not given */
            pData->inFIFO[ch][pData->FIFO_idx] = 0.0f;

        /* Pull output signals from outFIFO buffer */
        for(ch=0; ch<MIN(MIN(nOutputs, numOutputChannels),MAX_NUM_CHANNELS); ch++)
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
            if(pData->hMatrixConv != NULL && pData->ir_length>0)
                saf_matrixConv_apply(pData->hMatrixConv, FLATTEN2D(pData->inputFrameTD), FLATTEN2D(pData->outputFrameTD));
            /* if the matrix convolver handle has not been initialised yet (i.e. no filters have been loaded) then zero the output */
            else
                memset(FLATTEN2D(pData->outputFrameTD), 0, MAX_NUM_CHANNELS * (pData->hostBlockSize_clamped)*sizeof(float));

            /* copy signals to output buffer */
            for (i = 0; i < MIN(numOutputChannels, MAX_NUM_CHANNELS); i++)
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

void tvconv_refreshParams(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    pData->reInitFilters = 1;
}

void tvconv_checkReInit(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    
    /* reinitialise if needed */
    if ((pData->reInitFilters == 1) && (pData->irs != NULL)) {
        pData->reInitFilters = 2;
        saf_matrixConv_destroy(&(pData->hMatrixConv));
        pData->hMatrixConv = NULL;
        
        /* if length of the loaded sofa file was not divisable by the specified number of inputs, then the handle remains NULL,
         * and no convolution is applied */
        pData->hostBlockSize_clamped = CLAMP(pData->hostBlockSize, MIN_FRAME_SIZE, MAX_FRAME_SIZE);
        if(pData->ir_length>0){
            saf_matrixConv_create(&(pData->hMatrixConv),
                                  pData->hostBlockSize_clamped, /*pData->hostBlockSize,*/
                                  pData->irs[pData->position_idx],
                                  pData->ir_length,
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

void tvconv_setFiltersAndPositions
(
    void* const hTVCnv
)
{
    tvconv_data* pData = (tvconv_data*) hTVCnv;
    SAF_SOFA_ERROR_CODES error;
    saf_sofa_container sofa;
    int i;
    if(pData->sofa_filepath!=NULL){
        error = saf_sofa_open(&sofa, pData->sofa_filepath);
        
        if(error==SAF_SOFA_OK){
            pData->ir_fs = (int)sofa.DataSamplingRate;
            pData->ir_length = sofa.DataLengthIR;
            pData->nIrChannels = sofa.nReceivers;
            pData->nPositions = sofa.nListeners;
            pData->irs = (float**)realloc2d((void**)pData->irs, pData->nPositions, pData->nIrChannels*pData->ir_length, sizeof(float));
            int tmp_length = pData->nIrChannels * pData->ir_length;
            for(i=0; i<pData->nPositions; i++){
                memcpy(pData->irs[i], &(sofa.DataIR[i*tmp_length]), tmp_length*sizeof(float));
            }
            pData->positions = (vectorND*)realloc1d((void*)pData->positions, pData->nPositions*sizeof(vectorND));
            memcpy(pData->positions, sofa.ListenerPosition, pData->nPositions*sizeof(vectorND));
        }
    }
    
    pData->nOutputChannels = MIN(pData->nIrChannels, MAX_NUM_CHANNELS);
    saf_sofa_close(&sofa);
    tvconv_setMinMaxDimensions(hTVCnv);
    pData->reInitFilters = 1;
}

void tvconv_setEnablePart(void* const hTVCnv, int newState)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    if(pData->enablePartitionedConv!=newState){
        pData->enablePartitionedConv = newState;
        pData->reInitFilters = 1;
    }
}

void tvconv_setSofaFilePath(void* const hTVCnv, const char* path)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    
    pData->sofa_filepath = malloc1d(strlen(path) + 1);
    strcpy(pData->sofa_filepath, path);
    tvconv_setFiltersAndPositions(hTVCnv);
    pData->reInitFilters = 1;  // re-init and re-calc
}

void tvconv_setNumInputChannels(void* const hTVCnv, int newValue)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    pData->nInputChannels = 1;
//    pData->nInputChannels = CLAMP(newValue, 1, MAX_NUM_CHANNELS);
//    pData->nIrChannels = (pData->nOutputChannels) * (pData->nInputChannels);
//    if((pData->nOutputChannels > 0) && (pData->input_wav_length % pData->nInputChannels == 0))
//        pData->filter_length = (pData->input_wav_length) / (pData->nInputChannels);
//    else
//        pData->filter_length = 0;
    pData->reInitFilters = 1;
}

/*gets*/

int tvconv_getEnablePart(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->enablePartitionedConv;
}

int tvconv_getNumInputChannels(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->nInputChannels;
}

int tvconv_getNumOutputChannels(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->nOutputChannels;
}

int tvconv_getHostBlockSize(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->hostBlockSize;
}

int tvconv_getNIRs(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->nIrChannels;
}

int tvconv_getIRLength(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->ir_length;
}

int tvconv_getIRFs(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->ir_fs;
}

int tvconv_getHostFs(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->host_fs;
}

int tvconv_getProcessingDelay(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->hostBlockSize_clamped;
}

void tvconv_test(void* const hTVCnv)
{
    tvconv_setSofaFilePath(hTVCnv, "/Users/dauginr1/Documents/Special_Assignment/rir-interpolation-vst/rirs_unprocessed_M3_3D_rndLP.sofa");
    //tvconv_setFiltersAndPositions(hTVCnv);
    tvconv_data* pData = (tvconv_data*) hTVCnv;
    
//    pData->position = malloc1d(sizeof(vectorND));
    //pData->nPositions = 3;
    pData->position[0] = 0.04;
    pData->position[1] = 0.04;
    pData->position[2] = 0.05;
    
    tvconv_findNearestNeigbour(hTVCnv);
    printf("Nearest neighbour index: %i\n", pData->position_idx);
}
