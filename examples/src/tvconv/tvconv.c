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
 * @date 13.07.2021
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

    /* Default user parameters */
    pData->nInputChannels = 1;

    /* internal values */
    pData->hostBlockSize = -1; /* force initialisation */
    pData->inputFrameTD = NULL;
    pData->outputFrameTD = NULL;
    pData->hTVConv = NULL;
    pData->irs = NULL;
    pData->reInitFilters = 1;
    pData->nIrChannels = 0;
    pData->ir_length = 0;
    pData->ir_fs = 0;
    pData->nOutputChannels = 0;
    pData->sofa_filepath = NULL;

    /* set FIFO buffers */
    pData->FIFO_idx = 0;
    memset(pData->inFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
    memset(pData->outFIFO, 0, MAX_NUM_CHANNELS*MAX_FRAME_SIZE*sizeof(float));
    
    /* positions */
    pData->listenerPositions = NULL;
    pData->nListenerPositions = 0;
    pData->position_idx = 0;
    for (int d = 0; d < NUM_DIMENSIONS; d++){
        pData->sourcePosition[d] = 0.0f;
        pData->targetPosition[d] = 0.0f;
        pData->minDimensions[d] = 0.0f;
        pData->maxDimensions[d] = 0.0f;
    }

    /* flags/status */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

void tvconv_destroy
(
    void** const phTVCnv
)
{
    tvconv_data* pData = (tvconv_data*)(*phTVCnv);
    
    if (pData != NULL){
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            SAF_SLEEP(10);
        }
        
        free(pData->inputFrameTD);
        free(pData->outputFrameTD);
        free(pData->irs);
        free(pData->listenerPositions);
        saf_TVConv_destroy(&(pData->hTVConv));
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
        pData->hostBlockSize_clamped = SAF_CLAMP(pData->hostBlockSize, MIN_FRAME_SIZE, MAX_FRAME_SIZE);
        pData->reInitFilters = 1;
        tvconv_setCodecStatus(hTVCnv, CODEC_STATUS_NOT_INITIALISED);
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
    pData->procStatus = PROC_STATUS_ONGOING;
   
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
        if (pData->FIFO_idx >= pData->hostBlockSize_clamped && pData->reInitFilters == 0 &&
            pData->codecStatus == CODEC_STATUS_INITIALISED) {
            pData->FIFO_idx = 0;

            /* Load time-domain data */
            for(i=0; i < numInputChannels; i++)
                utility_svvcopy(pData->inFIFO[i], pData->hostBlockSize_clamped, pData->inputFrameTD[i]);

            if(pData->hTVConv != NULL && pData->ir_length>0){
             saf_TVConv_apply(pData->hTVConv,
                              FLATTEN2D(pData->inputFrameTD),
                              FLATTEN2D(pData->outputFrameTD),
                              pData->position_idx);
            }
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
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}


/*sets*/

void tvconv_refreshParams(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    pData->reInitFilters = 1;
    //tvconv_setCodecStatus(hTVCnv, CODEC_STATUS_NOT_INITIALISED);
}

void tvconv_checkReInit(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    
    while (pData->procStatus == CODEC_STATUS_INITIALISING){
        SAF_SLEEP(10);
    }
    /* reinitialise if needed */
    if ((pData->reInitFilters == 1) && (pData->irs != NULL)) {
        pData->reInitFilters = 2;
        saf_TVConv_destroy(&(pData->hTVConv));
        pData->hTVConv = NULL;
        /* if length of the loaded sofa file was not divisable by the specified number of inputs, then the handle remains NULL,
         * and no convolution is applied */
        pData->hostBlockSize_clamped = SAF_CLAMP(pData->hostBlockSize, MIN_FRAME_SIZE, MAX_FRAME_SIZE);
        if(pData->ir_length>0){
            saf_TVConv_create(&(pData->hTVConv),
                              pData->hostBlockSize_clamped,
                              pData->irs,
                              pData->ir_length,
                              pData->nListenerPositions,
                              pData->nOutputChannels,
                              pData->position_idx);
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
        pData->codecStatus = CODEC_STATUS_INITIALISED;
    }
}

void tvconv_setFiltersAndPositions
(
    void* const hTVCnv
)
{
    tvconv_data* pData = (tvconv_data*) hTVCnv;
    vectorND tmp;
    
    if (pData->codecStatus != CODEC_STATUS_NOT_INITIALISED)
        return; /* re-init not required, or already happening */
    while (pData->procStatus == PROC_STATUS_ONGOING){
        /* re-init required, but we need to wait for the current processing loop to end */
        pData->codecStatus = CODEC_STATUS_INITIALISING; /* indicate that we want to init */
        SAF_SLEEP(10);
    }
    
    /* for progress bar */
    pData->codecStatus = CODEC_STATUS_INITIALISING;
    strcpy(pData->progressBarText,"Initialising");
    pData->progressBar0_1 = 0.0f;
    
    SAF_SOFA_ERROR_CODES error;
    saf_sofa_container sofa;
    int i;
    if(pData->sofa_filepath!=NULL){
        strcpy(pData->progressBarText,"Opening SOFA file");
        pData->progressBar0_1 = 0.2f;
        error = saf_sofa_open(&sofa, pData->sofa_filepath, SAF_SOFA_READER_OPTION_NETCDF);
        
        if(error==SAF_SOFA_OK){
            strcpy(pData->progressBarText,"Loading IRs");
            pData->progressBar0_1 = 0.5f;
            
            pData->ir_fs = (int)sofa.DataSamplingRate;
            pData->ir_length = sofa.DataLengthIR;
            pData->nIrChannels = sofa.nReceivers;
            pData->nListenerPositions = sofa.nListeners;

            /* copy only the first source position, because number of source positions might be incorrect in sofa */
            if(!strcmp(sofa.SourcePositionType, "spherical")){
                memcpy(tmp, sofa.SourcePosition, sizeof(vectorND));
                unitSph2cart((float*)tmp, 1, 1, pData->sourcePosition);
            }
            else
                memcpy(pData->sourcePosition, sofa.SourcePosition, sizeof(vectorND));

            
            pData->irs = (float**)realloc2d((void**)pData->irs, pData->nListenerPositions, pData->nIrChannels*pData->ir_length, sizeof(float));
            int tmp_length = pData->nIrChannels * pData->ir_length;
            for(i=0; i<pData->nListenerPositions; i++){
                memcpy(pData->irs[i], &(sofa.DataIR[i*tmp_length]), tmp_length*sizeof(float));
            }
            
            strcpy(pData->progressBarText,"Loading positions");
            pData->progressBar0_1 = 0.8f;
            
            pData->listenerPositions = (vectorND*)realloc1d((void*)pData->listenerPositions, pData->nListenerPositions*sizeof(vectorND));
            memcpy(pData->listenerPositions, sofa.ListenerPosition, pData->nListenerPositions*sizeof(vectorND));
        }
    }
    
    pData->nOutputChannels = SAF_MIN(pData->nIrChannels, MAX_NUM_CHANNELS);
    saf_sofa_close(&sofa);
    tvconv_setMinMaxDimensions(hTVCnv);
    pData->position_idx = 0;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
    pData->reInitFilters = 1;
    
    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
}

void tvconv_setSofaFilePath(void* const hTVCnv, const char* path)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    
    pData->sofa_filepath = malloc1d(strlen(path) + 1);
    strcpy(pData->sofa_filepath, path);
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    tvconv_setFiltersAndPositions(hTVCnv);
}

void tvconv_setTargetPosition(void* const hTVCnv, float position, int dim){
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    saf_assert(dim >= 0 && dim < NUM_DIMENSIONS, "Dimension out of scope");
    pData->targetPosition[dim] = position;
    tvconv_findNearestNeigbour(hTVCnv);
}


/*gets*/

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

int tvconv_getNumIRs(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->nIrChannels;
}

int tvconv_getNumListenerPositions(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);

    return pData->codecStatus==CODEC_STATUS_INITIALISED ? pData->nListenerPositions : 0;
}

float tvconv_getListenerPosition(void* const hTVCnv, int index, int dim)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->codecStatus==CODEC_STATUS_INITIALISED ? pData->listenerPositions[index][dim] : 0.0f;
}

int tvconv_getListenerPositionIdx(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->position_idx;
}

float tvconv_getTargetPosition(void* const hTVCnv, int dim)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    saf_assert(dim >= 0 && dim < NUM_DIMENSIONS, "Dimension out of scope");
    return (float) pData->targetPosition[dim];
}

float tvconv_getSourcePosition(void* const hTVCnv, int dim)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    saf_assert(dim >= 0 && dim < NUM_DIMENSIONS, "Dimension out of scope");
    return (float) pData->sourcePosition[dim];
}

float tvconv_getMinDimension(void* const hTVCnv, int dim)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    saf_assert(dim >= 0 && dim < NUM_DIMENSIONS, "Dimension out of scope");
    return (float) pData->minDimensions[dim];
}

float tvconv_getMaxDimension(void* const hTVCnv, int dim)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    saf_assert(dim >= 0 && dim < NUM_DIMENSIONS, "Dimension out of scope");
    return (float) pData->maxDimensions[dim];
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

char* tvconv_getSofaFilePath(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    if(pData->sofa_filepath!=NULL)
        return pData->sofa_filepath;
    else
        return "no_file";
}

CODEC_STATUS tvconv_getCodecStatus(void* const hTVCnv)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    return pData->codecStatus;
}
