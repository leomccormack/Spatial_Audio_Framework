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
    //pData->nChannels = 1;
    //pData->enablePartitionedConv = 1;
    
    /* init values */
    pData->irs = NULL;
    pData->positions = NULL;
    
    matrixconv_create(&(pData->hMatrixconv_data));
    
}

void tvconv_destroy
(
    void** const phTVCnv
)
{
    tvconv_data* pData = (tvconv_data*)(*phTVCnv);
    if (pData != NULL){
        free(pData->irs);
        free(pData->positions);
        free(pData);
        if (&(pData->hMatrixconv_data) != NULL)
            matrixconv_destroy(&(pData->hMatrixconv_data));
        pData = NULL;
    }
}

void tvconv_init
(
    void* const hTVCnv,
    int samplerate,
    int hostBlockSize
)
{
    
}

void tvconv_setFiltersAndPositions(void* const hTVCnv)
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
    saf_sofa_close(&sofa);
    tvconv_setMinMaxDimensions(hTVCnv);
}

void tvconv_setSofaFilePath(void* const hTVCnv, const char* path)
{
    tvconv_data *pData = (tvconv_data*)(hTVCnv);
    
    pData->sofa_filepath = malloc1d(strlen(path) + 1);
    strcpy(pData->sofa_filepath, path);
    //binauraliser_refreshSettings(hBin);  // re-init and re-calc
}

void tvconv_test(void* const hTVCnv)
{
    tvconv_setSofaFilePath(hTVCnv, "/Users/dauginr1/Documents/Special_Assignment/rir-interpolation-vst/rirs_unprocessed_M3_3D_rndLP.sofa");
    tvconv_setFiltersAndPositions(hTVCnv);
    tvconv_data* pData = (tvconv_data*) hTVCnv;
    
//    pData->position = malloc1d(sizeof(vectorND));
    //pData->nPositions = 3;
    pData->position[0] = 0.04;
    pData->position[1] = 0.04;
    pData->position[2] = 0.05;
    
    tvconv_findNearestNeigbour(hTVCnv);
    printf("Nearest neighbour index: %i\n", pData->position_idx);
}
