/*
 * Copyright 2016-2018 Leo McCormack
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
 * @file ambi_enc.c
 * @brief A basic Ambisonic encoder
 *
 * @author Leo McCormack
 * @date 07.10.2016
 */

#include "ambi_enc.h"
#include "ambi_enc_internal.h"

void ambi_enc_create
(
    void ** const phAmbi
)
{
    ambi_enc_data* pData = (ambi_enc_data*)malloc1d(sizeof(ambi_enc_data));
    *phAmbi = (void*)pData;
    int i;
    pData->order = 1;
    
    /* default user parameters */
    loadSourceConfigPreset(SOURCE_CONFIG_PRESET_DEFAULT, pData->src_dirs_deg, &(pData->new_nSources));
    pData->nSources = pData->new_nSources;
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_SH_FLAG[i] = 1;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    pData->order = SH_ORDER_FIRST;
    pData->enablePostScaling = 1;
}

void ambi_enc_destroy
(
    void ** const phAmbi
)
{
    ambi_enc_data *pData = (ambi_enc_data*)(*phAmbi);
    
    if (pData != NULL) {
        free(pData);
        pData = NULL;
    }
}

void ambi_enc_init
(
    void * const hAmbi,
    int          sampleRate
)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    int i;
    
    pData->fs = (float)sampleRate;
    for(i=1; i<=FRAME_SIZE; i++)
        pData->interpolator[i-1] = (float)i*1.0f/(float)FRAME_SIZE;
    memset(pData->prev_Y, 0, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*sizeof(float));
    memset(pData->prev_inputFrameTD, 0, MAX_NUM_INPUTS*FRAME_SIZE*sizeof(float));
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_SH_FLAG[i] = 1;
}

void ambi_enc_process
(
    void  *  const hAmbi,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples
)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    int i, j, ch, nSources, nSH; 
    float src_dirs[MAX_NUM_INPUTS][2], azi_incl[2], scale;
    float* Y_src;

    /* local copies of user parameters */
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    int order;
    chOrdering = pData->chOrdering;
    norm = pData->norm;
    nSources = pData->nSources;
    memcpy(src_dirs, pData->src_dirs_deg, MAX_NUM_INPUTS*2*sizeof(float));
    order = MIN(pData->order, MAX_SH_ORDER);
    nSH = ORDER2NSH(order);

    /* Process frame */
    if (nSamples == FRAME_SIZE) {

        /* prep */
        Y_src = malloc1d(nSH*sizeof(float));

        /* Load time-domain data */
        for(i=0; i < MIN(nSources,nInputs); i++)
            utility_svvcopy(inputs[i], FRAME_SIZE, pData->inputFrameTD[i]);
        for(; i<MAX_NUM_INPUTS; i++)
            memset(pData->inputFrameTD[i], 0, FRAME_SIZE * sizeof(float));

        /* recalulate SHs */
        for(i=0; i<nSources; i++){
            if(pData->recalc_SH_FLAG[i]){
                azi_incl[0] = pData->src_dirs_deg[i][0]*M_PI/180.0f;
                azi_incl[1] =  M_PI/2.0f - pData->src_dirs_deg[i][1]*M_PI/180.0f;
                getSHreal_recur(order, azi_incl, 1, Y_src);
                for(j=0; j<nSH; j++)
                    pData->Y[j][i] = sqrtf(4.0f*M_PI)*Y_src[j];
                for(; j<MAX_NUM_SH_SIGNALS; j++)
                    pData->Y[j][i] = 0.0f;
                pData->recalc_SH_FLAG[i] = 0;
            }
            else{
                for(j=0; j<MAX_NUM_SH_SIGNALS; j++)
                    pData->Y[j][i] = pData->prev_Y[j][i];
            }
        }

        /* spatially encode the input signals into spherical harmonic signals */
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, FRAME_SIZE, nSources, 1.0f,
                    (float*)pData->prev_Y, MAX_NUM_INPUTS,
                    (float*)pData->prev_inputFrameTD, FRAME_SIZE, 0.0f,
                    (float*)pData->tempFrame, FRAME_SIZE);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, FRAME_SIZE, nSources, 1.0f,
                    (float*)pData->Y, MAX_NUM_INPUTS,
                    (float*)pData->prev_inputFrameTD, FRAME_SIZE, 0.0f,
                    (float*)pData->outputFrameTD, FRAME_SIZE);

        for (i=0; i < nSH; i++)
            for(j=0; j<FRAME_SIZE; j++)
                pData->outputFrameTD[i][j] = pData->interpolator[j] * pData->outputFrameTD[i][j] + (1.0f-pData->interpolator[j]) * pData->tempFrame[i][j];

        /* for next frame */
        utility_svvcopy((const float*)pData->inputFrameTD, nSources*FRAME_SIZE, (float*)pData->prev_inputFrameTD);
        utility_svvcopy((const float*)pData->Y, MAX_NUM_INPUTS*MAX_NUM_SH_SIGNALS, (float*)pData->prev_Y);

        /* scale by 1/sqrt(nSources) */
        if(pData->enablePostScaling){
            scale = 1.0f/sqrtf((float)nSources);
            utility_svsmul((float*)pData->outputFrameTD, &scale, nSH*FRAME_SIZE, (float*)pData->outputFrameTD);
        }

        /* account for output channel order */
        switch(chOrdering){
            case CH_ACN: /* already ACN */
                break;
            case CH_FUMA:
                convertHOAChannelConvention((float*)pData->outputFrameTD, order, FRAME_SIZE, HOA_CH_ORDER_ACN, HOA_CH_ORDER_FUMA);
                break;
        }

        /* account for normalisation scheme */
        switch(norm){
            case NORM_N3D: /* already N3D */
                break;
            case NORM_SN3D:
                convertHOANormConvention((float*)pData->outputFrameTD, order, FRAME_SIZE, HOA_NORM_N3D, HOA_NORM_SN3D);
                break;
            case NORM_FUMA:
                convertHOANormConvention((float*)pData->outputFrameTD, order, FRAME_SIZE, HOA_NORM_N3D, HOA_NORM_FUMA);
                break;
        }

        /* Copy to output */
        for(i = 0; i < MIN(nSH,nOutputs); i++)
            utility_svvcopy(pData->outputFrameTD[i], FRAME_SIZE, outputs[i]);
        for(; i < nOutputs; i++)
            memset(outputs[i], 0, FRAME_SIZE * sizeof(float));

        /* clean-up */
        free(Y_src);
    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, FRAME_SIZE*sizeof(float));
    }
}

/* Set Functions */

int ambi_enc_getFrameSize(void)
{
    return FRAME_SIZE;
}

void ambi_enc_refreshParams(void* const hAmbi)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    int i;
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_SH_FLAG[i] = 1;
}

void ambi_enc_setOutputOrder(void* const hAmbi, int newOrder)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    int i;
    if((SH_ORDERS)newOrder != pData->order){
        pData->order = (SH_ORDERS)newOrder;
        for(i=0; i<MAX_NUM_INPUTS; i++)
            pData->recalc_SH_FLAG[i] = 1;
        /* FUMA only supports 1st order */
        if(pData->order!=SH_ORDER_FIRST && pData->chOrdering == CH_FUMA)
            pData->chOrdering = CH_ACN;
        if(pData->order!=SH_ORDER_FIRST && pData->norm == NORM_FUMA)
            pData->norm = NORM_SN3D;
    }
}

void ambi_enc_setSourceAzi_deg(void* const hAmbi, int index, float newAzi_deg)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = MAX(newAzi_deg, -180.0f);
    newAzi_deg = MIN(newAzi_deg, 180.0f);
    pData->recalc_SH_FLAG[index] = 1;
    pData->src_dirs_deg[index][0] = newAzi_deg;
}

void ambi_enc_setSourceElev_deg(void* const hAmbi, int index, float newElev_deg)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    newElev_deg = MAX(newElev_deg, -90.0f);
    newElev_deg = MIN(newElev_deg, 90.0f);
    pData->recalc_SH_FLAG[index] = 1;
    pData->src_dirs_deg[index][1] = newElev_deg;
}

void ambi_enc_setNumSources(void* const hAmbi, int new_nSources)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    int i;
    pData->new_nSources = CLAMP(new_nSources, 1, MAX_NUM_INPUTS);
    pData->nSources = pData->new_nSources;
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_SH_FLAG[i] = 1;
}

void ambi_enc_setInputConfigPreset(void* const hAmbi, int newPresetID)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    int ch;
    loadSourceConfigPreset(newPresetID, pData->src_dirs_deg, &(pData->new_nSources));
    pData->nSources = pData->new_nSources;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_SH_FLAG[ch] = 1;
}

void ambi_enc_setChOrder(void* const hAmbi, int newOrder)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    if((CH_ORDER)newOrder != CH_FUMA || pData->order==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void ambi_enc_setNormType(void* const hAmbi, int newType)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    if((NORM_TYPES)newType != NORM_FUMA || pData->order==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
}

void ambi_enc_setEnablePostScaling(void* const hAmbi, int newStatus)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    pData->enablePostScaling = newStatus;
}


/* Get Functions */

int ambi_enc_getOutputOrder(void* const hAmbi)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    return (int)pData->order;
}

float ambi_enc_getSourceAzi_deg(void* const hAmbi, int index)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    return pData->src_dirs_deg[index][0];
}

float ambi_enc_getSourceElev_deg(void* const hAmbi, int index)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    return pData->src_dirs_deg[index][1];
}

int ambi_enc_getNumSources(void* const hAmbi)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    return pData->new_nSources;
}

int ambi_enc_getMaxNumSources()
{
    return MAX_NUM_INPUTS;
}

int ambi_enc_getNSHrequired(void* const hAmbi)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    return (pData->order+1)*(pData->order+1);
}

int ambi_enc_getChOrder(void* const hAmbi)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    return (int)pData->chOrdering;
}

int ambi_enc_getNormType(void* const hAmbi)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    return (int)pData->norm;
}

int ambi_enc_getEnablePostScaling(void* const hAmbi)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    return pData->enablePostScaling;
}

int ambi_enc_getProcessingDelay()
{
    return FRAME_SIZE;
}
