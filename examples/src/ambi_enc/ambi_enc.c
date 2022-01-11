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
 * @license ISC
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
    for(i=0; i<MAX_NUM_INPUTS; i++){
        pData->recalc_SH_FLAG[i] = 1;
        pData->src_gains[i] = 1.f;
    }
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
    for(i=1; i<=AMBI_ENC_FRAME_SIZE; i++){
        pData->interpolator_fadeIn[i-1]  = (float)i*1.0f/(float)AMBI_ENC_FRAME_SIZE;
        pData->interpolator_fadeOut[i-1] = 1.0f - pData->interpolator_fadeIn[i-1];
    }
    memset(pData->prev_Y, 0, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*sizeof(float));
    memset(pData->prev_inputFrameTD, 0, MAX_NUM_SH_SIGNALS*AMBI_ENC_FRAME_SIZE*sizeof(float));
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_SH_FLAG[i] = 1;
}

void ambi_enc_process
(
    void        *  const hAmbi,
    const float *const * inputs,
    float       ** const outputs,
    int                  nInputs,
    int                  nOutputs,
    int                  nSamples
)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    int i, j, ch, nSources, nSH, mixWithPreviousFLAG;
    float src_dirs[MAX_NUM_INPUTS][2], scale;
    float Y_src[MAX_NUM_SH_SIGNALS];

    /* local copies of user parameters */
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    int order;
    chOrdering = pData->chOrdering;
    norm = pData->norm;
    nSources = pData->nSources;
    memcpy(src_dirs, pData->src_dirs_deg, MAX_NUM_INPUTS*2*sizeof(float));
    order = SAF_MIN(pData->order, MAX_SH_ORDER);
    nSH = ORDER2NSH(order);

    /* Process frame */
    if (nSamples == AMBI_ENC_FRAME_SIZE) {
        /* Load time-domain data */
        for(i=0; i < SAF_MIN(nSources,nInputs); i++)
            utility_svvcopy(inputs[i], AMBI_ENC_FRAME_SIZE, pData->inputFrameTD[i]);
        for(; i<MAX_NUM_INPUTS; i++)
            memset(pData->inputFrameTD[i], 0, AMBI_ENC_FRAME_SIZE * sizeof(float));

        /* recalulate SHs (only if encoding direction has changed) */
        mixWithPreviousFLAG = 0;
        for(ch=0; ch<nSources; ch++){
            if(pData->recalc_SH_FLAG[ch]){
                getRSH_recur(order, pData->src_dirs_deg[ch], 1, (float*)Y_src);
                for(j=0; j<nSH; j++)
                    pData->Y[j][ch] = Y_src[j];
                for(; j<MAX_NUM_SH_SIGNALS; j++)
                    pData->Y[j][ch] = 0.0f;
                pData->recalc_SH_FLAG[ch] = 0;

                /* If encoding gains have changed, then we should also mix with and interpolate the previous gains */
                mixWithPreviousFLAG = 1;
            }
            /* Apply source gains */
            if(fabsf(pData->src_gains[ch] - 1.f) > 1e-6f)
                utility_svsmul(pData->inputFrameTD[ch], &(pData->src_gains[ch]), AMBI_ENC_FRAME_SIZE, NULL);
        }

        /* spatially encode the input signals into spherical harmonic signals */
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, AMBI_ENC_FRAME_SIZE, nSources, 1.0f,
                    (float*)pData->Y, MAX_NUM_INPUTS,
                    (float*)pData->prev_inputFrameTD, AMBI_ENC_FRAME_SIZE, 0.0f,
                    (float*)pData->outputFrameTD, AMBI_ENC_FRAME_SIZE);

        /* Fade between (linearly inerpolate) the new gains and the previous gains (only if the new gains are different) */
        if(mixWithPreviousFLAG){
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, AMBI_ENC_FRAME_SIZE, nSources, 1.0f,
                        (float*)pData->prev_Y, MAX_NUM_INPUTS,
                        (float*)pData->prev_inputFrameTD, AMBI_ENC_FRAME_SIZE, 0.0f,
                        (float*)pData->tempFrame, AMBI_ENC_FRAME_SIZE);

            /* Apply the linear interpolation */
            for (i=0; i < nSH; i++){
                utility_svvmul((float*)pData->interpolator_fadeIn, (float*)pData->outputFrameTD[i], AMBI_ENC_FRAME_SIZE, (float*)pData->outputFrameTD_fadeIn[i]);
                utility_svvmul((float*)pData->interpolator_fadeOut, (float*)pData->tempFrame[i], AMBI_ENC_FRAME_SIZE, (float*)pData->tempFrame_fadeOut[i]);
            }
            cblas_scopy(nSH*AMBI_ENC_FRAME_SIZE, (float*)pData->outputFrameTD_fadeIn, 1, (float*)pData->outputFrameTD, 1);
            cblas_saxpy(nSH*AMBI_ENC_FRAME_SIZE, 1.0f, (float*)pData->tempFrame_fadeOut, 1, (float*)pData->outputFrameTD, 1);

            /* for next frame */
            utility_svvcopy((const float*)pData->Y, MAX_NUM_INPUTS*MAX_NUM_SH_SIGNALS, (float*)pData->prev_Y);
        }

        /* for next frame */
        utility_svvcopy((const float*)pData->inputFrameTD, MAX_NUM_INPUTS*AMBI_ENC_FRAME_SIZE, (float*)pData->prev_inputFrameTD);

        /* scale by 1/sqrt(nSources) */
        if(pData->enablePostScaling){
            scale = 1.0f/sqrtf((float)nSources);
            cblas_sscal(nSH*AMBI_ENC_FRAME_SIZE, scale, (float*)pData->outputFrameTD, 1);
        }

        /* account for output channel order */
        switch(chOrdering){
            case CH_ACN:  /* already ACN, do nothing */  break;
            case CH_FUMA: convertHOAChannelConvention((float*)pData->outputFrameTD, order, AMBI_ENC_FRAME_SIZE, HOA_CH_ORDER_ACN, HOA_CH_ORDER_FUMA); break;
        }

        /* account for normalisation scheme */
        switch(norm){
            case NORM_N3D:  /* already N3D, do nothing */ break;  
            case NORM_SN3D: convertHOANormConvention((float*)pData->outputFrameTD, order, AMBI_ENC_FRAME_SIZE, HOA_NORM_N3D, HOA_NORM_SN3D); break;
            case NORM_FUMA: convertHOANormConvention((float*)pData->outputFrameTD, order, AMBI_ENC_FRAME_SIZE, HOA_NORM_N3D, HOA_NORM_FUMA); break;
        }

        /* Copy to output */
        for(i = 0; i < SAF_MIN(nSH,nOutputs); i++)
            utility_svvcopy(pData->outputFrameTD[i], AMBI_ENC_FRAME_SIZE, outputs[i]);
        for(; i < nOutputs; i++)
            memset(outputs[i], 0, AMBI_ENC_FRAME_SIZE * sizeof(float));
    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, AMBI_ENC_FRAME_SIZE*sizeof(float));
    }
}

/* Set Functions */

int ambi_enc_getFrameSize(void)
{
    return AMBI_ENC_FRAME_SIZE;
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
    newAzi_deg = SAF_MAX(newAzi_deg, -180.0f);
    newAzi_deg = SAF_MIN(newAzi_deg, 180.0f);
    pData->recalc_SH_FLAG[index] = 1;
    pData->src_dirs_deg[index][0] = newAzi_deg;
}

void ambi_enc_setSourceElev_deg(void* const hAmbi, int index, float newElev_deg)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    newElev_deg = SAF_MAX(newElev_deg, -90.0f);
    newElev_deg = SAF_MIN(newElev_deg, 90.0f);
    pData->recalc_SH_FLAG[index] = 1;
    pData->src_dirs_deg[index][1] = newElev_deg;
}

void ambi_enc_setNumSources(void* const hAmbi, int new_nSources)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    int i;
    pData->new_nSources = SAF_CLAMP(new_nSources, 1, MAX_NUM_INPUTS);
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

void ambi_enc_setSourceGain(void* const hAmbi, int srcIdx, float newGain)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    pData->src_gains[srcIdx] = newGain;
}

void ambi_enc_setSourceSolo(void* const hAmbi, int srcIdx)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    int i;
    for(i=0; i<pData->nSources; i++){
        if(i==srcIdx)
            pData->src_gains[i] = 1.f;
        else
            pData->src_gains[i] = 0.f;
    }
}

void ambi_enc_setUnSolo(void* const hAmbi)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    for(int i=0; i<pData->nSources; i++)
        pData->src_gains[i] = 1.f;
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
    return AMBI_ENC_FRAME_SIZE;
}
