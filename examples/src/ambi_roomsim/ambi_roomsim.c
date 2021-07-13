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
 * @file ambi_roomsim.c
 * @brief A simple shoebox room Ambisonic encoder.
 *
 * @author Leo McCormack
 * @date 10.08.2020
 * @license ISC
 */

#include "ambi_roomsim.h"
#include "ambi_roomsim_internal.h"

/** Default absorption coefficients per wall */
const float default_abs_wall[6] = { 0.341055000f, 0.431295000f, 0.351295000f, 0.344335000f, 0.401775000f, 0.482095000f};
/** Default room dimensions */
const float default_room_dims[3] = { 9.1f, 8.0f, 3.0f };

void ambi_roomsim_create
(
    void ** const phAmbi
)
{
    ambi_roomsim_data* pData = (ambi_roomsim_data*)malloc1d(sizeof(ambi_roomsim_data));
    *phAmbi = (void*)pData;
    
    /* default user parameters */
    pData->enableReflections = 1;
    pData->sh_order = 3;
    pData->refl_order = 3;
    pData->nSources = 1;
    pData->nReceivers = 1;
    memcpy(pData->abs_wall, default_abs_wall, 6*sizeof(float));
    memcpy(pData->room_dims, default_room_dims, 3*sizeof(float));
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    memset(pData->src_pos, 0, ROOM_SIM_MAX_NUM_SOURCES*3*sizeof(float));
    memset(pData->rec_pos, 0, ROOM_SIM_MAX_NUM_RECEIVERS*3*sizeof(float));

    /* Internal */
    pData->hIms = NULL;
    float src_pos[3]  = {5.2f, 1.5f, 1.4f};
    memcpy(pData->src_pos[0], src_pos, 3*sizeof(float));
    float src2_pos[3] = {2.1f, 1.0f, 1.3f};
    memcpy(pData->src_pos[1], src2_pos, 3*sizeof(float));
    float src3_pos[3] = {3.1f, 5.0f, 2.3f};
    memcpy(pData->src_pos[2], src3_pos, 3*sizeof(float));
    float src4_pos[3] = {7.1f, 2.0f, 1.4f};
    memcpy(pData->src_pos[3], src4_pos, 3*sizeof(float));
    float rec_pos[3]  = {5.2f, 3.5f, 1.4f};
    memcpy(pData->rec_pos[0], rec_pos, 3*sizeof(float));
    memcpy(pData->rec_pos[1], rec_pos, 3*sizeof(float));
    pData->new_sh_order = pData->sh_order;
    pData->new_nSources = pData->nSources;
    pData->new_nReceivers = pData->nReceivers;

    pData->src_sigs = (float**)malloc2d(MAX_NUM_CHANNELS, AMBI_ROOMSIM_FRAME_SIZE, sizeof(float));
    pData->rec_sh_outsigs = (float***)malloc3d(IMS_MAX_NUM_RECEIVERS, MAX_NUM_SH_SIGNALS, AMBI_ROOMSIM_FRAME_SIZE, sizeof(float));

    pData->reinit_room = 1;
}

void ambi_roomsim_destroy
(
    void ** const phAmbi
)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(*phAmbi);
    
    if (pData != NULL) {
        ims_shoebox_destroy(&(pData->hIms));
        free(pData->src_sigs);
        free(pData->rec_sh_outsigs);
        free(pData);
        pData = NULL;
    }
}

void ambi_roomsim_init
(
    void * const hAmbi,
    int          sampleRate
)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    pData->fs = (float)sampleRate;
}

void ambi_roomsim_process
(
    void        *  const hAmbi,
    const float *const * inputs,
    float       ** const outputs,
    int                  nInputs,
    int                  nOutputs,
    int                  nSamples
)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    int i, j, rec, ch, nSources, nReceivers, nSH, order;
    float maxTime_s;
    CH_ORDER chOrdering;
    NORM_TYPES norm;

    /* (ims_shoebox is actually much more flexible than this. So consider this as a minimal example, which will also make things easier when designing a GUI) */

    /* Reinitialise room if needed */
    if(pData->reinit_room){
        ims_shoebox_destroy(&(pData->hIms));
        ims_shoebox_create(&(pData->hIms), pData->room_dims, (float*)pData->abs_wall, 250.0f, 1, 343.0f, pData->fs);
        for(i=0; i<pData->new_nSources; i++) /* re-add source objects... */
            pData->sourceIDs[i] = ims_shoebox_addSource(pData->hIms, (float*)pData->src_pos[i], &(pData->src_sigs[i]));
        for(i=0; i<pData->new_nReceivers; i++) /* re-add receiver objects... */
            pData->receiverIDs[i] = ims_shoebox_addReceiverSH(pData->hIms, pData->new_sh_order, (float*)pData->rec_pos[i], &(pData->rec_sh_outsigs[i]));
        pData->nSources = pData->new_nSources;
        pData->nReceivers = pData->new_nReceivers;
        pData->sh_order = pData->new_sh_order;
        pData->reinit_room = 0;
    }

    /* Add/remove source objects */
    if(pData->new_nSources!=pData->nSources){
        if(pData->new_nSources>pData->nSources)
            for(i=pData->nSources; i<pData->new_nSources; i++)
                pData->sourceIDs[i] = ims_shoebox_addSource(pData->hIms, (float*)pData->src_pos[i], &(pData->src_sigs[i]));
        else
            for(i=pData->new_nSources; i<pData->nSources; i++)
                ims_shoebox_removeSource(pData->hIms, pData->sourceIDs[i]);
        pData->nSources = pData->new_nSources;
    }

    /* Add/remove receiver objects */
    if(pData->new_nReceivers!=pData->nReceivers){
        if(pData->new_nReceivers>pData->nReceivers)
            for(i=pData->nReceivers; i<pData->new_nReceivers; i++)
                pData->receiverIDs[i] = ims_shoebox_addReceiverSH(pData->hIms, pData->sh_order, (float*)pData->rec_pos[i], &(pData->rec_sh_outsigs[i]));
        else
            for(i=pData->new_nReceivers; i<pData->nReceivers; i++)
                ims_shoebox_removeReceiver(pData->hIms, pData->receiverIDs[i]);
        pData->nReceivers = pData->new_nReceivers;
    }

    /* local copies of user parameters */
    chOrdering = pData->chOrdering;
    norm = pData->norm;
    order = SAF_MIN(pData->sh_order, MAX_SH_ORDER);
    nSH = ORDER2NSH(order);
    nSources = pData->nSources;
    nReceivers = pData->nReceivers;
    maxTime_s = -0.05f; /* 50ms */

    /* Process frame */
    if (nSamples == AMBI_ROOMSIM_FRAME_SIZE) {
        /* Load time-domain data */
        for(i=0; i < SAF_MIN(nSources,nInputs); i++)
            memcpy(pData->src_sigs[i], inputs[i], AMBI_ROOMSIM_FRAME_SIZE * sizeof(float));
        for(; i < nInputs; i++)
            memset(pData->src_sigs[i], 0, AMBI_ROOMSIM_FRAME_SIZE * sizeof(float));

        /* Update source/receiver positions, room dims/coeffs and re-compute echgrams
         * (note, if nothing has changed since last frame then these calls will be bypassed internally) */
        for(i=0; i<nSources; i++)
            ims_shoebox_updateSource(pData->hIms, pData->sourceIDs[i], pData->src_pos[i]);
        for(i=0; i<nReceivers; i++)
            ims_shoebox_updateReceiver(pData->hIms, pData->receiverIDs[i], pData->rec_pos[i]);
        ims_shoebox_setRoomDimensions(pData->hIms, pData->room_dims);
        ims_shoebox_setWallAbsCoeffs(pData->hIms, (float*)pData->abs_wall);
        ims_shoebox_computeEchograms(pData->hIms, pData->enableReflections ? pData->refl_order : 0, maxTime_s);

        /* Render audio for each receiver */
        for(i=0; i<nReceivers; i++)
            ims_shoebox_applyEchogramTD(pData->hIms, pData->receiverIDs[i], nSamples, 0);

        /* Handle output */
        for(rec=0, i=0; rec<nReceivers; rec++){
            /* account for output channel order */
            switch(chOrdering){
                case CH_ACN: break;
                case CH_FUMA: convertHOAChannelConvention(FLATTEN2D(pData->rec_sh_outsigs[rec]), order, AMBI_ROOMSIM_FRAME_SIZE, HOA_CH_ORDER_ACN, HOA_CH_ORDER_FUMA); break;
            }

            /* account for normalisation scheme */
            switch(norm){
                case NORM_N3D: break;
                case NORM_SN3D: convertHOANormConvention(FLATTEN2D(pData->rec_sh_outsigs[rec]), order, AMBI_ROOMSIM_FRAME_SIZE, HOA_NORM_N3D, HOA_NORM_SN3D); break;
                case NORM_FUMA: convertHOANormConvention(FLATTEN2D(pData->rec_sh_outsigs[rec]), order, AMBI_ROOMSIM_FRAME_SIZE, HOA_NORM_N3D, HOA_NORM_FUMA); break;
            }

            /* Append this receiver's output channels to the master output buffer */
            for(j=0; (j<SAF_MIN(nSH, MAX_NUM_SH_SIGNALS) && i<SAF_MIN(nOutputs, MAX_NUM_CHANNELS)); j++, i++)
                memcpy(outputs[i], pData->rec_sh_outsigs[rec][j], AMBI_ROOMSIM_FRAME_SIZE * sizeof(float));
        }
        for(; i < nOutputs; i++)
            memset(outputs[i], 0, AMBI_ROOMSIM_FRAME_SIZE * sizeof(float));
    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, nSamples*sizeof(float));
    }
}

/* Set Functions */

void ambi_roomsim_refreshParams(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    pData->reinit_room = 1;
}

void ambi_roomsim_setEnableIMSflag(void* const hAmbi, int newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    pData->enableReflections = newValue;
}

void ambi_roomsim_setMaxReflectionOrder(void* const hAmbi, int newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    pData->refl_order = newValue;
}

void ambi_roomsim_setOutputOrder(void* const hAmbi, int newOrder)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    if(newOrder != pData->new_sh_order){
        pData->new_sh_order = newOrder;
        /* FUMA only supports 1st order */
        if(pData->new_sh_order!=SH_ORDER_FIRST && pData->chOrdering == CH_FUMA)
            pData->chOrdering = CH_ACN;
        if(pData->new_sh_order!=SH_ORDER_FIRST && pData->norm == NORM_FUMA)
            pData->norm = NORM_SN3D;

        pData->reinit_room = 1;
    }
}

void ambi_roomsim_setNumSources(void* const hAmbi, int new_nSources)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    pData->new_nSources = SAF_CLAMP(new_nSources, 1, ROOM_SIM_MAX_NUM_SOURCES);
}

void ambi_roomsim_setSourceX(void* const hAmbi, int index, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_SOURCES, "Index exceeds the maximum number of sources permitted");
    pData->src_pos[index][0] = newValue;
}

void ambi_roomsim_setSourceY(void* const hAmbi, int index, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_SOURCES, "Index exceeds the maximum number of sources permitted");
    pData->src_pos[index][1] = newValue;
}

void ambi_roomsim_setSourceZ(void* const hAmbi, int index, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_SOURCES, "Index exceeds the maximum number of sources permitted");
    pData->src_pos[index][2] = newValue;
}

void ambi_roomsim_setNumReceivers(void* const hAmbi, int new_nReceivers)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    pData->new_nReceivers = SAF_CLAMP(new_nReceivers, 1, ROOM_SIM_MAX_NUM_RECEIVERS);
}

void ambi_roomsim_setReceiverX(void* const hAmbi, int index, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_RECEIVERS, "Index exceeds the maximum number of receivers permitted");
    pData->rec_pos[index][0] = newValue;
}

void ambi_roomsim_setReceiverY(void* const hAmbi, int index, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_RECEIVERS, "Index exceeds the maximum number of receivers permitted");
    pData->rec_pos[index][1] = newValue;
}

void ambi_roomsim_setReceiverZ(void* const hAmbi, int index, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_RECEIVERS, "Index exceeds the maximum number of receivers permitted");
    pData->rec_pos[index][2] = newValue;
}

void ambi_roomsim_setRoomDimX(void* const hAmbi, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    if(newValue!=pData->room_dims[0]){
        pData->room_dims[0] = newValue;
    }
}

void ambi_roomsim_setRoomDimY(void* const hAmbi, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    if(newValue!=pData->room_dims[1]){
        pData->room_dims[1] = newValue;
    }
}

void ambi_roomsim_setRoomDimZ(void* const hAmbi, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    if(newValue!=pData->room_dims[2]){
        pData->room_dims[2] = newValue;
    }
}

void ambi_roomsim_setWallAbsCoeff(void* const hAmbi, int xyz_idx, int posNeg_idx, float new_value)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(xyz_idx<4, "xyz_idx indicates each spatial axis, so cannot exceed 4");
    saf_assert(posNeg_idx==0 || posNeg_idx==1, "posNeg_idx is a bool");
    if(new_value!=pData->abs_wall[2*xyz_idx+posNeg_idx]){
        pData->abs_wall[2*xyz_idx+posNeg_idx] = new_value;
    }
}

void ambi_roomsim_setChOrder(void* const hAmbi, int newOrder)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    if((CH_ORDER)newOrder != CH_FUMA || pData->new_sh_order==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void ambi_roomsim_setNormType(void* const hAmbi, int newType)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    if((NORM_TYPES)newType != NORM_FUMA || pData->new_sh_order==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
} 


/* Get Functions */

int ambi_roomsim_getFrameSize(void)
{
    return AMBI_ROOMSIM_FRAME_SIZE;
}

int ambi_roomsim_getEnableIMSflag(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return pData->enableReflections;
}

int ambi_roomsim_getMaxReflectionOrder(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return pData->refl_order;
}

int ambi_roomsim_getOutputOrder(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return pData->new_sh_order;
}

int ambi_roomsim_getNumSources(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return pData->new_nSources;
}

int ambi_roomsim_getMaxNumSources()
{
    return ROOM_SIM_MAX_NUM_SOURCES;
}

int ambi_roomsim_getNSHrequired(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return (pData->new_sh_order+1)*(pData->new_sh_order+1);
}

float ambi_roomsim_getSourceX(void* const hAmbi, int index)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_SOURCES, "Index exceeds the maximum number of sources permitted");
    return pData->src_pos[index][0];
}

float ambi_roomsim_getSourceY(void* const hAmbi, int index)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_SOURCES, "Index exceeds the maximum number of sources permitted");
    return pData->src_pos[index][1];
}

float ambi_roomsim_getSourceZ(void* const hAmbi, int index)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_SOURCES, "Index exceeds the maximum number of sources permitted");
    return pData->src_pos[index][2];
}

int ambi_roomsim_getNumReceivers(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return pData->new_nReceivers;
}

int ambi_roomsim_getMaxNumReceivers()
{
    return ROOM_SIM_MAX_NUM_RECEIVERS;
}

float ambi_roomsim_getReceiverX(void* const hAmbi, int index)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_RECEIVERS, "Index exceeds the maximum number of receivers permitted");
    return pData->rec_pos[index][0];
}

float ambi_roomsim_getReceiverY(void* const hAmbi, int index)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_RECEIVERS, "Index exceeds the maximum number of receivers permitted");
    return pData->rec_pos[index][1];
}

float ambi_roomsim_getReceiverZ(void* const hAmbi, int index)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    saf_assert(index<=ROOM_SIM_MAX_NUM_RECEIVERS, "Index exceeds the maximum number of receivers permitted");
    return pData->rec_pos[index][2];
}

float ambi_roomsim_getRoomDimX(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return pData->room_dims[0];
}

float ambi_roomsim_getRoomDimY(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return pData->room_dims[1];
}

float ambi_roomsim_getRoomDimZ(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return pData->room_dims[2];
}

float ambi_roomsim_getWallAbsCoeff(void* const hAmbi, int xyz_idx, int posNeg_idx)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return pData->abs_wall[2*xyz_idx+posNeg_idx];
}

int ambi_roomsim_getChOrder(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return (int)pData->chOrdering;
}

int ambi_roomsim_getNormType(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return (int)pData->norm;
}

int ambi_roomsim_getProcessingDelay()
{
    return AMBI_ROOMSIM_FRAME_SIZE;
}
