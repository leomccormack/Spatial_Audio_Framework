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
 */

#include "ambi_roomsim.h"
#include "ambi_roomsim_internal.h"

void ambi_roomsim_create
(
    void ** const phAmbi
)
{
    ambi_roomsim_data* pData = (ambi_roomsim_data*)malloc1d(sizeof(ambi_roomsim_data));
    *phAmbi = (void*)pData;
    int i;
    
    printf(SAF_VERSION_LICENSE_STRING);

    pData->order = 1;
    
    /* default user parameters */ 
    pData->nSources = pData->new_nSources;
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_SH_FLAG[i] = 1;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    pData->order = SH_ORDER_FIRST;

    /* sf */
    pData->hIms = NULL;
    pData->sh_order = 3;
    pData->nBands = 1;
    float abs_wall[5][6] = /* Absorption Coefficients per Octave band, and per wall */
      { {0.180791250f, 0.207307300f, 0.134990800f, 0.229002250f, 0.212128400f, 0.241055000f},
        {0.225971250f, 0.259113700f, 0.168725200f, 0.286230250f, 0.265139600f, 0.301295000f},
        {0.258251250f, 0.296128100f, 0.192827600f, 0.327118250f, 0.303014800f, 0.344335000f},
        {0.301331250f, 0.345526500f, 0.224994001f, 0.381686250f, 0.353562000f, 0.401775000f},
        {0.361571250f, 0.414601700f, 0.269973200f, 0.457990250f, 0.424243600f, 0.482095000f} };
    memcpy(pData->abs_wall,abs_wall, 5*6*sizeof(float));
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

    pData->src_sigs = (float**)malloc2d(MAX_NUM_CHANNELS, FRAME_SIZE, sizeof(float));
    pData->rec_sh_outsigs = (float**)malloc2d(MAX_NUM_CHANNELS, FRAME_SIZE, sizeof(float));

    pData->reinit_room = 1;
}

void ambi_roomsim_destroy
(
    void ** const phAmbi
)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(*phAmbi);
    
    if (pData != NULL) {
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
    int i;
    
    pData->fs = (float)sampleRate;
    for(i=1; i<=FRAME_SIZE; i++)
        pData->interpolator[i-1] = (float)i*1.0f/(float)FRAME_SIZE;
    memset(pData->prev_Y, 0, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*sizeof(float));
    memset(pData->prev_inputFrameTD, 0, MAX_NUM_INPUTS*FRAME_SIZE*sizeof(float));
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_SH_FLAG[i] = 1;


    /*  */
    if(pData->reinit_room){
        ims_shoebox_destroy(&(pData->hIms));
        ims_shoebox_create(&(pData->hIms), 10.5f, 7.0f, 3.0f, (float*)pData->abs_wall, 250.0f, pData->nBands, 343.0f, 48e3f);
        pData->sourceIDs[0] = ims_shoebox_addSource(pData->hIms, (float*)pData->src_pos[0], &(pData->src_sigs[0]));
        pData->sourceIDs[1] = ims_shoebox_addSource(pData->hIms, (float*)pData->src_pos[1], &(pData->src_sigs[1]));
        pData->sourceIDs[2] = ims_shoebox_addSource(pData->hIms, (float*)pData->src_pos[2], &pData->src_sigs[2]);
        pData->sourceIDs[3] = ims_shoebox_addSource(pData->hIms, (float*)pData->src_pos[3], &pData->src_sigs[3]);
    //    sourceIDs[1] = ims_shoebox_addSource(pData->hIms, (float*)pData->src2_pos, &src_sigs[1]);
    //    sourceIDs[2] = ims_shoebox_addSource(pData->hIms, (float*)pData->src3_pos, &src_sigs[2]);
    //    sourceIDs[3] = ims_shoebox_addSource(pData->hIms, (float*)pData->src4_pos, &src_sigs[3]);
        pData->receiverIDs[0] = ims_shoebox_addReceiverSH(pData->hIms, pData->sh_order, (float*)pData->rec_pos, &(pData->rec_sh_outsigs));
        pData->reinit_room = 0;
    }
}

void ambi_roomsim_process
(
    void  *  const hAmbi,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples
)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    int i, ch, nSources, nSH;

    /* local copies of user parameters */
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    int order;
    chOrdering = pData->chOrdering;
    norm = pData->norm;
    nSources = pData->nSources;
    order = MIN(pData->order, MAX_SH_ORDER);
    nSH = ORDER2NSH(order);


    float maxTime_s;

    /* Process frame */
    if (nSamples == FRAME_SIZE) {

        nSources = 4;

        /* Load time-domain data */
        for(i=0; i < MIN(nSources,nInputs); i++)
            memcpy(pData->src_sigs[i], inputs[i], FRAME_SIZE * sizeof(float));
        for(; i < nInputs; i++)
            memset(pData->src_sigs[i], 0, FRAME_SIZE * sizeof(float));

        maxTime_s = -0.05f; /* 50ms */

        //ims_shoebox_updateSource(pData->hIms, pData->sourceIDs[0], pData->mov_src_pos);
        //ims_shoebox_updateReceiver(pData->hIms, pData->receiverIDs[0], pData->mov_rec_pos);
        ims_shoebox_computeEchograms(pData->hIms, 6, maxTime_s);
        ims_shoebox_applyEchogramTD(pData->hIms, pData->receiverIDs[0], nSamples, 0);


        for(i=0; i<MIN(nOutputs, 16); i++)
            memcpy(outputs[i], pData->rec_sh_outsigs[i], FRAME_SIZE * sizeof(float));
        for(; i < nOutputs; i++)
            memset(outputs[i], 0, FRAME_SIZE * sizeof(float));
//
//        /* account for output channel order */
//        switch(chOrdering){
//            case CH_ACN: /* already ACN */
//                break;
//            case CH_FUMA:
//                convertHOAChannelConvention((float*)pData->outputFrameTD, order, FRAME_SIZE, HOA_CH_ORDER_ACN, HOA_CH_ORDER_FUMA);
//                break;
//        }
//
//        /* account for normalisation scheme */
//        switch(norm){
//            case NORM_N3D: /* already N3D */
//                break;
//            case NORM_SN3D:
//                convertHOANormConvention((float*)pData->outputFrameTD, order, FRAME_SIZE, HOA_NORM_N3D, HOA_NORM_SN3D);
//                break;
//            case NORM_FUMA:
//                convertHOANormConvention((float*)pData->outputFrameTD, order, FRAME_SIZE, HOA_NORM_N3D, HOA_NORM_FUMA);
//                break;
//        }
//
//        /* Copy to output */
//        for(i = 0; i < MIN(nSH,nOutputs); i++)
//            utility_svvcopy(pData->outputFrameTD[i], FRAME_SIZE, outputs[i]);
//        for(; i < nOutputs; i++)
//            memset(outputs[i], 0, FRAME_SIZE * sizeof(float));

    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, FRAME_SIZE*sizeof(float));
    }
}

/* Set Functions */

int ambi_roomsim_getFrameSize(void)
{
    return FRAME_SIZE;
}

void ambi_roomsim_refreshParams(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    int i;
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_SH_FLAG[i] = 1;
}

void ambi_roomsim_setOutputOrder(void* const hAmbi, int newOrder)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
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

void ambi_roomsim_setNumSources(void* const hAmbi, int new_nSources)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    int i;
    pData->new_nSources = CLAMP(new_nSources, 1, MAX_NUM_INPUTS);
    pData->nSources = pData->new_nSources;
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_SH_FLAG[i] = 1;
}

void ambi_roomsim_setSourceX(void* const hAmbi, int index, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    assert(index<=ROOM_SIM_MAX_NUM_SOURCES);
    pData->src_pos[index][0] = newValue;
    if(pData->hIms!=NULL)
        ims_shoebox_updateSource(pData->hIms, pData->sourceIDs[index], pData->src_pos[index]);
}

void ambi_roomsim_setSourceY(void* const hAmbi, int index, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    assert(index<=ROOM_SIM_MAX_NUM_SOURCES);
    pData->src_pos[index][1] = newValue;
    if(pData->hIms!=NULL)
        ims_shoebox_updateSource(pData->hIms, pData->sourceIDs[index], pData->src_pos[index]);
}

void ambi_roomsim_setSourceZ(void* const hAmbi, int index, float newValue)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    assert(index<=ROOM_SIM_MAX_NUM_SOURCES);
    pData->src_pos[index][2] = newValue;
    if(pData->hIms!=NULL)
        ims_shoebox_updateSource(pData->hIms, pData->sourceIDs[index], pData->src_pos[index]);
}

void ambi_roomsim_setInputConfigPreset(void* const hAmbi, int newPresetID)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    int ch; 
    pData->nSources = pData->new_nSources;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_SH_FLAG[ch] = 1;
}

void ambi_roomsim_setChOrder(void* const hAmbi, int newOrder)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    if((CH_ORDER)newOrder != CH_FUMA || pData->order==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void ambi_roomsim_setNormType(void* const hAmbi, int newType)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    if((NORM_TYPES)newType != NORM_FUMA || pData->order==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
} 


/* Get Functions */

int ambi_roomsim_getOutputOrder(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return (int)pData->order;
}

float ambi_roomsim_getSourceAzi_deg(void* const hAmbi, int index)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return 0.0f;// pData->src_dirs_deg[index][0];
}

float ambi_roomsim_getSourceElev_deg(void* const hAmbi, int index)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return 0.0f;// pData->src_dirs_deg[index][1];
}

int ambi_roomsim_getNumSources(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return pData->new_nSources;
}

int ambi_roomsim_getMaxNumSources()
{
    return MAX_NUM_INPUTS;
}

int ambi_roomsim_getNSHrequired(void* const hAmbi)
{
    ambi_roomsim_data *pData = (ambi_roomsim_data*)(hAmbi);
    return (pData->order+1)*(pData->order+1);
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
    return FRAME_SIZE;
}
