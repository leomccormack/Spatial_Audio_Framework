/*
 Copyright 2016-2018 Leo McCormack
 
 Permission to use, copy, modify, and/or distribute this software for any purpose with or
 without fee is hereby granted, provided that the above copyright notice and this permission
 notice appear in all copies.
 
 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 OR PERFORMANCE OF THIS SOFTWARE.
*/
/*
 * Filename:
 *     ambi_enc.c
 * Description:
 *     A simple, but flexible, Ambisonic encoder (aka: Ambisonic Panner).
 * Dependencies:
 *     saf_utilities, saf_sh
 * Author, date created:
 *     Leo McCormack, 07.10.2016
 */

#include "ambi_enc.h"
#include "ambi_enc_internal.h"

void ambi_enc_create
(
    void ** const phAmbi
)
{
    ambi_enc_data* pData = (ambi_enc_data*)malloc(sizeof(ambi_enc_data));
    if (pData == NULL) { return;/*error*/ }
    *phAmbi = (void*)pData;
    int i;
    pData->order = 1;
    
    /* user parameters */
    ambi_enc_loadPreset(PRESET_DEFAULT, pData->src_dirs_deg, &(pData->new_nSources)); /*check setStateInformation if you change default preset*/
    pData->nSources = pData->new_nSources;
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_SH_FLAG[i] = 1;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_N3D;
    pData->outputOrderPreset = OUTPUT_ORDER_FIRST;
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
    pData->fs = (float)sampleRate;
}

void ambi_enc_process
(
    void  *  const hAmbi,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples,
    int            isPlaying
)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    int i, j, ch, n, nSources, nSH;
    int o[MAX_ORDER+2];
    float src_dirs[MAX_NUM_INPUTS][2], scale;
    float* Y_src;
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    int order;
    
    if ( (nSamples == FRAME_SIZE) && (isPlaying == 1) ) {
        /* prep */
        for(n=0; n<MAX_ORDER+2; n++){  o[n] = n*n;  }
        chOrdering = pData->chOrdering;
        norm = pData->norm;
        nSources = pData->nSources;
        memcpy(src_dirs, pData->src_dirs_deg, MAX_NUM_INPUTS*2*sizeof(float));
        order = MIN(pData->order, MAX_ORDER);
        nSH = (order+1)*(order+1);
        Y_src = malloc(nSH*sizeof(float));
        
        /* Load time-domain data */
        for(i=0; i < MIN(nSources,nInputs); i++)
            memcpy(pData->inputFrameTD[i], inputs[i], FRAME_SIZE * sizeof(float));
        for(; i<MAX_NUM_INPUTS; i++)
            memset(pData->inputFrameTD[i], 0, FRAME_SIZE * sizeof(float));
        
        /* recalulate SHs */
        for(i=0; i<nSources; i++){
            if(pData->recalc_SH_FLAG[i]){
                getSHreal(order, pData->src_dirs_deg[i][0]*M_PI/180.0f, M_PI/2.0f - pData->src_dirs_deg[i][1]*M_PI/180.0f, Y_src);
                for(j=0; j<nSH; j++)
                    pData->Y[j][i] = sqrtf(4.0f*M_PI)*Y_src[j];
                for(; j<MAX_NUM_SH_SIGNALS; j++)
                    pData->Y[j][i] = 0.0f;
                pData->recalc_SH_FLAG[i] = 0;
            }
        }
        
        /* spatially encode the input signals into spherical harmonic signals */
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, FRAME_SIZE, nSources, 1.0,
                    (float*)pData->Y, MAX_NUM_INPUTS,
                    (float*)pData->inputFrameTD, FRAME_SIZE, 0.0,
                    (float*)pData->outputFrameTD, FRAME_SIZE);
        
        /* scale by 1/sqrt(nSources) */
        scale = 1.0f/sqrt(nSources);
        utility_svsmul((float*)pData->outputFrameTD, &scale, nSH*FRAME_SIZE, NULL);
        
        /* norm scheme */
        switch(norm){
            case NORM_N3D: /* already N3D */
                break;
            case NORM_SN3D:
                for (n = 0; n<order+1; n++)
                    for (ch = o[n]; ch<o[n+1]; ch++)
                        for(i = 0; i<FRAME_SIZE; i++)
                            pData->outputFrameTD[ch][i] /= sqrtf(2.0f*(float)n+1.0f);
                break;
        }
        
        /* save SH signals to output buffer */
        for(i = 0; i < MIN(nSH,nOutputs); i++)
            memcpy(outputs[i], pData->outputFrameTD[i], FRAME_SIZE * sizeof(float));
        for(; i < nOutputs; i++)
            memset(outputs[i], 0, FRAME_SIZE * sizeof(float));
        free((void*)Y_src);
    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, FRAME_SIZE*sizeof(float));
    }
}

/* Set Functions */

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
    if(newOrder != pData->outputOrderPreset){
        pData->outputOrderPreset = (OUTPUT_ORDERS)newOrder;
        switch(pData->outputOrderPreset){
            case OUTPUT_OMNI: pData->order = 0; break;
            case OUTPUT_ORDER_FIRST:   pData->order = 1; break;
            case OUTPUT_ORDER_SECOND:  pData->order = 2; break;
            case OUTPUT_ORDER_THIRD:   pData->order = 3; break;
            case OUTPUT_ORDER_FOURTH:  pData->order = 4; break;
            case OUTPUT_ORDER_FIFTH:   pData->order = 5; break;
            case OUTPUT_ORDER_SIXTH:   pData->order = 6; break;
            case OUTPUT_ORDER_SEVENTH: pData->order = 7; break;
        }
        for(i=0; i<MAX_NUM_INPUTS; i++)
            pData->recalc_SH_FLAG[i] = 1;
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
    pData->new_nSources = new_nSources > MAX_NUM_INPUTS ? MAX_NUM_INPUTS : new_nSources;
    pData->nSources = pData->new_nSources;
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_SH_FLAG[i] = 1;
}

void ambi_enc_setInputConfigPreset(void* const hAmbi, int newPresetID)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    int ch;
    ambi_enc_loadPreset(newPresetID, pData->src_dirs_deg, &(pData->new_nSources));
    pData->nSources = pData->new_nSources;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_SH_FLAG[ch] = 1;
}

void ambi_enc_setChOrder(void* const hAmbi, int newOrder)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    pData->chOrdering = (CH_ORDER)newOrder;
}

void ambi_enc_setNormType(void* const hAmbi, int newType)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    pData->norm = (NORM_TYPES)newType;
}

/* Get Functions */

int ambi_enc_getOutputOrder(void* const hAmbi)
{
    ambi_enc_data *pData = (ambi_enc_data*)(hAmbi);
    return (int)pData->outputOrderPreset;
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

