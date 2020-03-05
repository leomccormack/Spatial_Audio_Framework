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
 * @file beamformer.c
 * @brief Generates beamformers/virtual microphones in arbitrary directions
 *        with several different beam pattern to choose from
 *
 * @author Leo McCormack
 * @date 17.05.2019
 */
 
#include "beamformer_internal.h"

void beamformer_create
(
    void ** const phBeam
)
{
    beamformer_data* pData = (beamformer_data*)malloc1d(sizeof(beamformer_data));
    *phBeam = (void*)pData;
    int i, ch;

    /* default user parameters */
    pData->beamOrder = 1;
    pData->nSH = pData->new_nSH = (pData->beamOrder + 1)*(pData->beamOrder + 1);
    for(i=0; i<MAX_NUM_BEAMS; i++){
        pData->beam_dirs_deg[i][0] = default_LScoords64_rad[i][0]*180.0f/M_PI;
        pData->beam_dirs_deg[i][1] = (default_LScoords64_rad[i][1] - M_PI/2.0f) < -M_PI/2.0f ?
        (M_PI/2.0f + default_LScoords64_rad[i][1]) :  (default_LScoords64_rad[i][1] - M_PI/2.0f);
        pData->beam_dirs_deg[i][1] *= 180.0f/M_PI;
    }
    pData->nBeams = pData->new_nBeams = 1;
    pData->beamType = BEAM_TYPE_HYPERCARDIOID;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    
    /* internal parameters */
    pData->new_nSH = pData->new_nSH = (pData->beamOrder+1)*(pData->beamOrder+1);
    
    /* flags */
    pData->reInitTFT = 1;
    for(ch=0; ch<MAX_NUM_BEAMS; ch++)
        pData->recalc_beamWeights[ch] = 1;
}

void beamformer_destroy
(
    void ** const phBeam
)
{
    beamformer_data *pData = (beamformer_data*)(*phBeam);
    
    if (pData != NULL) {
        
        free(pData);
        pData = NULL;
    }
}

void beamformer_init
(
    void * const hBeam,
    int          sampleRate
)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int i, ch;
    
    /* define frequency vector */
    pData->fs = sampleRate;
   
    /* defaults */
    memset(pData->beamWeights, 0, MAX_NUM_BEAMS*MAX_NUM_SH_SIGNALS*sizeof(float));
    memset(pData->prev_beamWeights, 0, MAX_NUM_BEAMS*MAX_NUM_SH_SIGNALS*sizeof(float));
    memset(pData->prev_SHFrameTD, 0, MAX_NUM_SH_SIGNALS*FRAME_SIZE*sizeof(float));
    for(ch=0; ch<MAX_NUM_BEAMS; ch++)
        pData->recalc_beamWeights[ch] = 1;
    for(i=1; i<=FRAME_SIZE; i++)
        pData->interpolator[i-1] = (float)i*1.0f/(float)FRAME_SIZE;
}

void beamformer_process
(
    void  *  const hBeam,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples
)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int n, ch, i, j, bi;
    int o[MAX_SH_ORDER+2];

    /* local copies of user parameters */
    int nBeams, beamOrder;
    BEAMFORMER_NORM_TYPES norm;
    BEAMFORMER_CH_ORDER chOrdering;
    
    /* reinitialise if needed */
    if(pData->reInitTFT==1){
        pData->reInitTFT = 2;
        beamformer_initTFT(hBeam);
        pData->reInitTFT = 0;
    }

    /* decode audio to loudspeakers or headphones */
    if( (nSamples == FRAME_SIZE) && (pData->reInitTFT==0) ) {
        /* copy user parameters to local variables */
        for(n=0; n<MAX_SH_ORDER+2; n++){  o[n] = n*n;  }
        beamOrder = pData->beamOrder;
        nBeams = pData->nBeams;
        norm = pData->norm;
        chOrdering = pData->chOrdering;
        
        /* Load time-domain data */
        switch(chOrdering){
            case CH_ACN:
                for(i=0; i < MIN(pData->nSH, nInputs); i++)
                    utility_svvcopy(inputs[i], FRAME_SIZE, pData->SHFrameTD[i]);
                for(; i<pData->nSH; i++)
                    memset(pData->SHFrameTD[i], 0, FRAME_SIZE * sizeof(float)); /* fill remaining channels with zeros */
                break;
            case CH_FUMA:   /* only for first-order, convert to ACN */
                if(nInputs>=4){
                    utility_svvcopy(inputs[0], FRAME_SIZE, pData->SHFrameTD[0]);
                    utility_svvcopy(inputs[1], FRAME_SIZE, pData->SHFrameTD[3]);
                    utility_svvcopy(inputs[2], FRAME_SIZE, pData->SHFrameTD[1]);
                    utility_svvcopy(inputs[3], FRAME_SIZE, pData->SHFrameTD[2]);
                    for(i=4; i<pData->nSH; i++)
                        memset(pData->SHFrameTD[i], 0, FRAME_SIZE * sizeof(float)); /* fill remaining channels with zeros */
                }
                else
                    for(i=0; i<pData->nSH; i++)
                        memset(pData->SHFrameTD[i], 0, FRAME_SIZE * sizeof(float));
                break;
        }
        
        /* account for input normalisation scheme */
        switch(norm){
            case NORM_N3D:  /* already in N3D, do nothing */
                break;
            case NORM_SN3D: /* convert to N3D */
                for (n = 0; n<beamOrder+1; n++)
                    for (ch = o[n]; ch<o[n+1]; ch++)
                        for(i = 0; i<FRAME_SIZE; i++)
                            pData->SHFrameTD[ch][i] *= sqrtf(2.0f*(float)n+1.0f);
                break;
            case NORM_FUMA: /* only for first-order, convert to N3D */
                for(i = 0; i<FRAME_SIZE; i++)
                    pData->SHFrameTD[0][i] *= sqrtf(2.0f);
                for (ch = 1; ch<4; ch++)
                    for(i = 0; i<FRAME_SIZE; i++)
                        pData->SHFrameTD[ch][i] *= sqrtf(3.0f);
                break;
        }
        
        /* Main processing: */
        float* c_n;
        c_n = malloc1d((beamOrder+1)*sizeof(float));
            
        /* calculate beamforming coeffients */
        for(bi=0; bi<nBeams; bi++){
            if(pData->recalc_beamWeights[bi]){
                memset(pData->beamWeights[bi], 0, MAX_NUM_SH_SIGNALS*sizeof(float));
                switch(pData->beamType){
                    case BEAM_TYPE_CARDIOID: beamWeightsCardioid2Spherical(beamOrder, c_n); break;
                    case BEAM_TYPE_HYPERCARDIOID: beamWeightsHypercardioid2Spherical(beamOrder, c_n); break;
                    case BEAM_TYPE_MAX_EV: beamWeightsMaxEV(beamOrder, c_n); break;
                }
                rotateAxisCoeffsReal(beamOrder, c_n, M_PI/2.0f - pData->beam_dirs_deg[bi][1]*M_PI/180.0f,
                                        pData->beam_dirs_deg[bi][0]*M_PI/180.0f, (float*)pData->beamWeights[bi]);
                    
                pData->recalc_beamWeights[bi] = 0;
            }
            else
                memcpy(pData->beamWeights[bi], pData->prev_beamWeights[bi], pData->nSH*sizeof(float));
        }
        free(c_n);
            
        /* apply beam weights */
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nBeams, FRAME_SIZE, pData->nSH, 1.0f,
                    (const float*)pData->prev_beamWeights, MAX_NUM_SH_SIGNALS,
                    (const float*)pData->prev_SHFrameTD, FRAME_SIZE, 0.0f,
                    (float*)pData->tempFrame, FRAME_SIZE);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nBeams, FRAME_SIZE, pData->nSH, 1.0f,
                    (const float*)pData->beamWeights, MAX_NUM_SH_SIGNALS,
                    (const float*)pData->prev_SHFrameTD, FRAME_SIZE, 0.0f,
                    (float*)pData->outputFrameTD, FRAME_SIZE);
            
        for (i=0; i <nBeams; i++)
            for(j=0; j<FRAME_SIZE; j++)
                pData->outputFrameTD[i][j] =  pData->interpolator[j] * pData->outputFrameTD[i][j] + (1.0f-pData->interpolator[j]) * pData->tempFrame[i][j];
            
        /* for next frame */
        utility_svvcopy((const float*)pData->SHFrameTD, pData->nSH*FRAME_SIZE, (float*)pData->prev_SHFrameTD);
        utility_svvcopy((const float*)pData->beamWeights, MAX_NUM_BEAMS*MAX_NUM_SH_SIGNALS, (float*)pData->prev_beamWeights);
            
        /* copy to output buffer */
        for(ch = 0; ch < MIN(nBeams, nOutputs); ch++)
            utility_svvcopy(pData->outputFrameTD[ch], FRAME_SIZE, outputs[ch]);
        for (; ch < nOutputs; ch++)
            memset(outputs[ch], 0, FRAME_SIZE*sizeof(float));
    }
    else
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch], 0, FRAME_SIZE*sizeof(float));
}


/* Set Functions */

void beamformer_refreshSettings(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    pData->reInitTFT = 1;
}

void beamformer_setBeamOrder(void  * const hBeam, int newValue)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int ch;
    pData->beamOrder = MIN(MAX(newValue,1), MAX_SH_ORDER);
    pData->new_nSH = (pData->beamOrder+1)*(pData->beamOrder+1);
    pData->reInitTFT = 1;
    for(ch=0; ch<MAX_NUM_BEAMS; ch++)
        pData->recalc_beamWeights[ch] = 1;
    /* FUMA only supports 1st order */
    if(pData->beamOrder!=BEAM_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->beamOrder!=BEAM_ORDER_FIRST && pData->norm == NORM_FUMA)
        pData->norm = NORM_SN3D;
}

void beamformer_setBeamAzi_deg(void* const hBeam, int index, float newAzi_deg)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = MAX(newAzi_deg, -180.0f);
    newAzi_deg = MIN(newAzi_deg, 180.0f);
    pData->beam_dirs_deg[index][0] = newAzi_deg;
    pData->recalc_beamWeights[index] = 1;
}

void beamformer_setBeamElev_deg(void* const hBeam, int index, float newElev_deg)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    newElev_deg = MAX(newElev_deg, -90.0f);
    newElev_deg = MIN(newElev_deg, 90.0f);
    pData->beam_dirs_deg[index][1] = newElev_deg;
    pData->recalc_beamWeights[index] = 1;
}

void beamformer_setNumBeams(void* const hBeam, int new_nBeams)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int ch;
    pData->new_nBeams = CLAMP(new_nBeams, 1, MAX_NUM_BEAMS);
    if(pData->nBeams != pData->new_nBeams){
        pData->reInitTFT = 1;
        for(ch=0; ch<MAX_NUM_BEAMS; ch++)
            pData->recalc_beamWeights[ch] = 1;
    }
}

void beamformer_setChOrder(void* const hBeam, int newOrder)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    if((BEAMFORMER_CH_ORDER)newOrder != CH_FUMA || pData->beamOrder==BEAM_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (BEAMFORMER_CH_ORDER)newOrder;
}

void beamformer_setNormType(void* const hBeam, int newType)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    if((BEAMFORMER_NORM_TYPES)newType != NORM_FUMA || pData->beamOrder==BEAM_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (BEAMFORMER_NORM_TYPES)newType;
}

void beamformer_setBeamType(void* const hBeam, int newID)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int ch;
    pData->beamType = newID;
    for(ch=0; ch<MAX_NUM_BEAMS; ch++)
        pData->recalc_beamWeights[ch] = 1;
}

/* Get Functions */

int beamformer_getBeamOrder(void  * const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->beamOrder;
}

int beamformer_getNumberOfBands(void)
{
    return HYBRID_BANDS;
}

float beamformer_getBeamAzi_deg(void* const hBeam, int index)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->beam_dirs_deg[index][0];
}

float beamformer_getBeamElev_deg(void* const hBeam, int index)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->beam_dirs_deg[index][1];
}

int beamformer_getNumBeams(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->new_nBeams;
}

int beamformer_getMaxNumBeams()
{
    return MAX_NUM_BEAMS;
}

int  beamformer_getNSHrequired(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->nSH;
}

int beamformer_getChOrder(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return (int)pData->chOrdering;
}

int beamformer_getNormType(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return (int)pData->norm;
}

int beamformer_getBeamType(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->beamType;
}



