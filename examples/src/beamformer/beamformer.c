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
 *        with several different beam patterns to choose from
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
    for(i=0; i<MAX_NUM_BEAMS; i++){
        pData->beam_dirs_deg[i][0] = default_LScoords64_rad[i][0]*180.0f/M_PI;
        pData->beam_dirs_deg[i][1] = (default_LScoords64_rad[i][1] - M_PI/2.0f) < -M_PI/2.0f ?
        (M_PI/2.0f + default_LScoords64_rad[i][1]) :  (default_LScoords64_rad[i][1] - M_PI/2.0f);
        pData->beam_dirs_deg[i][1] *= 180.0f/M_PI;
    }
    pData->nBeams = 1;
    pData->beamType = STATIC_BEAM_TYPE_HYPERCARDIOID;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;

    /* flags */
    for(ch=0; ch<MAX_NUM_BEAMS; ch++)
        pData->recalc_beamWeights[ch] = 1;

    /* set FIFO buffers */
    pData->FIFO_idx = 0;
    memset(pData->inFIFO, 0, MAX_NUM_SH_SIGNALS*FRAME_SIZE*sizeof(float));
    memset(pData->outFIFO, 0, MAX_NUM_BEAMS*FRAME_SIZE*sizeof(float));
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
    int s, n, ch, i, j, bi, nSH;
    int o[MAX_SH_ORDER+2];

    /* local copies of user parameters */
    int nBeams, beamOrder;
    NORM_TYPES norm;
    CH_ORDER chOrdering;
    for(n=0; n<MAX_SH_ORDER+2; n++){  o[n] = n*n;  }
    beamOrder = pData->beamOrder;
    nSH = ORDER2NSH(beamOrder);
    nBeams = pData->nBeams;
    norm = pData->norm;
    chOrdering = pData->chOrdering;
     
    /* Loop over all samples */
    for(s=0; s<nSamples; s++){
        /* Load input signals into inFIFO buffer */
        for(ch=0; ch<MIN(nInputs,MAX_NUM_SH_SIGNALS); ch++)
            pData->inFIFO[ch][pData->FIFO_idx] = inputs[ch][s];
        for(; ch<nSH; ch++) /* Zero any channels that were not given */
            pData->inFIFO[ch][pData->FIFO_idx] = 0.0f;

        /* Pull output signals from outFIFO buffer */
        for(ch=0; ch<MIN(nOutputs, MAX_NUM_BEAMS); ch++)
            outputs[ch][s] = pData->outFIFO[ch][pData->FIFO_idx];
        for(; ch<nOutputs; ch++) /* Zero any extra channels */
            outputs[ch][s] = 0.0f;

        /* Increment buffer index */
        pData->FIFO_idx++;

        /* Process frame if inFIFO is full and codec is ready for it */
        if (pData->FIFO_idx >= FRAME_SIZE) {
            pData->FIFO_idx = 0;

            /* Load time-domain data */
            switch(chOrdering){
              case CH_ACN:
                  convertHOAChannelConvention((float*)pData->inFIFO, beamOrder, FRAME_SIZE, HOA_CH_ORDER_ACN, HOA_CH_ORDER_ACN, (float*)pData->SHFrameTD);
                  break;
              case CH_FUMA:
                  convertHOAChannelConvention((float*)pData->inFIFO, beamOrder, FRAME_SIZE, HOA_CH_ORDER_FUMA, HOA_CH_ORDER_ACN, (float*)pData->SHFrameTD);
                  break;
            }

            /* account for input normalisation scheme */
            switch(norm){
              case NORM_N3D:  /* already in N3D, do nothing */
                  break;
              case NORM_SN3D: /* convert to N3D */
                  convertHOANormConvention((float*)pData->SHFrameTD, beamOrder, FRAME_SIZE, HOA_NORM_SN3D, HOA_NORM_N3D);
                  break;
              case NORM_FUMA: /* only for first-order, convert to N3D */
                  convertHOANormConvention((float*)pData->SHFrameTD, beamOrder, FRAME_SIZE, HOA_NORM_FUMA, HOA_NORM_N3D);
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
                        case STATIC_BEAM_TYPE_CARDIOID: beamWeightsCardioid2Spherical(beamOrder, c_n); break;
                        case STATIC_BEAM_TYPE_HYPERCARDIOID: beamWeightsHypercardioid2Spherical(beamOrder, c_n); break;
                        case STATIC_BEAM_TYPE_MAX_EV: beamWeightsMaxEV(beamOrder, c_n); break;
                    }
                    rotateAxisCoeffsReal(beamOrder, c_n, M_PI/2.0f - pData->beam_dirs_deg[bi][1]*M_PI/180.0f,
                                            pData->beam_dirs_deg[bi][0]*M_PI/180.0f, (float*)pData->beamWeights[bi]);

                    pData->recalc_beamWeights[bi] = 0;
                }
                else
                    memcpy(pData->beamWeights[bi], pData->prev_beamWeights[bi], nSH*sizeof(float));
            }
            free(c_n);

            /* apply beam weights */
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nBeams, FRAME_SIZE, nSH, 1.0f,
                        (const float*)pData->prev_beamWeights, MAX_NUM_SH_SIGNALS,
                        (const float*)pData->prev_SHFrameTD, FRAME_SIZE, 0.0f,
                        (float*)pData->tempFrame, FRAME_SIZE);
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nBeams, FRAME_SIZE, nSH, 1.0f,
                        (const float*)pData->beamWeights, MAX_NUM_SH_SIGNALS,
                        (const float*)pData->prev_SHFrameTD, FRAME_SIZE, 0.0f,
                        (float*)pData->outputFrameTD, FRAME_SIZE);

            for (i=0; i <nBeams; i++)
                for(j=0; j<FRAME_SIZE; j++)
                    pData->outputFrameTD[i][j] =  pData->interpolator[j] * pData->outputFrameTD[i][j] + (1.0f-pData->interpolator[j]) * pData->tempFrame[i][j];

            /* for next frame */
            utility_svvcopy((const float*)pData->SHFrameTD, nSH*FRAME_SIZE, (float*)pData->prev_SHFrameTD);
            utility_svvcopy((const float*)pData->beamWeights, MAX_NUM_BEAMS*MAX_NUM_SH_SIGNALS, (float*)pData->prev_beamWeights);

            /* copy to output buffer */
            for(ch = 0; ch < MIN(nBeams, nOutputs); ch++)
                utility_svvcopy(pData->outputFrameTD[ch], FRAME_SIZE, pData->outFIFO[ch]);
        }
        else if(pData->FIFO_idx >= FRAME_SIZE){
            /* clear outFIFO if codec was not ready */
            pData->FIFO_idx = 0;
            memset(pData->outFIFO, 0, MAX_NUM_BEAMS*FRAME_SIZE*sizeof(float));
        }
    }
}


/* Set Functions */

void beamformer_refreshSettings(void* const hBeam)
{
}

void beamformer_setBeamOrder(void  * const hBeam, int newValue)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int ch;
    pData->beamOrder = MIN(MAX(newValue,1), MAX_SH_ORDER);
    for(ch=0; ch<MAX_NUM_BEAMS; ch++)
        pData->recalc_beamWeights[ch] = 1;
    /* FUMA only supports 1st order */
    if(pData->beamOrder!=SH_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->beamOrder!=SH_ORDER_FIRST && pData->norm == NORM_FUMA)
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
    if(pData->nBeams != new_nBeams){
        pData->nBeams = new_nBeams;
        for(ch=0; ch<MAX_NUM_BEAMS; ch++)
            pData->recalc_beamWeights[ch] = 1;
    }
}

void beamformer_setChOrder(void* const hBeam, int newOrder)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    if((CH_ORDER)newOrder != CH_FUMA || pData->beamOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void beamformer_setNormType(void* const hBeam, int newType)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    if((NORM_TYPES)newType != NORM_FUMA || pData->beamOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
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
    return pData->nBeams;
}

int beamformer_getMaxNumBeams()
{
    return MAX_NUM_BEAMS;
}

int  beamformer_getNSHrequired(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return ORDER2NSH(pData->beamOrder);
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

int beamformer_getProcessingDelay()
{
    return FRAME_SIZE;
}


