/*
 * Copyright 2017-2018 Leo McCormack
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
 * @file rotator.c
 * @brief  A simple spherical harmonic domain rotator, based on the recursive
 *         approach detailed in [1].
 *
 * @see [1] Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real
 *          Spherical Harmonics. Direct Determination by Recursion Page:
 *          Additions and Corrections. Journal of Physical Chemistry A, 102(45),
 *          9099?9100.
 *
 * @author Leo McCormack
 * @date 02.11.2017
 */

#include "rotator.h"
#include "rotator_internal.h"

void rotator_create
(
    void ** const phRot
)
{
    rotator_data* pData = (rotator_data*)malloc1d(sizeof(rotator_data));
    *phRot = (void*)pData;
    
    pData->recalc_M_rotFLAG = 1;
  
    /* Default user parameters */
    pData->yaw = 0.0f;
    pData->pitch = 0.0f;
    pData->roll = 0.0f;
    pData->bFlipYaw = 0;
    pData->bFlipPitch = 0;
    pData->bFlipRoll = 0;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    pData->useRollPitchYawFlag = 0;
    rotator_setOrder(*phRot, INPUT_ORDER_FIRST);

    /* set FIFO buffers */
    pData->FIFO_idx = 0;
    memset(pData->inFIFO, 0, MAX_NUM_SH_SIGNALS*FRAME_SIZE*sizeof(float));
    memset(pData->outFIFO, 0, MAX_NUM_SH_SIGNALS*FRAME_SIZE*sizeof(float));
}

void rotator_destroy
(
    void ** const phRot
)
{
    rotator_data *pData = (rotator_data*)(*phRot);

    if (pData != NULL) {
        free(pData);
        pData = NULL;
    }
}

void rotator_init
(
    void * const hRot,
    int          sampleRate
)
{
    rotator_data *pData = (rotator_data*)(hRot);
    int i;
    
    /* starting values */
    for(i=1; i<=FRAME_SIZE; i++)
        pData->interpolator[i-1] = (float)i*1.0f/(float)FRAME_SIZE;
    memset(pData->M_rot, 0, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*sizeof(float));
    memset(pData->prev_M_rot, 0, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*sizeof(float));
    memset(pData->prev_inputFrameTD, 0, MAX_NUM_SH_SIGNALS*FRAME_SIZE*sizeof(float));
    pData->recalc_M_rotFLAG = 1;
}

void rotator_process
(
    void  *  const hRot,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples
)
{
    rotator_data *pData = (rotator_data*)(hRot);
    int s, ch, i, j, order, nSH; 
    float Rxyz[3][3];
    float* M_rot_tmp;
    ROTATOR_CH_ORDER chOrdering;
    ROTATOR_NORM_TYPES norm;

    /* locals */
    chOrdering = pData->chOrdering;
    norm = pData->norm;
    order = (int)pData->inputOrder;
    nSH = ORDER2NSH(order);


   /* Loop over all samples */
    for(s=0; s<nSamples; s++){
        /* Load input signals into inFIFO buffer */
        for(ch=0; ch<MIN(nInputs,nSH); ch++)
            pData->inFIFO[ch][pData->FIFO_idx] = inputs[ch][s];
        for(; ch<nSH; ch++) /* Zero any channels that were not given */
            pData->inFIFO[ch][pData->FIFO_idx] = 0.0f;

        /* Pull output signals from outFIFO buffer */
        for(ch=0; ch<MIN(nOutputs, nSH); ch++)
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
                    convertHOAChannelConvention((float*)pData->inFIFO, order, FRAME_SIZE, HOA_CH_ORDER_ACN, HOA_CH_ORDER_ACN, (float*)pData->inputFrameTD);
                    break;
                case CH_FUMA:
                    convertHOAChannelConvention((float*)pData->inFIFO, order, FRAME_SIZE, HOA_CH_ORDER_FUMA, HOA_CH_ORDER_ACN, (float*)pData->inputFrameTD);
                    break;
            }

#if 0 /* actually doesn't matter, since only components of the same order are
       * used to rotate a given order of component; i.e, dipoles are used to
       * rotate dipoles, quadrapoles-qaudrapoles etc.. so this scaling doesn't
       * matter */
            /* account for input normalisation scheme */
            switch(norm){
                case NORM_N3D:  /* already in N3D, do nothing */
                    break;
                case NORM_SN3D: /* convert to N3D */
                    convertHOANormConvention((float*)pData->inputFrameTD, order, FRAME_SIZE, HOA_NORM_SN3D, HOA_NORM_N3D);
                    break;
                case NORM_FUMA: /* only for first-order, convert to N3D */
                    convertHOANormConvention((float*)pData->inputFrameTD, order, FRAME_SIZE, HOA_NORM_FUMA, HOA_NORM_N3D);
                    break;
            }
#endif

            if (order>0){
                /* calculate rotation matrix */
                if(pData->recalc_M_rotFLAG){
                    memset(pData->M_rot, 0, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*sizeof(float));
                    M_rot_tmp = malloc1d(nSH*nSH*sizeof(float));
                    yawPitchRoll2Rzyx (pData->yaw, pData->pitch, pData->roll, pData->useRollPitchYawFlag, Rxyz);
                    getSHrotMtxReal(Rxyz, M_rot_tmp, order);
                    for(i=0; i<nSH; i++)
                        for(j=0; j<nSH; j++)
                            pData->M_rot[i][j] = M_rot_tmp[i*nSH+j];
                    free(M_rot_tmp);
                    pData->recalc_M_rotFLAG = 0;
                }
                else
                    utility_svvcopy((const float*)pData->prev_M_rot, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS, (float*)pData->M_rot);

                /* apply rotation */
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, FRAME_SIZE, nSH, 1.0f,
                            (float*)(pData->prev_M_rot), MAX_NUM_SH_SIGNALS,
                            (float*)pData->prev_inputFrameTD, FRAME_SIZE, 0.0f,
                            (float*)pData->tempFrame, FRAME_SIZE);
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, FRAME_SIZE, nSH, 1.0f,
                            (float*)(pData->M_rot), MAX_NUM_SH_SIGNALS,
                            (float*)pData->prev_inputFrameTD, FRAME_SIZE, 0.0f,
                            (float*)pData->outputFrameTD, FRAME_SIZE);
                for (i=0; i < nSH; i++)
                    for(j=0; j<FRAME_SIZE; j++)
                        pData->outputFrameTD[i][j] = pData->interpolator[j] * pData->outputFrameTD[i][j] + (1.0f-pData->interpolator[j]) * pData->tempFrame[i][j];

                /* for next frame */
                utility_svvcopy((const float*)pData->inputFrameTD, nSH*FRAME_SIZE, (float*)pData->prev_inputFrameTD);
                utility_svvcopy((const float*)pData->M_rot, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS, (float*)pData->prev_M_rot);
            }
            else
                utility_svvcopy((const float*)pData->inputFrameTD[0], FRAME_SIZE, (float*)pData->outputFrameTD[0]);
        
            /* account for norm scheme */
            switch(norm){
                case NORM_N3D:
                    /* again, actually doesn't matter */
                    break;
                case NORM_SN3D:
                    /* again, actually doesn't matter */
                    break;
                case NORM_FUMA:
                    /* again, actually doesn't matter */
                    break;
            }

            /* copy rotated signals to output buffer */
            switch(chOrdering){
                case CH_ACN:
                    convertHOAChannelConvention((float*)pData->outputFrameTD, order, FRAME_SIZE, HOA_CH_ORDER_ACN, HOA_CH_ORDER_ACN, (float*)pData->outFIFO);
                    break;
                case CH_FUMA:
                    convertHOAChannelConvention((float*)pData->outputFrameTD, order, FRAME_SIZE, HOA_CH_ORDER_ACN,
                                                HOA_CH_ORDER_FUMA, (float*)pData->outFIFO);
                    break;
            }
        }
        else if(pData->FIFO_idx >= FRAME_SIZE){
            /* clear outFIFO if codec was not ready */
            pData->FIFO_idx = 0;
            memset(pData->outFIFO, 0, MAX_NUM_SH_SIGNALS*FRAME_SIZE*sizeof(float));
        }
    }
}

void rotator_setYaw(void  * const hRot, float newYaw)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->yaw = pData->bFlipYaw == 1 ? -DEG2RAD(newYaw) : DEG2RAD(newYaw);
    pData->recalc_M_rotFLAG = 1;
}

void rotator_setPitch(void* const hRot, float newPitch)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->pitch = pData->bFlipPitch == 1 ? -DEG2RAD(newPitch) : DEG2RAD(newPitch);
    pData->recalc_M_rotFLAG = 1;
}

void rotator_setRoll(void* const hRot, float newRoll)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->roll = pData->bFlipRoll == 1 ? -DEG2RAD(newRoll) : DEG2RAD(newRoll);
    pData->recalc_M_rotFLAG = 1;
}

void rotator_setFlipYaw(void* const hRot, int newState)
{
    rotator_data *pData = (rotator_data*)(hRot);
    if(newState !=pData->bFlipYaw ){
        pData->bFlipYaw = newState;
        rotator_setYaw(hRot, -rotator_getYaw(hRot));
    }
}

void rotator_setFlipPitch(void* const hRot, int newState)
{
    rotator_data *pData = (rotator_data*)(hRot);
    if(newState !=pData->bFlipPitch ){
        pData->bFlipPitch = newState;
        rotator_setPitch(hRot, -rotator_getPitch(hRot));
    }
}

void rotator_setFlipRoll(void* const hRot, int newState)
{
    rotator_data *pData = (rotator_data*)(hRot);
    if(newState !=pData->bFlipRoll ){
        pData->bFlipRoll = newState;
        rotator_setRoll(hRot, -rotator_getRoll(hRot));
    }
}

void rotator_setRPYflag(void* const hRot, int newState)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->useRollPitchYawFlag = newState;
}

void rotator_setChOrder(void* const hRot, int newOrder)
{
    rotator_data *pData = (rotator_data*)(hRot);
    if((ROTATOR_CH_ORDER)newOrder != CH_FUMA || pData->inputOrder==INPUT_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (ROTATOR_CH_ORDER)newOrder;
}

void rotator_setNormType(void* const hRot, int newType)
{
    rotator_data *pData = (rotator_data*)(hRot);
    if((ROTATOR_NORM_TYPES)newType != NORM_FUMA || pData->inputOrder==INPUT_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (ROTATOR_NORM_TYPES)newType;
}

void rotator_setOrder(void* const hRot, int newOrder)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->inputOrder = (ROTATOR_INPUT_ORDERS)newOrder;
    pData->recalc_M_rotFLAG = 1;
    /* FUMA only supports 1st order */
    if(pData->inputOrder!=INPUT_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->inputOrder!=INPUT_ORDER_FIRST && pData->norm == NORM_FUMA)
        pData->norm = NORM_SN3D;
}


/*gets*/

float rotator_getYaw(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->bFlipYaw == 1 ? -RAD2DEG(pData->yaw) : RAD2DEG(pData->yaw);
}

float rotator_getPitch(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->bFlipPitch == 1 ? -RAD2DEG(pData->pitch) : RAD2DEG(pData->pitch);
}

float rotator_getRoll(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->bFlipRoll == 1 ? -RAD2DEG(pData->roll) : RAD2DEG(pData->roll);
}

int rotator_getFlipYaw(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->bFlipYaw;
}

int rotator_getFlipPitch(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->bFlipPitch;
}

int rotator_getFlipRoll(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->bFlipRoll;
}

int rotator_getRPYflag(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->useRollPitchYawFlag;
}

int rotator_getChOrder(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return (int)pData->chOrdering;
}

int rotator_getNormType(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return (int)pData->norm;
}

int rotator_getOrder(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return (int)pData->inputOrder;
}

int rotator_getNSHrequired(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return (pData->inputOrder+1)*(pData->inputOrder+1);
}



