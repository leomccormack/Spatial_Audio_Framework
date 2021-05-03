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
 * @brief A basic spherical harmonic/ Ambisonic signals rotator, based on the
 *        recursive approach detailed in [1]
 *
 * @test test__saf_example_rotator()
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

    printf(SAF_VERSION_LICENSE_STRING);
    
    pData->M_rot_status = M_ROT_RECOMPUTE_QUATERNION;
  
    /* Default user parameters */
    pData->Q.w = 1.0f;
    pData->Q.x = 0.0f;
    pData->Q.y = 0.0f;
    pData->Q.z = 0.0f;
    pData->bFlipQuaternion = 0;
    pData->yaw = 0.0f;
    pData->pitch = 0.0f;
    pData->roll = 0.0f;
    pData->bFlipYaw = 0;
    pData->bFlipPitch = 0;
    pData->bFlipRoll = 0;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    pData->useRollPitchYawFlag = 0;
    rotator_setOrder(*phRot, SH_ORDER_FIRST);
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

    pData->fs = sampleRate;
    
    /* starting values */
    for(i=1; i<=FRAME_SIZE; i++)
        pData->interpolator[i-1] = (float)i*1.0f/(float)FRAME_SIZE;
    memset(pData->M_rot, 0, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*sizeof(float));
    memset(pData->prev_M_rot, 0, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*sizeof(float));
    memset(pData->prev_inputFrameTD, 0, MAX_NUM_SH_SIGNALS*FRAME_SIZE*sizeof(float));
    pData->M_rot_status = M_ROT_RECOMPUTE_QUATERNION;
}

void rotator_process
(
    void        *  const hRot,
    const float ** const inputs,
    float       ** const outputs,
    int                  nInputs,
    int                  nOutputs,
    int                  nSamples
)
{
    rotator_data *pData = (rotator_data*)(hRot);
    int i, j, order, nSH; 
    float Rxyz[3][3];
    float* M_rot_tmp;
    CH_ORDER chOrdering;

    /* locals */
    chOrdering = pData->chOrdering;
    order = (int)pData->inputOrder;
    nSH = ORDER2NSH(order);

    if (nSamples == FRAME_SIZE) {

        /* Load time-domain data */
        for(i=0; i < MIN(nSH, nInputs); i++)
            utility_svvcopy(inputs[i], FRAME_SIZE, pData->inputFrameTD[i]);
        for(; i<nSH; i++)
            memset(pData->inputFrameTD[i], 0, FRAME_SIZE * sizeof(float)); /* fill remaining channels with zeros */

        /* account for channel order */
        switch(chOrdering){
            case CH_ACN: /* already ACN */
                break;
            case CH_FUMA:
                convertHOAChannelConvention((float*)pData->inputFrameTD, order, FRAME_SIZE, HOA_CH_ORDER_FUMA, HOA_CH_ORDER_ACN);
                break;
        }

        if (order>0){
            /* calculate rotation matrix */
            if(pData->M_rot_status != M_ROT_READY){
                memset(pData->M_rot, 0, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*sizeof(float));
                M_rot_tmp = malloc1d(nSH*nSH*sizeof(float));
                if(pData->M_rot_status == M_ROT_RECOMPUTE_EULER){
                    yawPitchRoll2Rzyx (pData->yaw, pData->pitch, pData->roll, pData->useRollPitchYawFlag, Rxyz);
                    euler2Quaternion(pData->yaw, pData->pitch, pData->roll, 0,
                                     pData->useRollPitchYawFlag ? EULER_ROTATION_ROLL_PITCH_YAW : EULER_ROTATION_YAW_PITCH_ROLL, &(pData->Q));
                }
                else {/* M_ROT_RECOMPUTE_QUATERNION */
                    quaternion2rotationMatrix(&(pData->Q), Rxyz);
                    quaternion2euler(&(pData->Q), 0, pData->useRollPitchYawFlag ? EULER_ROTATION_ROLL_PITCH_YAW : EULER_ROTATION_YAW_PITCH_ROLL,
                                     &(pData->yaw), &(pData->pitch), &(pData->roll));
                }
                getSHrotMtxReal(Rxyz, M_rot_tmp, order);
                for(i=0; i<nSH; i++)
                    for(j=0; j<nSH; j++)
                        pData->M_rot[i][j] = M_rot_tmp[i*nSH+j];
                free(M_rot_tmp);
                pData->M_rot_status = M_ROT_READY;
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
  
        /* account for channel order */
        switch(chOrdering){
            case CH_ACN: /* already ACN */
                break;
            case CH_FUMA:
                convertHOAChannelConvention((float*)pData->outputFrameTD, order, FRAME_SIZE, HOA_CH_ORDER_ACN, HOA_CH_ORDER_FUMA);
                break;
        }

        /* Copy to output */
        for (i = 0; i < MIN(nSH, nOutputs); i++)
            utility_svvcopy(pData->outputFrameTD[i], FRAME_SIZE, outputs[i]);
        for (; i < nOutputs; i++)
            memset(outputs[i], 0, FRAME_SIZE*sizeof(float)); 
    }
    else{
        for (i = 0; i < nOutputs; i++)
            memset(outputs[i], 0, FRAME_SIZE*sizeof(float));
    }
}


/*sets*/

void rotator_setYaw(void* const hRot, float newYaw)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->yaw = pData->bFlipYaw == 1 ? -DEG2RAD(newYaw) : DEG2RAD(newYaw);
    pData->M_rot_status = M_ROT_RECOMPUTE_EULER;
}

void rotator_setPitch(void* const hRot, float newPitch)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->pitch = pData->bFlipPitch == 1 ? -DEG2RAD(newPitch) : DEG2RAD(newPitch);
    pData->M_rot_status = M_ROT_RECOMPUTE_EULER;
}

void rotator_setRoll(void* const hRot, float newRoll)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->roll = pData->bFlipRoll == 1 ? -DEG2RAD(newRoll) : DEG2RAD(newRoll);
    pData->M_rot_status = M_ROT_RECOMPUTE_EULER;
}

void rotator_setQuaternionW(void* const hRot, float newValue)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->Q.w = newValue;
    pData->M_rot_status = M_ROT_RECOMPUTE_QUATERNION;
}

void rotator_setQuaternionX(void* const hRot, float newValue)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->Q.x = pData->bFlipQuaternion == 1 ? -newValue : newValue;
    pData->M_rot_status = M_ROT_RECOMPUTE_QUATERNION;
}

void rotator_setQuaternionY(void* const hRot, float newValue)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->Q.y = pData->bFlipQuaternion == 1 ? -newValue : newValue;
    pData->M_rot_status = M_ROT_RECOMPUTE_QUATERNION;
}

void rotator_setQuaternionZ(void* const hRot, float newValue)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->Q.z = pData->bFlipQuaternion == 1 ? -newValue : newValue;
    pData->M_rot_status = M_ROT_RECOMPUTE_QUATERNION;
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

void rotator_setFlipQuaternion(void* const hRot, int newState)
{
    rotator_data *pData = (rotator_data*)(hRot);
    if(newState !=pData->bFlipQuaternion ){
        pData->bFlipQuaternion = newState;
        rotator_setQuaternionX(hRot, -rotator_getQuaternionX(hRot));
        rotator_setQuaternionY(hRot, -rotator_getQuaternionY(hRot));
        rotator_setQuaternionZ(hRot, -rotator_getQuaternionZ(hRot));
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
    if((CH_ORDER)newOrder != CH_FUMA || pData->inputOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void rotator_setNormType(void* const hRot, int newType)
{
    rotator_data *pData = (rotator_data*)(hRot);
    if((NORM_TYPES)newType != NORM_FUMA || pData->inputOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
}

void rotator_setOrder(void* const hRot, int newOrder)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->inputOrder = (SH_ORDERS)newOrder;
    pData->M_rot_status = M_ROT_RECOMPUTE_QUATERNION;
    /* FUMA only supports 1st order */
    if(pData->inputOrder!=SH_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->inputOrder!=SH_ORDER_FIRST && pData->norm == NORM_FUMA)
        pData->norm = NORM_SN3D;
}


/*gets*/

int rotator_getFrameSize(void)
{
    return FRAME_SIZE;
}

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

float rotator_getQuaternionW(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->Q.w;
}

float rotator_getQuaternionX(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->bFlipQuaternion == 1 ? -(pData->Q.x) : pData->Q.x;
}

float rotator_getQuaternionY(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->bFlipQuaternion == 1 ? -(pData->Q.y) : pData->Q.y;
}

float rotator_getQuaternionZ(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->bFlipQuaternion == 1 ? -(pData->Q.z) : pData->Q.z;
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

int rotator_getFlipQuaternion(void* const hRot)
{
    rotator_data *pData = (rotator_data*)(hRot);
    return pData->bFlipQuaternion;
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

int rotator_getProcessingDelay()
{
    return FRAME_SIZE;
}
