    /*
 Copyright 2017-2018 Leo McCormack
 
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
 *     rotator.c
 * Description:
 *     A simple spherical harmonic domain rotator.
 * Dependencies:
 *     saf_utilities, saf_sh
 * Author, date created:
 *     Leo McCormack, 02.11.2017
 */

#include "rotator.h"
#include "rotator_internal.h"


void rotator_create
(
    void ** const phRot
)
{
    rotator_data* pData = (rotator_data*)malloc(sizeof(rotator_data));
    if (pData == NULL) { return;/*error*/ }
    *phRot = (void*)pData;
  
    /* Default user parameters */
    pData->yaw = 0.0f;
    pData->pitch = 0.0f;
    pData->roll = 0.0f;
    pData->bFlipYaw = 0;
    pData->bFlipPitch = 0;
    pData->bFlipRoll = 0;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_N3D;
    rotator_setOrder(*phRot,  OUTPUT_ORDER_FIRST);
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
}


void rotator_process
(
    void  *  const hRot,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples,
    int            isPlaying
)
{
    rotator_data *pData = (rotator_data*)(hRot);
    int i, j, n, order, nSH;
    int o[MAX_SH_ORDER+2];
    float Rxyz[3][3];
    float* M_rot_tmp;
    CH_ORDER chOrdering;
    NORM_TYPES norm;
 
    if (nSamples == FRAME_SIZE && isPlaying) {
        /* prep */
        for(n=0; n<MAX_SH_ORDER+2; n++){  o[n] = n*n;  }
        chOrdering = pData->chOrdering;
        norm = pData->norm;
        order = pData->order;
        nSH = (order+1)*(order+1);
        for (i = 0; i < MIN(nSH, nInputs); i++)
            memcpy(pData->inputFrameTD[i], inputs[i], FRAME_SIZE * sizeof(float));
        for (; i < nSH; i++)
            memset(pData->inputFrameTD[i], 0, FRAME_SIZE * sizeof(float));
        
        /* account for norm scheme */
        switch(norm){
            case NORM_N3D: /* already N3D */
                break;
            case NORM_SN3D: /* convert to N3D before rotation */
#if 0 /* actually doesn't seem to matter */
                for (n = 0; n<order+1; n++)
                    for (ch = o[n]; ch<o[n+1]; ch++)
                        for(i = 0; i<FRAME_SIZE; i++)
                            pData->inputFrameTD[ch][i] *= sqrtf(2.0f*(float)n+1.0f);
#endif
                break;
        }
        
        if (order>0){
            /* calculate rotation matrix */
            M_rot_tmp = malloc(nSH*nSH*sizeof(float));
            yawPitchRoll2Rzyx (pData->yaw, pData->pitch, pData->roll, Rxyz);
            getSHrotMtxReal(Rxyz, M_rot_tmp, order);
            for(i=0; i<nSH; i++)
                for(j=0; j<nSH; j++)
                    pData->M_rot[i][j] = M_rot_tmp[i*nSH+j];
            free(M_rot_tmp);
            
            /* apply rotation (assumes ACN/N3D) */
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
            memcpy(pData->prev_inputFrameTD, pData->inputFrameTD, nSH*FRAME_SIZE*sizeof(float));
            for(i=0; i<nSH; i++)
                memcpy(pData->prev_M_rot[i], pData->M_rot[i], nSH*sizeof(float));
        }
        else
            memcpy(pData->outputFrameTD[0], pData->inputFrameTD[0], FRAME_SIZE*sizeof(float));
        
        /* account for norm scheme */
        switch(norm){
            case NORM_N3D: /* already N3D */
                break;
            case NORM_SN3D: /* convert back to SN3D after rotation */
#if 0 /* actually doesn't seem to matter */
                for (n = 0; n<order+1; n++)
                    for (ch = o[n]; ch<o[n+1]; ch++)
                        for(i = 0; i<FRAME_SIZE; i++)
                            pData->outputFrameTD[ch][i] /= sqrtf(2.0f*(float)n+1.0f);
#endif
                break;
        }
        for (i = 0; i < MIN(nSH, nOutputs); i++)
            memcpy(outputs[i], pData->outputFrameTD[i], FRAME_SIZE*sizeof(float));
        for (; i < nOutputs; i++)
            memset(outputs[i], 0, FRAME_SIZE*sizeof(float));
    }
    else{
        for (i = 0; i < nOutputs; i++)
            memset(outputs[i], 0, FRAME_SIZE*sizeof(float));
    }
}

void rotator_setYaw(void  * const hRot, float newYaw)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->yaw = pData->bFlipYaw == 1 ? -DEG2RAD(newYaw) : DEG2RAD(newYaw);
}

void rotator_setPitch(void* const hRot, float newPitch)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->pitch = pData->bFlipPitch == 1 ? -DEG2RAD(newPitch) : DEG2RAD(newPitch);
}

void rotator_setRoll(void* const hRot, float newRoll)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->roll = pData->bFlipRoll == 1 ? -DEG2RAD(newRoll) : DEG2RAD(newRoll);
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

void rotator_setChOrder(void* const hRot, int newOrder)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->chOrdering = (CH_ORDER)newOrder;
}

void rotator_setNormType(void* const hRot, int newType)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->norm = (NORM_TYPES)newType;
}

void rotator_setOrder(void* const hRot, int newOrder)
{
    rotator_data *pData = (rotator_data*)(hRot);
    pData->outputOrder = (OUTPUT_ORDERS)newOrder;
    switch(pData->outputOrder){
        case OUTPUT_OMNI:          pData->order = 0; break;
        default:
        case OUTPUT_ORDER_FIRST:   pData->order = 1; break;
        case OUTPUT_ORDER_SECOND:  pData->order = 2; break;
        case OUTPUT_ORDER_THIRD:   pData->order = 3; break;
        case OUTPUT_ORDER_FOURTH:  pData->order = 4; break;
        case OUTPUT_ORDER_FIFTH:   pData->order = 5; break;
        case OUTPUT_ORDER_SIXTH:   pData->order = 6; break;
        case OUTPUT_ORDER_SEVENTH: pData->order = 7; break;
    }
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
    return (int)pData->outputOrder;
}




