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
 * @file ambi_drc.c
 * @brief A frequency-dependent Ambisonic sound scene dynamic range compressor
 *        (DRC)
 *
 * The implementation can also keep track of the frequency-dependent gain
 * factors for the omnidirectional component over time (for optional plotting).
 * The design is based on the algorithm presented in [1].
 *
 * The DRC gain factors per band are determined based on the omnidirectional
 * component, which are then applied to all of the higher-order components;
 * thus, the spatial information of the Ambisonic sound scene is retained
 * (although, your perception of them may change due to the DRC).
 *
 * @see [1] McCormack, L., & Välimäki, V. (2017). "FFT-Based Dynamic Range
 *          Compression". in Proceedings of the 14th Sound and Music Computing
 *          Conference, July 5-8, Espoo, Finland.
 *
 * @author Leo McCormack
 * @date 07.01.2017
 * @license ISC
 */

#include "ambi_drc.h"
#include "ambi_drc_internal.h"

void ambi_drc_create
(
    void ** const phAmbi
)
{
    ambi_drc_data* pData = (ambi_drc_data*)malloc1d(sizeof(ambi_drc_data));
    *phAmbi = (void*)pData;
 
    /* afSTFT stuff and audio buffers*/
    pData->hSTFT = NULL;
    pData->frameTD = (float**)malloc2d(MAX_NUM_SH_SIGNALS, AMBI_DRC_FRAME_SIZE, sizeof(float));
    pData->inputFrameTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_SH_SIGNALS, TIME_SLOTS, sizeof(float_complex));
    pData->outputFrameTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_SH_SIGNALS, TIME_SLOTS, sizeof(float_complex));
    
    /* internal */
    pData->fs = 48000;
#ifdef ENABLE_TF_DISPLAY
    pData->gainsTF_bank0 = (float**)malloc2d(HYBRID_BANDS, AMBI_DRC_NUM_DISPLAY_TIME_SLOTS, sizeof(float));
    pData->gainsTF_bank1 = (float**)malloc2d(HYBRID_BANDS, AMBI_DRC_NUM_DISPLAY_TIME_SLOTS, sizeof(float));
#endif
  
    /* Default user parameters */
    pData->theshold = 0.0f;  
    pData->ratio = 8.0f;
    pData->knee = 0.0f;
    pData->inGain = 0.0f;
    pData->outGain = 0.0f;
    pData->attack_ms = 50.0f;
    pData->release_ms = 100.0f;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    pData->currentOrder = SH_ORDER_FIRST;
    
    /* for dynamically allocating the number of channels */
    ambi_drc_setInputOrder(pData->currentOrder, &(pData->new_nSH));
    pData->nSH = pData->new_nSH;
    pData->reInitTFT = 1;
}

void ambi_drc_destroy
(
    void ** const phAmbi
)
{
    ambi_drc_data *pData = (ambi_drc_data*)(*phAmbi);

    if (pData != NULL) {
        if (pData->hSTFT != NULL)
            afSTFT_destroy(&(pData->hSTFT));
        free(pData->frameTD);
        free(pData->inputFrameTF);
        free(pData->outputFrameTF);
#ifdef ENABLE_TF_DISPLAY
        free(pData->gainsTF_bank0);
        free(pData->gainsTF_bank1);
#endif
        free(pData);
        pData = NULL;
    }
}

void ambi_drc_init
(
    void * const hAmbi,
    int          sampleRate
)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    int band;

    pData->fs = (float)sampleRate;
    memset(pData->yL_z1, 0, HYBRID_BANDS * sizeof(float));
    afSTFT_getCentreFreqs(pData->hSTFT, (float)sampleRate, HYBRID_BANDS, pData->freqVector);
     
#ifdef ENABLE_TF_DISPLAY
    pData->rIdx = 0;
    pData->wIdx = 1;
    pData->storeIdx = 0;
    for (band = 0; band < HYBRID_BANDS; band++) {
        memset(pData->gainsTF_bank0[band], 0, AMBI_DRC_NUM_DISPLAY_TIME_SLOTS * sizeof(float));
        memset(pData->gainsTF_bank1[band], 0, AMBI_DRC_NUM_DISPLAY_TIME_SLOTS * sizeof(float));
    }
#endif

    /* reinitialise if needed */
    if (pData->reInitTFT == 1) {
        pData->reInitTFT = 2;
        ambi_drc_initTFT(hAmbi);
        pData->reInitTFT = 0;
    }
}

void ambi_drc_process
(
    void        *  const hAmbi,
    const float *const * inputs,
    float       ** const outputs,
    int                  nCh,
    int                  nSamples
)                                         
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    int i, t, ch, band;
    float xG, yG, xL, yL, cdB, alpha_a, alpha_r;
    float makeup, boost, theshold, ratio, knee;
    
    /* reinitialise if needed */
    if(pData->reInitTFT==1){
        pData->reInitTFT = 2;
        ambi_drc_initTFT(hAmbi);
        pData->reInitTFT = 0;
    }

    /* local copies of user parameters */
    alpha_a = expf(-1.0f / ( (pData->attack_ms  / ((float)AMBI_DRC_FRAME_SIZE / (float)TIME_SLOTS)) * pData->fs * 0.001f));
    alpha_r = expf(-1.0f / ( (pData->release_ms / ((float)AMBI_DRC_FRAME_SIZE / (float)TIME_SLOTS)) * pData->fs * 0.001f));
    boost = powf(10.0f, pData->inGain / 20.0f);
    makeup = powf(10.0f, pData->outGain / 20.0f);
    theshold = pData->theshold;
    ratio = pData->ratio;
    knee = pData->knee;

    /* Main processing loop */
    if (nSamples == AMBI_DRC_FRAME_SIZE && pData->reInitTFT == 0) {

        /* Load time-domain data */
        for(i=0; i < SAF_MIN(pData->nSH, nCh); i++)
            utility_svvcopy(inputs[i], AMBI_DRC_FRAME_SIZE, pData->frameTD[i]);
        for(; i<pData->nSH; i++)
            memset(pData->frameTD[i], 0, AMBI_DRC_FRAME_SIZE * sizeof(float));

        /* Apply time-frequency transform */
        afSTFT_forward_knownDimensions(pData->hSTFT, pData->frameTD, AMBI_DRC_FRAME_SIZE, MAX_NUM_SH_SIGNALS, TIME_SLOTS, pData->inputFrameTF);

        /* Main processing: */
        /* Calculate the dynamic range compression gain factors per frequency band based on the omnidirectional component.
            *     McCormack, L., & Välimäki, V. (2017). "FFT-Based Dynamic Range Compression". in Proceedings of the 14th
            *     Sound and Music Computing Conference, July 5-8, Espoo, Finland.*/
        for (t = 0; t < TIME_SLOTS; t++) {
            for (band = 0; band < HYBRID_BANDS; band++) {
                /* apply input boost */
                for (ch = 0; ch < pData->nSH; ch++)
                    pData->inputFrameTF[band][ch][t] = crmulf(pData->inputFrameTF[band][ch][t], boost);

                /* calculate gain factor for this frequency based on the omni component */
                xG = 10.0f*log10f(powf(cabsf(pData->inputFrameTF[band][0/* omni */][t]), 2.0f) + 2e-13f);
                yG = ambi_drc_gainComputer(xG, theshold, ratio, knee);
                xL = xG - yG;
                yL = ambi_drc_smoothPeakDetector(xL, pData->yL_z1[band], alpha_a, alpha_r);
                pData->yL_z1[band] = yL;
                cdB = -yL;
                cdB = SAF_MAX(AMBI_DRC_SPECTRAL_FLOOR, sqrtf(powf(10.0f, cdB / 20.0f)));

#ifdef ENABLE_TF_DISPLAY
                /* store gain factors in circular buffer for plotting */
                if(pData->storeIdx==0)
                    pData->gainsTF_bank0[band][pData->wIdx] = cdB;
                else
                    pData->gainsTF_bank1[band][pData->wIdx] = cdB;
#endif
                /* apply same gain factor to all SH components, the spatial characteristics will be preserved
                    * (although, ones perception of them may of course change) */
                for (ch = 0; ch < pData->nSH; ch++)
                    pData->outputFrameTF[band][ch][t] = crmulf(pData->inputFrameTF[band][ch][t], cdB*makeup);
            }
#ifdef ENABLE_TF_DISPLAY
            /* increment circular buffer indices */
            pData->wIdx++;
            pData->rIdx++;
            if (pData->wIdx >= AMBI_DRC_NUM_DISPLAY_TIME_SLOTS){
                pData->wIdx = 0;
                pData->storeIdx = pData->storeIdx == 0 ? 1 : 0;
            }
            if (pData->rIdx >= AMBI_DRC_NUM_DISPLAY_TIME_SLOTS)
                pData->rIdx = 0;
#endif
        }

        /* Inverse time-frequency transform */
        afSTFT_backward_knownDimensions(pData->hSTFT, pData->outputFrameTF, AMBI_DRC_FRAME_SIZE, MAX_NUM_SH_SIGNALS, TIME_SLOTS, pData->frameTD);

        /* Copy to output */
        for(ch = 0; ch < SAF_MIN(pData->nSH, nCh); ch++)
            utility_svvcopy(pData->frameTD[ch], AMBI_DRC_FRAME_SIZE, outputs[ch]);
        for (; ch < nCh; ch++)
            memset(outputs[ch], 0, AMBI_DRC_FRAME_SIZE*sizeof(float));
    }
    else {
        for (ch=0; ch < nCh; ch++)
            memset(outputs[ch], 0, AMBI_DRC_FRAME_SIZE*sizeof(float));
    }
}

/* SETS */

void ambi_drc_refreshSettings(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->reInitTFT = 1;
}

void ambi_drc_setThreshold(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->theshold = SAF_CLAMP(newValue, AMBI_DRC_THRESHOLD_MIN_VAL, AMBI_DRC_THRESHOLD_MAX_VAL);
}

void ambi_drc_setRatio(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->ratio = SAF_CLAMP(newValue, AMBI_DRC_RATIO_MIN_VAL, AMBI_DRC_RATIO_MAX_VAL);
}

void ambi_drc_setKnee(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->knee = SAF_CLAMP(newValue, AMBI_DRC_KNEE_MIN_VAL, AMBI_DRC_KNEE_MAX_VAL);
}

void ambi_drc_setInGain(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->inGain = SAF_CLAMP(newValue, AMBI_DRC_IN_GAIN_MIN_VAL, AMBI_DRC_IN_GAIN_MAX_VAL);
}

void ambi_drc_setOutGain(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->outGain = SAF_CLAMP(newValue, AMBI_DRC_OUT_GAIN_MIN_VAL, AMBI_DRC_OUT_GAIN_MAX_VAL);
}

void ambi_drc_setAttack(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->attack_ms = SAF_CLAMP(newValue, AMBI_DRC_ATTACK_MIN_VAL, AMBI_DRC_ATTACK_MAX_VAL);
}

void ambi_drc_setRelease(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->release_ms = SAF_CLAMP(newValue, AMBI_DRC_RELEASE_MIN_VAL, AMBI_DRC_RELEASE_MAX_VAL);
}

void ambi_drc_setChOrder(void* const hAmbi, int newOrder)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    if((CH_ORDER)newOrder != CH_FUMA || pData->currentOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void ambi_drc_setNormType(void* const hAmbi, int newType)
{
    ambi_drc_data *pData = (ambi_drc_data*)hAmbi;
    if((NORM_TYPES)newType != NORM_FUMA || pData->currentOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
}

void ambi_drc_setInputPreset(void* const hAmbi, SH_ORDERS newPreset)
{
    ambi_drc_data *pData = (ambi_drc_data*)hAmbi;
    ambi_drc_setInputOrder(newPreset, &(pData->new_nSH));
    pData->currentOrder = newPreset;
    if(pData->new_nSH!=pData->nSH)
        pData->reInitTFT = 1;
    /* FUMA only supports 1st order */
    if(pData->currentOrder!=SH_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->currentOrder!=SH_ORDER_FIRST && pData->norm == NORM_FUMA)
        pData->norm = NORM_SN3D;
}


/* GETS */

int ambi_drc_getFrameSize(void)
{
    return AMBI_DRC_FRAME_SIZE;
}

#ifdef ENABLE_TF_DISPLAY
float** ambi_drc_getGainTF(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    if(pData->storeIdx==0)
        return pData->gainsTF_bank0;
    else
        return pData->gainsTF_bank1;
}

int ambi_drc_getGainTFwIdx(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return pData->wIdx;
}

int ambi_drc_getGainTFrIdx(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return pData->rIdx;
}

float* ambi_drc_getFreqVector(void* const hAmbi, int* nFreqPoints)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    (*nFreqPoints) = HYBRID_BANDS;
    return pData->freqVector;
}
#endif

float ambi_drc_getThreshold(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return pData->theshold;
}

float ambi_drc_getRatio(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return pData->ratio;
}

float ambi_drc_getKnee(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return pData->knee;
}

float ambi_drc_getInGain(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return pData->inGain;
}

float ambi_drc_getOutGain(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return pData->outGain;
}

float ambi_drc_getAttack(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return pData->attack_ms;
}

float ambi_drc_getRelease(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return pData->release_ms;
}

int ambi_drc_getChOrder(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return (int)pData->chOrdering;
}

int ambi_drc_getNormType(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return (int)pData->norm;
}

SH_ORDERS ambi_drc_getInputPreset(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return pData->currentOrder;
}

int ambi_drc_getNSHrequired(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return pData->nSH;
}

int ambi_drc_getSamplerate(void* const hAmbi)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    return (int)(pData->fs+0.5f);
}

int ambi_drc_getProcessingDelay()
{
    return 12*HOP_SIZE;
}

