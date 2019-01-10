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
 *     ambi_drc.c
 * Description:
 *     A frequency-dependent spherical harmonic domain dynamic range compressor (DRC). The
 *     implementation can also keep track of the frequency-dependent gain factors for
 *     the omnidirectional component over time, for optional plotting. The design utilises
 *     a similar approach as in:
 *         McCormack, L., & Välimäki, V. (2017). "FFT-Based Dynamic Range Compression". in
 *         Proceedings of the 14th Sound and Music Computing Conference, July 5-8, Espoo,
 *         Finland.
 *     The DRC gain factors are determined based on analysing the omnidirectional component.
 *     These gain factors are then applied to the higher-order components, in a such a manner
 *     as to retain the spatial information within them.
 * Dependencies:
 *     saf_utilities, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 07.01.2017
 */

#include "ambi_drc.h"
#include "ambi_drc_internal.h"

void ambi_drc_create
(
    void ** const phAmbi
)
{
    ambi_drc_data* pData = (ambi_drc_data*)malloc(sizeof(ambi_drc_data));
    if (pData == NULL) { return;/*error*/ }
    *phAmbi = (void*)pData;
	int t, ch;
 
    /* afSTFT stuff */
    pData->hSTFT = NULL;
    pData->STFTInputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, MAX_NUM_SH_SIGNALS, sizeof(complexVector));
    for(t=0; t<TIME_SLOTS; t++) {
        for(ch=0; ch< MAX_NUM_SH_SIGNALS; ch++) {
            pData->STFTInputFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
            pData->STFTInputFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
        }
    }
    pData->tempHopFrameTD = (float**)malloc2d( MAX(MAX_NUM_SH_SIGNALS, MAX_NUM_SH_SIGNALS), HOP_SIZE, sizeof(float));
    pData->STFTOutputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, MAX_NUM_SH_SIGNALS, sizeof(complexVector));
    for(t=0; t<TIME_SLOTS; t++) {
        for(ch=0; ch< MAX_NUM_SH_SIGNALS; ch++) {
            pData->STFTOutputFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
            pData->STFTOutputFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
        }
    }
    
    /* internal */
    pData->fs = 48000;
     
#ifdef ENABLE_TF_DISPLAY
    pData->gainsTF_bank0 = (float**)malloc2d(HYBRID_BANDS, NUM_DISPLAY_TIME_SLOTS, sizeof(float));
    pData->gainsTF_bank1 = (float**)malloc2d(HYBRID_BANDS, NUM_DISPLAY_TIME_SLOTS, sizeof(float));
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
    pData->currentOrder = INPUT_ORDER_1;
    
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
    int t, ch;

    if (pData != NULL) {
        if (pData->hSTFT != NULL) {
            afSTFTfree(pData->hSTFT);
            for (t = 0; t < TIME_SLOTS; t++) {
                for (ch = 0; ch < MAX_NUM_SH_SIGNALS; ch++) {
                    free(pData->STFTInputFrameTF[t][ch].re);
                    free(pData->STFTInputFrameTF[t][ch].im);
                    free(pData->STFTOutputFrameTF[t][ch].re);
                    free(pData->STFTOutputFrameTF[t][ch].im);
                }
            }
            free2d((void**)pData->STFTInputFrameTF, TIME_SLOTS);
            free2d((void**)pData->STFTOutputFrameTF, TIME_SLOTS);
            free2d((void**)pData->tempHopFrameTD, MAX(MAX_NUM_SH_SIGNALS, MAX_NUM_SH_SIGNALS));
        }
     
#ifdef ENABLE_TF_DISPLAY
        free2d((void**)pData->gainsTF_bank0, HYBRID_BANDS);
        free2d((void**)pData->gainsTF_bank1, HYBRID_BANDS);
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
    for(band=0; band <HYBRID_BANDS; band++){
        if(sampleRate==44100)
            pData->freqVector[band] =  (float)__afCenterFreq44100[band];
        else /* assume 48e3 */
            pData->freqVector[band] =  (float)__afCenterFreq48e3[band];
    }

#ifdef ENABLE_TF_DISPLAY
    pData->rIdx = 0;
    pData->wIdx = 1;
    pData->storeIdx = 0;
    for (band = 0; band < HYBRID_BANDS; band++) {
        memset(pData->gainsTF_bank0[band], 0, NUM_DISPLAY_TIME_SLOTS * sizeof(float));
        memset(pData->gainsTF_bank1[band], 0, NUM_DISPLAY_TIME_SLOTS * sizeof(float));
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
    void*   const hAmbi,
    float** const inputs,
    float** const outputs,
    int nCh,
    int nSamples,
    int isPlaying
)                                         
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    int i, n, t, ch, band, sample;
    int o[MAX_ORDER+2];
    float xG, yG, xL, yL, cdB, alpha_a, alpha_r;
    float makeup, boost, theshold, ratio, knee;
    
    /* reinitialise if needed */
    if(pData->reInitTFT==1){
        pData->reInitTFT = 2;
        ambi_drc_initTFT(hAmbi);
        pData->reInitTFT = 0;
    }

    /* Main processing loop */
    if (nSamples == FRAME_SIZE && pData->reInitTFT == 0 && isPlaying) {
        /* prep */
        for(n=0; n<MAX_ORDER+2; n++){  o[n] = n*n;  }
        alpha_a = expf(-1.0f / ( (pData->attack_ms  / ((float)FRAME_SIZE / (float)TIME_SLOTS)) * pData->fs * 0.001f));
        alpha_r = expf(-1.0f / ( (pData->release_ms / ((float)FRAME_SIZE / (float)TIME_SLOTS)) * pData->fs * 0.001f));
        boost = powf(10.0f, pData->inGain / 20.0f);
        makeup = powf(10.0f, pData->outGain / 20.0f);
        theshold = pData->theshold;
        ratio = pData->ratio;
        knee = pData->knee;
        
        /* Load time-domain data */
        for(i=0; i < MIN(pData->nSH, nCh); i++)
            memcpy(pData->inputFrameTD[i], inputs[i], FRAME_SIZE * sizeof(float));
        for(; i<pData->nSH; i++)
            memset(pData->inputFrameTD[i], 0, FRAME_SIZE * sizeof(float));

        /* Apply time-frequency transform */
        for ( t=0; t< TIME_SLOTS; t++) {
            for( ch=0; ch < pData->nSH; ch++)
                for ( sample=0; sample < HOP_SIZE; sample++)
                    pData->tempHopFrameTD[ch][sample] = pData->inputFrameTD[ch][sample + t*HOP_SIZE];
            afSTFTforward(pData->hSTFT, (float**)pData->tempHopFrameTD, (complexVector*)pData->STFTInputFrameTF[t]);
        }
        for (band = 0; band < HYBRID_BANDS; band++)
            for (ch = 0; ch < pData->nSH; ch++)
                for (t = 0; t < TIME_SLOTS; t++)
                    pData->inputFrameTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[t][ch].re[band], pData->STFTInputFrameTF[t][ch].im[band]);
        
        /* Calculate the dynamic range compression gain factors per frequency band based on the omnidirectional component.
         * McCormack, L., & Välimäki, V. (2017). "FFT-Based Dynamic Range Compression". in Proceedings of the 14th
         * Sound and Music Computing Conference, July 5-8, Espoo, Finland.*/
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
                cdB = MAX(SPECTRAL_FLOOR, sqrtf(powf(10.0f, cdB / 20.0f)));

#ifdef ENABLE_TF_DISPLAY
                /* store gain factors in circular buffer for plotting */
                if(pData->storeIdx==0)
                    pData->gainsTF_bank0[band][pData->wIdx] = cdB;
                else
                    pData->gainsTF_bank1[band][pData->wIdx] = cdB;
#endif
                /* apply same gain factor to all SH components, the spatial characteristics will be preserved */
                for (ch = 0; ch < pData->nSH; ch++)
                    pData->outputFrameTF[band][ch][t] = crmulf(pData->inputFrameTF[band][ch][t], cdB*makeup);
            }
#ifdef ENABLE_TF_DISPLAY
            /* increment circular buffer indices */
            pData->wIdx++;
            pData->rIdx++;
            if (pData->wIdx >= NUM_DISPLAY_TIME_SLOTS){
                pData->wIdx = 0;
                pData->storeIdx = pData->storeIdx == 0 ? 1 : 0;
            }
            if (pData->rIdx >= NUM_DISPLAY_TIME_SLOTS)
                pData->rIdx = 0;
#endif
        }

        /* Inverse time-frequency transform */
        for (t = 0; t < TIME_SLOTS; t++) {
            for (ch = 0; ch < pData->nSH; ch++) {
                for (band = 0; band < HYBRID_BANDS; band++) {
                    pData->STFTOutputFrameTF[t][ch].re[band] = crealf(pData->outputFrameTF[band][ch][t]);
                    pData->STFTOutputFrameTF[t][ch].im[band] = cimagf(pData->outputFrameTF[band][ch][t]);
                }
            }
        }
        for (t = 0; t < TIME_SLOTS; t++) {
            afSTFTinverse(pData->hSTFT, pData->STFTOutputFrameTF[t], pData->tempHopFrameTD);
            for (ch = 0; ch < MIN(pData->nSH, nCh); ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    outputs[ch][sample + t* HOP_SIZE] = pData->tempHopFrameTD[ch][sample];
            for (; ch <  nCh; ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    outputs[ch][sample + t* HOP_SIZE] = 0.0f;
        }
        
    }
    else {
        for (ch=0; ch < nCh; ch++)
            memset(outputs[ch], 0, FRAME_SIZE*sizeof(float));
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
    pData->theshold = MAX(MIN(newValue, 0.0f), -60.0f);
}

void ambi_drc_setRatio(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->ratio = MAX(MIN(newValue, 30.0f), 1.0f);
}

void ambi_drc_setKnee(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->knee = MAX(MIN(newValue, 10.0f), 0.0f);
}

void ambi_drc_setInGain(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->inGain = MAX(MIN(newValue, 40.0f), -40.0f);
}

void ambi_drc_setOutGain(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->outGain = MAX(MIN(newValue, 40.0f), -20.0f);
}

void ambi_drc_setAttack(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->attack_ms = MAX(MIN(newValue, 200.0f), 10.0f);
}

void ambi_drc_setRelease(void* const hAmbi, float newValue)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->release_ms = MAX(MIN(newValue, 1000.0f), 50.0f);
}

void ambi_drc_setChOrder(void* const hAmbi, int newOrder)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);
    pData->chOrdering = (CH_ORDER)newOrder;
}

void ambi_drc_setNormType(void* const hAmbi, int newType)
{
    ambi_drc_data *pData = (ambi_drc_data*)hAmbi;
    pData->norm = (NORM_TYPES)newType;
}

void ambi_drc_setInputPreset(void* const hAmbi, INPUT_ORDER newPreset)
{
    ambi_drc_data *pData = (ambi_drc_data*)hAmbi;
    ambi_drc_setInputOrder(newPreset, &(pData->new_nSH));
    pData->currentOrder = newPreset;
    if(pData->new_nSH!=pData->nSH)
        pData->reInitTFT = 1;
}


/* GETS */

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

INPUT_ORDER ambi_drc_getInputPreset(void* const hAmbi)
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

