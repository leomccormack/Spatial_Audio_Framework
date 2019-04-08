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
 *     array2sh.c
 * Description:
 *     Spatially encodes spherical or cylindrical sensor array signals into spherical harmonic
 *     signals utilising theoretical encoding filters.
 *     The algorithms within array2sh were pieced together and developed in collaboration
 *     with Symeon Delikaris-Manias.
 *     A more detailed explanation of the algorithms in array2sh can be found in:
 *         McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki, V.,
 *         “Real-time conversion of sensor array signals into spherical harmonic signals with
 *         applications to spatially localised sub-band sound-field analysis,” in Audio
 *         Engineering Society Convention 144, Audio Engineering Society, 2018.
 *     Also included, is a diffuse-field equalisation option for frequencies past aliasing,
 *     developed in collaboration with Archontis Politis, 8.02.2019
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_sh
 * Author, date created:
 *     Leo McCormack, 13.09.2017
 */

#include "array2sh_internal.h" 

void array2sh_create
(
    void ** const phA2sh
)
{
    array2sh_data* pData = (array2sh_data*)malloc(sizeof(array2sh_data));
    if (pData == NULL) { return;/*error*/ }
    *phA2sh = (void*)pData;
    int ch;
     
    /* defualt parameters */
    array2sh_createArray(&(pData->arraySpecs)); 
    pData->filterType = FILTER_TIKHONOV;
    pData->regPar = 15.0f;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    pData->c = 343.0f;
    pData->gain_dB = 0.0f; /* post-gain */
    pData->maxFreq = 20e3f;
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    array2sh_initArray(arraySpecs, PRESET_DEFAULT, &(pData->order), 1);
    
    /* time-frequency transform + buffers */
    pData->hSTFT = NULL;
    pData->STFTInputFrameTF = (complexVector*)malloc(MAX_NUM_SENSORS * sizeof(complexVector));
    for(ch=0; ch< MAX_NUM_SENSORS; ch++) {
        pData->STFTInputFrameTF[ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
        pData->STFTInputFrameTF[ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
    }
    pData->STFTOutputFrameTF = (complexVector*)malloc(MAX_NUM_SH_SIGNALS * sizeof(complexVector));
    for(ch=0; ch< MAX_NUM_SH_SIGNALS; ch++) {
        pData->STFTOutputFrameTF[ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
        pData->STFTOutputFrameTF[ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
    }
    pData->tempHopFrameTD_in = (float**)malloc2d( MAX(MAX_NUM_SH_SIGNALS, MAX_NUM_SENSORS), HOP_SIZE, sizeof(float));
    pData->tempHopFrameTD_out = (float**)malloc2d( MAX(MAX_NUM_SH_SIGNALS, MAX_NUM_SENSORS), HOP_SIZE, sizeof(float));
    pData->reinitTFTFLAG = 1;
    pData->applyDiffEQFLAG = 1;
    
    /* internal */
    pData->reinitSHTmatrixFLAG = 1;
    pData->new_order = pData->order;
    pData->nSH = pData->new_nSH = (pData->order+1)*(pData->order+1);
    pData->bN = NULL;
    pData->evalReady = 0;
    
    /* display related stuff */
    pData->bN_modal_dB = (float**)malloc2d(HYBRID_BANDS, MAX_SH_ORDER + 1, sizeof(float));
    pData->bN_inv_dB = (float**)malloc2d(HYBRID_BANDS, MAX_SH_ORDER + 1, sizeof(float));
    pData->cSH = (float*)calloc((HYBRID_BANDS)*(MAX_SH_ORDER + 1),sizeof(float));
    pData->lSH = (float*)calloc((HYBRID_BANDS)*(MAX_SH_ORDER + 1),sizeof(float));
    
    pData->recalcEvalFLAG = 1;
}

void array2sh_destroy
(
    void ** const phM2sh
)
{
    array2sh_data *pData = (array2sh_data*)(*phM2sh);
    int ch;

    if (pData != NULL) {
        /* TFT stuff */
        if (pData->hSTFT != NULL)
            afSTFTfree(pData->hSTFT);
        for (ch = 0; ch< MAX_NUM_SENSORS; ch++) {
            free(pData->STFTInputFrameTF[ch].re);
            free(pData->STFTInputFrameTF[ch].im);
        }
        for(ch=0; ch< MAX_NUM_SH_SIGNALS; ch++) {
            free(pData->STFTOutputFrameTF[ch].re);
            free(pData->STFTOutputFrameTF[ch].im);
        }
        free(pData->STFTOutputFrameTF);
        free(pData->STFTInputFrameTF);
        free2d((void**)pData->tempHopFrameTD_in, MAX(MAX_NUM_SENSORS, MAX_NUM_SH_SIGNALS));
        free2d((void**)pData->tempHopFrameTD_out, MAX(MAX_NUM_SENSORS, MAX_NUM_SH_SIGNALS));
        array2sh_destroyArray(&(pData->arraySpecs));
        
        /* Display stuff */
        free2d((void**)pData->bN_modal_dB, HYBRID_BANDS-1);
        free2d((void**)pData->bN_inv_dB, HYBRID_BANDS-1);
        
        free(pData);
        pData = NULL;
    }
}

void array2sh_init
(
    void * const hA2sh,
    int          sampleRate
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    int band;
    
    pData->fs = sampleRate;
    for(band=0; band <HYBRID_BANDS; band++){
        if(sampleRate==44100)
            pData->freqVector[band] =  (float)__afCenterFreq44100[band];
        else /* assume 48e3 */
            pData->freqVector[band] =  (float)__afCenterFreq48e3[band];
    } 
    pData->freqVector[0] = pData->freqVector[1]/4.0f; /* avoids NaNs at DC */

    /* reinitialise if needed */
    array2sh_checkReInit(hA2sh);
}

void array2sh_process
(
    void  *  const hA2sh,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples,
    int            isPlaying
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    int n, t, sample, ch, i, band, Q, order, nSH;
    int o[MAX_SH_ORDER+2];
    const float_complex calpha = cmplxf(1.0f,0.0f), cbeta = cmplxf(0.0f, 0.0f);
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    float gain_lin, maxFreq;
    
    /* reinitialise if needed */
#ifdef __APPLE__
    array2sh_checkReInit(hA2sh);
#else
    if (pData->reinitTFTFLAG==1) {
        pData->reinitTFTFLAG = 2;
        array2sh_initTFT(hA2sh);
        pData->reinitTFTFLAG=0;
    }
    if (pData->reinitSHTmatrixFLAG==1) {
        pData->reinitSHTmatrixFLAG = 2;
        /* compute encoding matrix */
        array2sh_calculate_sht_matrix(hA2sh);
        /* calculate magnitude response curves */
        array2sh_calculate_mag_curves(hA2sh);
        pData->reinitSHTmatrixFLAG = 0;
    }
#endif
    
    if ((nSamples == FRAME_SIZE) && !(pData->recalcEvalFLAG) && !(pData->reinitSHTmatrixFLAG) && !(pData->reinitTFTFLAG)) {
        /* prep */
        for(n=0; n<MAX_SH_ORDER+2; n++){  o[n] = n*n;  }
        chOrdering = pData->chOrdering;
        norm = pData->norm;
        gain_lin = powf(10.0f, pData->gain_dB/20.0f);
        maxFreq = pData->maxFreq;
        Q = arraySpecs->Q;
        order = pData->order;
        nSH = pData->nSH;
        
        /* Load time-domain data */
        for(i=0; i < nInputs; i++)
            memcpy(pData->inputFrameTD[i], inputs[i], FRAME_SIZE * sizeof(float));
        for(; i<Q; i++)
            memset(pData->inputFrameTD[i], 0, FRAME_SIZE * sizeof(float));
        
        /* Apply time-frequency transform (TFT) */
        for(t=0; t< TIME_SLOTS; t++) {
            for(ch=0; ch < Q; ch++)
                for(sample=0; sample < HOP_SIZE; sample++)
                    pData->tempHopFrameTD_in[ch][sample] = pData->inputFrameTD[ch][sample + t*HOP_SIZE];
            afSTFTforward(pData->hSTFT, pData->tempHopFrameTD_in, pData->STFTInputFrameTF);
            for(band=0; band<HYBRID_BANDS; band++)
                for(ch=0; ch < Q; ch++)
                    pData->inputframeTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[ch].re[band], pData->STFTInputFrameTF[ch].im[band]);
        }
        
        /* Apply spherical harmonic transform (SHT) */
        if(isPlaying){
            for(band=0; band<HYBRID_BANDS; band++){
                cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, TIME_SLOTS, Q, &calpha,
                            pData->W[band], MAX_NUM_SENSORS,
                            pData->inputframeTF[band], TIME_SLOTS, &cbeta,
                            pData->SHframeTF[band], TIME_SLOTS);
            }
        }
        else
            memset(pData->SHframeTF, 0, HYBRID_BANDS*MAX_NUM_SH_SIGNALS*TIME_SLOTS*sizeof(float_complex));
        
        /* inverse-TFT */
        for(t = 0; t < TIME_SLOTS; t++) {
            for(band = 0; band < HYBRID_BANDS; band++) {
                if(pData->freqVector[band] < maxFreq){
                    for (ch = 0; ch < nSH; ch++) {
                        pData->STFTOutputFrameTF[ch].re[band] = crealf(pData->SHframeTF[band][ch][t]);
                        pData->STFTOutputFrameTF[ch].im[band] = cimagf(pData->SHframeTF[band][ch][t]);
                    }
                }
                else{
                    for (ch = 0; ch < nSH; ch++) {
                        pData->STFTOutputFrameTF[ch].re[band] = 0.0f;
                        pData->STFTOutputFrameTF[ch].im[band] = 0.0f;
                    }
                }
            }
            afSTFTinverse(pData->hSTFT, pData->STFTOutputFrameTF, pData->tempHopFrameTD_out);
            for (ch = 0; ch < MIN(nSH, nOutputs); ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    outputs[ch][sample + t* HOP_SIZE] = pData->tempHopFrameTD_out[ch][sample] * gain_lin;
            for (; ch < nOutputs; ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    outputs[ch][sample + t* HOP_SIZE] = 0.0f;
        }
        
        /* apply normalisation scheme */
        switch(norm){
            case NORM_N3D: /* already N3D */
                break;
            case NORM_SN3D:
                for (n = 0; n<order+1; n++)
                    for (ch = o[n]; ch<o[n+1]; ch++)
                        for(i = 0; i<FRAME_SIZE; i++)
                            outputs[ch][i] /= sqrtf(2.0f*(float)n+1.0f);
                break;
        }
    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, FRAME_SIZE*sizeof(float));
    }
}

/* Set Functions */

void array2sh_refreshSettings(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->reinitTFTFLAG = 1;
    pData->reinitSHTmatrixFLAG = 1; 
    pData->applyDiffEQFLAG = 1;
}

void array2sh_checkReInit(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    /* reinitialise if needed */
    if (pData->reinitTFTFLAG == 1) {
        pData->reinitTFTFLAG = 2;
        array2sh_initTFT(hA2sh);
        pData->reinitTFTFLAG = 0;
    }
    if (pData->reinitSHTmatrixFLAG == 1) {
        pData->reinitSHTmatrixFLAG = 2;
        /* compute encoding matrix */
        array2sh_calculate_sht_matrix(hA2sh);
        /* calculate magnitude response curves */
        array2sh_calculate_mag_curves(hA2sh);
        pData->reinitSHTmatrixFLAG = 0;
    }
    /* Too heavy to put in main loop: */
    if (pData->recalcEvalFLAG == 1) {
        pData->recalcEvalFLAG = 2;
        array2sh_evaluateSHTfilters(hA2sh);
        pData->recalcEvalFLAG = 0;
    }
    if (pData->applyDiffEQFLAG == 1) {
        pData->applyDiffEQFLAG = 2;
        array2sh_apply_diff_EQ(hA2sh);
        pData->applyDiffEQFLAG = 0;
    }
}

void array2sh_setEncodingOrder(void* const hA2sh, int newOrder)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh); 
    pData->new_order = newOrder;
    pData->new_nSH = (newOrder+1)*(newOrder+1);
    pData->reinitTFTFLAG = 1;
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_evaluateFilters(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->recalcEvalFLAG = 1;
}

void array2sh_applyDiffEQpastAliasing(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->applyDiffEQFLAG = 1;
}

void array2sh_setPreset(void* const hA2sh, int preset)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    
    array2sh_initArray(arraySpecs,(PRESETS)preset, &(pData->new_order), 0);
    pData->c = (PRESETS)preset == PRESET_AALTO_HYDROPHONE ? 1484.0f : 343.0f;
    pData->new_nSH = (pData->new_order+1)*(pData->new_order+1);
    pData->reinitTFTFLAG = 1;
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_setSensorAzi_rad(void* const hA2sh, int index, float newAzi_rad)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    
    arraySpecs->sensorCoords_rad[index][0] = newAzi_rad;
    arraySpecs->sensorCoords_deg[index][0] = newAzi_rad * (180.0f/M_PI);
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_setSensorElev_rad(void* const hA2sh, int index, float newElev_rad)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    
    arraySpecs->sensorCoords_rad[index][1] = newElev_rad;
    arraySpecs->sensorCoords_deg[index][1] = newElev_rad * (180.0f/M_PI);
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_setSensorAzi_deg(void* const hA2sh, int index, float newAzi_deg)

{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    
    arraySpecs->sensorCoords_rad[index][0] = newAzi_deg * (M_PI/180.0f);
    arraySpecs->sensorCoords_deg[index][0] = newAzi_deg;
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_setSensorElev_deg(void* const hA2sh, int index, float newElev_deg)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    
    arraySpecs->sensorCoords_rad[index][1] = newElev_deg * (M_PI/180.0f);
    arraySpecs->sensorCoords_deg[index][1] = newElev_deg;
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_setNumSensors(void* const hA2sh, int newQ)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    
    if (newQ < pData->nSH){
        pData->new_order = 1;
        pData->new_nSH = (pData->new_order+1)*(pData->new_order+1);
    }
    arraySpecs->newQ = newQ;
    //arraySpecs->newQ = newQ <= pData->new_nSH ? pData->new_nSH : newQ;
    pData->reinitTFTFLAG = arraySpecs->Q != arraySpecs->newQ ? 1 : 0;
    pData->reinitSHTmatrixFLAG = arraySpecs->Q != arraySpecs->newQ ? 1 : 0;
}

void array2sh_setr(void* const hA2sh, float newr)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    arraySpecs->r = newr;
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_setR(void* const hA2sh, float newR)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    arraySpecs->R = newR;
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_setArrayType(void* const hA2sh, int newType)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    arraySpecs->arrayType = (ARRAY_TYPES)newType;
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_setWeightType(void* const hA2sh, int newType)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    
    arraySpecs->weightType = (WEIGHT_TYPES)newType;
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_setFilterType(void* const hA2sh, int newType)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->filterType = (FILTER_TYPES)newType;
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_setRegPar(void* const hA2sh, float newVal)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->regPar = newVal;
    pData->reinitSHTmatrixFLAG = 1;
}

void array2sh_setChOrder(void* const hA2sh, int newOrder)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->chOrdering = (CH_ORDER)newOrder;
}

void array2sh_setNormType(void* const hA2sh, int newType)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->norm = (NORM_TYPES)newType;
}

void array2sh_setc(void* const hA2sh, float newc)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->c = newc;
    pData->reinitSHTmatrixFLAG = 1;
}


void array2sh_setGain(void* const hA2sh, float newGain)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->gain_dB = newGain;
}

void array2sh_setMaxFreq(void* const hA2sh, float newF)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->maxFreq = newF;
}


/* Get Functions */

int array2sh_getEvalReady(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    if (pData->evalReady) {
        pData->evalReady = 0;
        return 1;
    }
    else
        return 0;
}

int array2sh_getIsEvalValid(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->currentEvalIsValid;
}

int array2sh_getEncodingOrder(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->new_order;
}

float array2sh_getSensorAzi_rad(void* const hA2sh, int index)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    return arraySpecs->sensorCoords_rad[index][0];
}

float array2sh_getSensorElev_rad(void* const hA2sh, int index)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    return arraySpecs->sensorCoords_rad[index][1];
}

float array2sh_getSensorAzi_deg(void* const hA2sh, int index)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    return arraySpecs->sensorCoords_deg[index][0];
}

float array2sh_getSensorElev_deg(void* const hA2sh, int index)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    return arraySpecs->sensorCoords_deg[index][1];
}

int array2sh_getNumSensors(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
   // return arraySpecs->Q;
    return arraySpecs->newQ; /* return the new Q, incase the plug-in is still waiting for a refresh */
}

int array2sh_getMaxNumSensors(void)
{
    return MAX_NUM_SENSORS;
}

int array2sh_getMinNumSensors(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->new_nSH;
}

int array2sh_getNSHrequired(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->new_nSH;
}

float array2sh_getr(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    return arraySpecs->r;
}

float array2sh_getR(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    return arraySpecs->R;
} 

int array2sh_getArrayType(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    return (int)arraySpecs->arrayType;
}

int array2sh_getWeightType(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    return (int)arraySpecs->weightType;
}

int array2sh_getFilterType(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return (int)pData->filterType;
}

float array2sh_getRegPar(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->regPar;
}

int array2sh_getChOrder(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return (int)pData->chOrdering;
}

int array2sh_getNormType(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return (int)pData->norm;
}

float array2sh_getc(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->c;
}

float array2sh_getGain(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->gain_dB;
}

float array2sh_getMaxFreq(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->maxFreq;
}

float* array2sh_getFreqVector(void* const hA2sh, int* nFreqPoints)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    (*nFreqPoints) = HYBRID_BANDS;
    return &(pData->freqVector[0]);
}

float** array2sh_getbN_inv(void* const hA2sh, int* nCurves, int* nFreqPoints)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    (*nCurves) = pData->order+1;
    (*nFreqPoints) = HYBRID_BANDS;
    return pData->bN_inv_dB;
}

float** array2sh_getbN_modal(void* const hA2sh, int* nCurves, int* nFreqPoints)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    (*nCurves) = pData->order+1;
    (*nFreqPoints) = HYBRID_BANDS;
    return pData->bN_modal_dB;
}

float* array2sh_getSpatialCorrelation_Handle(void* const hA2sh, int* nCurves, int* nFreqPoints)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    (*nCurves) = pData->order+1;
    (*nFreqPoints) = HYBRID_BANDS;
    return pData->cSH;
}

float* array2sh_getLevelDifference_Handle(void* const hA2sh, int* nCurves, int* nFreqPoints)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    (*nCurves) = pData->order+1;
    (*nFreqPoints) = HYBRID_BANDS;
    return pData->lSH;
}

int array2sh_getSamplingRate(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->fs;
}

int array2sh_getProcessingDelay()
{
    return 12*HOP_SIZE;
}










    
    
