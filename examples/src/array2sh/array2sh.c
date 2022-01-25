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
 * @file array2sh.c
 * @brief Spatially encodes spherical microphone array signals into spherical
 *        harmonic signals (aka: Ambisonic signals) utilising theoretical
 *        encoding filters.
 *
 * The algorithms within array2sh were pieced together and developed in
 * collaboration with Symeon Delikaris-Manias and Angelo Farina.
 * A detailed explanation of the algorithms within array2sh can be found in [1].
 * Also included, is a diffuse-field equalisation option for frequencies past
 * aliasing, developed in collaboration with Archontis Politis, 8.02.2019
 *
 * @note Since the algorithms are based on theory, only array designs where
 *       there are analytical solutions available are supported. i.e. only
 *       spherical or cylindrical arrays, which have phase-matched sensors.
 *       For more information, the reader is referred to [2,3].
 * @test test__saf_example_array2sh()
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., "Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis," in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 * @see [2] Williams EG. Fourier acoustics: sound radiation and nearfield
 *          acoustical holography. Elsevier; 1999 Jun 10.
 * @see [3] Rafaely B. Fundamentals of spherical array processing. Berlin:
 *          Springer; 2015 Feb 18.
 *
 * @author Leo McCormack
 * @date 13.09.2017
 * @license ISC
 */

#include "array2sh_internal.h" 

void array2sh_create
(
    void ** const phA2sh
)
{
    array2sh_data* pData = (array2sh_data*)malloc1d(sizeof(array2sh_data));
    *phA2sh = (void*)pData;

    /* defualt parameters */
    array2sh_createArray(&(pData->arraySpecs)); 
    pData->filterType = FILTER_TIKHONOV;
    pData->regPar = 15.0f;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    pData->c = 343.0f;
    pData->gain_dB = 0.0f; /* post-gain */ 
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    array2sh_initArray(arraySpecs, MICROPHONE_ARRAY_PRESET_DEFAULT, &(pData->order), 1);
    pData->enableDiffEQpastAliasing = 1;
    
    /* time-frequency transform + buffers */
    pData->hSTFT = NULL;
    pData->inputFrameTD = (float**)malloc2d(MAX_NUM_SENSORS, ARRAY2SH_FRAME_SIZE, sizeof(float));
    pData->SHframeTD = (float**)malloc2d(MAX_NUM_SH_SIGNALS, ARRAY2SH_FRAME_SIZE, sizeof(float));
    pData->inputframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_SENSORS, TIME_SLOTS, sizeof(float_complex));
    pData->SHframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_SH_SIGNALS, TIME_SLOTS, sizeof(float_complex));

    /* internal */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    pData->evalStatus = EVAL_STATUS_NOT_EVALUATED;
    pData->evalRequestedFLAG = 0;
    pData->reinitSHTmatrixFLAG = 1;
    pData->new_order = pData->order;
    pData->bN = NULL;
    
    /* display related stuff */
    pData->bN_modal_dB = (float**)calloc2d(HYBRID_BANDS, MAX_SH_ORDER + 1, sizeof(float));
    pData->bN_inv_dB = (float**)calloc2d(HYBRID_BANDS, MAX_SH_ORDER + 1, sizeof(float));
    pData->cSH = (float*)calloc1d((HYBRID_BANDS)*(MAX_SH_ORDER + 1),sizeof(float));
    pData->lSH = (float*)calloc1d((HYBRID_BANDS)*(MAX_SH_ORDER + 1),sizeof(float));
}

void array2sh_destroy
(
    void ** const phM2sh
)
{
    array2sh_data *pData = (array2sh_data*)(*phM2sh);

    if (pData != NULL) {
        /* not safe to free memory during evaluation */
        while (pData->evalStatus == EVAL_STATUS_EVALUATING)
            SAF_SLEEP(10);
        
        /* free afSTFT and buffers */
        if (pData->hSTFT != NULL)
            afSTFT_destroy(&(pData->hSTFT));
        free(pData->inputFrameTD);
        free(pData->SHframeTD);
        free(pData->inputframeTF);
        free(pData->SHframeTF);
        array2sh_destroyArray(&(pData->arraySpecs));
        
        /* Display stuff */
        free(pData->bN_modal_dB);
        free(pData->bN_inv_dB);
        
        free(pData->progressBarText);

        free(pData->bN);
        free(pData->cSH);
        free(pData->lSH);
        
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
    
    pData->fs = sampleRate;
    afSTFT_getCentreFreqs(pData->hSTFT, (float)pData->fs, HYBRID_BANDS, pData->freqVector);
    pData->freqVector[0] = pData->freqVector[1]/4.0f; /* avoids NaNs at DC */
}

void array2sh_evalEncoder
(
    void* const hA2sh
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    
    if (pData->evalStatus != EVAL_STATUS_NOT_EVALUATED)
        return; /* eval not required */
    
    /* for progress bar */
    pData->evalStatus = EVAL_STATUS_EVALUATING;
    strcpy(pData->progressBarText,"Initialising evaluation");
    pData->progressBar0_1 = 0.0f;
    
    /* Evaluate Encoder */
    array2sh_evaluateSHTfilters(hA2sh);
    
    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f; 
    pData->evalStatus = EVAL_STATUS_RECENTLY_EVALUATED;
}

void array2sh_process
(
    void        *  const hA2sh,
    const float *const * inputs,
    float       ** const outputs,
    int                  nInputs,
    int                  nOutputs,
    int                  nSamples
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    int ch, i, band, Q, order, nSH;
    const float_complex calpha = cmplxf(1.0f,0.0f), cbeta = cmplxf(0.0f, 0.0f);
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    float gain_lin;
    
    /* reinit TFT if needed */
    array2sh_initTFT(hA2sh);
    
    /* compute encoding matrix if needed */
    if (pData->reinitSHTmatrixFLAG) {
        array2sh_calculate_sht_matrix(hA2sh); /* compute encoding matrix */
        array2sh_calculate_mag_curves(hA2sh); /* calculate magnitude response curves */
        pData->reinitSHTmatrixFLAG = 0;
    }

    /* local copy of user parameters */
    chOrdering = pData->chOrdering;
    norm = pData->norm;
    gain_lin = powf(10.0f, pData->gain_dB/20.0f);
    Q = arraySpecs->Q;
    order = pData->order;
    nSH = (order+1)*(order+1);

    /* processing loop */
    if ((nSamples == ARRAY2SH_FRAME_SIZE) && (pData->reinitSHTmatrixFLAG==0) ) {
        pData->procStatus = PROC_STATUS_ONGOING;

        /* Load time-domain data */
        for(i=0; i < nInputs; i++)
            utility_svvcopy(inputs[i], ARRAY2SH_FRAME_SIZE, pData->inputFrameTD[i]);
        for(; i<Q; i++)
            memset(pData->inputFrameTD[i], 0, ARRAY2SH_FRAME_SIZE * sizeof(float));

        /* Apply time-frequency transform (TFT) */
        afSTFT_forward_knownDimensions(pData->hSTFT, pData->inputFrameTD, ARRAY2SH_FRAME_SIZE, MAX_NUM_SENSORS, TIME_SLOTS, pData->inputframeTF);

        /* Apply spherical harmonic transform (SHT) */
        for(band=0; band<HYBRID_BANDS; band++){
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, TIME_SLOTS, Q, &calpha,
                        pData->W[band], MAX_NUM_SENSORS,
                        FLATTEN2D(pData->inputframeTF[band]), TIME_SLOTS, &cbeta,
                        FLATTEN2D(pData->SHframeTF[band]), TIME_SLOTS);
        }

        /* inverse-TFT */
        afSTFT_backward_knownDimensions(pData->hSTFT, pData->SHframeTF, ARRAY2SH_FRAME_SIZE, MAX_NUM_SH_SIGNALS, TIME_SLOTS, pData->SHframeTD);

        /* account for output channel order */
        switch(chOrdering){
            case CH_ACN:  /* already ACN, do nothing */ break;
            case CH_FUMA: convertHOAChannelConvention(FLATTEN2D(pData->SHframeTD), order, ARRAY2SH_FRAME_SIZE, HOA_CH_ORDER_ACN, HOA_CH_ORDER_FUMA); break;
        }

        /* account for normalisation scheme */
        switch(norm){
            case NORM_N3D:  /* already N3D, do nothing */ break;  
            case NORM_SN3D: convertHOANormConvention(FLATTEN2D(pData->SHframeTD), order, ARRAY2SH_FRAME_SIZE, HOA_NORM_N3D, HOA_NORM_SN3D); break;
            case NORM_FUMA: convertHOANormConvention(FLATTEN2D(pData->SHframeTD), order, ARRAY2SH_FRAME_SIZE, HOA_NORM_N3D, HOA_NORM_FUMA); break;
        }

        /* Apply post-gain */
        utility_svsmul(FLATTEN2D(pData->SHframeTD), &gain_lin, nSH*ARRAY2SH_FRAME_SIZE, NULL);

        /* Copy to output */
        for(i = 0; i < SAF_MIN(nSH,nOutputs); i++)
            utility_svvcopy(pData->SHframeTD[i], ARRAY2SH_FRAME_SIZE, outputs[i]);
        for(; i < nOutputs; i++)
            memset(outputs[i], 0, ARRAY2SH_FRAME_SIZE * sizeof(float));
    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, ARRAY2SH_FRAME_SIZE*sizeof(float));
    }

    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

/* Set Functions */

void array2sh_refreshSettings(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->reinitSHTmatrixFLAG = 1;
    array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
}

void array2sh_setEncodingOrder(void* const hA2sh, int newOrder)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    
    if(pData->new_order != newOrder){
        pData->new_order = newOrder;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
    /* FUMA only supports 1st order */
    if(pData->new_order!=SH_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->new_order!=SH_ORDER_FIRST && pData->norm == NORM_FUMA)
        pData->norm = NORM_SN3D;
}

void array2sh_setRequestEncoderEvalFLAG(void* const hA2sh, int newState)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->evalRequestedFLAG = newState;
}

void array2sh_setEvalStatus(void* const hA2sh, ARRAY2SH_EVAL_STATUS new_evalStatus)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    if(new_evalStatus==EVAL_STATUS_NOT_EVALUATED){
        /* Pause until current initialisation is complete */
        while(pData->evalStatus == EVAL_STATUS_EVALUATING)
            SAF_SLEEP(10);
    }
    pData->evalStatus = new_evalStatus;
}

void array2sh_setDiffEQpastAliasing(void* const hA2sh, int newState)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    if(pData->enableDiffEQpastAliasing != newState){
        pData->enableDiffEQpastAliasing = newState;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setPreset(void* const hA2sh, ARRAY2SH_MICROPHONE_ARRAY_PRESETS preset)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    
    array2sh_initArray(arraySpecs,(ARRAY2SH_MICROPHONE_ARRAY_PRESETS)preset, &(pData->new_order), 0);
    pData->c = (ARRAY2SH_MICROPHONE_ARRAY_PRESETS)preset == MICROPHONE_ARRAY_PRESET_AALTO_HYDROPHONE ? 1484.0f : 343.0f;
    pData->reinitSHTmatrixFLAG = 1;
    array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
}

void array2sh_setSensorAzi_rad(void* const hA2sh, int index, float newAzi_rad)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    
    if(arraySpecs->sensorCoords_rad[index][0] != newAzi_rad){
        arraySpecs->sensorCoords_rad[index][0] = newAzi_rad;
        arraySpecs->sensorCoords_deg[index][0] = newAzi_rad * (180.0f/SAF_PI);
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setSensorElev_rad(void* const hA2sh, int index, float newElev_rad)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    
    if(arraySpecs->sensorCoords_rad[index][1] != newElev_rad){
        arraySpecs->sensorCoords_rad[index][1] = newElev_rad;
        arraySpecs->sensorCoords_deg[index][1] = newElev_rad * (180.0f/SAF_PI);
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setSensorAzi_deg(void* const hA2sh, int index, float newAzi_deg)

{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    
    if(arraySpecs->sensorCoords_deg[index][0] != newAzi_deg){
        arraySpecs->sensorCoords_rad[index][0] = newAzi_deg * (SAF_PI/180.0f);
        arraySpecs->sensorCoords_deg[index][0] = newAzi_deg;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setSensorElev_deg(void* const hA2sh, int index, float newElev_deg)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    
    if(arraySpecs->sensorCoords_deg[index][1] != newElev_deg){
        arraySpecs->sensorCoords_rad[index][1] = newElev_deg * (SAF_PI/180.0f);
        arraySpecs->sensorCoords_deg[index][1] = newElev_deg;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setNumSensors(void* const hA2sh, int newQ)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    int nSH;
    
    nSH = (pData->new_order+1)*(pData->new_order+1);
    if (newQ < nSH){
        pData->new_order = 1;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
    if(arraySpecs->Q != newQ){
        arraySpecs->newQ = newQ;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setr(void* const hA2sh, float newr)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    
    newr = SAF_CLAMP(newr, ARRAY2SH_ARRAY_RADIUS_MIN_VALUE/1e3f, ARRAY2SH_ARRAY_RADIUS_MAX_VALUE/1e3f);
    if(arraySpecs->r!=newr){
        arraySpecs->r = newr;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setR(void* const hA2sh, float newR)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    
    newR = SAF_CLAMP(newR, ARRAY2SH_BAFFLE_RADIUS_MIN_VALUE/1e3f, ARRAY2SH_BAFFLE_RADIUS_MAX_VALUE/1e3f);
    if(arraySpecs->R!=newR){
        arraySpecs->R = newR;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setArrayType(void* const hA2sh, int newType)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    
    if(arraySpecs->arrayType != (ARRAY2SH_ARRAY_TYPES)newType){
        arraySpecs->arrayType = (ARRAY2SH_ARRAY_TYPES)newType;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setWeightType(void* const hA2sh, int newType)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    
    if(arraySpecs->weightType!=(ARRAY2SH_WEIGHT_TYPES)newType){
        arraySpecs->weightType = (ARRAY2SH_WEIGHT_TYPES)newType;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setFilterType(void* const hA2sh, int newType)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    
    if(pData->filterType!=(ARRAY2SH_FILTER_TYPES)newType){
        pData->filterType = (ARRAY2SH_FILTER_TYPES)newType;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setRegPar(void* const hA2sh, float newVal)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    newVal = SAF_CLAMP(newVal, ARRAY2SH_MAX_GAIN_MIN_VALUE, ARRAY2SH_MAX_GAIN_MAX_VALUE);
    if(pData->regPar!=newVal){
        pData->regPar = newVal;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setChOrder(void* const hA2sh, int newOrder)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    if((CH_ORDER)newOrder != CH_FUMA || pData->order==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void array2sh_setNormType(void* const hA2sh, int newType)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    if((NORM_TYPES)newType != NORM_FUMA || pData->order==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
}

void array2sh_setc(void* const hA2sh, float newc)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    newc = SAF_CLAMP(newc, ARRAY2SH_SPEED_OF_SOUND_MIN_VALUE, ARRAY2SH_SPEED_OF_SOUND_MAX_VALUE);
    if(newc!=pData->c){
        pData->c = newc;
        pData->reinitSHTmatrixFLAG = 1;
        array2sh_setEvalStatus(hA2sh, EVAL_STATUS_NOT_EVALUATED);
    }
}

void array2sh_setGain(void* const hA2sh, float newGain)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    pData->gain_dB = SAF_CLAMP(newGain, ARRAY2SH_POST_GAIN_MIN_VALUE, ARRAY2SH_POST_GAIN_MAX_VALUE);
}


/* Get Functions */

int array2sh_getFrameSize(void)
{
    return ARRAY2SH_FRAME_SIZE;
}

ARRAY2SH_EVAL_STATUS array2sh_getEvalStatus(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->evalStatus;
}

int array2sh_getReinitSHTmatrixFLAG(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->reinitSHTmatrixFLAG;
}

float array2sh_getProgressBar0_1(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->progressBar0_1;
}

void array2sh_getProgressBarText(void* const hA2sh, char* text)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    memcpy(text, pData->progressBarText, PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

int array2sh_getRequestEncoderEvalFLAG(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->evalRequestedFLAG;
}

int array2sh_getDiffEQpastAliasing(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->enableDiffEQpastAliasing;
}

int array2sh_getEncodingOrder(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return pData->new_order;
}

float array2sh_getSensorAzi_rad(void* const hA2sh, int index)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    return arraySpecs->sensorCoords_rad[index][0];
}

float array2sh_getSensorElev_rad(void* const hA2sh, int index)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    return arraySpecs->sensorCoords_rad[index][1];
}

float array2sh_getSensorAzi_deg(void* const hA2sh, int index)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    return arraySpecs->sensorCoords_deg[index][0];
}

float array2sh_getSensorElev_deg(void* const hA2sh, int index)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    return arraySpecs->sensorCoords_deg[index][1];
}

int array2sh_getNumSensors(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
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
    return (pData->new_order+1)*(pData->new_order+1);
}

int array2sh_getNSHrequired(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    return (pData->new_order+1)*(pData->new_order+1);
}

float array2sh_getr(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    return arraySpecs->r;
}

float array2sh_getR(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    return arraySpecs->R;
} 

int array2sh_getArrayType(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    return (int)arraySpecs->arrayType;
}

int array2sh_getWeightType(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
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
