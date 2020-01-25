/*
 * Copyright 2018 Leo McCormack
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

/*
 * Filename: ambi_bin.c
 * --------------------
 * A binaural Ambisonic decoder for reproducing ambisonic signals over
 * headphones. The decoder supports sound-field rotation for head-tracking and
 * may also accomodate custom HRIR sets via the SOFA standard.
 *
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_hrir, saf_vbap, saf_sh, saf_sofa_reader
 * Author, date created:
 *     Leo McCormack, 14.04.2018
 */
 
#include "ambi_bin_internal.h"

void ambi_bin_create
(
    void ** const phAmbi
)
{
    ambi_bin_data* pData = (ambi_bin_data*)malloc1d(sizeof(ambi_bin_data));
    *phAmbi = (void*)pData;
    int ch, band;

    /* default user parameters */
    for (band = 0; band<HYBRID_BANDS; band++)
        pData->EQ[band] = 1.0f;
    pData->useDefaultHRIRsFLAG = 1; /* pars->sofa_filepath must be valid to set this to 0 */
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    pData->enableMaxRE = 1;
    pData->enableDiffuseMatching = 1;
    pData->enablePhaseWarping = 0;
    pData->enableRotation = 0;
    pData->yaw = 0.0f;
    pData->pitch = 0.0f;
    pData->roll = 0.0f;
    pData->bFlipYaw = 0;
    pData->bFlipPitch = 0;
    pData->bFlipRoll = 0;
    pData->useRollPitchYawFlag = 0;
    pData->method = DECODING_METHOD_MAGLS;
    pData->order = pData->new_order = 1;
    pData->nSH =  (pData->order+1)*(pData->order+1);  
    
    /* afSTFT stuff */
    pData->hSTFT = NULL;
    pData->STFTOutputFrameTF = malloc1d(NUM_EARS * sizeof(complexVector));
    for(ch=0; ch< NUM_EARS; ch++) {
        pData->STFTOutputFrameTF[ch].re = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
        pData->STFTOutputFrameTF[ch].im = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
    }
    pData->tempHopFrameTD = (float**)malloc2d( MAX(MAX_NUM_SH_SIGNALS, NUM_EARS), HOP_SIZE, sizeof(float));
    pData->STFTInputFrameTF = malloc1d(MAX_NUM_SH_SIGNALS * sizeof(complexVector));
    for(ch=0; ch< MAX_NUM_SH_SIGNALS; ch++) {
        pData->STFTInputFrameTF[ch].re = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
        pData->STFTInputFrameTF[ch].im = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
    }

    /* codec data */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(AMBI_BIN_PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"HOSIRR");
    pData->pars = (codecPars*)malloc1d(sizeof(codecPars));
    codecPars* pars = pData->pars; 
    pars->sofa_filepath = NULL;
    pars->hrirs = NULL;
    pars->hrir_dirs_deg = NULL;
    pars->itds_s = NULL;
    pars->hrtf_fb = NULL;
    
    /* flags */
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->recalc_M_rotFLAG = 1;
}

void ambi_bin_destroy
(
    void ** const phAmbi
)
{
    ambi_bin_data *pData = (ambi_bin_data*)(*phAmbi);
    codecPars *pars = pData->pars;
    int ch;
    
    if (pData != NULL) {
        if(pData->hSTFT!=NULL)
            afSTFTfree(pData->hSTFT);
        if(pData->STFTInputFrameTF!=NULL){
            for (ch = 0; ch< MAX_NUM_SH_SIGNALS; ch++) {
                free(pData->STFTInputFrameTF[ch].re);
                free(pData->STFTInputFrameTF[ch].im);
            }
        }
        for (ch = 0; ch< NUM_EARS; ch++) {
            free(pData->STFTOutputFrameTF[ch].re);
            free(pData->STFTOutputFrameTF[ch].im);
        }
        free(pData->STFTInputFrameTF);
        free(pData->STFTOutputFrameTF);
        free(pData->tempHopFrameTD);
        free1d((void**)&(pars->hrtf_fb));
        free1d((void**)&(pars->itds_s));
        free1d((void**)&(pars->hrirs));
        free1d((void**)&(pars->hrir_dirs_deg));
        free(pars);
        free(pData->progressBarText);
        free(pData);
        pData = NULL;
    }
}

void ambi_bin_init
(
    void * const hAmbi,
    int          sampleRate
)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    int band;
    
    /* define frequency vector */
    pData->fs = sampleRate;
    for(band=0; band <HYBRID_BANDS; band++){
        if(sampleRate == 44100)
            pData->freqVector[band] =  (float)__afCenterFreq44100[band];
        else /* Assume 48kHz */
            pData->freqVector[band] =  (float)__afCenterFreq48e3[band];
    }
    
    /* reinitialise codec if needed */
    if (pData->codecStatus == CODEC_STATUS_NOT_INITIALISED)
        ambi_bin_initCodec(hAmbi);
 
    /* default starting values */
    memset(pData->M_rot, 0, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*sizeof(float_complex));
    pData->recalc_M_rotFLAG = 1;
}

void ambi_bin_initCodec
(
    void* const hAmbi
)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    codecPars* pars = pData->pars;
    int i, j, nSH, order, band;
    
    if (pData->codecStatus != CODEC_STATUS_NOT_INITIALISED)
        return; /* re-init not required */
    if (pData->procStatus == PROC_STATUS_ONGOING){
        /* re-init required, but need to wait for processing loop to end */
        pData->codecStatus = CODEC_STATUS_INITIALISING;
        ambi_bin_initCodec(hAmbi);
    }
    
    /* for progress bar */
    pData->codecStatus = CODEC_STATUS_INITIALISING;
    strcpy(pData->progressBarText,"Preparing HRIRs");
    pData->progressBar0_1 = 0.0f;
    
    /* (Re)Initialise afSTFT */
    order = pData->new_order;
    nSH = (order+1)*(order+1);
    if(pData->hSTFT==NULL)
        afSTFTinit(&(pData->hSTFT), HOP_SIZE, nSH, NUM_EARS, 0, 1);
    else if(pData->nSH != nSH) {/* Or change the number of channels */
        afSTFTchannelChange(pData->hSTFT, nSH, NUM_EARS);
        afSTFTclearBuffers(pData->hSTFT);
    }
    pData->nSH = nSH;
    
    /* load sofa file or default hrir data */
    strcpy(pData->progressBarText,"Preparing HRIRs");
    pData->progressBar0_1 = 0.15f;
    if(!pData->useDefaultHRIRsFLAG && pars->sofa_filepath!=NULL){
        loadSofaFile(pars->sofa_filepath,
                     &(pars->hrirs),
                     &(pars->hrir_dirs_deg),
                     &(pars->N_hrir_dirs),
                     &(pars->hrir_len),
                     &(pars->hrir_fs));
    }
    else{
        loadSofaFile(NULL, /* setting path to NULL loads default HRIR data */
                     &(pars->hrirs),
                     &(pars->hrir_dirs_deg),
                     &(pars->N_hrir_dirs),
                     &(pars->hrir_len),
                     &(pars->hrir_fs));
    }
    
    /* estimate the ITDs for each HRIR */
    pData->progressBar0_1 = 0.3f;
    pars->itds_s = realloc1d(pars->itds_s, pars->N_hrir_dirs*sizeof(float));
    estimateITDs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->hrir_fs, pars->itds_s);
    
    /* convert hrirs to filterbank coefficients */
    pData->progressBar0_1 = 0.9f;
    pars->hrtf_fb = realloc1d(pars->hrtf_fb, HYBRID_BANDS * NUM_EARS * (pars->N_hrir_dirs)*sizeof(float_complex));
    HRIRs2FilterbankHRTFs(pars->hrirs, pars->N_hrir_dirs, pars->hrir_len, pars->hrtf_fb);
    diffuseFieldEqualiseHRTFs(pars->N_hrir_dirs, pars->itds_s, pData->freqVector, HYBRID_BANDS, pars->hrtf_fb);
    
    /* get new decoder */
    strcpy(pData->progressBarText,"Computing Decoder");
    pData->progressBar0_1 = 0.95f;
    float_complex* decMtx;
    decMtx = calloc1d(HYBRID_BANDS*NUM_EARS*nSH, sizeof(float_complex));
    switch(pData->method){
        default:
        case DECODING_METHOD_LS:
            getBinauralAmbiDecoderMtx(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS,
                                      BINAURAL_DECODER_LS, order, pData->freqVector, pars->itds_s, NULL,
                                      pData->enableDiffuseMatching, pData->enableMaxRE, decMtx);
            break;
        case DECODING_METHOD_LSDIFFEQ:
            getBinauralAmbiDecoderMtx(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS,
                                      BINAURAL_DECODER_LSDIFFEQ, order, pData->freqVector, pars->itds_s, NULL,
                                      pData->enableDiffuseMatching, pData->enableMaxRE, decMtx);
            break;
        case DECODING_METHOD_SPR:
            getBinauralAmbiDecoderMtx(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS,
                                      BINAURAL_DECODER_SPR, order, pData->freqVector, pars->itds_s, NULL,
                                      pData->enableDiffuseMatching, pData->enableMaxRE, decMtx);
            break;
        case DECODING_METHOD_TA:
            getBinauralAmbiDecoderMtx(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS,
                                      BINAURAL_DECODER_TA, order, pData->freqVector, pars->itds_s, NULL,
                                      pData->enableDiffuseMatching, pData->enableMaxRE, decMtx);
            break;
        case DECODING_METHOD_MAGLS:
            getBinauralAmbiDecoderMtx(pars->hrtf_fb, pars->hrir_dirs_deg, pars->N_hrir_dirs, HYBRID_BANDS,
                                      BINAURAL_DECODER_MAGLS, order, pData->freqVector, pars->itds_s, NULL,
                                      pData->enableDiffuseMatching, pData->enableMaxRE, decMtx);
            break;
    }
    
    /* Apply Phase Warping */
    if(pData->enablePhaseWarping){
        // COMING SOON
    }
    
    /* replace current decoder */
    memset(pars->M_dec, 0, HYBRID_BANDS*NUM_EARS*MAX_NUM_SH_SIGNALS*sizeof(float_complex));
    for(band=0; band<HYBRID_BANDS; band++)
        for(i=0; i<NUM_EARS; i++)
            for(j=0; j<nSH; j++)
                pars->M_dec[band][i][j] = decMtx[band*2*nSH + i*nSH + j];
    free(decMtx);
    
    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
    pData->order = order;
}

void ambi_bin_process
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
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    codecPars* pars = pData->pars;
    int n, t, ch, i, j, band;
    int o[MAX_SH_ORDER+2];
    const float_complex calpha = cmplxf(1.0f,0.0f), cbeta = cmplxf(0.0f, 0.0f);
    float Rxyz[3][3];
    float* M_rot_tmp;
    
    /* local copies of user parameters */
    int order, nSH, enableRot;
    NORM_TYPES norm;
    CH_ORDER chOrdering;
  
    /* decode audio to loudspeakers or headphones */
    if ( (nSamples == FRAME_SIZE) && (pData->codecStatus == CODEC_STATUS_INITIALISED) ){
        pData->procStatus = PROC_STATUS_ONGOING;
        
        /* copy user parameters to local variables */
        for(n=0; n<MAX_SH_ORDER+2; n++){  o[n] = n*n;  }
        norm = pData->norm;
        chOrdering = pData->chOrdering;
        order = pData->order;
        nSH = (order+1)*(order+1);
        enableRot = pData->enableRotation;
        
        /* Load time-domain data */
        switch(chOrdering){
            case CH_ACN:
                for(i=0; i < MIN(nSH, nInputs); i++)
                    utility_svvcopy(inputs[i], FRAME_SIZE, pData->SHFrameTD[i]);
                for(; i<nSH; i++)
                    memset(pData->SHFrameTD[i], 0, FRAME_SIZE * sizeof(float)); /* fill remaining channels with zeros */
                break;
            case CH_FUMA:   /* only for first-order, convert to ACN */
                if(nInputs>=4){
                    utility_svvcopy(inputs[0], FRAME_SIZE, pData->SHFrameTD[0]);
                    utility_svvcopy(inputs[1], FRAME_SIZE, pData->SHFrameTD[3]);
                    utility_svvcopy(inputs[2], FRAME_SIZE, pData->SHFrameTD[1]);
                    utility_svvcopy(inputs[3], FRAME_SIZE, pData->SHFrameTD[2]);
                    for(i=4; i<nSH; i++)
                        memset(pData->SHFrameTD[i], 0, FRAME_SIZE * sizeof(float)); /* fill remaining channels with zeros */
                }
                else
                    for(i=0; i<nSH; i++)
                        memset(pData->SHFrameTD[i], 0, FRAME_SIZE * sizeof(float));  
                break;
        }
        
        /* account for input normalisation scheme */
        switch(norm){
            case NORM_N3D:  /* already in N3D, do nothing */
                break;
            case NORM_SN3D: /* convert to N3D */
                for (n = 0; n<order+1; n++)
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
        
        /* Apply time-frequency transform (TFT) */
        for(t=0; t< TIME_SLOTS; t++) {
            for(ch = 0; ch < nSH; ch++)
                utility_svvcopy(&(pData->SHFrameTD[ch][t*HOP_SIZE]), HOP_SIZE, pData->tempHopFrameTD[ch]);
            afSTFTforward(pData->hSTFT, pData->tempHopFrameTD, pData->STFTInputFrameTF);
            for(band=0; band<HYBRID_BANDS; band++)
                for(ch=0; ch < nSH; ch++)
                    pData->SHframeTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[ch].re[band], pData->STFTInputFrameTF[ch].im[band]);
        }
    
        /* Main processing: */
        if(isPlaying){
            /* Apply rotation */
            if(order > 0 && enableRot) {
                if(pData->recalc_M_rotFLAG){
                    memset(pData->M_rot, 0, MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*sizeof(float_complex));
                    M_rot_tmp = malloc1d(nSH*nSH * sizeof(float));
                    yawPitchRoll2Rzyx(pData->yaw, pData->pitch, pData->roll, pData->useRollPitchYawFlag, Rxyz);
                    getSHrotMtxReal(Rxyz, M_rot_tmp, order);
                    for (i = 0; i < nSH; i++)
                        for (j = 0; j < nSH; j++)
                            pData->M_rot[i][j] = cmplxf(M_rot_tmp[i*nSH + j], 0.0f);
                    free(M_rot_tmp);
                    pData->recalc_M_rotFLAG = 0;
                }
                for(band = 0; band < HYBRID_BANDS; band++) {
                    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, TIME_SLOTS, nSH, &calpha,
                                pData->M_rot, MAX_NUM_SH_SIGNALS,
                                pData->SHframeTF[band], TIME_SLOTS, &cbeta,
                                pData->SHframeTF_rot[band], TIME_SLOTS);
                }
            }
            else
                utility_cvvcopy((const float_complex*)pData->SHframeTF, HYBRID_BANDS*MAX_NUM_SH_SIGNALS*TIME_SLOTS, (float_complex*)pData->SHframeTF_rot);
            
            /* mix to headphones */
            for(band = 0; band < HYBRID_BANDS; band++) {
                cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, TIME_SLOTS, nSH, &calpha,
                            pars->M_dec[band], MAX_NUM_SH_SIGNALS,
                            pData->SHframeTF_rot[band], TIME_SLOTS, &cbeta,
                            pData->binframeTF[band], TIME_SLOTS);
            }
        }
        else
            memset(pData->binframeTF, 0, HYBRID_BANDS*NUM_EARS*TIME_SLOTS*sizeof(float_complex));
        
        /* inverse-TFT */
        //postGain = powf(10.0f, POST_GAIN/20.0f);
        for(t = 0; t < TIME_SLOTS; t++) {
            for(band = 0; band < HYBRID_BANDS; band++) {
                for(ch = 0; ch < NUM_EARS; ch++) {
                    pData->STFTOutputFrameTF[ch].re[band] = crealf(pData->binframeTF[band][ch][t]);
                    pData->STFTOutputFrameTF[ch].im[band] = cimagf(pData->binframeTF[band][ch][t]);
                }
            }
            afSTFTinverse(pData->hSTFT, pData->STFTOutputFrameTF, pData->tempHopFrameTD);
            for (ch = 0; ch < MIN(NUM_EARS, nOutputs); ch++)
                utility_svvcopy(pData->tempHopFrameTD[ch], HOP_SIZE, &(outputs[ch][t* HOP_SIZE]));
            for (; ch < nOutputs; ch++)
                memset(&(outputs[ch][t* HOP_SIZE]), 0, HOP_SIZE*sizeof(float));
        }
    }
    else
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, FRAME_SIZE*sizeof(float));
    
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}


/* Set Functions */

void ambi_bin_refreshParams(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
}

void ambi_bin_setUseDefaultHRIRsflag(void* const hAmbi, int newState)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    
    if((!pData->useDefaultHRIRsFLAG) && (newState)){
        pData->useDefaultHRIRsFLAG = newState;
        pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    }
}

void ambi_bin_setSofaFilePath(void* const hAmbi, const char* path)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    codecPars* pars = pData->pars;
    
    pars->sofa_filepath = malloc1d(strlen(path) + 1);
    strcpy(pars->sofa_filepath, path);
    pData->useDefaultHRIRsFLAG = 0;
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
}

void ambi_bin_setInputOrderPreset(void* const hAmbi, INPUT_ORDERS newOrder)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    if(pData->order != (int)newOrder){
        pData->new_order = (int)newOrder;
        pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
        
        /* FUMA only supports 1st order */
        if(pData->new_order!=INPUT_ORDER_FIRST && pData->chOrdering == CH_FUMA)
            pData->chOrdering = CH_ACN;
        if(pData->new_order!=INPUT_ORDER_FIRST && pData->norm == NORM_FUMA)
            pData->norm = NORM_SN3D;
    }
}

void ambi_bin_setDecodingMethod(void* const hAmbi, DECODING_METHODS newMethod)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    pData->method = newMethod;
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
}

void ambi_bin_setChOrder(void* const hAmbi, int newOrder)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    if((CH_ORDER)newOrder != CH_FUMA || pData->new_order==INPUT_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void ambi_bin_setNormType(void* const hAmbi, int newType)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    if((NORM_TYPES)newType != NORM_FUMA || pData->new_order==INPUT_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
}

void ambi_bin_setEnableMaxRE(void* const hAmbi, int newState)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    if(pData->enableMaxRE != newState){
        pData->enableMaxRE = newState;
        pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    }
}

void ambi_bin_setEnableDiffuseMatching(void* const hAmbi, int newState)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    if(pData->enableDiffuseMatching != newState){
        pData->enableDiffuseMatching = newState;
        pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    }
}

void ambi_bin_setEnablePhaseWarping(void* const hAmbi, int newState)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    if(pData->enablePhaseWarping != newState){
        pData->enablePhaseWarping = newState;
        pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    }
}

void ambi_bin_setEnableRotation(void* const hAmbi, int newState)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    pData->enableRotation = newState;
}

void ambi_bin_setYaw(void  * const hAmbi, float newYaw)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    pData->yaw = pData->bFlipYaw == 1 ? -DEG2RAD(newYaw) : DEG2RAD(newYaw);
    pData->recalc_M_rotFLAG = 1;
}

void ambi_bin_setPitch(void* const hAmbi, float newPitch)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    pData->pitch = pData->bFlipPitch == 1 ? -DEG2RAD(newPitch) : DEG2RAD(newPitch);
    pData->recalc_M_rotFLAG = 1;
}

void ambi_bin_setRoll(void* const hAmbi, float newRoll)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    pData->roll = pData->bFlipRoll == 1 ? -DEG2RAD(newRoll) : DEG2RAD(newRoll);
    pData->recalc_M_rotFLAG = 1;
}

void ambi_bin_setFlipYaw(void* const hAmbi, int newState)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    if(newState !=pData->bFlipYaw ){
        pData->bFlipYaw = newState;
        ambi_bin_setYaw(hAmbi, -ambi_bin_getYaw(hAmbi));
    }
}

void ambi_bin_setFlipPitch(void* const hAmbi, int newState)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    if(newState !=pData->bFlipPitch ){
        pData->bFlipPitch = newState;
        ambi_bin_setPitch(hAmbi, -ambi_bin_getPitch(hAmbi));
    }
}

void ambi_bin_setFlipRoll(void* const hAmbi, int newState)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    if(newState !=pData->bFlipRoll ){
        pData->bFlipRoll = newState;
        ambi_bin_setRoll(hAmbi, -ambi_bin_getRoll(hAmbi));
    }
}

void ambi_bin_setRPYflag(void* const hAmbi, int newState)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    pData->useRollPitchYawFlag = newState;
}


/* Get Functions */

CODEC_STATUS ambi_bin_getCodecStatus(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->codecStatus;
}

float ambi_bin_getProgressBar0_1(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->progressBar0_1;
}

void ambi_bin_getProgressBarText(void* const hAmbi, char* text)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    memcpy(text, pData->progressBarText, AMBI_BIN_PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

int ambi_bin_getUseDefaultHRIRsflag(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->useDefaultHRIRsFLAG;
}

int ambi_bin_getInputOrderPreset(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->order;
}

int ambi_bin_getDecodingMethod(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->method;
}

char* ambi_bin_getSofaFilePath(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    codecPars* pars = pData->pars;
    if(pars->sofa_filepath!=NULL)
        return pars->sofa_filepath;
    else
        return "no_file";
}

int ambi_bin_getChOrder(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return (int)pData->chOrdering;
}

int ambi_bin_getNormType(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return (int)pData->norm;
}

int ambi_bin_getEnableMaxRE(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->enableMaxRE;
}

int ambi_bin_getEnableDiffuseMatching(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->enableDiffuseMatching;
}

int ambi_bin_getEnablePhaseWarping(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->enablePhaseWarping;
}

int ambi_bin_getNumEars()
{ 
    return NUM_EARS;
}

int ambi_bin_getNSHrequired(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->nSH;
}

int ambi_bin_getEnableRotation(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->enableRotation;
}

float ambi_bin_getYaw(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->bFlipYaw == 1 ? -RAD2DEG(pData->yaw) : RAD2DEG(pData->yaw);
}

float ambi_bin_getPitch(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->bFlipPitch == 1 ? -RAD2DEG(pData->pitch) : RAD2DEG(pData->pitch);
}

float ambi_bin_getRoll(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->bFlipRoll == 1 ? -RAD2DEG(pData->roll) : RAD2DEG(pData->roll);
}

int ambi_bin_getFlipYaw(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->bFlipYaw;
}

int ambi_bin_getFlipPitch(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->bFlipPitch;
}

int ambi_bin_getFlipRoll(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->bFlipRoll;
}

int ambi_bin_getRPYflag(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->useRollPitchYawFlag;
}

int ambi_bin_getNDirs(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    codecPars* pars = pData->pars;
    return pars->N_hrir_dirs;
}

int ambi_bin_getHRIRlength(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    codecPars* pars = pData->pars;
    return pars->hrir_len;
}

int ambi_bin_getHRIRsamplerate(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    codecPars* pars = pData->pars;
    return pars->hrir_fs;
}

int ambi_bin_getDAWsamplerate(void* const hAmbi)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    return pData->fs;
}

int ambi_bin_getProcessingDelay()
{
    return 12*HOP_SIZE;
}
