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
 * @file: binauraliser.c
 * @brief Convolves input audio (up to 64 channels) with interpolated HRTFs in
 *        the time-frequency domain.
 *
 * The HRTFs are interpolated by applying amplitude-preserving VBAP gains to the
 * HRTF magnitude responses and inter-aural time differences (ITDs)
 * individually, before being re-combined. The example also allows the user to
 * specify an external SOFA file for the convolution, and rotations of the
 * source directions to accomodate head-tracking.
 *
 * @author Leo McCormack
 * @date 25.09.2017
 * @license ISC
 */

#include "binauraliser_internal.h"

void binauraliser_create
(
    void ** const phBin
)
{
    binauraliser_data* pData = (binauraliser_data*)malloc1d(sizeof(binauraliser_data));
    *phBin = (void*)pData;
    int ch, dummy;

    /* user parameters */
    binauraliser_loadPreset(SOURCE_CONFIG_PRESET_DEFAULT, pData->src_dirs_deg, &(pData->new_nSources), &(dummy)); /*check setStateInformation if you change default preset*/
    pData->useDefaultHRIRsFLAG = 1; /* pars->sofa_filepath must be valid to set this to 0 */
    pData->enableHRIRsDiffuseEQ = 1;
    pData->nSources = pData->new_nSources;
    pData->interpMode = INTERP_TRI;
    pData->yaw = 0.0f;
    pData->pitch = 0.0f;
    pData->roll = 0.0f;
    pData->bFlipYaw = 0;
    pData->bFlipPitch = 0;
    pData->bFlipRoll = 0;
    pData->useRollPitchYawFlag = 0;
    pData->enableRotation = 0;

    /* time-frequency transform + buffers */
    pData->hSTFT = NULL;
    pData->inputFrameTD = (float**)malloc2d(MAX_NUM_INPUTS, BINAURALISER_FRAME_SIZE, sizeof(float));
    pData->outframeTD = (float**)malloc2d(NUM_EARS, BINAURALISER_FRAME_SIZE, sizeof(float));
    pData->inputframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_INPUTS, TIME_SLOTS, sizeof(float_complex));
    pData->outputframeTF = (float_complex***)malloc3d(HYBRID_BANDS, NUM_EARS, TIME_SLOTS, sizeof(float_complex));
    
    /* hrir data */
    pData->hrirs = NULL;
    pData->hrir_dirs_deg = NULL;
    pData->sofa_filepath = NULL;
    pData->weights = NULL;
    pData->N_hrir_dirs = pData->hrir_loaded_len = pData->hrir_runtime_len = 0;
    pData->hrir_loaded_fs = pData->hrir_runtime_fs = -1; /* unknown */
    
    /* vbap (amplitude normalised) */
    pData->hrtf_vbap_gtableIdx = NULL;
    pData->hrtf_vbap_gtableComp = NULL;
    pData->nTriangles = pData->N_hrtf_vbap_gtable = 0;
    
    /* HRTF filterbank coefficients */
    pData->itds_s = NULL;
    pData->hrtf_fb = NULL;
    pData->hrtf_fb_mag = NULL;

    /* flags/status */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
    pData->reInitHRTFsAndGainTables = 1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++) {
        pData->recalc_hrtf_interpFLAG[ch] = 1;
        pData->src_gains[ch] = 1.f;
    }
    pData->recalc_M_rotFLAG = 1; 
}

void binauraliser_destroy
(
    void ** const phBin
)
{
    binauraliser_data *pData = (binauraliser_data*)(*phBin);

    if (pData != NULL) {
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            SAF_SLEEP(10);
        }
        
        /* free afSTFT and buffers */
        if(pData->hSTFT !=NULL)
            afSTFT_destroy(&(pData->hSTFT));
        free(pData->inputFrameTD);
        free(pData->outframeTD);
        free(pData->inputframeTF);
        free(pData->outputframeTF);
        free(pData->hrtf_vbap_gtableComp);
        free(pData->hrtf_vbap_gtableIdx);
        free(pData->hrtf_fb);
        free(pData->hrtf_fb_mag);
        free(pData->itds_s);
        free(pData->sofa_filepath);
        free(pData->hrirs);
        free(pData->hrir_dirs_deg);
        free(pData->weights);
        free(pData->progressBarText);
         
        free(pData);
        pData = NULL;
    }
}

void binauraliser_init
(
    void * const hBin,
    int          sampleRate
)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    
    /* define frequency vector */
    pData->fs = sampleRate;
    afSTFT_getCentreFreqs(pData->hSTFT, (float)sampleRate, HYBRID_BANDS, pData->freqVector);
    if(pData->hrir_runtime_fs!=pData->fs){
        pData->reInitHRTFsAndGainTables = 1;
        binauraliser_setCodecStatus(hBin, CODEC_STATUS_NOT_INITIALISED);
    }

    /* defaults */
    pData->recalc_M_rotFLAG = 1;
}

void binauraliser_initCodec
(
    void* const hBin
)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    
    if (pData->codecStatus != CODEC_STATUS_NOT_INITIALISED)
        return; /* re-init not required, or already happening */
    while (pData->procStatus == PROC_STATUS_ONGOING){
        /* re-init required, but we need to wait for the current processing loop to end */
        pData->codecStatus = CODEC_STATUS_INITIALISING; /* indicate that we want to init */
        SAF_SLEEP(10);
    }
    
    /* for progress bar */
    pData->codecStatus = CODEC_STATUS_INITIALISING;
    strcpy(pData->progressBarText,"Initialising");
    pData->progressBar0_1 = 0.0f;
    
    /* check if TFT needs to be reinitialised */
    binauraliser_initTFT(hBin);
    
    /* reinit HRTFs and interpolation tables */
    if(pData->reInitHRTFsAndGainTables){
        binauraliser_initHRTFsAndGainTables(hBin);
        pData->reInitHRTFsAndGainTables = 0;
    }
    
    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
    
}

void binauraliser_process
(
    void        *  const hBin,
    const float *const * inputs,
    float       ** const outputs,
    int                  nInputs,
    int                  nOutputs,
    int                  nSamples
)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    int ch, ear, i, band, nSources;
    float Rxyz[3][3], hypotxy;
    int enableRotation;

    /* copy user parameters to local variables */
    nSources = pData->nSources;
    enableRotation = pData->enableRotation;
    
    /* apply binaural panner */
    if ((nSamples == BINAURALISER_FRAME_SIZE) && (pData->hrtf_fb!=NULL) && (pData->codecStatus==CODEC_STATUS_INITIALISED) ){
        pData->procStatus = PROC_STATUS_ONGOING;

        /* Load time-domain data */
        for(i=0; i < SAF_MIN(nSources,nInputs); i++)
            utility_svvcopy(inputs[i], BINAURALISER_FRAME_SIZE, pData->inputFrameTD[i]);
        for(; i<nSources; i++)
            memset(pData->inputFrameTD[i], 0, BINAURALISER_FRAME_SIZE * sizeof(float));

        /* Apply source gains */
        for (ch = 0; ch < nSources; ch++) {
            if(fabsf(pData->src_gains[ch] - 1.f) > 1e-6f)
                utility_svsmul(pData->inputFrameTD[ch], &(pData->src_gains[ch]), BINAURALISER_FRAME_SIZE, NULL);
        }

        /* Apply time-frequency transform (TFT) */
        afSTFT_forward_knownDimensions(pData->hSTFT, pData->inputFrameTD, BINAURALISER_FRAME_SIZE, MAX_NUM_INPUTS, TIME_SLOTS, pData->inputframeTF);

        /* Rotate source directions */
        if(enableRotation && pData->recalc_M_rotFLAG){
            yawPitchRoll2Rzyx (pData->yaw, pData->pitch, pData->roll, pData->useRollPitchYawFlag, Rxyz);
            for(i=0; i<nSources; i++){
                pData->src_dirs_xyz[i][0] = cosf(DEG2RAD(pData->src_dirs_deg[i][1])) * cosf(DEG2RAD(pData->src_dirs_deg[i][0]));
                pData->src_dirs_xyz[i][1] = cosf(DEG2RAD(pData->src_dirs_deg[i][1])) * sinf(DEG2RAD(pData->src_dirs_deg[i][0]));
                pData->src_dirs_xyz[i][2] = sinf(DEG2RAD(pData->src_dirs_deg[i][1]));
                pData->recalc_hrtf_interpFLAG[i] = 1;
            }
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSources, 3, 3, 1.0f,
                        (float*)(pData->src_dirs_xyz), 3,
                        (float*)Rxyz, 3, 0.0f,
                        (float*)(pData->src_dirs_rot_xyz), 3);
            for(i=0; i<nSources; i++){
                hypotxy = sqrtf(powf(pData->src_dirs_rot_xyz[i][0], 2.0f) + powf(pData->src_dirs_rot_xyz[i][1], 2.0f));
                pData->src_dirs_rot_deg[i][0] = RAD2DEG(atan2f(pData->src_dirs_rot_xyz[i][1], pData->src_dirs_rot_xyz[i][0]));
                pData->src_dirs_rot_deg[i][1] = RAD2DEG(atan2f(pData->src_dirs_rot_xyz[i][2], hypotxy));
            }
            pData->recalc_M_rotFLAG = 0;
        }

        /* interpolate hrtfs and apply to each source */
        memset(FLATTEN3D(pData->outputframeTF), 0, HYBRID_BANDS*NUM_EARS*TIME_SLOTS * sizeof(float_complex));
        for (ch = 0; ch < nSources; ch++) {
            if(pData->recalc_hrtf_interpFLAG[ch]){
                if(enableRotation)
                    binauraliser_interpHRTFs(hBin, pData->interpMode, pData->src_dirs_rot_deg[ch][0], pData->src_dirs_rot_deg[ch][1], pData->hrtf_interp[ch]);
                else
                    binauraliser_interpHRTFs(hBin, pData->interpMode, pData->src_dirs_deg[ch][0], pData->src_dirs_deg[ch][1], pData->hrtf_interp[ch]);
                pData->recalc_hrtf_interpFLAG[ch] = 0;
            }

            /* Convolve this channel with the interpolated HRTF, and add it to the binaural buffer */
            for (band = 0; band < HYBRID_BANDS; band++)
                for (ear = 0; ear < NUM_EARS; ear++)
                    cblas_caxpy(TIME_SLOTS, &pData->hrtf_interp[ch][band][ear], pData->inputframeTF[band][ch], 1, pData->outputframeTF[band][ear], 1);
        }

        /* scale by number of sources */ 
        cblas_sscal(/*re+im*/2*HYBRID_BANDS*NUM_EARS*TIME_SLOTS, 1.0f/sqrtf((float)nSources), (float*)FLATTEN3D(pData->outputframeTF), 1);

        /* inverse-TFT */
        afSTFT_backward_knownDimensions(pData->hSTFT, pData->outputframeTF, BINAURALISER_FRAME_SIZE, NUM_EARS, TIME_SLOTS, pData->outframeTD);

        /* Copy to output buffer */
        for (ch = 0; ch < SAF_MIN(NUM_EARS, nOutputs); ch++)
            utility_svvcopy(pData->outframeTD[ch], BINAURALISER_FRAME_SIZE, outputs[ch]);
        for (; ch < nOutputs; ch++)
            memset(outputs[ch], 0, BINAURALISER_FRAME_SIZE*sizeof(float));
    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, BINAURALISER_FRAME_SIZE*sizeof(float));
    }

    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

/* Set Functions */

void binauraliser_refreshSettings(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    int ch;
    pData->reInitHRTFsAndGainTables = 1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_hrtf_interpFLAG[ch] = 1;
    binauraliser_setCodecStatus(hBin, CODEC_STATUS_NOT_INITIALISED);
}

void binauraliser_setSourceAzi_deg(void* const hBin, int index, float newAzi_deg)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = SAF_MAX(newAzi_deg, -180.0f);
    newAzi_deg = SAF_MIN(newAzi_deg, 180.0f);
    if(pData->src_dirs_deg[index][0]!=newAzi_deg){
        pData->src_dirs_deg[index][0] = newAzi_deg;
        pData->recalc_hrtf_interpFLAG[index] = 1;
        pData->recalc_M_rotFLAG = 1;
    }
}

void binauraliser_setSourceElev_deg(void* const hBin, int index, float newElev_deg)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    newElev_deg = SAF_MAX(newElev_deg, -90.0f);
    newElev_deg = SAF_MIN(newElev_deg, 90.0f);
    if(pData->src_dirs_deg[index][1] != newElev_deg){
        pData->src_dirs_deg[index][1] = newElev_deg;
        pData->recalc_hrtf_interpFLAG[index] = 1;
        pData->recalc_M_rotFLAG = 1;
    }
}

void binauraliser_setNumSources(void* const hBin, int new_nSources)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    pData->new_nSources = SAF_CLAMP(new_nSources, 1, MAX_NUM_INPUTS);
    pData->recalc_M_rotFLAG = 1;
    binauraliser_setCodecStatus(hBin, CODEC_STATUS_NOT_INITIALISED);
}

void binauraliser_setUseDefaultHRIRsflag(void* const hBin, int newState)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    if((!pData->useDefaultHRIRsFLAG) && (newState)){
        pData->useDefaultHRIRsFLAG = newState;
        binauraliser_refreshSettings(hBin);  // re-init and re-calc
    }
}

void binauraliser_setSofaFilePath(void* const hBin, const char* path)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    
    pData->sofa_filepath = realloc1d(pData->sofa_filepath, strlen(path) + 1);
    strcpy(pData->sofa_filepath, path);
    pData->useDefaultHRIRsFLAG = 0;
    binauraliser_refreshSettings(hBin);  // re-init and re-calc
}

void binauraliser_setEnableHRIRsDiffuseEQ(void* const hBin, int newState)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    if(newState!=pData->enableHRIRsDiffuseEQ){
        pData->enableHRIRsDiffuseEQ = newState;
        binauraliser_refreshSettings(hBin);  // re-init and re-calc
    }
}

void binauraliser_setInputConfigPreset(void* const hBin, int newPresetID)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    int ch, dummy;
    
    binauraliser_loadPreset(newPresetID, pData->src_dirs_deg, &(pData->new_nSources), &(dummy));
    if(pData->nSources != pData->new_nSources)
        binauraliser_setCodecStatus(hBin, CODEC_STATUS_NOT_INITIALISED);
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_hrtf_interpFLAG[ch] = 1;
}

void binauraliser_setEnableRotation(void* const hBin, int newState)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    int ch;

    pData->enableRotation = newState;
    if(!pData->enableRotation)
        for (ch = 0; ch<MAX_NUM_INPUTS; ch++) 
            pData->recalc_hrtf_interpFLAG[ch] = 1;
}

void binauraliser_setYaw(void  * const hBin, float newYaw)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    pData->yaw = pData->bFlipYaw == 1 ? -DEG2RAD(newYaw) : DEG2RAD(newYaw);
    pData->recalc_M_rotFLAG = 1;
}

void binauraliser_setPitch(void* const hBin, float newPitch)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    pData->pitch = pData->bFlipPitch == 1 ? -DEG2RAD(newPitch) : DEG2RAD(newPitch);
    pData->recalc_M_rotFLAG = 1;
}

void binauraliser_setRoll(void* const hBin, float newRoll)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    pData->roll = pData->bFlipRoll == 1 ? -DEG2RAD(newRoll) : DEG2RAD(newRoll);
    pData->recalc_M_rotFLAG = 1;
}

void binauraliser_setFlipYaw(void* const hBin, int newState)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    if(newState !=pData->bFlipYaw ){
        pData->bFlipYaw = newState;
        binauraliser_setYaw(hBin, -binauraliser_getYaw(hBin));
    }
}

void binauraliser_setFlipPitch(void* const hBin, int newState)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    if(newState !=pData->bFlipPitch ){
        pData->bFlipPitch = newState;
        binauraliser_setPitch(hBin, -binauraliser_getPitch(hBin));
    }
}

void binauraliser_setFlipRoll(void* const hBin, int newState)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    if(newState !=pData->bFlipRoll ){
        pData->bFlipRoll = newState;
        binauraliser_setRoll(hBin, -binauraliser_getRoll(hBin));
    }
}

void binauraliser_setRPYflag(void* const hBin, int newState)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    pData->useRollPitchYawFlag = newState;
}

void binauraliser_setInterpMode(void* const hBin, int newMode)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    int ch;
    pData->interpMode = newMode;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_hrtf_interpFLAG[ch] = 1;
}

void binauraliser_setSourceGain(void* const hAmbi, int srcIdx, float newGain)
{
    binauraliser_data *pData = (binauraliser_data*)(hAmbi);
    pData->src_gains[srcIdx] = newGain;
}

void binauraliser_setSourceSolo(void* const hAmbi, int srcIdx)
{
    binauraliser_data *pData = (binauraliser_data*)(hAmbi);
    int i;
    for(i=0; i<pData->nSources; i++){
        if(i==srcIdx)
            pData->src_gains[i] = 1.f;
        else
            pData->src_gains[i] = 0.f;
    }
}

void binauraliser_setUnSolo(void* const hAmbi)
{
    binauraliser_data *pData = (binauraliser_data*)(hAmbi);
    for(int i=0; i<pData->nSources; i++)
         pData->src_gains[i] = 1.f;
}


/* Get Functions */

int binauraliser_getFrameSize(void)
{
    return BINAURALISER_FRAME_SIZE;
}

CODEC_STATUS binauraliser_getCodecStatus(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->codecStatus;
}

float binauraliser_getProgressBar0_1(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->progressBar0_1;
}

void binauraliser_getProgressBarText(void* const hBin, char* text)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    memcpy(text, pData->progressBarText, PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

float binauraliser_getSourceAzi_deg(void* const hBin, int index)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->src_dirs_deg[index][0];
}

float binauraliser_getSourceElev_deg(void* const hBin, int index)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->src_dirs_deg[index][1];
}

int binauraliser_getNumSources(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->new_nSources;
}

int binauraliser_getMaxNumSources()
{
    return MAX_NUM_INPUTS;
}

int binauraliser_getNumEars(void)
{
    return NUM_EARS;
}

int binauraliser_getNDirs(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->N_hrir_dirs;
}

int binauraliser_getNTriangles(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->nTriangles;
}

float binauraliser_getHRIRAzi_deg(void* const hBin, int index)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    if(pData->hrir_dirs_deg!=NULL)
        return pData->hrir_dirs_deg[index*2+0];
    else
        return 0.0f;
}

float binauraliser_getHRIRElev_deg(void* const hBin, int index)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    if(pData->hrir_dirs_deg!=NULL)
        return pData->hrir_dirs_deg[index*2+1];
    else
        return 0.0f;
}

int binauraliser_getHRIRlength(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->hrir_loaded_len;
}

int binauraliser_getHRIRsamplerate(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->hrir_loaded_fs;
}

int binauraliser_getUseDefaultHRIRsflag(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->useDefaultHRIRsFLAG;
}

char* binauraliser_getSofaFilePath(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    if(pData->sofa_filepath!=NULL)
        return pData->sofa_filepath;
    else
        return "no_file";
}

int binauraliser_getEnableHRIRsDiffuseEQ(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->enableHRIRsDiffuseEQ;
}

int binauraliser_getDAWsamplerate(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->fs;
}

int binauraliser_getEnableRotation(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->enableRotation;
}

float binauraliser_getYaw(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->bFlipYaw == 1 ? -RAD2DEG(pData->yaw) : RAD2DEG(pData->yaw);
}

float binauraliser_getPitch(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->bFlipPitch == 1 ? -RAD2DEG(pData->pitch) : RAD2DEG(pData->pitch);
}

float binauraliser_getRoll(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->bFlipRoll == 1 ? -RAD2DEG(pData->roll) : RAD2DEG(pData->roll);
}

int binauraliser_getFlipYaw(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->bFlipYaw;
}

int binauraliser_getFlipPitch(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->bFlipPitch;
}

int binauraliser_getFlipRoll(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->bFlipRoll;
}

int binauraliser_getRPYflag(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->useRollPitchYawFlag;
}

int binauraliser_getInterpMode(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return (int)pData->interpMode;
}

int binauraliser_getProcessingDelay()
{
    return 12*HOP_SIZE;
}
