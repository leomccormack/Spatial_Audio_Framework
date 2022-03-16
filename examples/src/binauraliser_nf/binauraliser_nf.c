/*
 * Copyright 2022 Michael McCrea, Leo McCormack
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
 * @file: binauraliser_nf.c
 * @brief Convolves input audio (up to 64 channels) with interpolated HRTFs in
 *        the time-frequency domain, and applies optional near-field binaural
 *        filtering, as described in [1].
 *
 * The HRTFs are interpolated by applying amplitude-preserving VBAP gains to the
 * HRTF magnitude responses and inter-aural time differences (ITDs)
 * individually, before being re-combined. The example also allows the user to
 * specify an external SOFA file for the convolution, and rotations of the
 * source directions to accomodate head-tracking.
 *
 * @see [1] S. Spagnol, E. Tavazzi, and F. Avanzini, “Distance rendering and
 *          perception of nearby virtual sound sources with a near-field filter
 *          model,” Applied Acoustics, vol. 115, pp. 61–73, Jan. 2017,
 *          doi: 10.1016/j.apacoust.2016.08.015.
 *
 * @author Michael McCrea, Leo McCormack
 * @date 22.02.2022
 * @license ISC
 */

#include "binauraliser_nf_internal.h"

void binauraliserNF_create
(
    void ** const phBin
)
{
    binauraliserNF_data* pData = (binauraliserNF_data*)malloc1d(sizeof(binauraliserNF_data));
    *phBin = (void*)pData;
    int ch, nDim; // nDim not actually used

    /* user parameters */
    pData->useDefaultHRIRsFLAG  = 1; /* pars->sofa_filepath must be valid to set this to 0 */
    pData->enableHRIRsDiffuseEQ = 1;
    pData->nSources     = pData->new_nSources;
    pData->interpMode   = INTERP_TRI;
    pData->yaw          = 0.0f;
    pData->pitch        = 0.0f;
    pData->roll         = 0.0f;
    pData->bFlipYaw     = 0;
    pData->bFlipPitch   = 0;
    pData->bFlipRoll    = 0;
    pData->useRollPitchYawFlag = 0;
    pData->enableRotation = 0;

    /* Near field DVF settings
     * Head radius is set according to the linear combination of head width,
     * height and depth from:
     *      Algazi VR, Avendano C, Duda RO. Estimation of a spherical-head model
     *      from anthropometry. J Audio Eng Soc 2001; 49(6):472-9.
     * The far field threshold is set by rho (normalized distance) = 34, resulting
     * in a ~3 m far field, where the max DVF filter response is about +/-0.5 dB.
     * Near field limit set where filters are stable, in meters from head _center_.
     */
    pData->head_radius       = 0.09096; /* Should match a_head in in saf_utility_dvf.c */
    pData->head_radius_recip = 1.f / pData->head_radius;
    pData->farfield_thresh_m = pData->head_radius * 34.f;
    pData->farfield_headroom = 1.05f; /* 5% headroom above the far field threshold, for resetting to far field and UI range */
    pData->nearfield_limit_m = 0.15f;

    /* Set default source directions and distances */
    binauraliser_loadPreset(SOURCE_CONFIG_PRESET_DEFAULT, pData->src_dirs_deg, &(pData->new_nSources), &nDim); /* check setStateInformation if you change default preset */
    /* For now, any preset selected will reset sources to the far field */
    binauraliserNF_resetSourceDistances(pData); /* Must be called after pData->farfield_thresh_m is set */

    /* time-frequency transform + buffers */
    pData->hSTFT            = NULL;
    /* hrir data */
    pData->hrirs            = NULL;
    pData->hrir_dirs_deg    = NULL;
    pData->sofa_filepath    = NULL;
    pData->weights          = NULL;
    pData->N_hrir_dirs      = pData->hrir_loaded_len = pData->hrir_runtime_len = 0;
    pData->hrir_loaded_fs   = pData->hrir_runtime_fs = -1; /* unknown */
    /* time domain buffers */
    pData->inputFrameTD     = (float**)calloc2d(MAX_NUM_INPUTS, BINAURALISER_FRAME_SIZE, sizeof(float));
    pData->binsrcsTD        = (float**)calloc2d(MAX_NUM_INPUTS * NUM_EARS, BINAURALISER_FRAME_SIZE, sizeof(float));
    pData->outframeTD       = NULL; /* outframeTD inherited from _binauraliser but not used by binauraliserNF functions */
    /* frequency domain buffers */
    pData->inputframeTF     = (float_complex***)calloc3d(HYBRID_BANDS, MAX_NUM_INPUTS, TIME_SLOTS, sizeof(float_complex));
    pData->binauralTF       = (float_complex***)calloc3d(HYBRID_BANDS, MAX_NUM_INPUTS * NUM_EARS, TIME_SLOTS, sizeof(float_complex)); /* MAX_NUM_INPUTS*NUM_EARS dimensions are combined because afSTFT requires float_complex*** */
    pData->outputframeTF    = NULL; /* outputframeTF inherited from _binauraliser but not used by binauraliserNF functions */

    
    /* vbap (amplitude normalised) */
    pData->hrtf_vbap_gtableIdx  = NULL;
    pData->hrtf_vbap_gtableComp = NULL;
    pData->nTriangles = pData->N_hrtf_vbap_gtable = 0;

    /* HRTF filterbank coefficients */
    pData->itds_s      = NULL;
    pData->hrtf_fb     = NULL;
    pData->hrtf_fb_mag = NULL;

    /* flags/status */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PROGRESSBARTEXT_CHAR_LENGTH * sizeof(char));
    strcpy(pData->progressBarText, "");
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
    pData->reInitHRTFsAndGainTables = 1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++) {
        pData->recalc_hrtf_interpFLAG[ch] = 1;
        pData->src_gains[ch] = 1.f;
    }
    pData->recalc_M_rotFLAG = 1;
}

void binauraliserNF_destroy
(
    void ** const phBin
)
{
    binauraliserNF_data *pData = (binauraliserNF_data*)(*phBin);

    if (pData != NULL) {
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            SAF_SLEEP(10);
        }

        if(pData->hSTFT !=NULL)
            afSTFT_destroy(&(pData->hSTFT));
        free(pData->inputFrameTD);
        free(pData->outframeTD);
        free(pData->inputframeTF);
        free(pData->outputframeTF);
        free(pData->binsrcsTD);
        free(pData->binauralTF);
        free(pData->hrtf_vbap_gtableComp);
        free(pData->hrtf_vbap_gtableIdx);
        free(pData->hrtf_fb);
        free(pData->hrtf_fb_mag);
        free(pData->itds_s);
        free(pData->hrirs);
        free(pData->hrir_dirs_deg);
        free(pData->weights);
        free(pData->progressBarText);
        
        free(pData);
        pData = NULL;
    }
}

void binauraliserNF_init
(
      void * const hBin, int sampleRate
)
{
    binauraliser_init(hBin, sampleRate);
}

/* NOTE: This function is a copy of binauraliser_initCodec.
 * The only difference is that it calls binauraliserNF_initTFT, which needs
 * new_nSources * NUM_EARS output channels in the afSTFT. This function could
 * be omitted if binauraliser_initTFT could be refactored to set its number
 * of outputs differently for the regular and NF version of the binauraliser,
 * e.g. a member var afSTFT_nOuts, which changes when binauraliser_setNumSources
 * is called. */
void binauraliserNF_initCodec
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
    binauraliserNF_initTFT(hBin);
    
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

void binauraliserNF_process
(
    void        *  const hBin,
    const float *const * inputs,
    float       ** const outputs,
    int                  nInputs,
    int                  nOutputs,
    int                  nSamples
)
{
    binauraliserNF_data *pData = (binauraliserNF_data*)(hBin);
    int ch, ear, i, band, nSources, srci;
    float src_dirs[MAX_NUM_INPUTS][2], src_dists[MAX_NUM_INPUTS], Rxyz[3][3];
    float hypotxy, headRadiusRecip, fs, ffThresh;
    int enableRotation;
    float thetaLR[2] = { 0.0, 0.0 };
    float rho, wzL, wzR;

    /* copy user parameters to local variables */
    memcpy(src_dirs,  pData->src_dirs_deg, MAX_NUM_INPUTS * sizeof(float)*2);
    memcpy(src_dists, pData->src_dists_m,  MAX_NUM_INPUTS * sizeof(float));

    nSources        = pData->nSources;
    enableRotation  = pData->enableRotation;
    headRadiusRecip = pData->head_radius_recip;
    ffThresh        = pData->farfield_thresh_m;
    fs              = (float)pData->fs;

    /* apply binaural panner */
    if ((nSamples == BINAURALISER_FRAME_SIZE) && (pData->hrtf_fb!=NULL) && (pData->codecStatus==CODEC_STATUS_INITIALISED) ) {
        pData->procStatus = PROC_STATUS_ONGOING;

        /* Rotate source directions */
        if (enableRotation && pData->recalc_M_rotFLAG) {
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
            for (i = 0; i < nSources; i++) {
                hypotxy = sqrtf(powf(pData->src_dirs_rot_xyz[i][0], 2.0f) + powf(pData->src_dirs_rot_xyz[i][1], 2.0f));
                pData->src_dirs_rot_deg[i][0] = RAD2DEG(atan2f(pData->src_dirs_rot_xyz[i][1], pData->src_dirs_rot_xyz[i][0]));
                pData->src_dirs_rot_deg[i][1] = RAD2DEG(atan2f(pData->src_dirs_rot_xyz[i][2], hypotxy));
            }
            pData->recalc_M_rotFLAG = 0;
        }
        
        /* Load time-domain data */
        for (i = 0; i < SAF_MIN(nSources, nInputs); i++)
            utility_svvcopy(inputs[i], BINAURALISER_FRAME_SIZE, pData->inputFrameTD[i]);
        for (; i < nSources; i++)
            memset(pData->inputFrameTD[i], 0, BINAURALISER_FRAME_SIZE * sizeof(float));
        
        /* Apply source gains */
        for (ch = 0; ch < nSources; ch++) {
            if(fabsf(pData->src_gains[ch] - 1.f) > 1e-6f)
                utility_svsmul(pData->inputFrameTD[ch], &(pData->src_gains[ch]), BINAURALISER_FRAME_SIZE, NULL);
        }
        
        /* Zero out busses - not needed; overwritten, not summeed over */
        // memset(FLATTEN2D(pData->binsrcsTD), 0, NUM_EARS * BINAURALISER_FRAME_SIZE * sizeof(float));
        // memset(FLATTEN3D(pData->binauralTF), 0, HYBRID_BANDS * nSources * NUM_EARS * TIME_SLOTS * sizeof(float_complex));

        /* Apply time-frequency transform (TFT) */
        afSTFT_forward_knownDimensions(pData->hSTFT, pData->inputFrameTD, BINAURALISER_FRAME_SIZE, MAX_NUM_INPUTS, TIME_SLOTS, pData->inputframeTF);

        /* Interpolate and apply HRTFs, maintain individual binaural source streams */
        for (ch = 0; ch < nSources; ch++) {
            /* Interpolate hrtfs */
            if (pData->recalc_hrtf_interpFLAG[ch]) {
                if (enableRotation)
                    binauraliser_interpHRTFs(hBin, pData->interpMode, pData->src_dirs_rot_deg[ch][0], pData->src_dirs_rot_deg[ch][1], pData->hrtf_interp[ch]);
                else
                    binauraliser_interpHRTFs(hBin, pData->interpMode, pData->src_dirs_deg[ch][0], pData->src_dirs_deg[ch][1], pData->hrtf_interp[ch]);
                pData->recalc_hrtf_interpFLAG[ch] = 0;
            }
            
            /* Expand each source channel to binaural by applying HRTF filter (FD) */
            for (band = 0; band < HYBRID_BANDS; band++) {
                for (ear = 0; ear < NUM_EARS; ear++) {
                    utility_cvsmul(pData->inputframeTF[band][ch],
                                   &pData->hrtf_interp[ch][band][ear], TIME_SLOTS,
                                   &pData->binauralTF[band][NUM_EARS * ch + ear][0]);
                }
            }
        }

        /* Bring sources back to time domain.
         * Note, in the case of AFSTFT_TIME_CH_BANDS, args dataFD_nCH & dataFD_nHops are unused */
        afSTFT_backward_knownDimensions(pData->hSTFT, pData->binauralTF, BINAURALISER_FRAME_SIZE, MAX_NUM_INPUTS * NUM_EARS, TIME_SLOTS, pData->binsrcsTD);

        /* Iterate over sources, applying DVF to those in the near field */
        for (ch = 0; ch < nSources; ch++) {
            int l = NUM_EARS * ch;
            int r = l + 1;
            if (src_dists[ch] < ffThresh) {
                /* Get previous frame's last sample for DVF IIR filter */
                wzL = pData->binsrcsTD[l][BINAURALISER_FRAME_SIZE - 1];
                wzR = pData->binsrcsTD[r][BINAURALISER_FRAME_SIZE - 1];

                rho = src_dists[ch] * headRadiusRecip;
                convertFrontalDoAToIpsilateral(src_dirs[ch][0], &thetaLR[0]);
                applyDVF(thetaLR[0], rho, pData->binsrcsTD[l], BINAURALISER_FRAME_SIZE, fs, &wzL, pData->binsrcsTD[l]);
                applyDVF(thetaLR[1], rho, pData->binsrcsTD[r], BINAURALISER_FRAME_SIZE, fs, &wzR, pData->binsrcsTD[r]);
            }
        }

        /* Zero out all plugin output channels, cblas_saxpy below accumulates signals into the first 2 */
        for (ch = 0; ch < nOutputs; ch++)
            memset(outputs[ch], 0, BINAURALISER_FRAME_SIZE * sizeof(float));
                
        /* Sum binaural sources */
        // Note, it appears SAF doesn't support in-place vector addition, e.g. ippsAdd_32f_I
        for (srci = 0; srci < nSources; srci++) {
            for (ear = 0; ear < NUM_EARS; ear++) {
                cblas_saxpy(BINAURALISER_FRAME_SIZE, 1.f,
                            pData->binsrcsTD[srci * NUM_EARS + ear], 1, outputs[ear], 1);
            }
        }
        /* Scale binaural output by number of sources */
        for (ear = 0; ear < NUM_EARS; ear++) {
            cblas_sscal(BINAURALISER_FRAME_SIZE, 1.0f / sqrtf((float)nSources), outputs[ear], 1);
        }
    }
    else {
        for (ch = 0; ch < nOutputs; ch++)
            memset(outputs[ch], 0, BINAURALISER_FRAME_SIZE * sizeof(float));
    }

    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

/* Set Functions */

void binauraliserNF_setSourceDist_m(void* const hBin, int index, float newDist_m)
{
    binauraliserNF_data *pData = (binauraliserNF_data*)(hBin);
    newDist_m = SAF_MAX(newDist_m, pData->nearfield_limit_m);
    if(pData->src_dists_m[index] != newDist_m){
        pData->src_dists_m[index] = newDist_m;
    }
}

void binauraliserNF_setInputConfigPreset(void* const hBin, int newPresetID)
{
    binauraliserNF_data *pData = (binauraliserNF_data*)(hBin);
    int ch, nDim;

    binauraliser_loadPreset(newPresetID, pData->src_dirs_deg, &(pData->new_nSources), &nDim);
    // For now, any preset selected will reset sources to the far field
    binauraliserNF_resetSourceDistances(pData);
    
    if(pData->nSources != pData->new_nSources)
        binauraliser_setCodecStatus(hBin, CODEC_STATUS_NOT_INITIALISED);
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_hrtf_interpFLAG[ch] = 1;
}

float binauraliserNF_getSourceDist_m(void* const hBin, int index)
{
    binauraliserNF_data *pData = (binauraliserNF_data*)(hBin);
    
    return pData->src_dists_m[index];
}

float binauraliserNF_getFarfieldThresh_m(void* const hBin)
{
    binauraliserNF_data *pData = (binauraliserNF_data*)(hBin);
    
    return pData->farfield_thresh_m;
}

float binauraliserNF_getFarfieldHeadroom(void* const hBin)
{
    binauraliserNF_data *pData = (binauraliserNF_data*)(hBin);
    
    return pData->farfield_headroom;
}


float binauraliserNF_getNearfieldLimit_m(void* const hBin)
{
    binauraliserNF_data *pData = (binauraliserNF_data*)(hBin);
    
    return pData->nearfield_limit_m;
}
