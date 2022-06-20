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
 * @file: binauraliser_internal.c
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

void binauraliser_setCodecStatus(void* const hBin, CODEC_STATUS newStatus)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    if(newStatus==CODEC_STATUS_NOT_INITIALISED){
        /* Pause until current initialisation is complete */
        while(pData->codecStatus == CODEC_STATUS_INITIALISING)
            SAF_SLEEP(10);
    }
    pData->codecStatus = newStatus;
}

void binauraliser_interpHRTFs
(
    void* const hBin,
    INTERP_MODES mode,
    float azimuth_deg,
    float elevation_deg,
    float_complex h_intrp[HYBRID_BANDS][NUM_EARS]
)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    int i, band;
    int aziIndex, elevIndex, N_azi, idx3d;
    float_complex ipd;
    float_complex weights_cmplx[3], hrtf_fb3[NUM_EARS][3];
    float aziRes, elevRes, weights[3], itds3[3],  itdInterp;
    float magnitudes3[HYBRID_BANDS][3][NUM_EARS], magInterp[HYBRID_BANDS][NUM_EARS];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
     
    /* find closest pre-computed VBAP direction */
    aziRes = (float)pData->hrtf_vbapTableRes[0];
    elevRes = (float)pData->hrtf_vbapTableRes[1];
    N_azi = (int)(360.0f / aziRes + 0.5f) + 1;
    aziIndex = (int)(matlab_fmodf(azimuth_deg + 180.0f, 360.0f) / aziRes + 0.5f);
    elevIndex = (int)((elevation_deg + 90.0f) / elevRes + 0.5f);
    idx3d = elevIndex * N_azi + aziIndex;
    for (i = 0; i < 3; i++)
        weights[i] = pData->hrtf_vbap_gtableComp[idx3d*3 + i];

    switch(mode){
        case INTERP_TRI:
            for (i = 0; i < 3; i++)
                weights_cmplx[i] = cmplxf(weights[i], 0.0f);
            for (band = 0; band < HYBRID_BANDS; band++) {
                for (i = 0; i < 3; i++){
                    hrtf_fb3[0][i] = pData->hrtf_fb[band*NUM_EARS*(pData->N_hrir_dirs) + 0*(pData->N_hrir_dirs) + pData->hrtf_vbap_gtableIdx[idx3d*3+i]];
                    hrtf_fb3[1][i] = pData->hrtf_fb[band*NUM_EARS*(pData->N_hrir_dirs) + 1*(pData->N_hrir_dirs) + pData->hrtf_vbap_gtableIdx[idx3d*3+i]];
                } 
                cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, 1, 3, &calpha,
                            (float_complex*)hrtf_fb3, 3,
                            (float_complex*)weights_cmplx, 1, &cbeta,
                            (float_complex*)h_intrp[band], 1);
            }
            break;

        case INTERP_TRI_PS:
            /* retrieve the 3 itds and hrtf magnitudes */
            for (i = 0; i < 3; i++) {
                itds3[i] = pData->itds_s[pData->hrtf_vbap_gtableIdx[idx3d*3+i]];
                for (band = 0; band < HYBRID_BANDS; band++) {
                    magnitudes3[band][i][0] = pData->hrtf_fb_mag[band*NUM_EARS*(pData->N_hrir_dirs) + 0*(pData->N_hrir_dirs) + pData->hrtf_vbap_gtableIdx[idx3d*3+i]];
                    magnitudes3[band][i][1] = pData->hrtf_fb_mag[band*NUM_EARS*(pData->N_hrir_dirs) + 1*(pData->N_hrir_dirs) + pData->hrtf_vbap_gtableIdx[idx3d*3+i]];
                }
            }

            /* interpolate hrtf magnitudes and itd */
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, 1, 3, 1.0f,
                        (float*)weights, 3,
                        (float*)itds3, 1, 0.0f,
                        &itdInterp, 1);
            for (band = 0; band < HYBRID_BANDS; band++) {
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, 2, 3, 1.0f,
                            (float*)weights, 3,
                            (float*)magnitudes3[band], 2, 0.0f,
                            (float*)magInterp[band], 2);
            }

            /* introduce interaural phase difference */
            for (band = 0; band < HYBRID_BANDS; band++) {
                if(pData->freqVector[band]<1.5e3f)
                    ipd = cmplxf(0.0f, (matlab_fmodf(2.0f*SAF_PI*(pData->freqVector[band]) * itdInterp + SAF_PI, 2.0f*SAF_PI) - SAF_PI)/2.0f);
                else
                    ipd = cmplxf(0.0f, 0.0f);
                h_intrp[band][0] = crmulf(cexpf(ipd), magInterp[band][0]);
                h_intrp[band][1] = crmulf(conjf(cexpf(ipd)), magInterp[band][1]);
            }
            break;
    }
}

void binauraliser_initHRTFsAndGainTables(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    int i, new_len;
    float* hrtf_vbap_gtable, *hrirs_resampled;//, *hrir_dirs_rad;
#ifdef SAF_ENABLE_SOFA_READER_MODULE
    SAF_SOFA_ERROR_CODES error;
    saf_sofa_container sofa;
#endif
    
    strcpy(pData->progressBarText,"Loading HRIRs");
    pData->progressBar0_1 = 0.2f;
    
    /* load sofa file or load default hrir data */
#ifdef SAF_ENABLE_SOFA_READER_MODULE
    if(!pData->useDefaultHRIRsFLAG && pData->sofa_filepath!=NULL){
        /* Load SOFA file */
        error = saf_sofa_open(&sofa, pData->sofa_filepath, SAF_SOFA_READER_OPTION_DEFAULT);

        /* Load defaults instead */
        if(error!=SAF_SOFA_OK || sofa.nReceivers!=NUM_EARS){
            pData->useDefaultHRIRsFLAG = 1;
            saf_print_warning("Unable to load the specified SOFA file, or it contained something other than 2 channels. Using default HRIR data instead.");
        }
        else{
            /* Copy SOFA data */
            pData->hrir_loaded_fs = (int)sofa.DataSamplingRate;
            pData->hrir_loaded_len = sofa.DataLengthIR;
            pData->N_hrir_dirs = sofa.nSources;
            pData->hrirs = realloc1d(pData->hrirs, pData->N_hrir_dirs*NUM_EARS*(pData->hrir_loaded_len)*sizeof(float));
            memcpy(pData->hrirs, sofa.DataIR, pData->N_hrir_dirs*NUM_EARS*(pData->hrir_loaded_len)*sizeof(float));
            pData->hrir_dirs_deg = realloc1d(pData->hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));
            cblas_scopy(pData->N_hrir_dirs, sofa.SourcePosition, 3, pData->hrir_dirs_deg, 2); /* azi */
            cblas_scopy(pData->N_hrir_dirs, &sofa.SourcePosition[1], 3, &pData->hrir_dirs_deg[1], 2); /* elev */ 
        }

        /* Clean-up */
        saf_sofa_close(&sofa);
    }
#else
    pData->useDefaultHRIRsFLAG = 1; /* Can only load the default HRIR data */
#endif
    if(pData->useDefaultHRIRsFLAG){
        /* Copy default HRIR data */
        pData->hrir_loaded_fs = __default_hrir_fs;
        pData->hrir_loaded_len = __default_hrir_len;
        pData->N_hrir_dirs = __default_N_hrir_dirs;
        pData->hrirs = realloc1d(pData->hrirs, pData->N_hrir_dirs*NUM_EARS*(pData->hrir_loaded_len)*sizeof(float));
        memcpy(pData->hrirs, (float*)__default_hrirs, pData->N_hrir_dirs*NUM_EARS*(pData->hrir_loaded_len)*sizeof(float));
        pData->hrir_dirs_deg = realloc1d(pData->hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));
        memcpy(pData->hrir_dirs_deg, (float*)__default_hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));
    }

    /* Convert from the 0..360 convention, to -180..180 */
    convert_0_360To_m180_180(pData->hrir_dirs_deg, pData->N_hrir_dirs);

    /* estimate the ITDs for each HRIR */
    strcpy(pData->progressBarText,"Estimating ITDs");
    pData->progressBar0_1 = 0.4f;
    pData->itds_s = realloc1d(pData->itds_s, pData->N_hrir_dirs*sizeof(float));
    estimateITDs(pData->hrirs, pData->N_hrir_dirs, pData->hrir_loaded_len, pData->hrir_loaded_fs, pData->itds_s);

    /* Resample the HRIRs if needed */
    if(pData->hrir_loaded_fs!=pData->fs){
        strcpy(pData->progressBarText,"Resampling the HRIRs");
        pData->progressBar0_1 = 0.5f;
        hrirs_resampled = NULL;
        resampleHRIRs(pData->hrirs, pData->N_hrir_dirs, pData->hrir_loaded_len, pData->hrir_loaded_fs, pData->fs, 1, &hrirs_resampled, &new_len);
        pData->hrirs = realloc1d(pData->hrirs, pData->N_hrir_dirs*NUM_EARS*new_len*sizeof(float));
        cblas_scopy(pData->N_hrir_dirs*NUM_EARS*new_len, hrirs_resampled, 1, pData->hrirs, 1);
        free(hrirs_resampled);
        pData->hrir_runtime_fs = pData->fs;
        pData->hrir_runtime_len = new_len;
    }
    else{
        pData->hrir_runtime_fs = pData->hrir_loaded_fs;
        pData->hrir_runtime_len = pData->hrir_loaded_len;
    }
    
    /* generate VBAP gain table */
    strcpy(pData->progressBarText,"Generating interpolation table");
    pData->progressBar0_1 = 0.6f;
    hrtf_vbap_gtable = NULL;
    pData->hrtf_vbapTableRes[0] = 2;
    pData->hrtf_vbapTableRes[1] = 5;
    generateVBAPgainTable3D(pData->hrir_dirs_deg, pData->N_hrir_dirs, pData->hrtf_vbapTableRes[0], pData->hrtf_vbapTableRes[1], 1, 0, 0.0f,
                            &hrtf_vbap_gtable, &(pData->N_hrtf_vbap_gtable), &(pData->nTriangles));
    if(hrtf_vbap_gtable==NULL){
        /* if generating vbap gain tabled failed, re-calculate with default HRIR set */
        pData->useDefaultHRIRsFLAG = 1;
        binauraliser_initHRTFsAndGainTables(hBin);
    }
    
    /* compress VBAP table (i.e. remove the zero elements) */
    pData->hrtf_vbap_gtableComp = realloc1d(pData->hrtf_vbap_gtableComp, pData->N_hrtf_vbap_gtable * 3 * sizeof(float));
    pData->hrtf_vbap_gtableIdx  = realloc1d(pData->hrtf_vbap_gtableIdx,  pData->N_hrtf_vbap_gtable * 3 * sizeof(int));
    compressVBAPgainTable3D(hrtf_vbap_gtable, pData->N_hrtf_vbap_gtable, pData->N_hrir_dirs, pData->hrtf_vbap_gtableComp, pData->hrtf_vbap_gtableIdx);
    
    /* convert hrirs to filterbank coefficients */
    pData->progressBar0_1 = 0.6f;
    pData->hrtf_fb = realloc1d(pData->hrtf_fb, HYBRID_BANDS * NUM_EARS * (pData->N_hrir_dirs)*sizeof(float_complex));
    HRIRs2HRTFs_afSTFT(pData->hrirs, pData->N_hrir_dirs, pData->hrir_runtime_len, HOP_SIZE, 0, 1, pData->hrtf_fb);

    /* HRIR pre-processing */
    if(pData->enableHRIRsDiffuseEQ){
        /* get integration weights */
        strcpy(pData->progressBarText,"Applying HRIR diffuse-field EQ");
        pData->progressBar0_1 = 0.9f;
        if(pData->N_hrir_dirs<=1000){//3600
//            hrir_dirs_rad = malloc1d(pData->N_hrir_dirs*2*sizeof(float));
//            memcpy(hrir_dirs_rad, pData->hrir_dirs_deg, pData->N_hrir_dirs*2*sizeof(float));
//            cblas_sscal(pData->N_hrir_dirs*2, SAF_PI/180.0f, hrir_dirs_rad, 1);
            pData->weights = realloc1d(pData->weights, pData->N_hrir_dirs*sizeof(float));
//            calculateGridWeights(hrir_dirs_rad, pData->N_hrir_dirs, -1, pData->weights);
//            free(hrir_dirs_rad);
            getVoronoiWeights(pData->hrir_dirs_deg, pData->N_hrir_dirs, 0, pData->weights);
        }
        else{
            pData->weights = realloc1d(pData->weights, pData->N_hrir_dirs*sizeof(float));
            for(int idx=0; idx < pData->N_hrir_dirs; idx++)
                pData->weights[idx] = 4.f*SAF_PI / (float)pData->N_hrir_dirs;
        }
        diffuseFieldEqualiseHRTFs(pData->N_hrir_dirs, pData->itds_s, pData->freqVector, HYBRID_BANDS, pData->weights, 1, 0, pData->hrtf_fb);
    }

    /* calculate magnitude responses */
    pData->hrtf_fb_mag = realloc1d(pData->hrtf_fb_mag, HYBRID_BANDS*NUM_EARS*(pData->N_hrir_dirs)*sizeof(float)); 
    for(i=0; i<HYBRID_BANDS*NUM_EARS* (pData->N_hrir_dirs); i++)
        pData->hrtf_fb_mag[i] = cabsf(pData->hrtf_fb[i]);

    /* The HRTFs should be re-interpolated */
    for(i=0; i<MAX_NUM_INPUTS; i++)
        pData->recalc_hrtf_interpFLAG[i] = 1;
    
    /* clean-up */
    free(hrtf_vbap_gtable);
}

void binauraliser_initTFT
(
    void* const hBin
)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
 
    if(pData->hSTFT==NULL)
        afSTFT_create(&(pData->hSTFT), pData->new_nSources, NUM_EARS, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    else if(pData->new_nSources!=pData->nSources){
        afSTFT_channelChange(pData->hSTFT, pData->new_nSources, NUM_EARS);
        afSTFT_clearBuffers(pData->hSTFT);
    }
    pData->nSources = pData->new_nSources;
}

void binauraliser_loadPreset
(
    SOURCE_CONFIG_PRESETS preset,
    float dirs_deg[MAX_NUM_INPUTS][2],
    int* newNCH,
    int* nDims
)
{
    float sum_elev;
    int ch, i, nCH;
    
    switch(preset){
        default:
        case SOURCE_CONFIG_PRESET_DEFAULT:
            nCH = 1;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = 0.0f;
            break;
        case SOURCE_CONFIG_PRESET_MONO:
            nCH = 1;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __mono_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_STEREO:
            nCH = 2;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __stereo_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_5PX:
            nCH = 5;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __5pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_7PX:
            nCH = 7;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __7pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_8PX:
            nCH = 8;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __8pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_9PX:
            nCH = 9;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __9pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_10PX:
            nCH = 10;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __10pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_11PX:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_11PX_7_4:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_7_4_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_13PX:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __13pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_22PX:
            nCH = 22;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __22pX_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_22P2_9_10_3:
            nCH = 24;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __9_10_3p2_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_AALTO_MCC:
            nCH = 45;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_MCC_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_AALTO_MCC_SUBSET:
            nCH = 37;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_MCCsubset_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_AALTO_APAJA:
            nCH = 29;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_Apaja_dirs_deg[ch][i];
            break; 
        case SOURCE_CONFIG_PRESET_AALTO_LR:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_LR_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_DTU_AVIL:
            nCH = 64;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __DTU_AVIL_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_ZYLIA_LAB:
            nCH = 22;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Zylia_Lab_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_4:
            nCH = 4;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_2_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_12:
            nCH = 12;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_4_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_24:
            nCH = 24;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_6_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_36:
            nCH = 36;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_8_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_48:
            nCH = 48;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_9_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_T_DESIGN_60:
            nCH = 60;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_10_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_9:
            nCH = 9;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_9_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_16:
            nCH = 16;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_16_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_25:
            nCH = 25;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_25_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_49:
            nCH = 49;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_49_dirs_deg[ch][i];
            break;
        case SOURCE_CONFIG_PRESET_SPH_COV_64:
            nCH = 64;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __SphCovering_64_dirs_deg[ch][i];
            break;
    }
    
    /* Fill remaining slots with default coords */
    for(; ch<MAX_NUM_INPUTS; ch++){
        for(i=0; i<2; i++){
            dirs_deg[ch][i] = __default_LScoords64_rad[ch][i]* (180.0f/SAF_PI);
        }
    }
    
    /* For dynamically changing the number of TFT channels */
    (*newNCH) = nCH;
    
    /* estimate number of dimensions. (Obviously fails if using 2D setups thare are on an angle.
       However, in these cases, triangulation should fail and revert to 2D anyway) */
    sum_elev = 0.0f;
    for(i=0; i<nCH; i++){
        sum_elev += dirs_deg[i][1];
    }
    if(sum_elev < 0.01f)
        (*nDims) = 2;
    else
        (*nDims) = 3;
}
 










 
