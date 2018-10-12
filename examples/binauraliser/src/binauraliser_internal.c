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
 *     binauraliser_internal.h
 * Description:
 *     Convolves input audio (up to 64 channels) with interpolated HRTFs in the time-frequency
 *     domain. The HRTFs are interpolated by applying amplitude-preserving VBAP gains to the
 *     HRTF magnitude responses and inter-aural time differences (ITDs) individually, before
 *     being re-combined. The example allows the user to specify an external SOFA file for the
 *     convolution.
 * Dependencies:
 *     saf_utilities, saf_hrir, saf_vbap, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 25.09.2017
 */

#include "binauraliser_internal.h"
#ifdef __APPLE__
  #define SAF_ENABLE_SOFA_READER
  #include "saf_sofa_reader.h"
#endif

static inline float matlab_fmodf(float x, float y) {
    float tmp = fmodf(x, y);
    return tmp >= 0 ? tmp : tmp + y; 
}

void binauraliser_interpHRTFs
(
    void* const hBin,
    float azimuth_deg,
    float elevation_deg,
    float_complex h_intrp[HYBRID_BANDS][NUM_EARS]
)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    int i, band;
    int aziIndex, elevIndex, N_azi, idx3d;
    float_complex ipd;
    float aziRes, elevRes, weights[3], itds3[3],  itdInterp;
    float magnitudes3[HYBRID_BANDS][3][NUM_EARS], magInterp[HYBRID_BANDS][NUM_EARS];
     
    /* find closest pre-computed VBAP direction */
    aziRes = (float)pData->hrtf_vbapTableRes[0];
    elevRes = (float)pData->hrtf_vbapTableRes[1];
    N_azi = (int)(360.0f / aziRes + 0.5f) + 1;
    aziIndex = (int)(matlab_fmodf(azimuth_deg + 180.0f, 360.0f) / aziRes + 0.5f);
    elevIndex = (int)((elevation_deg + 90.0f) / elevRes + 0.5f);
    idx3d = elevIndex * N_azi + aziIndex;
    for (i = 0; i < 3; i++)
        weights[i] = pData->hrtf_vbap_gtableComp[idx3d*3 + i];
    
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
            ipd = cmplxf(0.0f, 1.3f*(matlab_fmodf(2.0f*PI*(pData->freqVector[band]) * itdInterp + PI, 2.0f*PI) - PI)/2.0f);
        else
            ipd = cmplxf(0.0f, 0.0f);
        h_intrp[band][0] = crmulf(cexpf(ipd), magInterp[band][0]);
        h_intrp[band][1] = crmulf(conjf(cexpf(ipd)), magInterp[band][1]);
    }
}

void binauraliser_initHRTFsAndGainTables(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    int i, j, k;
    float* hrtf_vbap_gtable;
    
    /* load sofa file or load default hrir data */
    if(!pData->useDefaultHRIRsFLAG){
#ifdef SAF_ENABLE_SOFA_READER
        loadSofaFile(pData->sofa_filepath,
                     &(pData->hrirs),
                     &(pData->hrir_dirs_deg),
                     &(pData->N_hrir_dirs),
                     &(pData->hrir_len),
                     &(pData->hrir_fs));
        if((pData->hrirs==NULL) || (pData->hrir_dirs_deg == NULL))
            pData->useDefaultHRIRsFLAG = 1;
#else
        pData->useDefaultHRIRsFLAG = 1;
#endif
    }
    if(pData->useDefaultHRIRsFLAG){
        /* load defaults */
        pData->N_hrir_dirs = __default_N_hrir_dirs;
        pData->hrir_len = __default_hrir_len;
        pData->hrir_fs = __default_hrir_fs;
        if(pData->hrirs != NULL)
            free(pData->hrirs);
        pData->hrirs = malloc(pData->N_hrir_dirs * 2 * (pData->hrir_len)*sizeof(float));
        for(i=0; i<pData->N_hrir_dirs; i++)
            for(j=0; j<2; j++)
                for(k=0; k< pData->hrir_len; k++)
                    pData->hrirs[i*2*(pData->hrir_len) + j*(pData->hrir_len) + k] = (float)__default_hrirs[i][j][k];
        if(pData->hrir_dirs_deg != NULL)
            free(pData->hrir_dirs_deg);
        pData->hrir_dirs_deg = malloc(pData->N_hrir_dirs * 2 * sizeof(float));
        for(i=0; i<pData->N_hrir_dirs; i++)
            for(j=0; j<2; j++)
                pData->hrir_dirs_deg[i*2+j] = (float)__default_hrir_dirs_deg[i][j];
    }
    
    /* estimate the ITDs for each HRIR */
    if(pData->itds_s!= NULL){
        free(pData->itds_s);
        pData->itds_s = NULL;
    }
    estimateITDs(pData->hrirs, pData->N_hrir_dirs, pData->hrir_len, pData->hrir_fs, &(pData->itds_s));
    
    /* generate VBAP gain table */
    hrtf_vbap_gtable = NULL;
    pData->hrtf_vbapTableRes[0] = 2;
    pData->hrtf_vbapTableRes[1] = 5;
    generateVBAPgainTable3D(pData->hrir_dirs_deg, pData->N_hrir_dirs, pData->hrtf_vbapTableRes[0], pData->hrtf_vbapTableRes[1], 1, 0,
                            &hrtf_vbap_gtable, &(pData->N_hrtf_vbap_gtable), &(pData->nTriangles));
    if(hrtf_vbap_gtable==NULL){
        /* if generating vbap gain tabled failed, re-calculate with default HRIR set */
        pData->useDefaultHRIRsFLAG = 1;
        binauraliser_initHRTFsAndGainTables(hBin);
    }
    
    /* compress VBAP table */
    if(pData->hrtf_vbap_gtableComp!= NULL){
        free(pData->hrtf_vbap_gtableComp);
        pData->hrtf_vbap_gtableComp = NULL;
    }
    if(pData->hrtf_vbap_gtableIdx!= NULL){
        free(pData->hrtf_vbap_gtableIdx);
        pData->hrtf_vbap_gtableIdx = NULL;
    }
    compressVBAPgainTable3D(hrtf_vbap_gtable, pData->N_hrtf_vbap_gtable, pData->N_hrir_dirs, &(pData->hrtf_vbap_gtableComp), &(pData->hrtf_vbap_gtableIdx));
    
    /* convert hrirs to filterbank coefficients */
    if(pData->hrtf_fb!= NULL){
        free(pData->hrtf_fb);
        pData->hrtf_fb = NULL;
    }
    HRIRs2FilterbankHRTFs(pData->hrirs, pData->N_hrir_dirs, pData->hrir_len, pData->itds_s, pData->freqVector, HYBRID_BANDS, 0, &(pData->hrtf_fb));
    
    /* calculate magnitude responses */
    if(pData->hrtf_fb_mag!= NULL)
        free(pData->hrtf_fb_mag);
    pData->hrtf_fb_mag = malloc(HYBRID_BANDS*NUM_EARS* (pData->N_hrir_dirs)*sizeof(float));
    for(i=0; i<HYBRID_BANDS*NUM_EARS* (pData->N_hrir_dirs); i++)
        pData->hrtf_fb_mag[i] = cabsf(pData->hrtf_fb[i]);
    
    /* clean-up */
    if(hrtf_vbap_gtable!=NULL)
        free(hrtf_vbap_gtable);
}

void binauraliser_initTFT
(
    void* const hBin
)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
 
    if(pData->hSTFT==NULL)
        afSTFTinit(&(pData->hSTFT), HOP_SIZE, pData->new_nSources, NUM_EARS, 0, 1);
    else 
        afSTFTchannelChange(pData->hSTFT, pData->new_nSources, NUM_EARS);
    pData->nSources = pData->new_nSources;
}

void binauraliser_loadPreset(PRESETS preset, float dirs_deg[MAX_NUM_INPUTS][2], int* newNCH, int* nDims)
{
    float sum_elev;
    int ch, i, nCH;
    
    switch(preset){
        default:
        case PRESET_DEFAULT:
            nCH = 1;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = 0.0f;
            break;
#ifdef ENABLE_MONO_PRESET
        case PRESET_MONO:
            nCH = 1;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __mono_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_STEREO_PRESET
        case PRESET_STEREO:
            nCH = 2;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __stereo_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_5PX_PRESET
        case PRESET_5PX:
            nCH = 5;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __5pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_7PX_PRESET
        case PRESET_7PX:
            nCH = 7;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __7pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_8PX_PRESET
        case PRESET_8PX:
            nCH = 8;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __8pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_9PX_PRESET
        case PRESET_9PX:
            nCH = 9;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __9pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_10PX_PRESET
        case PRESET_10PX:
            nCH = 10;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __10pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_11PX_PRESET
        case PRESET_11PX:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_11PX_7_4_PRESET
        case PRESET_11PX_7_4:
            nCH = 11;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __11pX_7_4_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_13PX_PRESET
        case PRESET_13PX:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __13pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_22PX_PRESET
        case PRESET_22PX:
            nCH = 22;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __22pX_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_MCC_PRESET
        case PRESET_AALTO_MCC:
            nCH = 29;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_MCC_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_APAJA_PRESET
        case PRESET_AALTO_APAJA:
            nCH = 29;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_Apaja_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_APAJA2_PRESET
        case PRESET_AALTO_APAJA2:
            nCH = 39;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_Apaja2_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_AALTO_LR_PRESET
        case PRESET_AALTO_LR:
            nCH = 13;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Aalto_LR_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_DTU_AVIL_PRESET
        case PRESET_DTU_AVIL:
            nCH = 64;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __DTU_AVIL_dirs_deg[ch][i];
            break;
#endif
#ifdef  ENABLE_ZYLIA_LAB_PRESET
        case PRESET_ZYLIA_LAB:
            nCH = 22;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Zylia_Lab_dirs_deg[ch][i];
            break; 
#endif
#ifdef ENABLE_T_DESIGN_4_PRESET
        case PRESET_T_DESIGN_4:
            nCH = 4;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_2_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_12_PRESET
        case PRESET_T_DESIGN_12:
            nCH = 12;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_4_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_24_PRESET
        case PRESET_T_DESIGN_24:
            nCH = 24;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_6_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_36_PRESET
        case PRESET_T_DESIGN_36:
            nCH = 36;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_8_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_48_PRESET
        case PRESET_T_DESIGN_48:
            nCH = 48;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_9_dirs_deg[ch][i];
            break;
#endif
#ifdef ENABLE_T_DESIGN_60_PRESET
        case PRESET_T_DESIGN_60:
            nCH = 60;
            for(ch=0; ch<nCH; ch++)
                for(i=0; i<2; i++)
                    dirs_deg[ch][i] = __Tdesign_degree_10_dirs_deg[ch][i];
            break;
#endif
    }
    
    /* Fill remaining slots with default coords */
    for(; ch<MAX_NUM_INPUTS; ch++){
        for(i=0; i<2; i++){
            dirs_deg[ch][i] = default_LScoords64_rad[ch][i]* (180.0f/M_PI);
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
 










 
