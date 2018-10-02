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
 *     binauraliser.c
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

void binauraliser_create
(
    void ** const phBin
)
{
    binauraliser_data* pData = (binauraliser_data*)malloc(sizeof(binauraliser_data));
    if (pData == NULL) { return;/*error*/ }
    *phBin = (void*)pData;
    
    /* time-frequency transform + buffers */
    pData->hSTFT = NULL;
    pData->STFTInputFrameTF = NULL;
    int t, ch;
    pData->STFTOutputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, NUM_EARS, sizeof(complexVector));
    for(t=0; t<TIME_SLOTS; t++) {
        for(ch=0; ch< NUM_EARS; ch++) {
            pData->STFTOutputFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
            pData->STFTOutputFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
        }
    }
    pData->tempHopFrameTD = NULL;
    
    /* hrir data */
    pData->useDefaultHRIRsFLAG=1;
    pData->hrirs = NULL;
    pData->hrir_dirs_deg = NULL;
    pData->sofa_filepath = "/Users/mccorml1/Documents/SourceTree/AkisPlugins/database/AALTO/leo_aalto2016.sofa";
    //pData->sofa_filename = "/Users/mccorml1/Documents/SourceTree/AkisPlugins/database/CIPIC/subject_003.sofa";
    /* vbap */
    pData->hrtf_vbap_gtableIdx = NULL;
    pData->hrtf_vbap_gtableComp = NULL;
    
    /* HRTF filterbank coefficients */
    pData->itds_s = NULL;
    pData->hrtf_fb = NULL;
    pData->hrtf_fb_mag = NULL;
    
    /* flags */
    pData->reInitHRTFsAndGainTables = 1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_hrtf_interpFLAG[ch] = 1;
    pData->reInitTFT = 1;
    
    /* user parameters */
    binauraliser_loadPreset(PRESET_DEFAULT, pData->src_dirs_deg, &(pData->new_nSources), &(pData->input_nDims)); /*check setStateInformation if you change default preset*/
    pData->nSources = pData->new_nSources;
}


void binauraliser_destroy
(
    void ** const phBin
)
{
    binauraliser_data *pData = (binauraliser_data*)(*phBin);
    int t, ch;

    if (pData != NULL) {
        if(pData->hSTFT !=NULL)
            afSTFTfree(pData->hSTFT);
        for (t = 0; t<TIME_SLOTS; t++) {
            if(pData->STFTInputFrameTF!=NULL){
                for(ch=0; ch< pData->nSources; ch++) {
                    free(pData->STFTInputFrameTF[t][ch].re);
                    free(pData->STFTInputFrameTF[t][ch].im);
                }
            }
            for (ch = 0; ch< NUM_EARS; ch++) {
                free(pData->STFTOutputFrameTF[t][ch].re);
                free(pData->STFTOutputFrameTF[t][ch].im);
            }
        }
        if(pData->STFTInputFrameTF!=NULL)
            free2d((void**)pData->STFTInputFrameTF, TIME_SLOTS);
        if(pData->STFTOutputFrameTF!=NULL)
            free2d((void**)pData->STFTOutputFrameTF, TIME_SLOTS);
        if(pData->tempHopFrameTD!=NULL)
            free2d((void**)pData->tempHopFrameTD, MAX(pData->nSources, NUM_EARS));
        
        if(pData->hrtf_vbap_gtableComp!= NULL)
            free(pData->hrtf_vbap_gtableComp);
        if(pData->hrtf_vbap_gtableIdx!= NULL)
            free(pData->hrtf_vbap_gtableIdx);
        if(pData->hrtf_fb!= NULL)
            free(pData->hrtf_fb);
        if(pData->hrtf_fb_mag!= NULL)
            free(pData->hrtf_fb_mag);
        if(pData->itds_s!= NULL)
            free(pData->itds_s);
        if(pData->hrirs!= NULL)
            free(pData->hrirs);
        if(pData->hrir_dirs_deg!= NULL)
            free(pData->hrir_dirs_deg);
         
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
    int band;
    
    /* define frequency vector */
    pData->fs = sampleRate;
    for(band=0; band <HYBRID_BANDS; band++){
        if(sampleRate == 44100)
            pData->freqVector[band] =  (float)__afCenterFreq44100[band];
        else
            pData->freqVector[band] =  (float)__afCenterFreq48e3[band];
    }
}

void binauraliser_process
(
    void  *  const hBin,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples,
    int            isPlaying
)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    int t, sample, ch, ear, i, band, nSources; 
    float src_dirs[MAX_NUM_INPUTS][2];
    
#ifdef ENABLE_FADE_IN_OUT
    int applyFadeIn;
    if(pData->reInitTFT || pData->reInitHRTFsAndGainTables)
        applyFadeIn = 1;
    else
        applyFadeIn = 0;
#endif
    
    /* reinitialise if needed */
    if(pData->reInitTFT){
        binauraliser_initTFT(hBin);
        pData->reInitTFT = 0;
        nSources = pData->nSources;
    }
    else
        nSources = pData->nSources;
    if(pData->reInitHRTFsAndGainTables){
        binauraliser_initHRTFsAndGainTables(hBin);
        pData->reInitHRTFsAndGainTables = 0;
    }
    
    /* apply binaural panner */
    if ((nSamples == FRAME_SIZE) && (isPlaying == 1) && (pData->hrtf_fb!=NULL)) {
        nSources = pData->nSources;  
        memcpy(src_dirs, pData->src_dirs_deg, MAX_NUM_INPUTS*2*sizeof(float));
        /* Load time-domain data */
        for(i=0; i < MIN(nSources,nInputs); i++)
            memcpy(pData->inputFrameTD[i], inputs[i], FRAME_SIZE * sizeof(float));
        for(; i<MAX_NUM_INPUTS; i++)
            memset(pData->inputFrameTD[i], 0, FRAME_SIZE * sizeof(float));
#ifdef ENABLE_FADE_IN_OUT
        if(applyFadeIn)
            for(ch=0; ch < nSources;ch++)
                for(i=0; i<FRAME_SIZE; i++)
                    pData->inputFrameTD[ch][i] *= (float)i/(float)FRAME_SIZE;
#endif
        
        /* Apply time-frequency transform (TFT) */
        for ( t=0; t< TIME_SLOTS; t++) {
            for( ch=0; ch < nSources; ch++)
                for ( sample=0; sample < HOP_SIZE; sample++)
                    pData->tempHopFrameTD[ch][sample] = pData->inputFrameTD[ch][sample + t*HOP_SIZE];
            afSTFTforward(pData->hSTFT, (float**)pData->tempHopFrameTD, (complexVector*)pData->STFTInputFrameTF[t]);
        }
        for(band=0; band<HYBRID_BANDS; band++)
            for( ch=0; ch < nSources; ch++)
                for ( t=0; t<TIME_SLOTS; t++)
                    pData->inputframeTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[t][ch].re[band], pData->STFTInputFrameTF[t][ch].im[band]);
     
        /* interpolate hrtfs and apply to each source */
        memset(pData->outputframeTF, 0, HYBRID_BANDS*NUM_EARS*TIME_SLOTS * sizeof(float_complex));
        for (ch = 0; ch < nSources; ch++) {
            if(pData->recalc_hrtf_interpFLAG[ch]){
                binauraliser_interpHRTFs(hBin, pData->src_dirs_deg[ch][0], pData->src_dirs_deg[ch][1], pData->hrtf_interp[ch]);
                pData->recalc_hrtf_interpFLAG[ch] = 0;
            }
            for (band = 0; band < HYBRID_BANDS; band++)
                for (ear = 0; ear < NUM_EARS; ear++)
                    for (t = 0; t < TIME_SLOTS; t++)
                        pData->outputframeTF[band][ear][t] = ccaddf(pData->outputframeTF[band][ear][t], ccmulf(pData->inputframeTF[band][ch][t], pData->hrtf_interp[ch][band][ear]));
        }
        
        /* scale by number of sources */
        for (band = 0; band < HYBRID_BANDS; band++)
            for (ear = 0; ear < NUM_EARS; ear++)
                for (t = 0; t < TIME_SLOTS; t++)
                    pData->outputframeTF[band][ear][t] = crmulf(pData->outputframeTF[band][ear][t], 1.0f/sqrtf((float)nSources));
        
        /* inverse-TFT */
        for (band = 0; band < HYBRID_BANDS; band++) {
            for (ch = 0; ch < NUM_EARS; ch++) {
                for (t = 0; t < TIME_SLOTS; t++) {
                    pData->STFTOutputFrameTF[t][ch].re[band] = crealf(pData->outputframeTF[band][ch][t]);
                    pData->STFTOutputFrameTF[t][ch].im[band] = cimagf(pData->outputframeTF[band][ch][t]);
                }
            }
        }
        for (t = 0; t < TIME_SLOTS; t++) {
            afSTFTinverse(pData->hSTFT, pData->STFTOutputFrameTF[t], pData->tempHopFrameTD);
            for (ch = 0; ch < MIN(NUM_EARS, nOutputs); ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    outputs[ch][sample + t* HOP_SIZE] = pData->tempHopFrameTD[ch][sample];
            for (; ch < nOutputs; ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    outputs[ch][sample + t* HOP_SIZE] = 0.0f;
        }
#ifdef ENABLE_FADE_IN_OUT
        if(pData->reInitTFT || pData->reInitHRTFsAndGainTables)
            for(ch=0; ch < NUM_EARS;ch++)
                for(i=0; i<FRAME_SIZE; i++)
                    outputs[ch][i] *= (1.0f - (float)(i+1)/(float)FRAME_SIZE);
#endif
    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, FRAME_SIZE*sizeof(float));
    }
}

/* Set Functions */

void binauraliser_refreshSettings(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    pData->reInitHRTFsAndGainTables = 1;
    pData->reInitTFT = 1;
}

void binauraliser_setSourceAzi_deg(void* const hBin, int index, float newAzi_deg)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    newAzi_deg = MAX(newAzi_deg, -180.0f);
    newAzi_deg = MIN(newAzi_deg, 180.0f);
    pData->recalc_hrtf_interpFLAG[index] = 1;
    pData->src_dirs_deg[index][0] = newAzi_deg;
}

void binauraliser_setSourceElev_deg(void* const hBin, int index, float newElev_deg)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    newElev_deg = MAX(newElev_deg, -90.0f);
    newElev_deg = MIN(newElev_deg, 90.0f);
    pData->recalc_hrtf_interpFLAG[index] = 1;
    pData->src_dirs_deg[index][1] = newElev_deg;
}

void binauraliser_setNumSources(void* const hBin, int new_nSources)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    pData->new_nSources = new_nSources > MAX_NUM_INPUTS ? MAX_NUM_INPUTS : new_nSources;
    if(pData->nSources != pData->new_nSources)
        pData->reInitTFT = 1;
}

void binauraliser_setUseDefaultHRIRsflag(void* const hBin, int newState)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    if((!pData->useDefaultHRIRsFLAG) && (newState)){
        pData->useDefaultHRIRsFLAG = newState;
        pData->reInitHRTFsAndGainTables = 1;
    }
}

void binauraliser_setSofaFilePath(void* const hBin, const char* path)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    
    pData->sofa_filepath = malloc(strlen(path) + 1);
    strcpy(pData->sofa_filepath, path);
    pData->useDefaultHRIRsFLAG = 0;
    pData->reInitHRTFsAndGainTables = 1; 
}

void binauraliser_setInputConfigPreset(void* const hBin, int newPresetID)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    int ch;
    
    binauraliser_loadPreset(newPresetID, pData->src_dirs_deg, &(pData->new_nSources), &(pData->input_nDims));
    if(pData->nSources != pData->new_nSources)
        pData->reInitTFT = 1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++){
        pData->recalc_hrtf_interpFLAG[ch] = 1;
    }
}


/* Get Functions */

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
    return pData->hrir_len;
}

int binauraliser_getHRIRsamplerate(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->hrir_fs;
}

int binauraliser_getUseDefaultHRIRsflag(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->useDefaultHRIRsFLAG;
}

char* binauraliser_getSofaFilePath(void* const hCmp)
{
    binauraliser_data *pData = (binauraliser_data*)(hCmp);
    if(pData->sofa_filepath!=NULL)
        return pData->sofa_filepath;
    else
        return "no_file";
}

int binauraliser_getDAWsamplerate(void* const hBin)
{
    binauraliser_data *pData = (binauraliser_data*)(hBin);
    return pData->fs;
} 


    
    
