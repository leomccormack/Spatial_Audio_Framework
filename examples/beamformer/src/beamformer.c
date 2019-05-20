/*
 Copyright 2019 Leo McCormack
 
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
 *     beamformer.c
 * Description:
 *     Generates beamformers/virtual microphones in arbitrary directions. Several
 *     different beam pattern types are included.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_sh
 * Author, date created:
 *     Leo McCormack, 17.05.2019
 */
 
#include "beamformer_internal.h"

void beamformer_create
(
    void ** const phBeam
)
{
    beamformer_data* pData = (beamformer_data*)malloc(sizeof(beamformer_data));
    if (pData == NULL) { return;/*error*/ }
    *phBeam = (void*)pData;
    int i, j, ch, band;

    /* default user parameters */
    pData->beamOrder = pData->new_beamOrder = 1;
    pData->nSH = pData->new_nSH = (pData->beamOrder + 1)*(pData->beamOrder + 1);
    for(i=0; i<MAX_NUM_BEAMS; i++){
        pData->beam_dirs_deg[i][0] = default_LScoords64_rad[i][0]*180.0f/M_PI;
        pData->beam_dirs_deg[i][1] = default_LScoords64_rad[i][0]*180.0f/M_PI;
    }
    
    
    beamformer_loadPreset(PRESET_T_DESIGN_24, pData->loudpkrs_dirs_deg, &(pData->new_nLoudpkrs), &(pData->loudpkrs_nDims));
    pData->nLoudpkrs = pData->new_nLoudpkrs;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D; 
    pData->rE_WEIGHT[0] = 0;
    pData->rE_WEIGHT[1] = 1;
    pData->diffEQmode[0] = AMPLITUDE_PRESERVING;
    pData->diffEQmode[1] = ENERGY_PRESERVING;
    pData->transitionFreq = 1000.0f;
    
    /* afSTFT stuff */
    pData->hSTFT = NULL;
    pData->STFTInputFrameTF = malloc(MAX_NUM_SH_SIGNALS * sizeof(complexVector));
    for(ch=0; ch< MAX_NUM_SH_SIGNALS; ch++) {
        pData->STFTInputFrameTF[ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
        pData->STFTInputFrameTF[ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
    }
    pData->tempHopFrameTD = (float**)malloc2d( MAX(MAX_NUM_SH_SIGNALS, MAX_NUM_LOUDSPEAKERS), HOP_SIZE, sizeof(float));
    pData->STFTOutputFrameTF = malloc(MAX_NUM_LOUDSPEAKERS * sizeof(complexVector));
    for(ch=0; ch< MAX_NUM_LOUDSPEAKERS; ch++) {
        pData->STFTOutputFrameTF[ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
        pData->STFTOutputFrameTF[ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
    }
    
    /* codec data */
    pData->pars = (codecPars*)malloc(sizeof(codecPars));
    codecPars* pars = pData->pars;
    for (i=0; i<NUM_DECODERS; i++){
        for(j=0; j<MAX_SH_ORDER; j++){
            pars->M_dec[i][j] = NULL;
            pars->M_dec_cmplx[i][j] = NULL;
            pars->M_dec_maxrE[i][j] = NULL;
            pars->M_dec_cmplx_maxrE[i][j] = NULL;
        }
    }
    pars->sofa_filepath = NULL;
    pars->hrirs = NULL;
    pars->hrir_dirs_deg = NULL;
    pars->hrtf_vbap_gtableIdx = NULL;
    pars->hrtf_vbap_gtableComp = NULL;
    pars->itds_s = NULL;
    pars->hrtf_fb = NULL;
    pars->hrtf_fb_mag = NULL;
    
    /* internal parameters */
    pData->new_binauraliseLS = 0;
    pData->binauraliseLS = pData->new_binauraliseLS;
    
    /* flags */
    pData->reInitCodec = 1;
    pData->reInitTFT = 1;
    pData->reInitHRTFs = 1;
    for(ch=0; ch<MAX_NUM_LOUDSPEAKERS; ch++)
        pData->recalc_hrtf_interpFLAG[ch] = 1; 
}

void beamformer_destroy
(
    void ** const phBeam
)
{
    beamformer_data *pData = (beamformer_data*)(*phBeam);
    codecPars *pars = pData->pars;
    int i, j, ch;
    
    if (pData != NULL) {
        if(pData->hSTFT!=NULL)
            afSTFTfree(pData->hSTFT);
        if(pData->STFTInputFrameTF!=NULL){
            for (ch = 0; ch< MAX_NUM_SH_SIGNALS; ch++) {
                free(pData->STFTInputFrameTF[ch].re);
                free(pData->STFTInputFrameTF[ch].im);
            }
        }
        if((pData->STFTOutputFrameTF!=NULL)){
            for (ch = 0; ch < MAX_NUM_LOUDSPEAKERS; ch++) {
                free(pData->STFTOutputFrameTF[ch].re);
                free(pData->STFTOutputFrameTF[ch].im);
            }
        }
        free(pData->STFTInputFrameTF);
        free(pData->STFTOutputFrameTF);
        if(pData->binauraliseLS && (pData->tempHopFrameTD!=NULL) )
            free2d((void**)pData->tempHopFrameTD, MAX(NUM_EARS, MAX_NUM_SH_SIGNALS));
        else if(pData->tempHopFrameTD!=NULL)
            free2d((void**)pData->tempHopFrameTD, MAX(MAX_NUM_LOUDSPEAKERS, MAX_NUM_SH_SIGNALS));
        
        if(pars->hrtf_vbap_gtableComp!= NULL)
            free(pars->hrtf_vbap_gtableComp);
        if(pars->hrtf_vbap_gtableIdx!= NULL)
            free(pars->hrtf_vbap_gtableIdx);
        if(pars->hrtf_fb!= NULL)
            free(pars->hrtf_fb);
        if(pars->hrtf_fb_mag!= NULL)
            free(pars->hrtf_fb_mag);
        if(pars->itds_s!= NULL)
            free(pars->itds_s);
        if(pars->hrirs!= NULL)
            free(pars->hrirs);
        if(pars->hrir_dirs_deg!= NULL)
            free(pars->hrir_dirs_deg);
        for (i=0; i<NUM_DECODERS; i++){
            for(j=0; j<MAX_SH_ORDER; j++){
                free(pars->M_dec[i][j]);
                free(pars->M_dec_cmplx[i][j]);
                free(pars->M_dec_maxrE[i][j]);
                free(pars->M_dec_cmplx_maxrE[i][j]);
            }
        }

        free(pData);
        pData = NULL;
    }
}

void beamformer_init
(
    void * const hBeam,
    int          sampleRate
)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int band;
    
    /* define frequency vector */
    pData->fs = sampleRate;
    for(band=0; band <HYBRID_BANDS; band++){
        if(sampleRate == 44100)
            pData->freqVector[band] =  (float)__afCenterFreq44100[band];
        else /* Assume 48kHz */
            pData->freqVector[band] =  (float)__afCenterFreq48e3[band];
    } 

    /* reinitialise if needed */
    beamformer_checkReInit(hBeam);
}

void beamformer_process
(
    void  *  const hBeam,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples,
    int            isPlaying
)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    codecPars* pars = pData->pars;
    int n, t, ch, ear, i, band, orderBand, nSH_band, decIdx;
    int o[MAX_SH_ORDER+2];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    /* local copies of user parameters */
    int nLoudspeakers, binauraliseLS, masterOrder;
    int orderPerBand[HYBRID_BANDS], rE_WEIGHT[NUM_DECODERS];
    float transitionFreq;
    DIFFUSE_FIELD_EQ_APPROACH diffEQmode[NUM_DECODERS];
    NORM_TYPES norm;
    
    /* reinitialise if needed */
#ifdef __APPLE__
    beamformer_checkReInit(hBeam);
#else
    if(pData->reInitTFT==1){
        pData->reInitTFT = 2;
        beamformer_initTFT(hBeam); /* always init before codec or hrtfs  */
        pData->reInitTFT = 0;
    }  
    if(pData->reInitCodec==1){
        pData->reInitCodec = 2;
        beamformer_initCodec(hBeam);
        pData->reInitCodec = 0;
    }
#endif

    /* decode audio to loudspeakers or headphones */
    if( (nSamples == FRAME_SIZE) && (pData->reInitCodec==0) && (pData->reInitTFT==0) && (pData->reInitHRTFs==0) ) {
        /* copy user parameters to local variables */
        for(n=0; n<MAX_SH_ORDER+2; n++){  o[n] = n*n;  }
        masterOrder = pData->masterOrder;
        nLoudspeakers = pData->nLoudpkrs;
        memcpy(orderPerBand, pData->orderPerBand, HYBRID_BANDS*sizeof(int));
        transitionFreq = pData->transitionFreq;
        memcpy(diffEQmode, pData->diffEQmode, NUM_DECODERS*sizeof(int));
        binauraliseLS = pData->binauraliseLS;
        norm = pData->norm;
        memcpy(rE_WEIGHT, pData->rE_WEIGHT, NUM_DECODERS*sizeof(int));
        
        /* Load time-domain data */
        for(i=0; i < MIN(pData->nSH, nInputs); i++)
            utility_svvcopy(inputs[i], FRAME_SIZE, pData->SHFrameTD[i]);
        for(; i<pData->nSH; i++)
            memset(pData->SHFrameTD[i], 0, FRAME_SIZE * sizeof(float));
        
        /* account for input normalisation scheme */
        switch(norm){
            case NORM_N3D:  /* already in N3D, do nothing */
                break;
            case NORM_SN3D: /* convert to N3D */
                for (n = 0; n<masterOrder+1; n++)
                    for (ch = o[n]; ch<o[n+1]; ch++)
                        for(i = 0; i<FRAME_SIZE; i++)
                            pData->SHFrameTD[ch][i] *= sqrtf(2.0f*(float)n+1.0f);
                break;
        }
        
        /* Apply time-frequency transform (TFT) */
        for(t=0; t< TIME_SLOTS; t++) {
            for(ch = 0; ch < pData->nSH; ch++)
                utility_svvcopy(&(pData->SHFrameTD[ch][t*HOP_SIZE]), HOP_SIZE, pData->tempHopFrameTD[ch]);
            afSTFTforward(pData->hSTFT, (float**)pData->tempHopFrameTD, (complexVector*)pData->STFTInputFrameTF);
            for(band=0; band<HYBRID_BANDS; band++)
                for(ch=0; ch < pData->nSH; ch++)
                    pData->SHframeTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[ch].re[band], pData->STFTInputFrameTF[ch].im[band]);
        }
        
        /* Main processing: */
        if(isPlaying){
            /* Decode to loudspeaker set-up */
            memset(pData->outputframeTF, 0, HYBRID_BANDS*MAX_NUM_LOUDSPEAKERS*TIME_SLOTS*sizeof(float_complex));
            for(band=0; band<HYBRID_BANDS; band++){
                orderBand = MAX(MIN(orderPerBand[band], masterOrder),1);
                nSH_band = (orderBand+1)*(orderBand+1);
                decIdx = pData->freqVector[band] < transitionFreq ? 0 : 1; /* different decoder for low (0) and high (1) frequencies */
                if(rE_WEIGHT[decIdx]){
                    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nLoudspeakers, TIME_SLOTS, nSH_band, &calpha,
                                pars->M_dec_cmplx_maxrE[decIdx][orderBand-1], nSH_band,
                                pData->SHframeTF[band], TIME_SLOTS, &cbeta,
                                pData->outputframeTF[band], TIME_SLOTS);
                }
                else{
                    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nLoudspeakers, TIME_SLOTS, nSH_band, &calpha,
                                pars->M_dec_cmplx[decIdx][orderBand-1], nSH_band,
                                pData->SHframeTF[band], TIME_SLOTS, &cbeta,
                                pData->outputframeTF[band], TIME_SLOTS);
                }
                for(i=0; i<nLoudspeakers; i++){
                    for(t=0; t<TIME_SLOTS; t++){
                        if(diffEQmode[decIdx]==AMPLITUDE_PRESERVING)
                            pData->outputframeTF[band][i][t] = crmulf(pData->outputframeTF[band][i][t], pars->M_norm[decIdx][orderBand-1][0]);
                        else
                            pData->outputframeTF[band][i][t] = crmulf(pData->outputframeTF[band][i][t], pars->M_norm[decIdx][orderBand-1][1]);
                    }
                }
            }
            
            /* binauralise the loudspeaker signals */
            if(binauraliseLS){
                memset(pData->binframeTF, 0, HYBRID_BANDS*NUM_EARS*TIME_SLOTS * sizeof(float_complex));
                /* interpolate hrtfs and apply to each source */
                for (ch = 0; ch < nLoudspeakers; ch++) {
                    if(pData->recalc_hrtf_interpFLAG[ch]){
                        beamformer_interpHRTFs(hBeam, pData->loudpkrs_dirs_deg[ch][0], pData->loudpkrs_dirs_deg[ch][1], pars->hrtf_interp[ch]);
                        pData->recalc_hrtf_interpFLAG[ch] = 0;
                    }
                    for (band = 0; band < HYBRID_BANDS; band++)
                        for (ear = 0; ear < NUM_EARS; ear++)
                            for (t = 0; t < TIME_SLOTS; t++)
                                pData->binframeTF[band][ear][t] = ccaddf(pData->binframeTF[band][ear][t], ccmulf(pData->outputframeTF[band][ch][t], pars->hrtf_interp[ch][band][ear]));
                }
                
                /* scale by sqrt(number of loudspeakers) */
                for (band = 0; band < HYBRID_BANDS; band++)
                    for (ear = 0; ear < NUM_EARS; ear++)
                        for (t = 0; t < TIME_SLOTS; t++)
                            pData->binframeTF[band][ear][t] = crmulf(pData->binframeTF[band][ear][t], 1.0f/sqrtf((float)nLoudspeakers));
            }
        }
        else{
            memset(pData->outputframeTF, 0, HYBRID_BANDS*MAX_NUM_LOUDSPEAKERS*TIME_SLOTS*sizeof(float_complex));
            memset(pData->binframeTF, 0, HYBRID_BANDS*NUM_EARS*TIME_SLOTS*sizeof(float_complex));
        }
        
        /* inverse-TFT */
        for(t = 0; t < TIME_SLOTS; t++) {
            for(band = 0; band < HYBRID_BANDS; band++) {
                if(binauraliseLS){
                    for (ch = 0; ch < NUM_EARS; ch++) {
                        pData->STFTOutputFrameTF[ch].re[band] = crealf(pData->binframeTF[band][ch][t]);
                        pData->STFTOutputFrameTF[ch].im[band] = cimagf(pData->binframeTF[band][ch][t]);
                    }
                }
                else {
                    for(ch = 0; ch < nLoudspeakers; ch++) {
                        pData->STFTOutputFrameTF[ch].re[band] = crealf(pData->outputframeTF[band][ch][t]);
                        pData->STFTOutputFrameTF[ch].im[band] = cimagf(pData->outputframeTF[band][ch][t]);
                    }
                }
            }
            afSTFTinverse(pData->hSTFT, pData->STFTOutputFrameTF, pData->tempHopFrameTD);
            for(ch = 0; ch < MIN(binauraliseLS==1 ? NUM_EARS : nLoudspeakers, nOutputs); ch++)
                utility_svvcopy(pData->tempHopFrameTD[ch], HOP_SIZE, &(outputs[ch][t* HOP_SIZE]));
            for (; ch < nOutputs; ch++)
                memset(&(outputs[ch][t* HOP_SIZE]), 0, HOP_SIZE*sizeof(float));
        }
    }
    else
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch], 0, FRAME_SIZE*sizeof(float));
}


/* Set Functions */

void beamformer_refreshSettings(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    pData->reInitCodec = 1;
    pData->reInitTFT = 1;
    pData->reInitHRTFs = 1;
}

void beamformer_checkReInit(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    /* reinitialise if needed */
    if (pData->reInitTFT == 1) {
        pData->reInitTFT = 2;
        beamformer_initTFT(hBeam); /* always init before codec or hrtfs  */
        pData->reInitTFT = 0;
    }
    if (pData->reInitCodec == 1) {
        pData->reInitCodec = 2;
        beamformer_initCodec(hBeam);
        pData->reInitCodec = 0;
    }
    if (pData->reInitHRTFs == 1) {
        pData->reInitHRTFs = 2;
        beamformer_initHRTFs(hBeam);
        pData->reInitHRTFs = 0;
    }
}

void beamformer_setMasterDecOrder(void  * const hBeam, int newValue)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    pData->new_masterOrder = MIN(MAX(newValue,1), MAX_SH_ORDER);
    pData->new_nSH = (pData->new_masterOrder+1)*(pData->new_masterOrder+1);
    pData->reInitTFT = 1;
    pData->reInitCodec = 1;
}

void beamformer_setDecOrder(void  * const hBeam, int newValue, int bandIdx)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    pData->orderPerBand[bandIdx] = MIN(MAX(newValue,1), pData->new_masterOrder);
}

void beamformer_setDecOrderAllBands(void  * const hBeam, int newValue)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int band;
    
    for(band=0; band<HYBRID_BANDS; band++)
        pData->orderPerBand[band] = MIN(MAX(newValue,1), pData->new_masterOrder);
}

void beamformer_setLoudspeakerAzi_deg(void* const hBeam, int index, float newAzi_deg)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = MAX(newAzi_deg, -180.0f);
    newAzi_deg = MIN(newAzi_deg, 180.0f);
    pData->loudpkrs_dirs_deg[index][0] = newAzi_deg;
    pData->reInitCodec=1;
    pData->recalc_hrtf_interpFLAG[index] = 1;
}

void beamformer_setLoudspeakerElev_deg(void* const hBeam, int index, float newElev_deg)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    newElev_deg = MAX(newElev_deg, -90.0f);
    newElev_deg = MIN(newElev_deg, 90.0f);
    pData->loudpkrs_dirs_deg[index][1] = newElev_deg;
    pData->reInitCodec=1;
    pData->recalc_hrtf_interpFLAG[index] = 1;
}

void beamformer_setNumLoudspeakers(void* const hBeam, int new_nLoudspeakers)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int ch; 
    pData->new_nLoudpkrs = new_nLoudspeakers > MAX_NUM_LOUDSPEAKERS ? MAX_NUM_LOUDSPEAKERS : new_nLoudspeakers;
    pData->new_nLoudpkrs = MAX(MIN_NUM_LOUDSPEAKERS, pData->new_nLoudpkrs);
    if(pData->nLoudpkrs != pData->new_nLoudpkrs){
        pData->reInitTFT = 1;
        pData->reInitCodec=1;
        for(ch=0; ch<MAX_NUM_LOUDSPEAKERS; ch++)
            pData->recalc_hrtf_interpFLAG[ch] = 1;
    }
}

void beamformer_setBinauraliseLSflag(void* const hBeam, int newState)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    
    pData->new_binauraliseLS = newState;
    if(pData->new_binauraliseLS != pData->binauraliseLS)
        pData->reInitTFT = 1;
    if(pData->new_binauraliseLS)
        pData->reInitHRTFs = 1;
}

void beamformer_setUseDefaultHRIRsflag(void* const hBeam, int newState)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    
    if((!pData->useDefaultHRIRsFLAG) && (newState)){
        pData->useDefaultHRIRsFLAG = newState;
        pData->reInitHRTFs = 1;
    }
}

void beamformer_setSofaFilePath(void* const hBeam, const char* path)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    codecPars* pars = pData->pars;
    
    pars->sofa_filepath = malloc(strlen(path) + 1);
    strcpy(pars->sofa_filepath, path);
    pData->useDefaultHRIRsFLAG = 0;
    pData->reInitHRTFs = 1;
}

void beamformer_setOutputConfigPreset(void* const hBeam, int newPresetID)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int ch;
    
    beamformer_loadPreset(newPresetID, pData->loudpkrs_dirs_deg, &(pData->new_nLoudpkrs), &(pData->loudpkrs_nDims));
    if(pData->nLoudpkrs != pData->new_nLoudpkrs)
        pData->reInitTFT = 1;
    pData->reInitCodec = 1;
    for(ch=0; ch<MAX_NUM_LOUDSPEAKERS; ch++)
        pData->recalc_hrtf_interpFLAG[ch] = 1;
}

void beamformer_setSourcePreset(void* const hBeam, int newPresetID)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int band, rangeIdx, curOrder, reverse;
    
    rangeIdx = 0;
    curOrder = 1;
    reverse = 0;
    switch(newPresetID){
        /* Ideal spherical harmonics will have SH_ORDER at all frequencies */
        case MIC_PRESET_IDEAL:
            for(band=0; band<HYBRID_BANDS; band++)
                pData->orderPerBand[band] = pData->masterOrder;
            break;
            
        /* For real microphone arrays, the maximum usable spherical harmonic order will depend on frequency  */
#ifdef ENABLE_ZYLIA_MIC_PRESET
        case MIC_PRESET_ZYLIA:
            for(band=0; band<HYBRID_BANDS; band++){
                if(rangeIdx<2*(__Zylia_maxOrder-1)){
                    if(pData->freqVector[band]>__Zylia_freqRange[rangeIdx]){
                        if(!reverse)
                            curOrder++;
                        else
                            curOrder--;
                        reverse = (curOrder == __Zylia_maxOrder) || (reverse) ? 1 : 0;
                        rangeIdx++;
                    }
                }
                pData->orderPerBand[band] = MIN(pData->masterOrder,curOrder);
            }
            break;
#endif
            
#ifdef ENABLE_EIGENMIKE32_MIC_PRESET
        case MIC_PRESET_EIGENMIKE32:
            for(band=0; band<HYBRID_BANDS; band++){
                if(rangeIdx<2*(__Eigenmike32_maxOrder-1)){
                    if(pData->freqVector[band]>__Eigenmike32_freqRange[rangeIdx]){
                        if(!reverse)
                            curOrder++;
                        else
                            curOrder--;
                        reverse = (curOrder == __Eigenmike32_maxOrder) || (reverse) ? 1 : 0;
                        rangeIdx++;
                    }
                }
                pData->orderPerBand[band] = MIN(pData->masterOrder,curOrder);
            }
            break;
#endif
            
#ifdef ENABLE_DTU_MIC_MIC_PRESET
        case MIC_PRESET_DTU_MIC:
            for(band=0; band<HYBRID_BANDS; band++){
                if(rangeIdx<2*(__DTU_mic_maxOrder-1)){
                    if(pData->freqVector[band]>__DTU_mic_freqRange[rangeIdx]){
                        if(!reverse)
                            curOrder++;
                        else
                            curOrder--;
                        reverse = (curOrder == __DTU_mic_maxOrder) || (reverse) ? 1 : 0;
                        rangeIdx++;
                    }
                }
                pData->orderPerBand[band] = MIN(pData->masterOrder,curOrder);
            }
            break;
#endif
    }
}

void beamformer_setChOrder(void* const hBeam, int newOrder)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    pData->chOrdering = (CH_ORDER)newOrder;
}

void beamformer_setNormType(void* const hBeam, int newType)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    pData->norm = (NORM_TYPES)newType;
}

void beamformer_setDecMethod(void* const hBeam, int index, int newID)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    pData->dec_method[index] = newID;
    pData->reInitCodec = 1;
}

void beamformer_setDecEnableMaxrE(void* const hBeam, int index, int newID)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    pData->rE_WEIGHT[index] = newID;
}

void beamformer_setDecNormType(void* const hBeam, int index, int newID)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    pData->diffEQmode[index] = newID;
}

void beamformer_setTransitionFreq(void* const hBeam, float newValue)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    pData->transitionFreq = newValue;
}


/* Get Functions */

int beamformer_getMasterDecOrder(void  * const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->new_masterOrder;
}

int beamformer_getDecOrder(void  * const hBeam, int bandIdx)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->orderPerBand[bandIdx];
}

int beamformer_getDecOrderAllBands(void  * const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->orderPerBand[0];
}

void beamformer_getDecOrderHandle
(
    void* const hBeam,
    float** pX_vector,
    int** pY_values,
    int* pNpoints
)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    (*pX_vector) = &pData->freqVector[0];
    (*pY_values) = &pData->orderPerBand[0];
    (*pNpoints) = HYBRID_BANDS;
}

int beamformer_getNumberOfBands(void)
{
    return HYBRID_BANDS;
}

float beamformer_getLoudspeakerAzi_deg(void* const hBeam, int index)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->loudpkrs_dirs_deg[index][0];
}

float beamformer_getLoudspeakerElev_deg(void* const hBeam, int index)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->loudpkrs_dirs_deg[index][1];
}

int beamformer_getNumLoudspeakers(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->new_nLoudpkrs;
}

int beamformer_getMaxNumLoudspeakers()
{
    return MAX_NUM_LOUDSPEAKERS;
}

int  beamformer_getNSHrequired(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->nSH;
}

int beamformer_getBinauraliseLSflag(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->new_binauraliseLS;
}

int beamformer_getUseDefaultHRIRsflag(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->useDefaultHRIRsFLAG;
}

char* beamformer_getSofaFilePath(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    codecPars* pars = pData->pars;
    if(pars->sofa_filepath!=NULL)
        return pars->sofa_filepath;
    else
        return "no_file";
}

int beamformer_getChOrder(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return (int)pData->chOrdering;
}

int beamformer_getNormType(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return (int)pData->norm;
}

int beamformer_getDecMethod(void* const hBeam, int index)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->dec_method[index];
}

int beamformer_getDecEnableMaxrE(void* const hBeam, int index)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->rE_WEIGHT[index];
}

int beamformer_getDecNormType(void* const hBeam, int index)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->diffEQmode[index];
}

float beamformer_getTransitionFreq(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->transitionFreq;
}

int beamformer_getHRIRsamplerate(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    codecPars* pars = pData->pars;
    return pars->hrir_fs;
}

int beamformer_getDAWsamplerate(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->fs;
}

int beamformer_getProcessingDelay()
{
    return 12*HOP_SIZE;
}



