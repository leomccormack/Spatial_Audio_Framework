    /*
 Copyright 2016-2018 Leo McCormack
 
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
 *     powermap.c
 * Description:
 *     A powermap-based sound-field visualiser, which utilises spherical harmonic
 *     signals as input.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 26.04.2016
 */

#include "powermap.h"
#include "powermap_internal.h"
#include "powermap_database.h"

void powermap_create
(
    void ** const phPm
)
{
    powermap_data* pData = (powermap_data*)malloc(sizeof(powermap_data));
    if (pData == NULL) { return;/*error*/ }
    *phPm = (void*)pData;
    int n, i, t, ch, band;
    
    afSTFTinit(&(pData->hSTFT), HOP_SIZE, MAX_NUM_SH_SIGNALS, 0, 0, 1);
    pData->STFTInputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, MAX_NUM_SH_SIGNALS, sizeof(complexVector));
    for(t=0; t<TIME_SLOTS; t++) {
        for(ch=0; ch< MAX_NUM_SH_SIGNALS; ch++) {
            pData->STFTInputFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
            pData->STFTInputFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
        }
    }
    pData->tempHopFrameTD = (float**)malloc2d(MAX_NUM_SH_SIGNALS, HOP_SIZE, sizeof(float));
    
    /* codec data */
    pData->pars = (codecPars*)malloc(sizeof(codecPars));
    codecPars* pars = pData->pars;
    pars->interp_dirs_deg = NULL;
    for(n=0; n<MAX_SH_ORDER; n++){
        pars->Y_grid[n] = NULL;
        pars->Y_grid_cmplx[n] = NULL;
    }
    pars->interp_table = NULL;
    
    /* internal */
	pData->reInitTFT = 1;
	pData->reInitAna = 1;
    pData->dispWidth = 140;

    /* display */
    pData->pmap = NULL;
    pData->prev_pmap = NULL;
    for(i=0; i<NUM_DISP_SLOTS; i++)
        pData->pmap_grid[i] = NULL;
    pData->pmapReady = 0;
    pData->recalcPmap = 1;
    
    /* Default user parameters */
    pData->masterOrder = pData->new_masterOrder = MASTER_ORDER_FIRST;
    pData->nSH = pData->new_nSH = (pData->masterOrder+1)*(pData->masterOrder+1);
    for(band=0; band<HYBRID_BANDS; band++){
        pData->analysisOrderPerBand[band] = pData->masterOrder;
        pData->pmapEQ[band] = 1.0f;
    }
    pData->covAvgCoeff = 0.0f;
    pData->pmapAvgCoeff = 0.666f;
    pData->nSources = 1;
    pData->pmap_mode = PM_MODE_MUSIC;
    pData->HFOVoption = HFOV_360;
    pData->aspectRatioOption = ASPECT_RATIO_2_1;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
}

void powermap_destroy
(
    void ** const phPm
)
{
    powermap_data *pData = (powermap_data*)(*phPm);
    codecPars* pars = pData->pars;
    int i, t, ch;
    
    if (pData != NULL) {
        afSTFTfree(pData->hSTFT);
        for (t = 0; t<TIME_SLOTS; t++) {
            for (ch = 0; ch< MAX_NUM_SH_SIGNALS; ch++) {
                free(pData->STFTInputFrameTF[t][ch].re);
                free(pData->STFTInputFrameTF[t][ch].im);
            }
        }
        free2d((void**)pData->STFTInputFrameTF, TIME_SLOTS);
        free2d((void**)pData->tempHopFrameTD, MAX_NUM_SH_SIGNALS);
        
        if(pData->pmap!=NULL)
            free(pData->pmap);
        if(pData->prev_pmap!=NULL)
            free(pData->prev_pmap);
        for(i=0; i<NUM_DISP_SLOTS; i++)
            if(pData->pmap_grid[i] !=NULL)
                free(pData->pmap_grid[i]);
        if(pars->interp_dirs_deg!=NULL)
            free(pars->interp_dirs_deg);
        for(i=0; i<MAX_SH_ORDER; i++){
            if(pars->Y_grid[i] !=NULL)
                free(pars->Y_grid[i]);
            if(pars->Y_grid_cmplx[i]!=NULL)
                free(pars->Y_grid_cmplx[i]);
        }
        if(pars->interp_table!=NULL)
            free(pars->interp_table);
        free(pData->pars);
        free(pData);
        pData = NULL;
    }
}

void powermap_init
(
    void * const hPm,
    float        sampleRate
)
{
    powermap_data *pData = (powermap_data*)(hPm);
    codecPars* pars = pData->pars;
    int band;
    
    pData->fs = sampleRate;
    
    /* specify frequency vector and determine the number of bands */
    switch ((int)(sampleRate+0.5f)){
        case 44100:
            for(band=0; band<HYBRID_BANDS; band++)
                pData->freqVector[band] = (float)__afCenterFreq44100[band];
            break;
        default:
        case 48000:
            for(band=0; band<HYBRID_BANDS; band++)
                pData->freqVector[band] = (float)__afCenterFreq48e3[band];
            break;
    }
    
    /* intialise parameters */
    memset(pData->Cx, 0 , MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*HYBRID_BANDS*sizeof(float_complex));
    if(pData->prev_pmap!=NULL)
        memset(pData->prev_pmap, 0, pars->grid_nDirs*sizeof(float));
    pData->pmapReady = 0;
    pData->dispSlotIdx = 0;

	/* reinitialise if needed */
	powermap_checkReInit(hPm);
}


void powermap_analysis
(
    void  *  const hPm,
    float ** const inputs,
    int            nInputs,
    int            nSamples,
    int            isPlaying
)
{
    powermap_data *pData = (powermap_data*)(hPm);
    codecPars* pars = pData->pars;
    int i, j, t, n, ch, sample, band, nSH_order, order_band, nSH_maxOrder, maxOrder;
    float C_grp_trace, covScale, pmapEQ_band;
    int o[MAX_SH_ORDER+2];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    float_complex new_Cx[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];
    float_complex* C_grp;
    
    /* local parameters */
    int analysisOrderPerBand[HYBRID_BANDS];
    int nSources, masterOrder, nSH;
    float covAvgCoeff, pmapAvgCoeff;
    float pmapEQ[HYBRID_BANDS];
    NORM_TYPES norm;
    POWERMAP_MODES pmap_mode;
    
    /* reinitialise if needed */
#ifdef __APPLE__
    powermap_checkReInit(hPm);
#else
    if(pData->reInitTFT==1){
        pData->reInitTFT = 2;
        powermap_initTFT(hPm);
        pData->reInitTFT = 0;
    }
#endif
    
    /* The main processing: */
    if (nSamples == FRAME_SIZE && (pData->reInitAna == 0) && (pData->reInitTFT == 0) && isPlaying ) {
        /* copy current parameters to be thread safe */
        memcpy(analysisOrderPerBand, pData->analysisOrderPerBand, HYBRID_BANDS*sizeof(int));
        memcpy(pmapEQ, pData->pmapEQ, HYBRID_BANDS*sizeof(float));
        norm = pData->norm;
        nSources = pData->nSources;
        covAvgCoeff = MIN(pData->covAvgCoeff, MAX_COV_AVG_COEFF);
        pmapAvgCoeff = pData->pmapAvgCoeff;
        pmap_mode = pData->pmap_mode;
        masterOrder = pData->masterOrder;
        nSH = pData->nSH;
        
        /* load intput time-domain data */
        for (i = 0; i < MIN(nSH, nInputs); i++)
            memcpy(pData->SHframeTD[i], inputs[i], FRAME_SIZE*sizeof(float));
        for (; i < nSH; i++)
            memset(pData->SHframeTD[i], 0, FRAME_SIZE*sizeof(float));
        
        /* account for input normalisation scheme */
        switch(norm){
            case NORM_N3D:  /* already in N3D, do nothing */
                break;
            case NORM_SN3D: /* convert to N3D */
                for(n=0; n<masterOrder+2; n++){  o[n] = n*n;  };
                for (n = 0; n<masterOrder+1; n++)
                    for (ch = o[n]; ch<o[n+1]; ch++)
                        for(i = 0; i<FRAME_SIZE; i++)
                            pData->SHframeTD[ch][i] *= sqrtf(2.0f*(float)n+1.0f);
                break;
        }
        
        /* apply the time-frequency transform */
        for (t = 0; t < TIME_SLOTS; t++) {
            for (ch = 0; ch < nSH; ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    pData->tempHopFrameTD[ch][sample] = pData->SHframeTD[ch][sample + t*HOP_SIZE];
            afSTFTforward(pData->hSTFT, (float**)pData->tempHopFrameTD, (complexVector*)pData->STFTInputFrameTF[t]);
        }
        for (band = 0; band < HYBRID_BANDS; band++)
            for (ch = 0; ch < nSH; ch++)
                for (t = 0; t < TIME_SLOTS; t++)
                    pData->SHframeTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[t][ch].re[band], pData->STFTInputFrameTF[t][ch].im[band]);

        /* Update covarience matrix per band */
        covScale = 1.0f/(float)(nSH);
        for(band=0; band<HYBRID_BANDS; band++){
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, nSH, TIME_SLOTS, &calpha,
                        pData->SHframeTF[band], TIME_SLOTS,
                        pData->SHframeTF[band], TIME_SLOTS, &cbeta,
                        new_Cx, MAX_NUM_SH_SIGNALS);

            /* scale with nSH */
            for(i=0; i<nSH; i++)
                for(j=0; j<nSH; j++)
                    new_Cx[i][j] = crmulf(new_Cx[i][j], covScale);

            /* average over time */
            for(i=0; i<nSH; i++)
                for(j=0; j<nSH; j++)
                    pData->Cx[band][i][j] = ccaddf( crmulf(new_Cx[i][j], 1.0f-covAvgCoeff), crmulf(pData->Cx[band][i][j], covAvgCoeff));
        }
        
        /* update the powermap */
        if(pData->recalcPmap==1){
            pData->recalcPmap = 0;
            pData->pmapReady = 0;

            /* determine maximum analysis order */
            maxOrder = 1;
            for(i=0; i<HYBRID_BANDS; i++)
                maxOrder = MAX(maxOrder, MIN(analysisOrderPerBand[i], masterOrder));
            nSH_maxOrder = (maxOrder+1)*(maxOrder+1);

            /* group covarience matrices */
            C_grp = calloc(nSH_maxOrder*nSH_maxOrder, sizeof(float_complex));
            for (band=0; band<HYBRID_BANDS; band++){
                order_band = MAX(MIN(pData->analysisOrderPerBand[band], masterOrder),1);
                nSH_order = (order_band+1)*(order_band+1);
                pmapEQ_band = MIN(MAX(pmapEQ[band], 0.0f), 2.0f);
                for(i=0; i<nSH_order; i++)
                    for(j=0; j<nSH_order; j++)
                        C_grp[i*nSH_maxOrder+j] = ccaddf(C_grp[i*nSH_maxOrder+j], crmulf(pData->Cx[band][i][j], 1e10f*pmapEQ_band));
            }

            /* generate powermap */
            C_grp_trace = 0.0f;
            for(i=0; i<nSH_maxOrder; i++)
                C_grp_trace+=crealf(C_grp[i*nSH_maxOrder+ i]);
            switch(pmap_mode){
                default:
                case PM_MODE_PWD:
                    generatePWDmap(maxOrder, C_grp, pars->Y_grid_cmplx[maxOrder-1], pars->grid_nDirs, pData->pmap);
                    break;

                case PM_MODE_MVDR:
                    if(C_grp_trace>1e-8f)
                        generateMVDRmap(maxOrder, C_grp, pars->Y_grid_cmplx[maxOrder-1], pars->grid_nDirs, 8.0f, pData->pmap, NULL);
                    else
                        memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                    break;

                case PM_MODE_CROPAC_LCMV:
                    if(C_grp_trace>1e-8f)
                        generateCroPaCLCMVmap(maxOrder, C_grp, pars->Y_grid_cmplx[maxOrder-1], pars->grid_nDirs, 8.0f, 0.0f, pData->pmap);
                    else
                        memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                    break;

                case PM_MODE_MUSIC:
                    if(C_grp_trace>1e-8f)
                        generateMUSICmap(maxOrder, C_grp, pars->Y_grid_cmplx[maxOrder-1], nSources, pars->grid_nDirs, 0, pData->pmap);
                    else
                        memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                    break;
                    
                case PM_MODE_MUSIC_LOG:
                    if(C_grp_trace>1e-8f)
                        generateMUSICmap(maxOrder, C_grp, pars->Y_grid_cmplx[maxOrder-1], nSources, pars->grid_nDirs, 1, pData->pmap);
                    else
                        memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                    break; 
                    
                case PM_MODE_MINNORM:
                    if(C_grp_trace>1e-8f)
                        generateMinNormMap(maxOrder, C_grp, pars->Y_grid_cmplx[maxOrder-1], nSources, pars->grid_nDirs, 0, pData->pmap);
                    else
                        memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                    break;
                    
                case PM_MODE_MINNORM_LOG:
                    if(C_grp_trace>1e-8f)
                        generateMinNormMap(maxOrder, C_grp, pars->Y_grid_cmplx[maxOrder-1], nSources, pars->grid_nDirs, 1, pData->pmap);
                    else
                        memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                    break;
            }
            free(C_grp);
            
            /* average powermap over time */
            for(i=0; i<pars->grid_nDirs; i++)
                pData->pmap[i] =  (1.0f-pmapAvgCoeff) * (pData->pmap[i] )+ pmapAvgCoeff * (pData->prev_pmap[i]);
            memcpy(pData->prev_pmap,  pData->pmap , pars->grid_nDirs*sizeof(float));

            /* interpolate powermap */
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->interp_nDirs, 1, pars->grid_nDirs, 1.0f,
                        pars->interp_table, pars->grid_nDirs,
                        pData->pmap, 1, 0.0f,
                        pData->pmap_grid[pData->dispSlotIdx], 1);

            /* ascertain minimum and maximum values for powermap colour scaling */
            pData->pmap_grid_minVal = FLT_MAX;
            pData->pmap_grid_maxVal = FLT_MIN;
            for(i=0; i<pars->interp_nDirs; i++){
                pData->pmap_grid_minVal = pData->pmap_grid[pData->dispSlotIdx][i] < pData->pmap_grid_minVal ? pData->pmap_grid[pData->dispSlotIdx][i] : pData->pmap_grid_minVal;
                pData->pmap_grid_maxVal = pData->pmap_grid[pData->dispSlotIdx][i] > pData->pmap_grid_maxVal ? pData->pmap_grid[pData->dispSlotIdx][i] : pData->pmap_grid_maxVal;
            }

            /* normalise the powermap to 0..1 */
            for(i=0; i<pars->interp_nDirs; i++)
                pData->pmap_grid[pData->dispSlotIdx][i] = (pData->pmap_grid[pData->dispSlotIdx][i]-pData->pmap_grid_minVal)/(pData->pmap_grid_maxVal-pData->pmap_grid_minVal+1e-11f);

            /* signify that the powermap in current slot is ready for plotting */
            pData->dispSlotIdx++;
            if(pData->dispSlotIdx>=NUM_DISP_SLOTS)
                pData->dispSlotIdx = 0;
            pData->pmapReady = 1;
        }
    }
}

/* SETS */
 
void powermap_refreshSettings(void* const hPm)
{
	powermap_data *pData = (powermap_data*)(hPm);
	pData->reInitTFT = 1;
	pData->reInitAna = 1; 
}
 
void powermap_checkReInit(void* const hPm)
{
	powermap_data *pData = (powermap_data*)(hPm); 
	/* reinitialise if needed */
	if (pData->reInitTFT == 1) {
		pData->reInitTFT = 2;
		powermap_initTFT(hPm);
		pData->reInitTFT = 0;
	}
	if (pData->reInitAna == 1) {
		pData->reInitAna = 2;  /* indicate init in progress */
		pData->pmapReady = 0;  /* avoid trying to draw pmap during reinit */
		powermap_initAna(hPm);
		pData->reInitAna = 0;  /* indicate init complete */
		pData->recalcPmap = 1; /* recalculate powermap with new configuration */
	}
}

void powermap_setPowermapMode(void* const hPm, int newMode)
{
    powermap_data *pData = (powermap_data*)(hPm);
    codecPars* pars = pData->pars;
    pData->pmap_mode = (POWERMAP_MODES)newMode;
    if(pData->prev_pmap!=NULL)
        memset(pData->prev_pmap, 0, pars->grid_nDirs*sizeof(float));
}

void powermap_setMasterOrder(void* const hPm,  int newValue)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->new_masterOrder = newValue;
    pData->new_nSH = (newValue+1)*(newValue+1);
    pData->reInitTFT = 1;
    pData->reInitAna = 1;
}

void powermap_setCovAvgCoeff(void* const hPm, float newAvg)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->covAvgCoeff = MIN(MAX(0.0f, newAvg), 0.99999999f);
}

void powermap_setNumSources(void* const hPm, int newValue)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->nSources = newValue;
}

void powermap_setSourcePreset(void* const hPm, int newPresetID)
{
    powermap_data *pData = (powermap_data*)(hPm);
    int band, rangeIdx, curOrder, reverse;

    rangeIdx = 0;
    curOrder = 1;
    reverse = 0;
    switch(newPresetID){
        case MIC_PRESET_IDEAL:
            /* Ideal SH should have maximum order per frequency */
            for(band=0; band<HYBRID_BANDS; band++)
                pData->analysisOrderPerBand[band] = pData->new_masterOrder;
            break;
            
            /* In the case of real microphone arrays, the analysis order should be frequency dependent
            *  and the frequencies above the spatial-aliasing limit should be EQ's out. */
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
                pData->analysisOrderPerBand[band] = MIN(pData->new_masterOrder,curOrder);
                if(pData->freqVector[band] > __Zylia_freqRange[(__Zylia_maxOrder-1)*2 - 1])
                    pData->pmapEQ[band] = 0.0f;
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
                pData->analysisOrderPerBand[band] = MIN(pData->new_masterOrder,curOrder);
                if(pData->freqVector[band] > __Eigenmike32_freqRange[(__Eigenmike32_maxOrder-1)*2 - 1])
                    pData->pmapEQ[band] = 0.0f;
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
                pData->analysisOrderPerBand[band] = MIN(pData->new_masterOrder,curOrder);
                if(pData->freqVector[band] > __DTU_mic_freqRange[(__DTU_mic_maxOrder-1)*2 - 1])
                    pData->pmapEQ[band] = 0.0f;
            }
            break;
#endif
    }
}

void powermap_setAnaOrder(void  * const hPm, int newValue, int bandIdx)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->analysisOrderPerBand[bandIdx] = MIN(MAX(newValue,1), pData->new_masterOrder);
}

void powermap_setAnaOrderAllBands(void  * const hPm, int newValue)
{
    powermap_data *pData = (powermap_data*)(hPm);
    int band;

    for(band=0; band<HYBRID_BANDS; band++)
        pData->analysisOrderPerBand[band] = MIN(MAX(newValue,1), pData->new_masterOrder);
}

void powermap_setPowermapEQ(void  * const hPm, float newValue, int bandIdx)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->pmapEQ[bandIdx] = newValue;
}

void powermap_setPowermapEQAllBands(void  * const hPm, float newValue)
{
    powermap_data *pData = (powermap_data*)(hPm);
    int band;
    
    for(band=0; band<HYBRID_BANDS; band++)
        pData->pmapEQ[band] = newValue;
}

void powermap_setChOrder(void* const hPm, int newOrder)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->chOrdering = (CH_ORDER)newOrder;
}

void powermap_setNormType(void* const hPm, int newType)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->norm = (NORM_TYPES)newType;
}

void powermap_setDispFOV(void* const hPm, int newOption)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->HFOVoption = (HFOV_OPTIONS)newOption;
}

void powermap_setAspectRatio(void* const hPm, int newOption)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->aspectRatioOption = (ASPECT_RATIO_OPTIONS)newOption;
}

void powermap_setPowermapAvgCoeff(void* const hPm, float newValue)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->pmapAvgCoeff = MIN(MAX(0.0f, newValue), 0.99999999f);
}

void powermap_requestPmapUpdate(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->recalcPmap = 1;
}


/* GETS */

int powermap_getMasterOrder(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return pData->new_masterOrder;
}

int powermap_getPowermapMode(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return (int)pData->pmap_mode;
}

int powermap_getSamplingRate(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return (int)(pData->fs+0.5f);
}

float powermap_getCovAvgCoeff(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return pData->covAvgCoeff;
}

float powermap_getPowermapEQ(void  * const hPm, int bandIdx)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return pData->pmapEQ[bandIdx];
}

float powermap_getPowermapEQAllBands(void  * const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return pData->pmapEQ[0];
}

void powermap_getPowermapEQHandle
(
    void* const hPm,
    float** pX_vector,
    float** pY_values,
    int* pNpoints
)
{
    powermap_data *pData = (powermap_data*)(hPm);
    (*pX_vector) = &(pData->freqVector[0]);
    (*pY_values) = &(pData->pmapEQ[0]);
    (*pNpoints) = HYBRID_BANDS;
}

int powermap_getAnaOrder(void  * const hPm, int bandIdx)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return pData->analysisOrderPerBand[bandIdx];
}

int powermap_getAnaOrderAllBands(void  * const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return pData->analysisOrderPerBand[0];
}

void powermap_getAnaOrderHandle
(
    void* const hPm,
    float** pX_vector,
    int** pY_values,
    int* pNpoints
)
{
    powermap_data *pData = (powermap_data*)(hPm);
    (*pX_vector) = &(pData->freqVector[0]);
    (*pY_values) = &(pData->analysisOrderPerBand[0]);
    (*pNpoints) = HYBRID_BANDS;
}

int powermap_getNumberOfBands(void)
{
    return HYBRID_BANDS;
}

int powermap_getNSHrequired(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return pData->new_nSH;
}

int powermap_getChOrder(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return (int)pData->chOrdering;
}

int powermap_getNormType(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return (int)pData->norm;
}

int powermap_getNumSources(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return pData->nSources;
}

int powermap_getDispFOV(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return (int)pData->HFOVoption;
}

int powermap_getAspectRatio(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return (int)pData->aspectRatioOption;
}

float powermap_getPowermapAvgCoeff(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return pData->pmapAvgCoeff;
}

int powermap_getPmap(void* const hPm, float** grid_dirs, float** pmap, int* nDirs,int* pmapWidth, int* hfov, int* aspectRatio) //TODO: hfov and aspectRatio should be float, if 16:9 etc options are added
{
    powermap_data *pData = (powermap_data*)(hPm);
    codecPars* pars = pData->pars;
    if((pData->reInitAna == 0) && pData->pmapReady){
        (*grid_dirs) = pars->interp_dirs_deg;
        (*pmap) = pData->pmap_grid[pData->dispSlotIdx-1 < 0 ? NUM_DISP_SLOTS-1 : pData->dispSlotIdx-1];
        (*nDirs) = pars->interp_nDirs;
        (*pmapWidth) = pData->dispWidth;
        switch(pData->HFOVoption){
            default:
            case HFOV_360:
                (*hfov) = 360;
                break;
        }
        switch(pData->aspectRatioOption){
            default:
            case ASPECT_RATIO_2_1:
                (*aspectRatio) = 2;
                break;
        }
    }
    return pData->pmapReady;
}




