/*
 * Copyright 2016-2018 Leo McCormack
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
 * @file powermap.c
 * @brief A sound-field visualiser, which utilises spherical harmonic signals as
 *        input; note this code is a remnant from the work conducted in [1]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S. and Pulkki, V., 2017. Parametric
 *          acoustic camera for real-time sound capture, analysis and tracking.
 *          In Proceedings of the 20th International Conference on Digital Audio
 *          Effects (DAFx-17) (pp. 412-419)
 *
 * @author Leo McCormack
 * @date 26.04.2016
 * @license ISC
 */

#include "powermap.h"
#include "powermap_internal.h" 

void powermap_create
(
    void ** const phPm
)
{
    powermap_data* pData = (powermap_data*)malloc1d(sizeof(powermap_data));
    *phPm = (void*)pData;
    int n, i, band;

    /* Default user parameters */
    pData->masterOrder = pData->new_masterOrder = SH_ORDER_FIRST;
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
    
    afSTFT_create(&(pData->hSTFT), MAX_NUM_SH_SIGNALS, 0, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    pData->SHframeTD = (float**)malloc2d(MAX_NUM_SH_SIGNALS, POWERMAP_FRAME_SIZE, sizeof(float));
    pData->SHframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_SH_SIGNALS, TIME_SLOTS, sizeof(float_complex));

    /* codec data */
    pData->pars = (powermap_codecPars*)malloc1d(sizeof(powermap_codecPars));
    powermap_codecPars* pars = pData->pars;
    pars->interp_dirs_deg = NULL;
    for(n=0; n<MAX_SH_ORDER; n++){
        pars->Y_grid[n] = NULL;
        pars->Y_grid_cmplx[n] = NULL;
    }
    pars->interp_table = NULL;
    
    /* internal */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
    pData->dispWidth = 140;

    /* display */
    pData->pmap = NULL;
    pData->prev_pmap = NULL;
    for(i=0; i<NUM_DISP_SLOTS; i++)
        pData->pmap_grid[i] = NULL;
    pData->pmapReady = 0;
    pData->recalcPmap = 1;

    /* set FIFO buffer */
    pData->FIFO_idx = 0;
    memset(pData->inFIFO, 0, MAX_NUM_SH_SIGNALS*POWERMAP_FRAME_SIZE*sizeof(float));
}

void powermap_destroy
(
    void ** const phPm
)
{
    powermap_data *pData = (powermap_data*)(*phPm);
    powermap_codecPars* pars;
    int i;
    
    if (pData != NULL) {
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            SAF_SLEEP(10);
        }

        pars = pData->pars;

        /* free afSTFT and buffers */
        afSTFT_destroy(&(pData->hSTFT));
        free(pData->SHframeTD);
        free(pData->SHframeTF);
        
        free(pData->pmap);
        free(pData->prev_pmap);
        for(i=0; i<NUM_DISP_SLOTS; i++)
            free(pData->pmap_grid[i]);
        free(pars->interp_dirs_deg);
        for(i=0; i<MAX_SH_ORDER; i++){
            free(pars->Y_grid[i]);
            free(pars->Y_grid_cmplx[i]);
        }
        free(pars->interp_table);
        free(pData->pars);
        free(pData->progressBarText);
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
    powermap_codecPars* pars = pData->pars;
    
    pData->fs = sampleRate;
    
    /* specify frequency vector and determine the number of bands */
    afSTFT_getCentreFreqs(pData->hSTFT, sampleRate, HYBRID_BANDS, pData->freqVector);
    
    /* intialise parameters */
    memset(pData->Cx, 0 , MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*HYBRID_BANDS*sizeof(float_complex));
    if(pData->prev_pmap!=NULL)
        memset(pData->prev_pmap, 0, pars->grid_nDirs*sizeof(float));
    pData->pmapReady = 0;
    pData->dispSlotIdx = 0;
}

void powermap_initCodec
(
    void* const hPm
)
{
    powermap_data *pData = (powermap_data*)(hPm);
    
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
    
    powermap_initTFT(hPm);
    powermap_initAna(hPm);
    
    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
}

void powermap_analysis
(
    void        *  const hPm,
    const float *const * inputs,
    int                  nInputs,
    int                  nSamples,
    int                  isPlaying
)
{
    powermap_data *pData = (powermap_data*)(hPm);
    powermap_codecPars* pars = pData->pars;
    int s, i, j, ch, band, nSH_order, order_band, nSH_maxOrder, maxOrder;
    float C_grp_trace, pmapEQ_band;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    float_complex new_Cx[MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS];
    float_complex C_grp[MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS];
    
    /* local parameters */
    int analysisOrderPerBand[HYBRID_BANDS];
    int nSources, masterOrder, nSH;
    float covAvgCoeff, pmapAvgCoeff;
    float pmapEQ[HYBRID_BANDS];
    NORM_TYPES norm;
    CH_ORDER chOrdering;
    POWERMAP_MODES pmap_mode;
    memcpy(analysisOrderPerBand, pData->analysisOrderPerBand, HYBRID_BANDS*sizeof(int));
    memcpy(pmapEQ, pData->pmapEQ, HYBRID_BANDS*sizeof(float));
    norm = pData->norm;
    chOrdering = pData->chOrdering;
    nSources = pData->nSources;
    covAvgCoeff = SAF_MIN(pData->covAvgCoeff, MAX_COV_AVG_COEFF);
    pmapAvgCoeff = pData->pmapAvgCoeff;
    pmap_mode = pData->pmap_mode;
    masterOrder = pData->masterOrder;
    nSH = (masterOrder+1)*(masterOrder+1);

    /* Loop over all samples */
    for(s=0; s<nSamples; s++){
        /* Load input signals into inFIFO buffer */
        for(ch=0; ch<SAF_MIN(nInputs,nSH); ch++)
            pData->inFIFO[ch][pData->FIFO_idx] = inputs[ch][s];
        for(; ch<nSH; ch++) /* Zero any channels that were not given */
            pData->inFIFO[ch][pData->FIFO_idx] = 0.0f;

        /* Increment buffer index */
        pData->FIFO_idx++;

        /* Process frame if inFIFO is full and codec is ready for it */
        if (pData->FIFO_idx >= POWERMAP_FRAME_SIZE && (pData->codecStatus == CODEC_STATUS_INITIALISED) && isPlaying ) {
            pData->FIFO_idx = 0;
            pData->procStatus = PROC_STATUS_ONGOING;

            /* Load time-domain data */
            for(ch=0; ch<nSH; ch++)
                memcpy(pData->SHframeTD[ch], pData->inFIFO[ch], POWERMAP_FRAME_SIZE*sizeof(float));

            /* account for input channel order */
            switch(chOrdering){
                case CH_ACN:  /* already ACN */ break; /* Otherwise, convert to ACN... */
                case CH_FUMA: convertHOAChannelConvention(FLATTEN2D(pData->SHframeTD), masterOrder, POWERMAP_FRAME_SIZE, HOA_CH_ORDER_FUMA, HOA_CH_ORDER_ACN); break;
            }

            /* account for input normalisation scheme */
            switch(norm){
                case NORM_N3D:  /* already in N3D, do nothing */ break; /* Otherwise, convert to N3D... */
                case NORM_SN3D: convertHOANormConvention(FLATTEN2D(pData->SHframeTD), masterOrder, POWERMAP_FRAME_SIZE, HOA_NORM_SN3D, HOA_NORM_N3D); break;
                case NORM_FUMA: convertHOANormConvention(FLATTEN2D(pData->SHframeTD), masterOrder, POWERMAP_FRAME_SIZE, HOA_NORM_FUMA, HOA_NORM_N3D); break;
            }

            /* apply the time-frequency transform */
            afSTFT_forward_knownDimensions(pData->hSTFT, pData->SHframeTD, POWERMAP_FRAME_SIZE, MAX_NUM_SH_SIGNALS, TIME_SLOTS, pData->SHframeTF);

            /* Update covarience matrix per band */
            for(band=0; band<HYBRID_BANDS; band++){
                cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, nSH, TIME_SLOTS, &calpha,
                            FLATTEN2D(pData->SHframeTF[band]), TIME_SLOTS,
                            FLATTEN2D(pData->SHframeTF[band]), TIME_SLOTS, &cbeta,
                            new_Cx, nSH);

                /* average over time */
                cblas_sscal(nSH*nSH*2, covAvgCoeff, (float*)pData->Cx[band], 1);
                cblas_saxpy(nSH*nSH*2, 1.0f-covAvgCoeff, (float*)new_Cx, 1, (float*)pData->Cx[band], 1);
            }

            /* update the powermap */
            if(pData->recalcPmap==1){
                pData->recalcPmap = 0;
                pData->pmapReady = 0;

                /* determine maximum analysis order */
                maxOrder = 1;
                for(i=0; i<HYBRID_BANDS; i++)
                    maxOrder = SAF_MAX(maxOrder, SAF_MIN(analysisOrderPerBand[i], masterOrder));
                nSH_maxOrder = (maxOrder+1)*(maxOrder+1);

                /* group covarience matrices */
                memset(C_grp, 0, nSH_maxOrder*nSH_maxOrder*sizeof(float_complex));
                for (band=0; band<HYBRID_BANDS; band++){
                    order_band = SAF_MAX(SAF_MIN(pData->analysisOrderPerBand[band], masterOrder),1);
                    nSH_order = (order_band+1)*(order_band+1);
                    pmapEQ_band = SAF_MIN(SAF_MAX(pmapEQ[band], 0.0f), 2.0f);
                    for(i=0; i<nSH_order; i++)
                        for(j=0; j<nSH_order; j++)
                            C_grp[i*nSH_maxOrder+j] = ccaddf(C_grp[i*nSH_maxOrder+j], crmulf(pData->Cx[band][i*nSH+j], 1e3f*pmapEQ_band));
                }

                /* generate powermap */
                C_grp_trace = 0.0f;
                for(i=0; i<nSH_maxOrder; i++)
                    C_grp_trace+=crealf(C_grp[i*nSH_maxOrder+ i]);
                switch(pmap_mode){
                    default:
                    case PM_MODE_PWD:
                        generatePWDmap(maxOrder, (float_complex*)C_grp, pars->Y_grid_cmplx[maxOrder-1], pars->grid_nDirs, pData->pmap);
                        break;

                    case PM_MODE_MVDR:
                        if(C_grp_trace>1e-8f)
                            generateMVDRmap(maxOrder, (float_complex*)C_grp, pars->Y_grid_cmplx[maxOrder-1], pars->grid_nDirs, 8.0f, pData->pmap, NULL);
                        else
                            memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                        break;

                    case PM_MODE_CROPAC_LCMV:
                        if(C_grp_trace>1e-8f)
                            generateCroPaCLCMVmap(maxOrder, (float_complex*)C_grp, pars->Y_grid_cmplx[maxOrder-1], pars->grid_nDirs, 8.0f, 0.0f, pData->pmap);
                        else
                            memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                        break;

                    case PM_MODE_MUSIC:
                        if(C_grp_trace>1e-8f)
                            generateMUSICmap(maxOrder, (float_complex*)C_grp, pars->Y_grid_cmplx[maxOrder-1], nSources, pars->grid_nDirs, 0, pData->pmap);
                        else
                            memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                        break;

                    case PM_MODE_MUSIC_LOG:
                        if(C_grp_trace>1e-8f)
                            generateMUSICmap(maxOrder, (float_complex*)C_grp, pars->Y_grid_cmplx[maxOrder-1], nSources, pars->grid_nDirs, 1, pData->pmap);
                        else
                            memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                        break;

                    case PM_MODE_MINNORM:
                        if(C_grp_trace>1e-8f)
                            generateMinNormMap(maxOrder, (float_complex*)C_grp, pars->Y_grid_cmplx[maxOrder-1], nSources, pars->grid_nDirs, 0, pData->pmap);
                        else
                            memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                        break;

                    case PM_MODE_MINNORM_LOG:
                        if(C_grp_trace>1e-8f)
                            generateMinNormMap(maxOrder, (float_complex*)C_grp, pars->Y_grid_cmplx[maxOrder-1], nSources, pars->grid_nDirs, 1, pData->pmap);
                        else
                            memset(pData->pmap, 0, pars->grid_nDirs*sizeof(float));
                        break;
                }

                /* average powermap over time */
                for(i=0; i<pars->grid_nDirs; i++)
                    pData->pmap[i] =  (1.0f-pmapAvgCoeff) * (pData->pmap[i] )+ pmapAvgCoeff * (pData->prev_pmap[i]);
                utility_svvcopy(pData->pmap, pars->grid_nDirs, pData->prev_pmap);

                /* interpolate powermap */
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->interp_nDirs, 1, pars->grid_nDirs, 1.0f,
                            pars->interp_table, pars->grid_nDirs,
                            pData->pmap, 1, 0.0f,
                            pData->pmap_grid[pData->dispSlotIdx], 1);

                /* ascertain minimum and maximum values for powermap colour scaling */
                int ind;
                utility_siminv(pData->pmap_grid[pData->dispSlotIdx], pars->interp_nDirs, &ind);
                pData->pmap_grid_minVal = pData->pmap_grid[pData->dispSlotIdx][ind];
                utility_simaxv(pData->pmap_grid[pData->dispSlotIdx], pars->interp_nDirs, &ind);
                pData->pmap_grid_maxVal = pData->pmap_grid[pData->dispSlotIdx][ind];

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
        else if(pData->FIFO_idx >= POWERMAP_FRAME_SIZE){
            /* reset FIFO_idx if codec was not ready */
            pData->FIFO_idx = 0;
        }
    }

    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

/* SETS */
 
void powermap_refreshSettings(void* const hPm)
{
    powermap_setCodecStatus(hPm, CODEC_STATUS_NOT_INITIALISED);
}

void powermap_setPowermapMode(void* const hPm, int newMode)
{
    powermap_data *pData = (powermap_data*)(hPm);
    powermap_codecPars* pars = pData->pars;
    pData->pmap_mode = (POWERMAP_MODES)newMode;
    if(pData->prev_pmap!=NULL)
        memset(pData->prev_pmap, 0, pars->grid_nDirs*sizeof(float));
}

void powermap_setMasterOrder(void* const hPm,  int newValue)
{
    powermap_data *pData = (powermap_data*)(hPm);
    if(pData->new_masterOrder != newValue){
        pData->new_masterOrder = newValue;
        powermap_setCodecStatus(hPm, CODEC_STATUS_NOT_INITIALISED);
    }
    /* FUMA only supports 1st order */
    if(pData->new_masterOrder!=SH_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->new_masterOrder!=SH_ORDER_FIRST && pData->norm == NORM_FUMA)
        pData->norm = NORM_SN3D;
}

void powermap_setCovAvgCoeff(void* const hPm, float newAvg)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->covAvgCoeff = SAF_MIN(SAF_MAX(0.0f, newAvg), 0.99999999f);
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
                pData->analysisOrderPerBand[band] = SAF_MIN(pData->new_masterOrder,curOrder);
                if(pData->freqVector[band] > __Zylia_freqRange[(__Zylia_maxOrder-1)*2 - 1])
                    pData->pmapEQ[band] = 0.0f;
            }
            break;

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
                pData->analysisOrderPerBand[band] = SAF_MIN(pData->new_masterOrder,curOrder);
                if(pData->freqVector[band] > __Eigenmike32_freqRange[(__Eigenmike32_maxOrder-1)*2 - 1])
                    pData->pmapEQ[band] = 0.0f;
            }
            break;

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
                pData->analysisOrderPerBand[band] = SAF_MIN(pData->new_masterOrder,curOrder);
                if(pData->freqVector[band] > __DTU_mic_freqRange[(__DTU_mic_maxOrder-1)*2 - 1])
                    pData->pmapEQ[band] = 0.0f;
            }
            break;
    }
}

void powermap_setAnaOrder(void  * const hPm, int newValue, int bandIdx)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->analysisOrderPerBand[bandIdx] = SAF_MIN(SAF_MAX(newValue,1), pData->new_masterOrder);
}

void powermap_setAnaOrderAllBands(void  * const hPm, int newValue)
{
    powermap_data *pData = (powermap_data*)(hPm);
    int band;

    for(band=0; band<HYBRID_BANDS; band++)
        pData->analysisOrderPerBand[band] = SAF_MIN(SAF_MAX(newValue,1), pData->new_masterOrder);
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
    if((CH_ORDER)newOrder != CH_FUMA || pData->new_masterOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void powermap_setNormType(void* const hPm, int newType)
{
    powermap_data *pData = (powermap_data*)(hPm);
    if((NORM_TYPES)newType != NORM_FUMA || pData->new_masterOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
}

void powermap_setDispFOV(void* const hPm, int newOption)
{
    powermap_data *pData = (powermap_data*)(hPm);
    if(pData->HFOVoption != (HFOV_OPTIONS)newOption){
        pData->HFOVoption = (HFOV_OPTIONS)newOption;
        powermap_setCodecStatus(hPm, CODEC_STATUS_NOT_INITIALISED);
    }
}

void powermap_setAspectRatio(void* const hPm, int newOption)
{
    powermap_data *pData = (powermap_data*)(hPm);
    if(pData->aspectRatioOption != (ASPECT_RATIO_OPTIONS)newOption){
        pData->aspectRatioOption = (ASPECT_RATIO_OPTIONS)newOption;
        powermap_setCodecStatus(hPm, CODEC_STATUS_NOT_INITIALISED);
    }
}

void powermap_setPowermapAvgCoeff(void* const hPm, float newValue)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->pmapAvgCoeff = SAF_MIN(SAF_MAX(0.0f, newValue), 0.99999999f);
}

void powermap_requestPmapUpdate(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    pData->recalcPmap = 1;
}

/* GETS */

int powermap_getFrameSize(void)
{
    return POWERMAP_FRAME_SIZE;
}

CODEC_STATUS powermap_getCodecStatus(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return pData->codecStatus;
}

float powermap_getProgressBar0_1(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    return pData->progressBar0_1;
}

void powermap_getProgressBarText(void* const hPm, char* text)
{
    powermap_data *pData = (powermap_data*)(hPm);
    memcpy(text, pData->progressBarText, PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

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
    return (pData->new_masterOrder+1)*(pData->new_masterOrder+1);
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
    powermap_codecPars* pars = pData->pars;
    if((pData->codecStatus == CODEC_STATUS_INITIALISED) && pData->pmapReady){
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

int powermap_getProcessingDelay()
{
    return POWERMAP_FRAME_SIZE + 12*HOP_SIZE;
}

