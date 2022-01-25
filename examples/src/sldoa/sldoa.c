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
 * @file sldoa.c
 * @brief A spatially-localised active-intensity (SLAI) based direction-of-
 *        arrival estimator (SLDoA)
 *
 * VBAP gain patterns are imposed on the spherical harmonic signals, such that
 * the DoA can be estimated in a spatially-constrained region; thus mitigating
 * the effect of interferes and reflections arriving from other directions.
 * The DoA is estimated per sector for each frequency band.
 *
 * The algorithms within sldoa were developed in collaboration with Symeon
 * Delikaris-Manias and Angelo Farina, and are explained in more detail in [1,2]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 * @see [2] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., "Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis," in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 *
 * @author Leo McCormack
 * @date 18.10.2017
 * @license ISC
 */

#include "sldoa.h"
#include "sldoa_internal.h"
#include "sldoa_database.h" 

void sldoa_create
(
    void ** const phSld
)
{
    sldoa_data* pData = (sldoa_data*)malloc1d(sizeof(sldoa_data));
    *phSld = (void*)pData;
    int i, j, band;

    /* Default user parameters */
    pData->new_masterOrder = pData->masterOrder = 1;
    for(band=0; band<HYBRID_BANDS; band++){
        pData->analysisOrderPerBand[band] = pData->masterOrder;
        pData->nSectorsPerBand[band] = ORDER2NUMSECTORS(pData->analysisOrderPerBand[band]);
    }
    pData->minFreq = 500.0f;
    pData->maxFreq = 5e3f;
    pData->avg_ms = 500.0f;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;

    /* TFT */
    afSTFT_create(&(pData->hSTFT), MAX_NUM_SH_SIGNALS, 0, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    pData->SHframeTD = (float**)malloc2d(MAX_NUM_SH_SIGNALS, SLDOA_FRAME_SIZE, sizeof(float));
    pData->SHframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_SH_SIGNALS, TIME_SLOTS, sizeof(float_complex));

    /* internal */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
    for(i=0; i<MAX_SH_ORDER-1; i++)
        pData->secCoeffs[i] = NULL;
    for(i=0; i<64; i++)
        for(j=0; j<NUM_GRID_DIRS; j++)
            pData->grid_Y[i][j] = (float)__grid_Y[i][j] * sqrtf(4.0f*SAF_PI);
    for(i=0; i<3; i++)
        for(j=0; j<NUM_GRID_DIRS; j++)
            pData->grid_Y_dipoles_norm[i][j] = pData->grid_Y[i+1][j]/sqrtf(3); /* scale to [0..1] */
    for(i=0; i<NUM_GRID_DIRS; i++)
        for(j=0; j<2; j++)
            pData->grid_dirs_deg[i][j] = (float)__grid_dirs_deg[i][j];
    
    /* display */
    for(i=0; i<NUM_DISP_SLOTS; i++){
        pData->azi_deg[i] = malloc1d(HYBRID_BANDS*MAX_NUM_SECTORS * sizeof(float));
        pData->elev_deg[i] = malloc1d(HYBRID_BANDS*MAX_NUM_SECTORS * sizeof(float));
        pData->colourScale[i] = malloc1d(HYBRID_BANDS*MAX_NUM_SECTORS * sizeof(float));
        pData->alphaScale[i] = malloc1d(HYBRID_BANDS*MAX_NUM_SECTORS * sizeof(float));
    }

    /* set FIFO buffer */
    pData->FIFO_idx = 0;
    memset(pData->inFIFO, 0, MAX_NUM_SH_SIGNALS*SLDOA_FRAME_SIZE*sizeof(float));
}

void sldoa_destroy
(
    void ** const phSld
)
{
    sldoa_data *pData = (sldoa_data*)(*phSld);
    int i;

    if (pData != NULL) {
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            SAF_SLEEP(10);
        }
        
        /* free afSTFT and buffers */
        afSTFT_destroy(&(pData->hSTFT));
        free(pData->SHframeTD);
        free(pData->SHframeTF);

        for(i=0; i<NUM_DISP_SLOTS; i++){
            free(pData->azi_deg[i]);
            free(pData->elev_deg[i]);
            free(pData->colourScale[i]);
            free(pData->alphaScale[i]);
        }
        free(pData->progressBarText);
        free(pData);
        pData = NULL;
    }
}

void sldoa_init
(
    void * const hSld,
    float        sampleRate
)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    int i;
    
    pData->fs = sampleRate;
    
    /* specify frequency vector and determine the number of bands */
    afSTFT_getCentreFreqs(pData->hSTFT, sampleRate, HYBRID_BANDS, pData->freqVector);

    /* intialise display parameters */
    pData->current_disp_idx = 0;
    memset(pData->doa_rad, 0, HYBRID_BANDS*MAX_NUM_SECTORS*2* sizeof(float));
    memset(pData->energy, 0, HYBRID_BANDS*MAX_NUM_SECTORS* sizeof(float));
    for(i=0; i<NUM_DISP_SLOTS; i++){
        memset(pData->azi_deg[i], 0, HYBRID_BANDS*MAX_NUM_SECTORS* sizeof(float));
        memset(pData->elev_deg[i], 0, HYBRID_BANDS*MAX_NUM_SECTORS * sizeof(float));
        memset(pData->colourScale[i], 0, HYBRID_BANDS*MAX_NUM_SECTORS * sizeof(float));
        memset(pData->alphaScale[i], 0, HYBRID_BANDS*MAX_NUM_SECTORS * sizeof(float));
    }
}

void sldoa_initCodec
(
    void* const hSld
)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    
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
    
    sldoa_initTFT(hSld);
    sldoa_initAna(hSld);
    
    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
}

void sldoa_analysis
(
    void        *  const hSld,
    const float *const * inputs,
    int                  nInputs,
    int                  nSamples,
    int                  isPlaying 
)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    int s, i, j, t, ch, band, nSectors, min_band, numAnalysisBands, current_disp_idx;
    float avgCoeff, max_en[HYBRID_BANDS], min_en[HYBRID_BANDS];
    float new_doa[MAX_NUM_SECTORS][TIME_SLOTS][2], new_doa_xyz[3], doa_xyz[3], avg_xyz[3];
    float new_energy[MAX_NUM_SECTORS][TIME_SLOTS];
    
    /* local parameters */
    int nSH, masterOrder;
    int analysisOrderPerBand[HYBRID_BANDS];
    int nSectorsPerBand[HYBRID_BANDS];
    float minFreq, maxFreq, avg_ms;
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    memcpy(analysisOrderPerBand, pData->analysisOrderPerBand, HYBRID_BANDS*sizeof(int));
    memcpy(nSectorsPerBand, pData->nSectorsPerBand, HYBRID_BANDS*sizeof(int));
    minFreq = pData->minFreq;
    maxFreq = pData->maxFreq;
    avg_ms = pData->avg_ms;
    chOrdering = pData->chOrdering;
    norm = pData->norm;
    masterOrder = pData->masterOrder;
    nSH = ORDER2NSH(masterOrder);

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
        if (pData->FIFO_idx >= SLDOA_FRAME_SIZE && (pData->codecStatus == CODEC_STATUS_INITIALISED) && isPlaying) {
            pData->FIFO_idx = 0;
            pData->procStatus = PROC_STATUS_ONGOING;
            current_disp_idx = pData->current_disp_idx;

            /* Load time-domain data */
            for(ch=0; ch<nSH; ch++)
                memcpy(pData->SHframeTD[ch],pData->inFIFO[ch], SLDOA_FRAME_SIZE*sizeof(float));

            /* account for input channel order */
            switch(chOrdering){
                case CH_ACN:  /* already ACN */ break; /* Otherwise, convert to ACN... */
                case CH_FUMA: convertHOAChannelConvention(FLATTEN2D(pData->SHframeTD), masterOrder, SLDOA_FRAME_SIZE, HOA_CH_ORDER_FUMA, HOA_CH_ORDER_ACN); break;
            }

            /* account for input normalisation scheme */
            switch(norm){
                case NORM_N3D:  /* already in N3D, do nothing */ break; /* Otherwise, convert to N3D... */
                case NORM_SN3D: convertHOANormConvention(FLATTEN2D(pData->SHframeTD), masterOrder, SLDOA_FRAME_SIZE, HOA_NORM_SN3D, HOA_NORM_N3D); break;
                case NORM_FUMA: convertHOANormConvention(FLATTEN2D(pData->SHframeTD), masterOrder, SLDOA_FRAME_SIZE, HOA_NORM_FUMA, HOA_NORM_N3D); break;
            }
        
            /* apply the time-frequency transform */
            afSTFT_forward_knownDimensions(pData->hSTFT, pData->SHframeTD, SLDOA_FRAME_SIZE, MAX_NUM_SH_SIGNALS, TIME_SLOTS, pData->SHframeTF);

            /* apply sector-based, frequency-dependent DOA analysis */
            numAnalysisBands = 0;
            min_band = 0;
            for(band=1/* ignore DC */; band<HYBRID_BANDS; band++){
                if(pData->freqVector[band] <= minFreq)
                    min_band = band;
                if(pData->freqVector[band] >= minFreq && pData->freqVector[band]<=maxFreq){
                    nSectors = nSectorsPerBand[band];
                    avgCoeff = avg_ms < 10.0f ? 1.0f : 1.0f / ((avg_ms/1e3f) / (1.0f/(float)HOP_SIZE) + 2.23e-9f);
                    avgCoeff = SAF_MAX(SAF_MIN(avgCoeff, 0.99999f), 0.0f); /* ensures stability */
                    sldoa_estimateDoA(pData->SHframeTF[band],
                                      analysisOrderPerBand[band],
                                      pData->secCoeffs[analysisOrderPerBand[band]-2], /* -2, as first order is skipped */
                                      new_doa,
                                      new_energy);

                    /* average the raw data over time */
                    for(i=0; i<nSectors; i++){
                        for( t = 0; t<TIME_SLOTS; t++){
                            /* avg doa estimate */
                            unitSph2cart((float*)new_doa[i][t], 1, 0, (float*)new_doa_xyz);
                            unitSph2cart((float*)pData->doa_rad[band][i], 1, 0, (float*)doa_xyz);
                            for(j=0; j<3; j++)
                                avg_xyz[j] = new_doa_xyz[j]*avgCoeff + doa_xyz[j] * (1.0f-avgCoeff); 
                            unitCart2sph((float*)avg_xyz, 1, 0, (float*)pData->doa_rad[band][i]);

                            /* avg energy */
                            pData->energy[band][i] = new_energy[i][t]*avgCoeff + pData->energy[band][i] * (1.0f-avgCoeff);
                        }
                    }
                    numAnalysisBands++;
                }
            }

            /* determine the minimum and maximum sector energies per frequency (to scale them 0..1) */
            for(band=1/* ignore DC */; band<HYBRID_BANDS; band++){
                if(pData->freqVector[band] >= minFreq && pData->freqVector[band]<=maxFreq){
                    nSectors = nSectorsPerBand[band];
                    max_en[band] = 2.3e-13f; min_en[band] = 2.3e13f; /* starting values */
                    for(i=0; i<nSectors; i++){
                        max_en[band] = pData->energy[band][i] > max_en[band] ? pData->energy[band][i] : max_en[band];
                        min_en[band] = pData->energy[band][i] < min_en[band] ? pData->energy[band][i] : min_en[band];
                    }
                }
            }

            /* prep data for plotting */
            for(band=1/* ignore DC */; band<HYBRID_BANDS; band++){
                if(pData->freqVector[band] >= minFreq && pData->freqVector[band]<=maxFreq){
                    nSectors = nSectorsPerBand[band];
                    /* store averaged values */
                    for(i=0; i<nSectors; i++){
                        pData->azi_deg [current_disp_idx][band*MAX_NUM_SECTORS + i] = pData->doa_rad[band][i][0]*180.0f/SAF_PI;
                        pData->elev_deg[current_disp_idx][band*MAX_NUM_SECTORS + i] = pData->doa_rad[band][i][1]*180.0f/SAF_PI;

                        /* colour should indicate the different frequencies */
                        pData->colourScale[current_disp_idx][band*MAX_NUM_SECTORS + i] = (float)(band-min_band)/(float)(numAnalysisBands+1);

                        /* transparancy should indicate the energy of the sector for each DoA estimate, for each frequency */
                        if( analysisOrderPerBand[band]==1  )
                            pData->alphaScale[current_disp_idx][band*MAX_NUM_SECTORS + i] = 1.0f;
                        else
                            pData->alphaScale[current_disp_idx][band*MAX_NUM_SECTORS + i] = SAF_MIN(SAF_MAX((pData->energy[band][i]-min_en[band])/(max_en[band]-min_en[band]+2.3e-10f), 0.05f),1.0f);
                    }
                }
                else{
                    memset(&(pData->azi_deg [current_disp_idx][band*MAX_NUM_SECTORS]), 0, MAX_NUM_SECTORS*sizeof(float));
                    memset(&(pData->elev_deg [current_disp_idx][band*MAX_NUM_SECTORS]), 0, MAX_NUM_SECTORS*sizeof(float));
                    memset(&(pData->colourScale [current_disp_idx][band*MAX_NUM_SECTORS]), 0, MAX_NUM_SECTORS*sizeof(float));
                    memset(&(pData->alphaScale [current_disp_idx][band*MAX_NUM_SECTORS]), 0, MAX_NUM_SECTORS*sizeof(float));
                }
            }
        }
        else if(pData->FIFO_idx >= SLDOA_FRAME_SIZE){
            /* reset FIFO_idx if codec was not ready */
            pData->FIFO_idx = 0;
        }
    }

    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

/* SETS */

void sldoa_setMasterOrder(void* const hSld,  int newValue)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    if(pData->new_masterOrder != newValue){
        pData->new_masterOrder = newValue;
        sldoa_setCodecStatus(hSld, CODEC_STATUS_NOT_INITIALISED);
    }
    /* FUMA only supports 1st order */
    if(pData->new_masterOrder!=SH_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->new_masterOrder!=SH_ORDER_FIRST && pData->norm == NORM_FUMA)
        pData->norm = NORM_SN3D;
}

void sldoa_refreshSettings(void* const hSld)
{
    sldoa_setCodecStatus(hSld, CODEC_STATUS_NOT_INITIALISED);
}

void sldoa_setMaxFreq(void* const hSld, float newFreq)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    newFreq = SAF_MAX(SAF_MIN(newFreq, pData->fs/2.0f),0.0f);
    if(newFreq < pData->minFreq )
        pData->minFreq = newFreq;
    pData->maxFreq = newFreq;
}

void sldoa_setMinFreq(void* const hSld, float newFreq)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    newFreq = SAF_MAX(SAF_MIN(newFreq, pData->fs/2.0f),0.0f);
    if(newFreq > pData->maxFreq )
        pData->maxFreq = newFreq;
    pData->minFreq = newFreq;
}

void sldoa_setAvg(void* const hSld, float newAvg)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    pData->avg_ms = newAvg;
}

void sldoa_setSourcePreset(void* const hSld, int newPresetID)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    int band, rangeIdx, curOrder, reverse;
    
    rangeIdx = 0;
    curOrder = 1;
    reverse = 0;
    switch(newPresetID){
        case MIC_PRESET_IDEAL:
            for(band=0; band<HYBRID_BANDS; band++)
                pData->analysisOrderPerBand[band] = pData->new_masterOrder;
            break;

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
            }
            pData->maxFreq = __Zylia_freqRange[(__Zylia_maxOrder-1)*2-1];
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
            }
            pData->maxFreq = __Eigenmike32_freqRange[(__Eigenmike32_maxOrder-1)*2-1];
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
            }
            pData->maxFreq = __DTU_mic_freqRange[(__DTU_mic_maxOrder-1)*2-1];
            break;
    }
    for(band=0; band<HYBRID_BANDS; band++)
        pData->nSectorsPerBand[band] = ORDER2NUMSECTORS(pData->analysisOrderPerBand[band]);
}

void sldoa_setAnaOrder(void * const hSld, int newValue, int bandIdx)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    pData->analysisOrderPerBand[bandIdx] = SAF_MIN(SAF_MAX(newValue,1), pData->new_masterOrder);
    pData->nSectorsPerBand[bandIdx] = ORDER2NUMSECTORS(pData->analysisOrderPerBand[bandIdx]);
}

void sldoa_setAnaOrderAllBands(void * const hSld, int newValue)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    int band;
    
    for(band=0; band<HYBRID_BANDS; band++){
        pData->analysisOrderPerBand[band] = SAF_MIN(SAF_MAX(newValue,1), pData->new_masterOrder);
        pData->nSectorsPerBand[band] = ORDER2NUMSECTORS(pData->analysisOrderPerBand[band]);
    }
}

void sldoa_setChOrder(void* const hSld, int newOrder)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    if((CH_ORDER)newOrder != CH_FUMA || pData->new_masterOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void sldoa_setNormType(void* const hSld, int newType)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    if((NORM_TYPES)newType != NORM_FUMA || pData->new_masterOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
}


/* GETS */

int sldoa_getFrameSize(void)
{
    return SLDOA_FRAME_SIZE;
}

CODEC_STATUS sldoa_getCodecStatus(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return pData->codecStatus;
}

float sldoa_getProgressBar0_1(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return pData->progressBar0_1;
}

void sldoa_getProgressBarText(void* const hSld, char* text)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    memcpy(text, pData->progressBarText, PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

int sldoa_getMasterOrder(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return pData->new_masterOrder;
}

int sldoa_getSamplingRate(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return (int)(pData->fs+0.5f);
}

float sldoa_getMaxFreq(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return pData->maxFreq;
}

float sldoa_getMinFreq(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return pData->minFreq;
}

float sldoa_getAvg(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return pData->avg_ms;
}

/* Not very elegent, but does the job */
void sldoa_getDisplayData
(
    void  *  const hSld,
    float** azi_deg,
    float** elev_deg,
    float** colourScale,
    float** alphaScale,
    int** pNsectorsPerBand,
    int* maxNumSectors,
    int* startBand,
    int* endBand
)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    int i;
    
    (*azi_deg) = pData->azi_deg[pData->current_disp_idx];
    (*elev_deg) = pData->elev_deg[pData->current_disp_idx];
    (*colourScale) = pData->colourScale[pData->current_disp_idx];
    (*alphaScale) = pData->alphaScale[pData->current_disp_idx];
    (*pNsectorsPerBand) = pData->nSectorsPerBand;
    (*maxNumSectors) = MAX_NUM_SECTORS;
    (*startBand) =1;
    (*endBand) =1;
    for(i=1/*ignore DC*/; i<HYBRID_BANDS; i++){
        if(pData->freqVector[i]<pData->minFreq)
            (*startBand) = i+1;
        if(pData->freqVector[i]<pData->maxFreq)
            (*endBand) = i;
    }
    /* read the next buffer for the next call */
    pData->current_disp_idx++;
    if(pData->current_disp_idx >= NUM_DISP_SLOTS)
        pData->current_disp_idx = 0;
}

int sldoa_getAnaOrder(void  * const hSld, int bandIdx)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return pData->analysisOrderPerBand[bandIdx];
}

int sldoa_getAnaOrderAllBands(void  * const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return pData->analysisOrderPerBand[0];
}

void sldoa_getAnaOrderHandle
(
    void* const hSld,
    float** pX_vector,
    int** pY_values,
    int* pNpoints
)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    (*pX_vector) = &(pData->freqVector[0]);
    (*pY_values) = &(pData->analysisOrderPerBand[0]);
    (*pNpoints) = HYBRID_BANDS;
}

int sldoa_getNumberOfBands(void)
{
    return HYBRID_BANDS;
}

int sldoa_getNSHrequired(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return (pData->new_masterOrder+1)*(pData->new_masterOrder+1);
}

int sldoa_getChOrder(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return (int)pData->chOrdering;
}

int sldoa_getNormType(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return (int)pData->norm;
}

int sldoa_getProcessingDelay()
{
    return SLDOA_FRAME_SIZE + 12*HOP_SIZE;
}
