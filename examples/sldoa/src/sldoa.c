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
 *     sldoa.c
 * Description:
 *     A spatially-localised pressure-intensity based direction-of-arrival estimator (SLDoA).
 *     VBAP gain patterns are imposed on the spherical harmonic signals, such that the DoA
 *     can be estimated in a spatially-constrained region; thus mitigating interferes and
 *     reflections arriving from other directions. The DoA is estimated per sector for
 *     each frequency band.
 *     The algorithms within sldoa were developed in collaboration with Symeon Delikaris-
 *     Manias, and are explained in more detail in:
 *     McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki, V.,
 *     “Real-time conversion of sensor array signals into spherical harmonic signals with
 *     applications to spatially localised sub-band sound-field analysis,” in Audio
 *     Engineering Society Convention 144, Audio Engineering Society, 2018.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 18.10.2017
 */

#include "sldoa.h"
#include "sldoa_internal.h"
#include "sldoa_database.h" 

void sldoa_create
(
    void ** const phSld
)
{
    sldoa_data* pData = (sldoa_data*)malloc(sizeof(sldoa_data));
    if (pData == NULL) { return;/*error*/ }
    *phSld = (void*)pData;
    int i, j, t, ch, band;
    
    afSTFTinit(&(pData->hSTFT), HOP_SIZE, NUM_SH_SIGNALS, 0, 0, 1);
    pData->STFTInputFrameTF = (complexVector**)malloc2d(TIME_SLOTS, NUM_SH_SIGNALS, sizeof(complexVector));
    for(t=0; t<TIME_SLOTS; t++) {
        for(ch=0; ch< NUM_SH_SIGNALS; ch++) {
            pData->STFTInputFrameTF[t][ch].re = (float*)calloc(HYBRID_BANDS, sizeof(float));
            pData->STFTInputFrameTF[t][ch].im = (float*)calloc(HYBRID_BANDS, sizeof(float));
        }
    }
    pData->tempHopFrameTD = (float**)malloc2d(NUM_SH_SIGNALS, HOP_SIZE, sizeof(float));
    
    /* internal */
    pData->reInitAna = 1;
#if SH_ORDER > 1
    for(i=0; i<SH_ORDER; i++)
        pData->secCoeffs[i] = NULL;
#endif
    for(i=0; i<64; i++)
        for(j=0; j<NUM_GRID_DIRS; j++)
            pData->grid_Y[i][j] = (float)__grid_Y[i][j];
    for(i=0; i<NUM_GRID_DIRS; i++)
        for(j=0; j<2; j++)
            pData->grid_dirs_deg[i][j] = (float)__grid_dirs_deg[i][j];
    
    /* display */
    for(i=0; i<NUM_DISP_SLOTS; i++){
        pData->azi_deg[i] = malloc(HYBRID_BANDS*MAX_NUM_SECTORS * sizeof(float));
        pData->elev_deg[i] = malloc(HYBRID_BANDS*MAX_NUM_SECTORS * sizeof(float));
        pData->colourScale[i] = malloc(HYBRID_BANDS*MAX_NUM_SECTORS * sizeof(float));
        pData->alphaScale[i] = malloc(HYBRID_BANDS*MAX_NUM_SECTORS * sizeof(float));
    }
    
    /* Default user parameters */
    for(band=0; band<HYBRID_BANDS; band++){
        pData->analysisOrderPerBand[band] = SH_ORDER;
        pData->nSectorsPerBand[band] = ORDER2NUMSECTORS(pData->analysisOrderPerBand[band]);
    }
    pData->minFreq = 500.0f;
    pData->maxFreq = 5e3f;
    pData->avg_ms = 500.0f;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_N3D;
}

void sldoa_destroy
(
    void ** const phSld
)
{
    sldoa_data *pData = (sldoa_data*)(*phSld);
    int i, t, ch;

    if (pData != NULL) {
        afSTFTfree(pData->hSTFT);
        for (t = 0; t<TIME_SLOTS; t++) {
            for (ch = 0; ch< NUM_SH_SIGNALS; ch++) {
                free(pData->STFTInputFrameTF[t][ch].re);
                free(pData->STFTInputFrameTF[t][ch].im);
            }
        }
        free2d((void**)pData->STFTInputFrameTF, TIME_SLOTS);
        free2d((void**)pData->tempHopFrameTD, NUM_SH_SIGNALS);
        for(i=0; i<NUM_DISP_SLOTS; i++){
            free(pData->azi_deg[i]);
            free(pData->elev_deg[i]);
            free(pData->colourScale[i]);
            free(pData->alphaScale[i]);
        }
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
    int i, band;
    
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


void sldoa_analysis
(
    void  *  const hSld,
    float ** const inputs,
    int            nInputs,
    int            nSamples,
    int            isPlaying 
)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    int i, j, t, n, ch, sample, band, nSectors, min_band, numAnalysisBands, current_disp_idx;
    float avgCoeff, max_en[HYBRID_BANDS], min_en[HYBRID_BANDS];
    float new_doa[MAX_NUM_SECTORS][TIME_SLOTS][2], new_doa_xyz[3], doa_xyz[3], avg_xyz[3];
    float new_energy[MAX_NUM_SECTORS][TIME_SLOTS];
    int o[SH_ORDER+2];
    
    /* local parameters */
    int analysisOrderPerBand[HYBRID_BANDS];
    int nSectorsPerBand[HYBRID_BANDS];
    float minFreq, maxFreq, avg_ms;
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    
    /* reinitialise if needed */
    if(pData->reInitAna == 1){
        pData->reInitAna = 2; /* indicate init in progress */
        sldoa_initAna(hSld);
        pData->reInitAna = 0; /* indicate init complete */
    }
    if (nSamples == FRAME_SIZE && (pData->reInitAna == 0) && isPlaying) {
        current_disp_idx = pData->current_disp_idx;
        /* copy current parameters to be thread safe */
        memcpy(analysisOrderPerBand, pData->analysisOrderPerBand, HYBRID_BANDS*sizeof(int));
        memcpy(nSectorsPerBand, pData->nSectorsPerBand, HYBRID_BANDS*sizeof(int));
        minFreq = pData->minFreq;
        maxFreq = pData->maxFreq;
        avg_ms = pData->avg_ms;
        chOrdering = pData->chOrdering;
        norm = pData->norm;
        
        /* load intput time-domain data */
        for (i = 0; i < MIN(NUM_SH_SIGNALS, nInputs); i++)
            memcpy(pData->SHframeTD[i], inputs[i], FRAME_SIZE*sizeof(float));
        for (; i < NUM_SH_SIGNALS; i++)
            memset(pData->SHframeTD[i], 0, FRAME_SIZE*sizeof(float));
        
        /* account for input normalisation scheme */
        switch(norm){
            case NORM_N3D:  /* already in N3D, do nothing */
                break;
            case NORM_SN3D: /* convert to N3D */
                for(n=0; n<SH_ORDER+2; n++){  o[n] = n*n;  };
                for (n = 0; n<SH_ORDER+1; n++)
                    for (ch = o[n]; ch<o[n+1]; ch++)
                        for(i = 0; i<FRAME_SIZE; i++)
                            pData->SHframeTD[ch][i] *= sqrtf(2.0f*(float)n+1.0f);
                break;
        }
        
        /* apply the time-frequency transform */
        for (t = 0; t < TIME_SLOTS; t++) {
            for (ch = 0; ch < NUM_SH_SIGNALS; ch++)
                for (sample = 0; sample < HOP_SIZE; sample++)
                    pData->tempHopFrameTD[ch][sample] = pData->SHframeTD[ch][sample + t*HOP_SIZE];
            afSTFTforward(pData->hSTFT, (float**)pData->tempHopFrameTD, (complexVector*)pData->STFTInputFrameTF[t]);
        }
        for (band = 0; band < HYBRID_BANDS; band++)
            for (ch = 0; ch < NUM_SH_SIGNALS; ch++)
                for (t = 0; t < TIME_SLOTS; t++)
                    pData->SHframeTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[t][ch].re[band], pData->STFTInputFrameTF[t][ch].im[band]);
        
        /* apply sector-based, frequency-dependent DOA analysis */
        numAnalysisBands = 0;
        min_band = 0;
        for(band=1/* ignore DC */; band<HYBRID_BANDS; band++){
            if(pData->freqVector[band] <= minFreq)
                min_band = band;
            if(pData->freqVector[band] >= minFreq && pData->freqVector[band]<=maxFreq){
                nSectors = nSectorsPerBand[band];
                avgCoeff = avg_ms < 10.0f ? 1.0f : 1.0f / ((avg_ms/1e3f) / (1.0f/(float)HOP_SIZE) + 2.23e-9f);
                avgCoeff = MAX(MIN(avgCoeff, 0.99999f), 0.0f); /* ensures stability */
                sldoa_estimateDoA(pData->SHframeTF[band],
                                  analysisOrderPerBand[band],
#if SH_ORDER > 1
                                  pData->secCoeffs[analysisOrderPerBand[band]-2],
#else
                                  NULL,
#endif
                                  new_doa,
                                  new_energy);
                
                /* average the raw data over time */
                for(i=0; i<nSectors; i++){
                    for( t = 0; t<TIME_SLOTS; t++){
                        /* avg doa estimate */
                        unitSph2Cart(new_doa[i][t][0], new_doa[i][t][1], new_doa_xyz);
                        unitSph2Cart(pData->doa_rad[band][i][0],
                                     pData->doa_rad[band][i][1],
                                     doa_xyz);
                        for(j=0; j<3; j++)
                            avg_xyz[j] = new_doa_xyz[j]*avgCoeff + doa_xyz[j] * (1.0f-avgCoeff);
                        unitCart2Sph_aziElev(avg_xyz, &(pData->doa_rad[band][i][0]), &(pData->doa_rad[band][i][1]));
                        
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
                    pData->azi_deg [current_disp_idx][band*MAX_NUM_SECTORS + i] = pData->doa_rad[band][i][0]*180.0f/M_PI;
                    pData->elev_deg[current_disp_idx][band*MAX_NUM_SECTORS + i] = pData->doa_rad[band][i][1]*180.0f/M_PI;
                    
                    /* colour should indicate the different frequencies */
                    pData->colourScale[current_disp_idx][band*MAX_NUM_SECTORS + i] = (float)(band-min_band)/(float)(numAnalysisBands+1);
                    
                    /* transparancy should indicate the energy of the sector for each DoA estimate, for each frequency */
                    if( analysisOrderPerBand[band]==1  )
                        pData->alphaScale[current_disp_idx][band*MAX_NUM_SECTORS + i] = 1.0f;
                    else
                        pData->alphaScale[current_disp_idx][band*MAX_NUM_SECTORS + i] = MIN(MAX((pData->energy[band][i]-min_en[band])/(max_en[band]-min_en[band]+2.3e-10f), 0.11f),1.0f);
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
}

/* SETS */

void sldoa_refreshSettings(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    pData->reInitAna = 1;
}

void sldoa_setMaxFreq(void* const hSld, float newFreq)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    newFreq = MAX(MIN(newFreq, pData->fs/2.0f),0.0f);
    if(newFreq < pData->minFreq )
        pData->minFreq = newFreq;
    pData->maxFreq = newFreq;
}

void sldoa_setMinFreq(void* const hSld, float newFreq)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    newFreq = MAX(MIN(newFreq, pData->fs/2.0f),0.0f);
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
                pData->analysisOrderPerBand[band] = SH_ORDER;
            break;
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
                pData->analysisOrderPerBand[band] = MIN(SH_ORDER,curOrder);
            }
            pData->maxFreq = __Zylia_freqRange[(__Zylia_maxOrder-1)*2-1];
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
                pData->analysisOrderPerBand[band] = MIN(SH_ORDER,curOrder);
            }
            pData->maxFreq = __Eigenmike32_freqRange[(__Eigenmike32_maxOrder-1)*2-1];
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
                pData->analysisOrderPerBand[band] = MIN(SH_ORDER,curOrder);
            }
            pData->maxFreq = __DTU_mic_freqRange[(__DTU_mic_maxOrder-1)*2-1];
            break;
#endif
    }
    for(band=0; band<HYBRID_BANDS; band++)
        pData->nSectorsPerBand[band] = ORDER2NUMSECTORS(pData->analysisOrderPerBand[band]);
}

void sldoa_setAnaOrder(void * const hSld, int newValue, int bandIdx)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    pData->analysisOrderPerBand[bandIdx] = MIN(MAX(newValue,1), SH_ORDER);
    pData->nSectorsPerBand[bandIdx] = ORDER2NUMSECTORS(pData->analysisOrderPerBand[bandIdx]);
}

void sldoa_setAnaOrderAllBands(void * const hSld, int newValue)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    int band;
    
    for(band=0; band<HYBRID_BANDS; band++){
        pData->analysisOrderPerBand[band] = MIN(MAX(newValue,1), SH_ORDER);
        pData->nSectorsPerBand[band] = ORDER2NUMSECTORS(pData->analysisOrderPerBand[band]);
    }
}

void sldoa_setChOrder(void* const hSld, int newOrder)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    pData->chOrdering = (CH_ORDER)newOrder;
}

void sldoa_setNormType(void* const hSld, int newType)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    pData->norm = (NORM_TYPES)newType;
}

/* GETS */

float sldoa_getSamplingRate(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    return pData->fs;
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

int sldoa_getNSHrequired(void)
{
    return NUM_SH_SIGNALS;
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

