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
 *     dirass.c
 * Description:
 *     A sound-field visualiser based on the directional re-assignment of beamformer energy,
 *     utilising the DoA estimates extracted from spatially-localised active-intensity
 *     (SLAI) vectors; which correspond to the scanning grid directions.
 *     For more information on the method, refer to:
 *         McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of angular
 *         spectra based on a directional re-assignment approach for ambisonic sound-field
 *         visualisation". IEEE International Conference on Acoustics, Speech and Signal
 *         Processing (ICASSP).
 *
 * Dependencies:
 *     saf_utilities, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 21.02.2019
 */

#include "dirass.h"
#include "dirass_internal.h"

void dirass_create
(
    void ** const phDir
)
{
    dirass_data* pData = (dirass_data*)malloc(sizeof(dirass_data));
    if (pData == NULL) { return;/*error*/ }
    *phDir = (void*)pData;
    int i;
    
    /* codec data */
    pData->pars = (codecPars*)malloc(sizeof(codecPars));
    codecPars* pars = pData->pars;
    pars->interp_dirs_deg = NULL;
    pars->interp_dirs_rad = NULL;
    pars->Y_up = NULL;
    pars->interp_table = NULL;
    pars->w = NULL;
    pars->Cw = NULL;
    pars->Uw = NULL;
    pars->Cxyz = NULL;
    pars->ss = NULL;
    pars->ssxyz = NULL;
    pars->est_dirs = NULL;
    pars->est_dirs_idx = NULL;
    pars->prev_intensity = NULL;
    
    /* internal */
	pData->reInitAna = 1; 

    /* display */
    pData->pmap = NULL;
    pData->prev_pmap = NULL;
    for(i=0; i<NUM_DISP_SLOTS; i++)
        pData->pmap_grid[i] = NULL;
    pData->pmapReady = 0;
    pData->recalcPmap = 1;
    
    /* Default user parameters */
    pData->inputOrder = pData->new_inputOrder = INPUT_ORDER_FIRST;
    pData->beamType = BEAM_TYPE_HYPERCARD;
    pData->DirAssMode = REASS_UPSCALE;
    pData->upscaleOrder = pData->new_upscaleOrder = UPSCALE_ORDER_TENTH;
    pData->gridOption = GRID_GEOSPHERE_8;
    pData->pmapAvgCoeff = 0.666f;
    pData->minFreq_hz = 100.0f;
    pData->maxFreq_hz = 8e3f;
    pData->dispWidth = 120;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    pData->HFOVoption = HFOV_360;
    pData->aspectRatioOption = ASPECT_RATIO_2_1;
}

void dirass_destroy
(
    void ** const phDir
)
{
    dirass_data *pData = (dirass_data*)(*phDir);
    codecPars* pars = pData->pars;
    int i;
    
    if (pData != NULL) {
        if(pData->pmap!=NULL)
            free(pData->pmap);
        if(pData->prev_pmap!=NULL)
            free(pData->prev_pmap);
        for(i=0; i<NUM_DISP_SLOTS; i++)
            if(pData->pmap_grid[i] !=NULL)
                free(pData->pmap_grid[i]);
        
        if(pars->interp_dirs_deg!=NULL)
            free(pars->interp_dirs_deg);
        if(pars->Y_up !=NULL)
            free(pars->Y_up);
        if(pars->interp_table!=NULL)
            free(pars->interp_table);
        if(pars->ss !=NULL)
            free(pars->ss);
        if(pars->ssxyz!=NULL)
            free(pars->ssxyz);
        if(pars->Cxyz !=NULL)
            free(pars->Cxyz);
        if(pars->w!=NULL)
            free(pars->w);
        if(pars->Cw!=NULL)
            free(pars->Cw);
        if(pars->Uw!=NULL)
            free(pars->Uw);
        if(pars->est_dirs!=NULL)
            free(pars->est_dirs);
        
        free(pData->pars);
        free(pData);
        pData = NULL;
    }
}

void dirass_init
(
    void * const hDir,
    float        sampleRate
)
{
    dirass_data *pData = (dirass_data*)(hDir);
    codecPars* pars = pData->pars;

    pData->fs = sampleRate;
    
    /* intialise parameters */
    if(pData->prev_pmap!=NULL)
        memset(pData->prev_pmap, 0, pars->grid_nDirs*sizeof(float));
    if(pars->prev_intensity!=NULL)
        memset(pars->prev_intensity, 0, pars->grid_nDirs*3*sizeof(float));
    memset(pData->Wz12_hpf, 0, MAX_NUM_INPUT_SH_SIGNALS*2*sizeof(float));
    memset(pData->Wz12_lpf, 0, MAX_NUM_INPUT_SH_SIGNALS*2*sizeof(float));
    pData->pmapReady = 0;
    pData->dispSlotIdx = 0;

	/* reinitialise if needed */
	dirass_checkReInit(hDir);
}


void dirass_analysis
(
    void  *  const hDir,
    float ** const inputs,
    int            nInputs,
    int            nSamples,
    int            isPlaying
)
{
    dirass_data *pData = (dirass_data*)(hDir);
    codecPars* pars = pData->pars;
    int i, j, k, n, ch, sec_nSH, secOrder, nSH, up_nSH;
    int o[MAX_INPUT_SH_ORDER+2];
    float intensity[3];
    
    /* local parameters */
    int inputOrder, DirAssMode, upscaleOrder;
    float pmapAvgCoeff, minFreq_hz, maxFreq_hz;
    NORM_TYPES norm;
    
    /* reinitialise if needed */
#ifdef __APPLE__
    dirass_checkReInit(hDir);
#endif
    
    /* The main processing: */
    if (nSamples == FRAME_SIZE && (pData->reInitAna == 0)  && isPlaying ) {
        for(n=0; n<MAX_INPUT_SH_ORDER+2; n++){  o[n] = n*n;  }
        
        /* copy current parameters to be thread safe */
        norm = pData->norm;
        pmapAvgCoeff = pData->pmapAvgCoeff;
        DirAssMode = pData->DirAssMode;
        upscaleOrder = pData->upscaleOrder;
        minFreq_hz = pData->minFreq_hz;
        maxFreq_hz = pData->maxFreq_hz;
        inputOrder = pData->inputOrder;
        secOrder = inputOrder-1;
        nSH = (inputOrder+1)*(inputOrder+1);
        sec_nSH = (secOrder+1)*(secOrder+1);
        up_nSH = (upscaleOrder+1)*(upscaleOrder+1);
        
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
                for(n=0; n<inputOrder+2; n++){  o[n] = n*n;  };
                for (n = 0; n<inputOrder+1; n++)
                    for (ch = o[n]; ch<o[n+1]; ch++)
                        for(i = 0; i<FRAME_SIZE; i++)
                            pData->SHframeTD[ch][i] *= sqrtf(2.0f*(float)n+1.0f);
                break;
        }
        
        /* update the dirass powermap */
        if(pData->recalcPmap==1){
            pData->recalcPmap = 0;
            pData->pmapReady = 0;
            
            /* filter input signals */
            float b[3], a[3];
            biQuadCoeffs(BIQUAD_FILTER_HPF, minFreq_hz, pData->fs, 0.7071f, 0.0f, b, a);
            for(i=0; i<nSH; i++)
                applyBiQuadFilter(b, a, pData->Wz12_hpf[i], pData->SHframeTD[i], FRAME_SIZE);
            biQuadCoeffs(BIQUAD_FILTER_LPF, maxFreq_hz, pData->fs, 0.7071f, 0.0f, b, a);
            for(i=0; i<nSH; i++)
                applyBiQuadFilter(b, a, pData->Wz12_lpf[i], pData->SHframeTD[i], FRAME_SIZE);
            
            /* DoA estimation for each spatially-localised sector */
            if(DirAssMode==REASS_UPSCALE || DirAssMode==REASS_NEAREST){
                /* Beamform using the sector patterns */
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->grid_nDirs, FRAME_SIZE, sec_nSH, 1.0f,
                            pars->Cw, sec_nSH,
                            (const float*)pData->SHframeTD, FRAME_SIZE, 0.0f,
                            pars->ss, FRAME_SIZE);
                
                for(i=0; i<pars->grid_nDirs; i++){
                    /* beamforming to get velocity patterns */
                    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, FRAME_SIZE, nSH, 1.0f,
                                &(pars->Cxyz[i*nSH*3]), 3,
                                (const float*)pData->SHframeTD, FRAME_SIZE, 0.0f,
                                pars->ssxyz, FRAME_SIZE);
                    
                    /* take the sum or mean ss.*ssxyz, to get intensity vector */
                    memset(intensity, 0, 3*sizeof(float));
                    for(k=0; k<3; k++){ 
                        for(j=0; j<FRAME_SIZE; j++)
                            intensity[k] += pars->ssxyz[k*FRAME_SIZE + j] * pars->ss[i*FRAME_SIZE+j];
                        intensity[k] /= (float)FRAME_SIZE;
                        
                        /* average over time */
                        intensity[k] = pmapAvgCoeff * (pars->prev_intensity[i*3+k]) + (1.0f-pmapAvgCoeff) * intensity[k];
                        pars->prev_intensity[i*3+k] = intensity[k];
                    }

                    /* extract DoA [azi elev] convention */
                    pars->est_dirs[i*2] = atan2f(intensity[1], intensity[0]);
                    pars->est_dirs[i*2+1] = atan2f(intensity[2], sqrtf(powf(intensity[0], 2.0f) + powf(intensity[1], 2.0f)));
                    if(DirAssMode==REASS_UPSCALE)
                        pars->est_dirs[i*2+1] = M_PI/2.0f - pars->est_dirs[i*2+1]; /* convert to inclination */
                }
            }
            
            /* Obtain pmap/upscaled pmap in the case of REASS_MODE_OFF and REASS_UPSCALE modes, respectively.
             * OR find the nearest display grid indices, corresponding to the DoA estimates, for the REASS_NEAREST mode */
            switch(DirAssMode) {
                default:
                case REASS_MODE_OFF:
                    /* Standard beamformer-based pmap */
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->grid_nDirs, FRAME_SIZE, nSH, 1.0f,
                                pars->w, nSH,
                                (const float*)pData->SHframeTD, FRAME_SIZE, 0.0f,
                                pars->ss, FRAME_SIZE);
                    
                    /* sum energy over the length of the frame to obtain the pmap */
                    memset(pData->pmap, 0, pars->grid_nDirs *sizeof(float));
                    for(i=0; i<pars->grid_nDirs; i++)
                        for(j=0; j<FRAME_SIZE; j++)
                            pData->pmap[i] += (pars->ss[i*FRAME_SIZE+j])*(pars->ss[i*FRAME_SIZE+j]);
                    
                    /* average the actual pmap over time (averaging is achieved for the reassignment modes via averaging the intensity) */
                    for(i=0; i<pars->grid_nDirs; i++)
                        pData->pmap[i] =  (1.0f-pmapAvgCoeff) * (pData->pmap[i] )+ pmapAvgCoeff * (pData->prev_pmap[i]);
                    memcpy(pData->prev_pmap,  pData->pmap , pars->grid_nDirs*sizeof(float));
                    
                    /* interpolate the pmap */
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->interp_nDirs, 1, pars->grid_nDirs, 1.0f,
                                pars->interp_table, pars->grid_nDirs,
                                pData->pmap, 1, 0.0f,
                                pData->pmap_grid[pData->dispSlotIdx], 1);
                    break;
                    
                case REASS_UPSCALE:
                    /* upscale */
                    getSHreal_recur(upscaleOrder, pars->est_dirs, pars->grid_nDirs, pars->Y_up);
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, up_nSH, FRAME_SIZE, pars->grid_nDirs, 1.0f,
                                pars->Y_up, pars->grid_nDirs,
                                pars->ss, FRAME_SIZE, 0.0f,
                                (float*)pData->SHframe_upTD, FRAME_SIZE);

                    /* Beamform using the new spatially upscaled frame */
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->grid_nDirs, FRAME_SIZE, up_nSH, 1.0f,
                                pars->Uw, up_nSH,
                                (float*)pData->SHframe_upTD, FRAME_SIZE, 0.0f,
                                pars->ss, FRAME_SIZE);
                    
                    /* sum energy over the length of the frame to obtain the pmap */
                    memset(pData->pmap, 0, pars->grid_nDirs *sizeof(float));
                    for(i=0; i<pars->grid_nDirs; i++)
                        for(j=0; j<FRAME_SIZE; j++)
                            pData->pmap[i] += (pars->ss[i*FRAME_SIZE+j])*(pars->ss[i*FRAME_SIZE+j]);
                    
                    /* interpolate the pmap */
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->interp_nDirs, 1, pars->grid_nDirs, 1.0f,
                                pars->interp_table, pars->grid_nDirs,
                                pData->pmap, 1, 0.0f,
                                pData->pmap_grid[pData->dispSlotIdx], 1);
                    break;
                 
                case REASS_NEAREST:
                    /* Assign the sector energies to the nearest display grid point */
                    findClosestGridPoints(pars->interp_dirs_rad, pars->interp_nDirs, pars->est_dirs, pars->grid_nDirs, 0, pars->est_dirs_idx, NULL, NULL);
                    memset(pData->pmap_grid[pData->dispSlotIdx], 0, pars->interp_nDirs * sizeof(float));
                    for(i=0; i< pars->grid_nDirs; i++)
                        for(j=0; j<FRAME_SIZE; j++)
                            pData->pmap_grid[pData->dispSlotIdx][pars->est_dirs_idx[i]] = (pars->ss[i*FRAME_SIZE+j])*(pars->ss[i*FRAME_SIZE+j]);
                    break;
            }
             
            /* ascertain the minimum and maximum values for pmap colour scaling */
            pData->pmap_grid_minVal = FLT_MAX;
            pData->pmap_grid_maxVal = FLT_MIN;
            for(i=0; i<pars->interp_nDirs; i++){
                pData->pmap_grid_minVal = pData->pmap_grid[pData->dispSlotIdx][i] < pData->pmap_grid_minVal ? pData->pmap_grid[pData->dispSlotIdx][i] : pData->pmap_grid_minVal;
                pData->pmap_grid_maxVal = pData->pmap_grid[pData->dispSlotIdx][i] > pData->pmap_grid_maxVal ? pData->pmap_grid[pData->dispSlotIdx][i] : pData->pmap_grid_maxVal;
            }

            /* normalise the pmap to 0..1 */
            for(i=0; i<pars->interp_nDirs; i++)
                pData->pmap_grid[pData->dispSlotIdx][i] = (pData->pmap_grid[pData->dispSlotIdx][i]-pData->pmap_grid_minVal)/(pData->pmap_grid_maxVal-pData->pmap_grid_minVal+1e-11f);

            /* signify that the pmap in the current slot is ready for plotting */
            pData->dispSlotIdx++;
            if(pData->dispSlotIdx>=NUM_DISP_SLOTS)
                pData->dispSlotIdx = 0;
            pData->pmapReady = 1;
        }
    }
}

/* SETS */
 
void dirass_refreshSettings(void* const hDir)
{
	dirass_data *pData = (dirass_data*)(hDir);
	pData->reInitAna = 1; 
}
 
void dirass_checkReInit(void* const hDir)
{
	dirass_data *pData = (dirass_data*)(hDir); 
	/* reinitialise if needed */
	if (pData->reInitAna == 1) {
		pData->reInitAna = 2;  /* indicate init in progress */
		pData->pmapReady = 0;  /* avoid trying to draw pmap during reinit */
		dirass_initAna(hDir);
		pData->reInitAna = 0;  /* indicate init complete */
		pData->recalcPmap = 1; /* recalculate dirass with new configuration */
	}
}

void dirass_setBeamType(void* const hDir, int newType)
{
    dirass_data *pData = (dirass_data*)(hDir);
    codecPars* pars = pData->pars;
    pData->beamType = (BEAM_TYPES)newType;
    if(pData->prev_pmap!=NULL)
        memset(pData->prev_pmap, 0, pars->grid_nDirs*sizeof(float));
    pData->reInitAna = 1;
}

void dirass_setInputOrder(void* const hDir,  int newValue)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->new_inputOrder = newValue;
    pData->reInitAna = 1;
}

void dirass_setDisplayGridOption(void* const hDir,  int newState)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->gridOption = newState;
    pData->reInitAna = 1;
}

void dirass_setDispWidth(void* const hDir,  int newValue)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->dispWidth = newValue;
    pData->reInitAna = 1;
}

void dirass_setUpscaleOrder(void* const hDir,  int newValue)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->new_upscaleOrder = newValue;
    pData->reInitAna = 1;
}

void dirass_setDiRAssMode(void* const hDir,  int newMode)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->DirAssMode = newMode;
}

void dirass_setMinFreq(void* const hDir,  float newValue)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->minFreq_hz = newValue;
}

void dirass_setMaxFreq(void* const hDir,  float newValue)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->maxFreq_hz = newValue;
}

void dirass_setChOrder(void* const hDir, int newOrder)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->chOrdering = (CH_ORDER)newOrder;
}

void dirass_setNormType(void* const hDir, int newType)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->norm = (NORM_TYPES)newType;
}

void dirass_setDispFOV(void* const hDir, int newOption)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->HFOVoption = (HFOV_OPTIONS)newOption;
}

void dirass_setAspectRatio(void* const hDir, int newOption)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->aspectRatioOption = (ASPECT_RATIO_OPTIONS)newOption;
}

void dirass_setMapAvgCoeff(void* const hDir, float newValue)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->pmapAvgCoeff = MIN(MAX(0.0f, newValue), 0.999f);
}

void dirass_requestPmapUpdate(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->recalcPmap = 1;
}


/* GETS */

int dirass_getInputOrder(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return pData->new_inputOrder;
}

int dirass_getBeamType(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (int)pData->beamType;
}

int dirass_getDisplayGridOption(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (int)pData->gridOption;
}

int dirass_getDispWidth(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return pData->dispWidth;
}

int dirass_getUpscaleOrder(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (int)pData->gridOption;
}

int dirass_getDiRAssMode(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return pData->DirAssMode;
}

float dirass_getMinFreq(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return pData->minFreq_hz;
}

float dirass_getMaxFreq(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return pData->maxFreq_hz;
}

int dirass_getSamplingRate(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (int)(pData->fs+0.5f);
}

int dirass_getNSHrequired(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (pData->inputOrder+1)*(pData->inputOrder+1);
}

int dirass_getChOrder(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (int)pData->chOrdering;
}

int dirass_getNormType(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (int)pData->norm;
}

int dirass_getDispFOV(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (int)pData->HFOVoption;
}

int dirass_getAspectRatio(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (int)pData->aspectRatioOption;
}

float dirass_getMapAvgCoeff(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return pData->pmapAvgCoeff;
}

int dirass_getPmap(void* const hDir, float** grid_dirs, float** pmap, int* nDirs,int* pmapWidth, int* hfov, int* aspectRatio)
//TODO: hfov and aspectRatio should be float, if 16:9 etc options are added
{
    dirass_data *pData = (dirass_data*)(hDir);
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




