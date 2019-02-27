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
    pars->Y_grid = NULL;
    pars->interp_table = NULL;
    pars->Cw = NULL;
    pars->Cxyz = NULL;
    
    /* internal */
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
    pData->inputOrder = pData->new_inputOrder = INPUT_ORDER_THIRD;
    pData->pmap_mode = PM_MODE_PWD;
    pData->enableDirAss = 1;
    pData->upscaleOrder = pData->new_upscaleOrder = UPSCALE_ORDER_SECOND;
    pData->gridOption = GRID_GEOSPHERE_9;
    pData->pmapAvgCoeff = 0.666f;
    pData->minFreq_hz = 100.0f;
    pData->maxFreq_hz = 20e3f;
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
        if(pars->Y_grid !=NULL)
            free(pars->Y_grid);
        if(pars->interp_table!=NULL)
            free(pars->interp_table);
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
    int i, j, t, n, ch, sample, band, nSH_order, order_band;
    int o[MAX_INPUT_SH_ORDER+2];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    float_complex* C_grp;
    
    /* local parameters */
    int nSources, inputOrder, nSH;
    float covAvgCoeff, pmapAvgCoeff;
    NORM_TYPES norm;
    MAP_MODES pmap_mode;
    
    /* reinitialise if needed */
#ifdef __APPLE__
    dirass_checkReInit(hDir);
#endif
    
    /* The main processing: */
    if (nSamples == FRAME_SIZE && (pData->reInitAna == 0)  && isPlaying ) {
        /* copy current parameters to be thread safe */
        norm = pData->norm;
        pmapAvgCoeff = pData->pmapAvgCoeff;
        pmap_mode = pData->pmap_mode;
        inputOrder = pData->inputOrder;
        nSH = (inputOrder+1)*(inputOrder+1);
        
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
        
        
        return;
        /* update the dirass powermap */
        if(pData->recalcPmap==1){
            pData->recalcPmap = 0;
            pData->pmapReady = 0;
 
            
            
            pData->pmap; 
            
            
            /* average dirass over time */
            for(i=0; i<pars->grid_nDirs; i++)
                pData->pmap[i] =  (1.0f-pmapAvgCoeff) * (pData->pmap[i] )+ pmapAvgCoeff * (pData->prev_pmap[i]);
            memcpy(pData->prev_pmap,  pData->pmap , pars->grid_nDirs*sizeof(float));

            /* interpolate dirass */
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->interp_nDirs, 1, pars->grid_nDirs, 1.0f,
                        pars->interp_table, pars->grid_nDirs,
                        pData->pmap, 1, 0.0f,
                        pData->pmap_grid[pData->dispSlotIdx], 1);

            /* ascertain minimum and maximum values for dirass colour scaling */
            pData->pmap_grid_minVal = FLT_MAX;
            pData->pmap_grid_maxVal = FLT_MIN;
            for(i=0; i<pars->interp_nDirs; i++){
                pData->pmap_grid_minVal = pData->pmap_grid[pData->dispSlotIdx][i] < pData->pmap_grid_minVal ? pData->pmap_grid[pData->dispSlotIdx][i] : pData->pmap_grid_minVal;
                pData->pmap_grid_maxVal = pData->pmap_grid[pData->dispSlotIdx][i] > pData->pmap_grid_maxVal ? pData->pmap_grid[pData->dispSlotIdx][i] : pData->pmap_grid_maxVal;
            }

            /* normalise the dirass to 0..1 */
            for(i=0; i<pars->interp_nDirs; i++)
                pData->pmap_grid[pData->dispSlotIdx][i] = (pData->pmap_grid[pData->dispSlotIdx][i]-pData->pmap_grid_minVal)/(pData->pmap_grid_maxVal-pData->pmap_grid_minVal+1e-11f);

            /* signify that the dirass in current slot is ready for plotting */
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

void dirass_setMapMode(void* const hDir, int newMode)
{
    dirass_data *pData = (dirass_data*)(hDir);
    codecPars* pars = pData->pars;
    pData->pmap_mode = (MAP_MODES)newMode;
    if(pData->prev_pmap!=NULL)
        memset(pData->prev_pmap, 0, pars->grid_nDirs*sizeof(float));
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

void dirass_setUpscaleOrder(void* const hDir,  int newState)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->new_upscaleOrder = newState;
    pData->reInitAna = 1;
}

void dirass_setEnableDiRAss(void* const hDir,  int newState)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->enableDirAss = newState;
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
    pData->pmapAvgCoeff = MIN(MAX(0.0f, newValue), 0.99999999f);
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

int dirass_getMapMode(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (int)pData->pmap_mode;
}

int dirass_getDisplayGridOption(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (int)pData->gridOption;
}

int dirass_getUpscaleOrder(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return (int)pData->gridOption;
}

int dirass_getEnableDiRAss(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return pData->enableDirAss;
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




