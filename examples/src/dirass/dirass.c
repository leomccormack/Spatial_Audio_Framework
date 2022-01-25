/*
 * Copyright 2019 Leo McCormack
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
 * @file dirass.c
 * @brief A sound-field visualiser based on the directional re-assignment of
 *        beamformer energy based on local DoA estimates [1,2]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 * @see [2] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of
 *          angular spectra based on a directional re-assignment approach for
 *          ambisonic sound-field visualisation". IEEE International Conference
 *          on Acoustics, Speech and Signal Processing (ICASSP).
 *
 * @author Leo McCormack
 * @date 21.02.2019
 * @license ISC
 */

#include "dirass.h"
#include "dirass_internal.h"

void dirass_create
(
    void ** const phDir
)
{
    dirass_data* pData = (dirass_data*)malloc1d(sizeof(dirass_data));
    *phDir = (void*)pData;
    int i;

    /* Default user parameters */
    pData->inputOrder = pData->new_inputOrder = SH_ORDER_FIRST;
    pData->beamType = STATIC_BEAM_TYPE_HYPERCARDIOID;
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

    /* codec data */
    pData->pars = (dirass_codecPars*)malloc1d(sizeof(dirass_codecPars));
    dirass_codecPars* pars = pData->pars;
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
    pars->prev_energy = NULL;
    
    /* internal */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->procStatus = PROC_STATUS_NOT_ONGOING;

    /* display */
    pData->pmap = NULL;
    for(i=0; i<NUM_DISP_SLOTS; i++)
        pData->pmap_grid[i] = NULL;
    pData->pmapReady = 0;
    pData->recalcPmap = 1;

    /* set FIFO buffers */
    pData->FIFO_idx = 0;
    memset(pData->inFIFO, 0, MAX_NUM_INPUT_SH_SIGNALS*DIRASS_FRAME_SIZE*sizeof(float));
}

void dirass_destroy
(
    void ** const phDir
)
{
    dirass_data *pData = (dirass_data*)(*phDir);
    dirass_codecPars* pars;
    int i;
    
    if (pData != NULL) {
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            SAF_SLEEP(10);
        }
        
        if(pData->pmap!=NULL)
            free(pData->pmap);
        for(i=0; i<NUM_DISP_SLOTS; i++)
            free(pData->pmap_grid[i]);

        pars = pData->pars;
        free(pars->interp_dirs_deg);
        free(pars->interp_dirs_rad);
        free(pars->Y_up);
        free(pars->interp_table);
        free(pars->ss);
        free(pars->ssxyz);
        free(pars->Cxyz);
        free(pars->w);
        free(pars->Cw);
        free(pars->Uw);
        free(pars->est_dirs);
        free(pars->est_dirs_idx);
        free(pars->prev_intensity);
        free(pars->prev_energy);
        
        free(pData->pars);
        free(pData->progressBarText);
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
    dirass_codecPars* pars = pData->pars;

    pData->fs = sampleRate;
    
    /* intialise parameters */
    if(pars->prev_intensity!=NULL)
        memset(pars->prev_intensity, 0, pars->grid_nDirs*3*sizeof(float));
    if(pars->prev_energy!=NULL)
        memset(pars->prev_energy, 0, pars->grid_nDirs*sizeof(float));
    memset(pData->Wz12_hpf, 0, MAX_NUM_INPUT_SH_SIGNALS*2*sizeof(float));
    memset(pData->Wz12_lpf, 0, MAX_NUM_INPUT_SH_SIGNALS*2*sizeof(float));
    pData->pmapReady = 0;
    pData->dispSlotIdx = 0;
}

void dirass_initCodec
(
    void* const hDir
)
{
    dirass_data *pData = (dirass_data*)(hDir);
    
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
    
    dirass_initAna(hDir);
    
    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
}

void dirass_analysis
(
    void        *  const hDir,
    const float *const * inputs,
    int                  nInputs,
    int                  nSamples,
    int                  isPlaying
)
{
    dirass_data *pData = (dirass_data*)(hDir);
    dirass_codecPars* pars = pData->pars;
    int s, i, j, k, ch, sec_nSH, secOrder, nSH, up_nSH;
    float intensity[3];
    
    /* local copy of user parameters */
    int inputOrder, DirAssMode, upscaleOrder;
    float pmapAvgCoeff, minFreq_hz, maxFreq_hz;
    NORM_TYPES norm;
    CH_ORDER chOrdering;
    norm = pData->norm;
    chOrdering = pData->chOrdering;
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
        if (pData->FIFO_idx >= DIRASS_FRAME_SIZE && (pData->codecStatus == CODEC_STATUS_INITIALISED) && isPlaying) {
            pData->FIFO_idx = 0;
            pData->procStatus = PROC_STATUS_ONGOING;

            /* Load time-domain data */
            for(ch=0; ch<nSH; ch++)
                memcpy(pData->SHframeTD[ch],pData->inFIFO[ch], DIRASS_FRAME_SIZE*sizeof(float));

            /* account for input channel order */
            switch(chOrdering){
                case CH_ACN:  /* already ACN */ break; /* Otherwise, convert to ACN... */
                case CH_FUMA: convertHOAChannelConvention((float*)pData->SHframeTD, inputOrder, DIRASS_FRAME_SIZE, HOA_CH_ORDER_FUMA, HOA_CH_ORDER_ACN); break;
            }

            /* account for input normalisation scheme */
            switch(norm){
                case NORM_N3D:  /* already in N3D, do nothing */ break; /* Otherwise, convert to N3D... */
                case NORM_SN3D: convertHOANormConvention((float*)pData->SHframeTD, inputOrder, DIRASS_FRAME_SIZE, HOA_NORM_SN3D, HOA_NORM_N3D); break;
                case NORM_FUMA: convertHOANormConvention((float*)pData->SHframeTD, inputOrder, DIRASS_FRAME_SIZE, HOA_NORM_FUMA, HOA_NORM_N3D); break;
            }

            /* update the dirass powermap */
            if(pData->recalcPmap==1){
                pData->recalcPmap = 0;
                pData->pmapReady = 0;

                /* filter input signals */
                float b[3], a[3];
                biQuadCoeffs(BIQUAD_FILTER_HPF, minFreq_hz, pData->fs, 0.7071f, 0.0f, b, a);
                for(i=0; i<nSH; i++)
                    applyBiQuadFilter(b, a, pData->Wz12_hpf[i], pData->SHframeTD[i], DIRASS_FRAME_SIZE);
                biQuadCoeffs(BIQUAD_FILTER_LPF, maxFreq_hz, pData->fs, 0.7071f, 0.0f, b, a);
                for(i=0; i<nSH; i++)
                    applyBiQuadFilter(b, a, pData->Wz12_lpf[i], pData->SHframeTD[i], DIRASS_FRAME_SIZE);

                /* DoA estimation for each spatially-localised sector */
                if(DirAssMode==REASS_UPSCALE || DirAssMode==REASS_NEAREST){
                    /* Beamform using the sector patterns */
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->grid_nDirs, DIRASS_FRAME_SIZE, sec_nSH, 1.0f,
                                pars->Cw, sec_nSH,
                                (const float*)pData->SHframeTD, DIRASS_FRAME_SIZE, 0.0f,
                                pars->ss, DIRASS_FRAME_SIZE);

                    for(i=0; i<pars->grid_nDirs; i++){
                        /* beamforming to get velocity patterns */
                        cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, DIRASS_FRAME_SIZE, nSH, 1.0f,
                                    &(pars->Cxyz[i*nSH*3]), 3,
                                    (const float*)pData->SHframeTD, DIRASS_FRAME_SIZE, 0.0f,
                                    pars->ssxyz, DIRASS_FRAME_SIZE);

                        /* take the sum or mean ss.*ssxyz, to get intensity vector */
                        memset(intensity, 0, 3*sizeof(float));
                        for(k=0; k<3; k++){
                            for(j=0; j<DIRASS_FRAME_SIZE; j++)
                                intensity[k] += pars->ssxyz[k*DIRASS_FRAME_SIZE + j] * pars->ss[i*DIRASS_FRAME_SIZE+j];
                            intensity[k] /= (float)DIRASS_FRAME_SIZE;

                            /* average over time */
                            intensity[k] = pmapAvgCoeff * (pars->prev_intensity[i*3+k]) + (1.0f-pmapAvgCoeff) * intensity[k];
                            pars->prev_intensity[i*3+k] = intensity[k];
                        }

                        /* extract DoA [azi elev] convention */
                        pars->est_dirs[i*2] = atan2f(intensity[1], intensity[0]);
                        pars->est_dirs[i*2+1] = atan2f(intensity[2], sqrtf(powf(intensity[0], 2.0f) + powf(intensity[1], 2.0f)));
                        if(DirAssMode==REASS_UPSCALE)
                            pars->est_dirs[i*2+1] = SAF_PI/2.0f - pars->est_dirs[i*2+1]; /* convert to inclination */
                    }
                }

                /* Obtain pmap/upscaled pmap in the case of REASS_MODE_OFF and REASS_UPSCALE modes, respectively.
                 * OR find the nearest display grid indices, corresponding to the DoA estimates, for the REASS_NEAREST mode */
                switch(DirAssMode) {
                    default:
                    case REASS_MODE_OFF:
                        /* Standard beamformer-based pmap */
                        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->grid_nDirs, DIRASS_FRAME_SIZE, nSH, 1.0f,
                                    pars->w, nSH,
                                    (const float*)pData->SHframeTD, DIRASS_FRAME_SIZE, 0.0f,
                                    pars->ss, DIRASS_FRAME_SIZE);

                        /* sum energy over the length of the frame to obtain the pmap */
                        memset(pData->pmap, 0, pars->grid_nDirs *sizeof(float));
                        for(i=0; i<pars->grid_nDirs; i++)
                            for(j=0; j<DIRASS_FRAME_SIZE; j++)
                                pData->pmap[i] += (pars->ss[i*DIRASS_FRAME_SIZE+j])*(pars->ss[i*DIRASS_FRAME_SIZE+j]);

                        /* average energy over time */
                        for(i=0; i<pars->grid_nDirs; i++){
                            pData->pmap[i] = pmapAvgCoeff * (pars->prev_energy[i]) + (1.0f-pmapAvgCoeff) * (pData->pmap[i]);
                            pars->prev_energy[i] = pData->pmap[i];
                        }

                        /* interpolate the pmap */
                        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->interp_nDirs, 1, pars->grid_nDirs, 1.0f,
                                    pars->interp_table, pars->grid_nDirs,
                                    pData->pmap, 1, 0.0f,
                                    pData->pmap_grid[pData->dispSlotIdx], 1);
                        break;

                    case REASS_UPSCALE:
                        /* upscale */
                        getSHreal_recur(upscaleOrder, pars->est_dirs, pars->grid_nDirs, pars->Y_up);
                        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, up_nSH, DIRASS_FRAME_SIZE, pars->grid_nDirs, 1.0f,
                                    pars->Y_up, pars->grid_nDirs,
                                    pars->ss, DIRASS_FRAME_SIZE, 0.0f,
                                    (float*)pData->SHframe_upTD, DIRASS_FRAME_SIZE);

                        /* Beamform using the new spatially upscaled frame */
                        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pars->grid_nDirs, DIRASS_FRAME_SIZE, up_nSH, 1.0f,
                                    pars->Uw, up_nSH,
                                    (float*)pData->SHframe_upTD, DIRASS_FRAME_SIZE, 0.0f,
                                    pars->ss, DIRASS_FRAME_SIZE);

                        /* sum energy over the length of the frame to obtain the pmap */
                        memset(pData->pmap, 0, pars->grid_nDirs *sizeof(float));
                        for(i=0; i<pars->grid_nDirs; i++)
                            for(j=0; j<DIRASS_FRAME_SIZE; j++)
                                pData->pmap[i] += (pars->ss[i*DIRASS_FRAME_SIZE+j])*(pars->ss[i*DIRASS_FRAME_SIZE+j]);

                        /* average energy over time */
                        for(i=0; i<pars->grid_nDirs; i++){
                            pData->pmap[i] = pmapAvgCoeff * (pars->prev_energy[i]) + (1.0f-pmapAvgCoeff) * (pData->pmap[i]);
                            pars->prev_energy[i] = pData->pmap[i];
                        }

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
                            for(j=0; j<DIRASS_FRAME_SIZE; j++)
                                pData->pmap[i] = (pars->ss[i*DIRASS_FRAME_SIZE+j])*(pars->ss[i*DIRASS_FRAME_SIZE+j]);

                        /* average energy over time, and assign to nearest grid direction */
                        for(i=0; i<pars->grid_nDirs; i++){
                            pData->pmap[i] = pmapAvgCoeff * (pars->prev_energy[i]) + (1.0f-pmapAvgCoeff) * (pData->pmap[i]);
                            pars->prev_energy[i] = pData->pmap[i];
                            pData->pmap_grid[pData->dispSlotIdx][pars->est_dirs_idx[i]] += pData->pmap[i];
                        }
                        break;
                }

                /* ascertain the minimum and maximum values for pmap colour scaling */
                int ind;
                utility_siminv(pData->pmap_grid[pData->dispSlotIdx], pars->interp_nDirs, &ind);
                pData->pmap_grid_minVal = pData->pmap_grid[pData->dispSlotIdx][ind];
                utility_simaxv(pData->pmap_grid[pData->dispSlotIdx], pars->interp_nDirs, &ind);
                pData->pmap_grid_maxVal = pData->pmap_grid[pData->dispSlotIdx][ind];

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
        else if(pData->FIFO_idx >= DIRASS_FRAME_SIZE){
            /* reset FIFO_idx index if codec was not ready */
            pData->FIFO_idx = 0;
        }
    }
    
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

/* SETS */
 
void dirass_refreshSettings(void* const hDir)
{
    dirass_setCodecStatus(hDir, CODEC_STATUS_NOT_INITIALISED);
}

void dirass_setBeamType(void* const hDir, int newType)
{
    dirass_data *pData = (dirass_data*)(hDir);
    if(pData->beamType != (STATIC_BEAM_TYPES)newType){
        pData->beamType = (STATIC_BEAM_TYPES)newType;
        dirass_setCodecStatus(hDir, CODEC_STATUS_NOT_INITIALISED);
    }
}

void dirass_setInputOrder(void* const hDir,  int newValue)
{
    dirass_data *pData = (dirass_data*)(hDir);
    if(pData->new_inputOrder != newValue){
        pData->new_inputOrder = newValue;
        dirass_setCodecStatus(hDir, CODEC_STATUS_NOT_INITIALISED);
    }
    /* FUMA only supports 1st order */
    if(pData->new_inputOrder!=SH_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->new_inputOrder!=SH_ORDER_FIRST && pData->norm == NORM_FUMA)
        pData->norm = NORM_SN3D;
}

void dirass_setDisplayGridOption(void* const hDir,  int newState)
{
    dirass_data *pData = (dirass_data*)(hDir);
    if(pData->gridOption != (DIRASS_GRID_OPTIONS)newState){
        pData->gridOption = (DIRASS_GRID_OPTIONS)newState;
        dirass_setCodecStatus(hDir, CODEC_STATUS_NOT_INITIALISED);
    }
}

void dirass_setDispWidth(void* const hDir,  int newValue)
{
    dirass_data *pData = (dirass_data*)(hDir);
    if(pData->dispWidth != newValue){
        pData->dispWidth = newValue;
        dirass_setCodecStatus(hDir, CODEC_STATUS_NOT_INITIALISED);
    }
}

void dirass_setUpscaleOrder(void* const hDir,  int newValue)
{
    dirass_data *pData = (dirass_data*)(hDir);
    if(pData->new_upscaleOrder != newValue){
        pData->new_upscaleOrder = newValue;
        dirass_setCodecStatus(hDir, CODEC_STATUS_NOT_INITIALISED);
    }
}

void dirass_setDiRAssMode(void* const hDir,  int newMode)
{
    dirass_data *pData = (dirass_data*)(hDir);
    dirass_codecPars* pars = pData->pars;
    if(pData->DirAssMode!=(DIRASS_REASS_MODES)newMode){
        pData->DirAssMode = (DIRASS_REASS_MODES)newMode;
        if(pars->prev_intensity!=NULL)
            memset(pars->prev_intensity, 0, pars->grid_nDirs*3*sizeof(float));
        memset(pars->prev_energy, 0, pars->grid_nDirs*sizeof(float));
    }
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
    if((CH_ORDER)newOrder != CH_FUMA || pData->new_inputOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void dirass_setNormType(void* const hDir, int newType)
{
    dirass_data *pData = (dirass_data*)(hDir);
    if((NORM_TYPES)newType != NORM_FUMA || pData->new_inputOrder==SH_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
}

void dirass_setDispFOV(void* const hDir, int newOption)
{
    dirass_data *pData = (dirass_data*)(hDir);
    if(pData->HFOVoption != (HFOV_OPTIONS)newOption){
        pData->HFOVoption = (HFOV_OPTIONS)newOption;
        dirass_setCodecStatus(hDir, CODEC_STATUS_NOT_INITIALISED);
    }
}

void dirass_setAspectRatio(void* const hDir, int newOption)
{
    dirass_data *pData = (dirass_data*)(hDir);
    if(pData->aspectRatioOption != (ASPECT_RATIO_OPTIONS)newOption){
        pData->aspectRatioOption = (ASPECT_RATIO_OPTIONS)newOption;
        dirass_setCodecStatus(hDir, CODEC_STATUS_NOT_INITIALISED);
    }
}

void dirass_setMapAvgCoeff(void* const hDir, float newValue)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->pmapAvgCoeff = SAF_MIN(SAF_MAX(0.0f, newValue), 0.999f);
}

void dirass_requestPmapUpdate(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    pData->recalcPmap = 1;
}


/* GETS */

int dirass_getFrameSize(void)
{
    return DIRASS_FRAME_SIZE;
}

CODEC_STATUS dirass_getCodecStatus(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return pData->codecStatus;
}

float dirass_getProgressBar0_1(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    return pData->progressBar0_1;
}

void dirass_getProgressBarText(void* const hDir, char* text)
{
    dirass_data *pData = (dirass_data*)(hDir);
    memcpy(text, pData->progressBarText, PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

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

int dirass_getPmap(void* const hDir, float** grid_dirs, float** pmap, int* nDirs,int* pmapWidth, int* hfov, float* aspectRatio) 
{
    dirass_data *pData = (dirass_data*)(hDir);
    dirass_codecPars* pars = pData->pars;
    if((pData->codecStatus == CODEC_STATUS_INITIALISED) && pData->pmapReady){
        (*grid_dirs) = pars->interp_dirs_deg;
        (*pmap) = pData->pmap_grid[pData->dispSlotIdx-1 < 0 ? NUM_DISP_SLOTS-1 : pData->dispSlotIdx-1];
        (*nDirs) = pars->interp_nDirs;
        (*pmapWidth) = pData->dispWidth;
        switch(pData->HFOVoption){
            default:
            case HFOV_360: (*hfov) = 360; break;
            case HFOV_180: (*hfov) = 180; break;
            case HFOV_90:  (*hfov) = 90;  break;
            case HFOV_60:  (*hfov) = 60;  break;
        }
        switch(pData->aspectRatioOption){
            default:
            case ASPECT_RATIO_2_1:  (*aspectRatio) = 2.0f; break;
            case ASPECT_RATIO_16_9: (*aspectRatio) = 16.0f/9.0f; break;
            case ASPECT_RATIO_4_3:  (*aspectRatio) = 4.0f/3.0f; break;
        }
    }
    return pData->pmapReady;
}

int dirass_getProcessingDelay()
{
    return 2*DIRASS_FRAME_SIZE;
}
