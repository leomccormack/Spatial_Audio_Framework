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

/*
 * Filename: panner.c
 * ------------------
 * A frequency-dependent 3D panner, based on the Vector-base Amplitude Panning
 * (VBAP) method. Depending on the room, it may be beneficial to employ
 * amplitude-normalised gains for low frequencies, and energy-normalised gains
 * for high frequencies. Therefore, this VBAP implementation uses the method
 * described in [1], to do just that.
 *
 * Dependencies:
 *     saf_utilities, saf_vbap, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 25.09.2017
 *
 * [1] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V. (2014).
 *     Gain normalisation in amplitude panning as a function of frequency and
 *     room reverberance. 55th International Conference of the AES. Helsinki,
 *     Finland.
 */
 
#include "panner_internal.h"

void panner_create
(
    void ** const phPan
)
{
    panner_data* pData = (panner_data*)malloc1d(sizeof(panner_data));
    *phPan = (void*)pData;
    int ch, dummy;
    
    /* time-frequency transform + buffers */
    pData->hSTFT = NULL;
    pData->STFTInputFrameTF = malloc1d(MAX_NUM_INPUTS*sizeof(complexVector));
    pData->STFTOutputFrameTF = malloc1d(MAX_NUM_OUTPUTS*sizeof(complexVector));
    for(ch=0; ch< MAX_NUM_INPUTS; ch++) {
        pData->STFTInputFrameTF[ch].re = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
        pData->STFTInputFrameTF[ch].im = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
    }
    for(ch=0; ch< MAX_NUM_OUTPUTS; ch++) {
        pData->STFTOutputFrameTF[ch].re = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
        pData->STFTOutputFrameTF[ch].im = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
    }
    pData->tempHopFrameTD = (float**)malloc2d( MAX(MAX_NUM_INPUTS, MAX_NUM_OUTPUTS), HOP_SIZE, sizeof(float));
    
    /* flags and gain table */
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
    pData->reInitGainTables = 1;
    pData->vbap_gtable = NULL;
    pData->reInitTFT = 1;
    pData->recalc_M_rotFLAG = 1;
    
    /* default user parameters */
    panner_loadPreset(PRESET_DEFAULT, pData->src_dirs_deg, &(pData->new_nSources), &(dummy)); /*check setStateInformation if you change default preset*/
    pData->nSources = pData->new_nSources;
    pData->DTT = 0.0f;
    pData->spread_deg = 0.0f;
    panner_loadPreset(PRESET_5PX, pData->loudpkrs_dirs_deg, &(pData->new_nLoudpkrs), &(pData->output_nDims)); /*check setStateInformation if you change default preset*/
    pData->nLoudpkrs = pData->new_nLoudpkrs;
    pData->yaw = 0.0f;
    pData->pitch = 0.0f;
    pData->roll = 0.0f;
    pData->bFlipYaw = 0;
    pData->bFlipPitch = 0;
    pData->bFlipRoll = 0;
}

void panner_destroy
(
    void ** const phPan
)
{
    panner_data *pData = (panner_data*)(*phPan);
    int ch;

    if (pData != NULL) {
        if(pData->hSTFT !=NULL)
            afSTFTfree(pData->hSTFT);
        for(ch=0; ch< MAX_NUM_INPUTS; ch++) {
            free(pData->STFTInputFrameTF[ch].re);
            free(pData->STFTInputFrameTF[ch].im);
        }
    
        for (ch = 0; ch< MAX_NUM_OUTPUTS; ch++) {
            free(pData->STFTOutputFrameTF[ch].re);
            free(pData->STFTOutputFrameTF[ch].im);
        }
        free(pData->STFTInputFrameTF);
        free(pData->STFTOutputFrameTF);
        free(pData->tempHopFrameTD);
        free1d((void**)&(pData->vbap_gtable));
         
        free(pData);
        pData = NULL;
    }
}

void panner_init
(
    void * const hPan,
    int          sampleRate
)
{
    panner_data *pData = (panner_data*)(hPan);
    int band;
    
    /* define frequency vector */
    pData->fs = sampleRate;
    for(band=0; band <HYBRID_BANDS; band++){
        if(sampleRate == 44100)
            pData->freqVector[band] =  (float)__afCenterFreq44100[band];
        else
            pData->freqVector[band] =  (float)__afCenterFreq48e3[band];
    }
    
    /* calculate pValue per frequency */
    getPvalues(pData->DTT, pData->freqVector, HYBRID_BANDS, pData->pValue);

    /* reinitialise if needed */
    panner_checkReInit(hPan);
    pData->recalc_M_rotFLAG = 1;
}

void panner_process
(
    void  *  const hPan,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples,
    int            isPlaying
)
{
    panner_data *pData = (panner_data*)(hPan);
    int t, ch, ls, i, band, nSources, nLoudspeakers, N_azi, aziIndex, elevIndex, idx3d, idx2D;
    float aziRes, elevRes, pv_f, gains3D_sum_pvf, gains2D_sum_pvf, Rxyz[3][3], hypotxy;
    float src_dirs[MAX_NUM_INPUTS][2], pValue[HYBRID_BANDS], gains3D[MAX_NUM_OUTPUTS], gains2D[MAX_NUM_OUTPUTS];
	const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
	float_complex outputTemp[MAX_NUM_OUTPUTS][TIME_SLOTS];
    
    /* reinitialise if needed */
#ifdef __APPLE__
    panner_checkReInit(hPan);
#else
    if(pData->reInitTFT==1){
        pData->reInitTFT = 2;
        panner_initTFT(hPan);
        pData->reInitTFT = 0;
    }
#endif

    /* apply panner */
    if ((nSamples == FRAME_SIZE) && (pData->vbap_gtable != NULL)
        && (pData->reInitTFT == 0) && (pData->reInitGainTables == 0)) {
        memcpy(src_dirs, pData->src_dirs_deg, MAX_NUM_INPUTS*2*sizeof(float));
        memcpy(pValue, pData->pValue, HYBRID_BANDS*sizeof(float));
        nSources = pData->nSources;
        nLoudspeakers = pData->nLoudpkrs;
        
        /* Load time-domain data */
        for(i=0; i < MIN(nSources,nInputs); i++)
            utility_svvcopy(inputs[i], FRAME_SIZE, pData->inputFrameTD[i]);
        for(; i<MAX_NUM_INPUTS; i++)
            memset(pData->inputFrameTD[i], 0, FRAME_SIZE * sizeof(float));
        
        /* Apply time-frequency transform (TFT) */
        for(t=0; t< TIME_SLOTS; t++) {
            for(ch = 0; ch < nSources; ch++)
                utility_svvcopy(&(pData->inputFrameTD[ch][t*HOP_SIZE]), HOP_SIZE, pData->tempHopFrameTD[ch]);
            afSTFTforward(pData->hSTFT, (float**)pData->tempHopFrameTD, (complexVector*)pData->STFTInputFrameTF);
            for(band=0; band<HYBRID_BANDS; band++)
                for(ch=0; ch < nSources; ch++)
                    pData->inputframeTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[ch].re[band], pData->STFTInputFrameTF[ch].im[band]);
        }
        memset(pData->outputframeTF, 0, HYBRID_BANDS*MAX_NUM_OUTPUTS*TIME_SLOTS * sizeof(float_complex));
		memset(outputTemp, 0, MAX_NUM_OUTPUTS*TIME_SLOTS * sizeof(float_complex));
        
        /* Main processing: */
        if(isPlaying){
            /* Rotate source directions */
            if(pData->recalc_M_rotFLAG){
                yawPitchRoll2Rzyx (pData->yaw, pData->pitch, pData->roll, 0, Rxyz);
                for(i=0; i<nSources; i++){
                    pData->src_dirs_xyz[i][0] = cosf(DEG2RAD(pData->src_dirs_deg[i][1])) * cosf(DEG2RAD(pData->src_dirs_deg[i][0]));
                    pData->src_dirs_xyz[i][1] = cosf(DEG2RAD(pData->src_dirs_deg[i][1])) * sinf(DEG2RAD(pData->src_dirs_deg[i][0]));
                    pData->src_dirs_xyz[i][2] = sinf(DEG2RAD(pData->src_dirs_deg[i][1]));
                    pData->recalc_gainsFLAG[i] = 1;
                }
                cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSources, 3, 3, 1.0f,
                            (float*)(pData->src_dirs_xyz), 3,
                            (float*)Rxyz, 3, 0.0f,
                            (float*)(pData->src_dirs_rot_xyz), 3);
                for(i=0; i<nSources; i++){
                    hypotxy = sqrtf(powf(pData->src_dirs_rot_xyz[i][0], 2.0f) + powf(pData->src_dirs_rot_xyz[i][1], 2.0f));
                    pData->src_dirs_rot_deg[i][0] = RAD2DEG(atan2f(pData->src_dirs_rot_xyz[i][1], pData->src_dirs_rot_xyz[i][0]));
                    pData->src_dirs_rot_deg[i][1] = RAD2DEG(atan2f(pData->src_dirs_rot_xyz[i][2], hypotxy));
                }
                pData->recalc_M_rotFLAG = 0;
            }
            
            /* Apply VBAP Panning */
            if(pData->output_nDims == 3){/* 3-D case */
                aziRes = (float)pData->vbapTableRes[0];
                elevRes = (float)pData->vbapTableRes[1];
                N_azi = (int)(360.0f / aziRes + 0.5f) + 1;
                for (ch = 0; ch < nSources; ch++) {
                    /* recalculate frequency dependent panning gains */
                    if(pData->recalc_gainsFLAG[ch]){
                        //aziIndex = (int)(matlab_fmodf(pData->src_dirs_deg[ch][0] + 180.0f, 360.0f) / aziRes + 0.5f);
                        //elevIndex = (int)((pData->src_dirs_deg[ch][1] + 90.0f) / elevRes + 0.5f);
                        aziIndex = (int)(matlab_fmodf(pData->src_dirs_rot_deg[ch][0] + 180.0f, 360.0f) / aziRes + 0.5f);
                        elevIndex = (int)((pData->src_dirs_rot_deg[ch][1] + 90.0f) / elevRes + 0.5f);
                        idx3d = elevIndex * N_azi + aziIndex;
                        for (ls = 0; ls < nLoudspeakers; ls++)
                            gains3D[ls] =  pData->vbap_gtable[idx3d*nLoudspeakers+ls];
                        for (band = 0; band < HYBRID_BANDS; band++){
                            /* apply pValue per frequency */
                            pv_f = pData->pValue[band];
                            if(pv_f != 2.0f){
                                gains3D_sum_pvf = 0.0f;
                                for (ls = 0; ls < nLoudspeakers; ls++)
                                    gains3D_sum_pvf += powf(MAX(gains3D[ls], 0.0f), pv_f);
                                gains3D_sum_pvf = powf(gains3D_sum_pvf, 1.0f/(pv_f+2.23e-9f));
                                for (ls = 0; ls < nLoudspeakers; ls++)
                                    pData->G_src[band][ch][ls] = cmplxf(gains3D[ls] / (gains3D_sum_pvf+2.23e-9f), 0.0f);
                            }
                            else
                                for (ls = 0; ls < nLoudspeakers; ls++)
                                    pData->G_src[band][ch][ls] = cmplxf(gains3D[ls], 0.0f);
                        }
                        pData->recalc_gainsFLAG[ch] = 0;
                    } 
                }
				/* apply panning gains */
				for (band = 0; band < HYBRID_BANDS; band++) {
					cblas_cgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nLoudspeakers, TIME_SLOTS, nSources, &calpha,
						pData->G_src[band], MAX_NUM_OUTPUTS,
						pData->inputframeTF[band], TIME_SLOTS, &cbeta,
						outputTemp, TIME_SLOTS);
					for (i = 0; i < nLoudspeakers; i++)
						for (t = 0; t < TIME_SLOTS; t++)
							pData->outputframeTF[band][i][t] = ccaddf(pData->outputframeTF[band][i][t], outputTemp[i][t]);
				}
            }
            else{/* 2-D case */
                aziRes = (float)pData->vbapTableRes[0];
                for (ch = 0; ch < nSources; ch++) {
                    /* recalculate frequency dependent panning gains */
                    if(pData->recalc_gainsFLAG[ch]){
                        //idx2D = (int)((matlab_fmodf(pData->src_dirs_deg[ch][0]+180.0f,360.0f)/aziRes)+0.5f);
                        idx2D = (int)((matlab_fmodf(pData->src_dirs_rot_deg[ch][0]+180.0f,360.0f)/aziRes)+0.5f);
                        for (ls = 0; ls < nLoudspeakers; ls++)
                            gains2D[ls] = pData->vbap_gtable[idx2D*nLoudspeakers+ls];
                        for (band = 0; band < HYBRID_BANDS; band++){
                            /* apply pValue per frequency */
                            pv_f = pData->pValue[band];
                            if(pv_f != 2.0f){
                                gains2D_sum_pvf = 0.0f;
                                for (ls = 0; ls < nLoudspeakers; ls++)
                                    gains2D_sum_pvf += powf(MAX(gains2D[ls], 0.0f), pv_f);
                                gains2D_sum_pvf = powf(gains2D_sum_pvf, 1.0f/(pv_f+2.23e-9f));
                                for (ls = 0; ls < nLoudspeakers; ls++)
                                    pData->G_src[band][ch][ls] = cmplxf(gains2D[ls] / (gains2D_sum_pvf+2.23e-9f), 0.0f);
                            }
                            else
                                for (ls = 0; ls < nLoudspeakers; ls++)
                                    pData->G_src[band][ch][ls] = cmplxf(gains2D[ls], 0.0f);
                        }
                        pData->recalc_gainsFLAG[ch] = 0;
                    }
                    /* apply panning gains */
                    for (band = 0; band < HYBRID_BANDS; band++){
                        for (ls = 0; ls < nLoudspeakers; ls++)
                            for (t = 0; t < TIME_SLOTS; t++)
                                pData->outputframeTF[band][ls][t] = ccaddf(pData->outputframeTF[band][ls][t], ccmulf(pData->inputframeTF[band][ch][t], pData->G_src[band][ch][ls]));
                    }
                }
            }
            /* scale by sqrt(number of sources) */
            for (band = 0; band < HYBRID_BANDS; band++)
                for (ls = 0; ls < nLoudspeakers; ls++)
                    for (t = 0; t < TIME_SLOTS; t++)
                        pData->outputframeTF[band][ls][t] = crmulf(pData->outputframeTF[band][ls][t], 1.0f/sqrtf((float)nSources));
        }
        else
            memset(pData->outputframeTF, 0, HYBRID_BANDS*MAX_NUM_OUTPUTS*TIME_SLOTS*sizeof(float_complex));
            
        /* inverse-TFT */
        for(t = 0; t < TIME_SLOTS; t++) {
            for(band = 0; band < HYBRID_BANDS; band++) {
                for(ch = 0; ch < nLoudspeakers; ch++) {
                    pData->STFTOutputFrameTF[ch].re[band] = crealf(pData->outputframeTF[band][ch][t]);
                    pData->STFTOutputFrameTF[ch].im[band] = cimagf(pData->outputframeTF[band][ch][t]);
                }
            }
            afSTFTinverse(pData->hSTFT, pData->STFTOutputFrameTF, pData->tempHopFrameTD);
            for (ch = 0; ch < MIN(nLoudspeakers, nOutputs); ch++)
                utility_svvcopy(pData->tempHopFrameTD[ch], HOP_SIZE, &(outputs[ch][t* HOP_SIZE]));
            for (; ch < nOutputs; ch++)
                memset(&(outputs[ch][t* HOP_SIZE]), 0, HOP_SIZE*sizeof(float));
        }
    }
    else 
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, FRAME_SIZE*sizeof(float));
}


/* Set Functions */

void panner_refreshSettings(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    pData->reInitGainTables = 1;
    pData->reInitTFT = 1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
}

void panner_checkReInit(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);

    /* reinitialise if needed */
    if (pData->reInitTFT==1) {
        pData->reInitTFT = 2;
        panner_initTFT(hPan);
        pData->reInitTFT = 0;
    }
    if (pData->reInitGainTables==1) {
        pData->reInitGainTables = 2;
        panner_initGainTables(hPan);
        pData->reInitGainTables = 0;
    }
}

void panner_setSourceAzi_deg(void* const hPan, int index, float newAzi_deg)
{
    panner_data *pData = (panner_data*)(hPan);
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = MAX(newAzi_deg, -180.0f);
    newAzi_deg = MIN(newAzi_deg, 180.0f);
    pData->src_dirs_deg[index][0] = newAzi_deg;
    pData->recalc_gainsFLAG[index] = 1;
    pData->recalc_M_rotFLAG = 1;
}

void panner_setSourceElev_deg(void* const hPan, int index, float newElev_deg)
{
    panner_data *pData = (panner_data*)(hPan);
    newElev_deg = MAX(newElev_deg, -90.0f);
    newElev_deg = MIN(newElev_deg, 90.0f);
    pData->src_dirs_deg[index][1] = newElev_deg;
    pData->recalc_gainsFLAG[index] = 1;
    pData->recalc_M_rotFLAG = 1;
}

void panner_setNumSources(void* const hPan, int new_nSources)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    /* determine if TFT must be reinitialised */
    pData->new_nSources = new_nSources > MAX_NUM_INPUTS ? MAX_NUM_INPUTS : new_nSources;
    if(pData->nSources != pData->new_nSources){
        pData->reInitTFT = 1;
        for(ch=pData->nSources; ch<pData->new_nSources; ch++)
            pData->recalc_gainsFLAG[ch] = 1;
        pData->recalc_M_rotFLAG = 1;
    }
}

void panner_setLoudspeakerAzi_deg(void* const hPan, int index, float newAzi_deg)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = MAX(newAzi_deg, -180.0f);
    newAzi_deg = MIN(newAzi_deg, 180.0f);
    pData->loudpkrs_dirs_deg[index][0] = newAzi_deg;
    pData->reInitGainTables=1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
    pData->recalc_M_rotFLAG = 1;
}

void panner_setLoudspeakerElev_deg(void* const hPan, int index, float newElev_deg)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    newElev_deg = MAX(newElev_deg, -90.0f);
    newElev_deg = MIN(newElev_deg, 90.0f);
    pData->loudpkrs_dirs_deg[index][1] = newElev_deg;
    pData->reInitGainTables=1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
    pData->recalc_M_rotFLAG = 1;
}

void panner_setNumLoudspeakers(void* const hPan, int new_nLoudspeakers)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    pData->new_nLoudpkrs = new_nLoudspeakers > MAX_NUM_OUTPUTS ? MAX_NUM_OUTPUTS : new_nLoudspeakers;
    if(pData->nLoudpkrs != pData->new_nLoudpkrs)
        pData->reInitTFT = 1; 
    pData->reInitGainTables=1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
    pData->recalc_M_rotFLAG = 1;
}

void panner_setOutputConfigPreset(void* const hPan, int newPresetID)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch, dummy;
    panner_loadPreset(newPresetID, pData->loudpkrs_dirs_deg, &(pData->new_nLoudpkrs), &dummy);
    if(pData->nLoudpkrs != pData->new_nLoudpkrs)
        pData->reInitTFT = 1;
    pData->reInitGainTables=1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
    pData->recalc_M_rotFLAG = 1;
}

void panner_setInputConfigPreset(void* const hPan, int newPresetID)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch, dummy;
    panner_loadPreset(newPresetID, pData->src_dirs_deg, &(pData->new_nSources), &dummy);
    if(pData->nSources != pData->new_nSources)
        pData->reInitTFT = 1;
    for(ch=0; ch<pData->new_nSources; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
    pData->recalc_M_rotFLAG = 1;
}

void panner_setDTT(void* const hPan, float newValue)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    pData->DTT = newValue;
    getPvalues(pData->DTT, pData->freqVector, HYBRID_BANDS, pData->pValue);
    for(ch=0; ch<pData->new_nSources; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
    pData->recalc_M_rotFLAG = 1;
}

void panner_setSpread(void* const hPan, float newValue)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    if(pData->spread_deg!=newValue){
        pData->spread_deg = CLAMP(newValue, PANNER_SPREAD_MIN_VALUE, PANNER_SPREAD_MAX_VALUE);
        pData->reInitGainTables=1;
        for(ch=0; ch<MAX_NUM_INPUTS; ch++)
            pData->recalc_gainsFLAG[ch] = 1;
        pData->recalc_M_rotFLAG = 1;
    }
}

void panner_setYaw(void  * const hBin, float newYaw)
{
    panner_data *pData = (panner_data*)(hBin);
    pData->yaw = pData->bFlipYaw == 1 ? -DEG2RAD(newYaw) : DEG2RAD(newYaw);
    pData->recalc_M_rotFLAG = 1;
}

void panner_setPitch(void* const hBin, float newPitch)
{
    panner_data *pData = (panner_data*)(hBin);
    pData->pitch = pData->bFlipPitch == 1 ? -DEG2RAD(newPitch) : DEG2RAD(newPitch);
    pData->recalc_M_rotFLAG = 1;
}

void panner_setRoll(void* const hBin, float newRoll)
{
    panner_data *pData = (panner_data*)(hBin);
    pData->roll = pData->bFlipRoll == 1 ? -DEG2RAD(newRoll) : DEG2RAD(newRoll);
    pData->recalc_M_rotFLAG = 1;
}

void panner_setFlipYaw(void* const hBin, int newState)
{
    panner_data *pData = (panner_data*)(hBin);
    if(newState !=pData->bFlipYaw ){
        pData->bFlipYaw = newState;
        panner_setYaw(hBin, -panner_getYaw(hBin));
    }
}

void panner_setFlipPitch(void* const hBin, int newState)
{
    panner_data *pData = (panner_data*)(hBin);
    if(newState !=pData->bFlipPitch ){
        pData->bFlipPitch = newState;
        panner_setPitch(hBin, -panner_getPitch(hBin));
    }
}

void panner_setFlipRoll(void* const hBin, int newState)
{
    panner_data *pData = (panner_data*)(hBin);
    if(newState !=pData->bFlipRoll ){
        pData->bFlipRoll = newState;
        panner_setRoll(hBin, -panner_getRoll(hBin));
    }
}


/* Get Functions */

float panner_getSourceAzi_deg(void* const hPan, int index)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->src_dirs_deg[index][0];
}

float panner_getSourceElev_deg(void* const hPan, int index)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->src_dirs_deg[index][1];
}

int panner_getNumSources(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->new_nSources;
}

int panner_getMaxNumSources()
{
    return MAX_NUM_INPUTS;
}

float panner_getLoudspeakerAzi_deg(void* const hPan, int index)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->loudpkrs_dirs_deg[index][0];
}

float panner_getLoudspeakerElev_deg(void* const hPan, int index)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->loudpkrs_dirs_deg[index][1];
}

int panner_getNumLoudspeakers(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->new_nLoudpkrs;
}

int panner_getMaxNumLoudspeakers()
{
    return MAX_NUM_OUTPUTS;
}

int panner_getDAWsamplerate(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->fs;
}

float panner_getDTT(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->DTT;
}

float panner_getSpread(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->spread_deg;
}

float panner_getYaw(void* const hBin)
{
    panner_data *pData = (panner_data*)(hBin);
    return pData->bFlipYaw == 1 ? -RAD2DEG(pData->yaw) : RAD2DEG(pData->yaw);
}

float panner_getPitch(void* const hBin)
{
    panner_data *pData = (panner_data*)(hBin);
    return pData->bFlipPitch == 1 ? -RAD2DEG(pData->pitch) : RAD2DEG(pData->pitch);
}

float panner_getRoll(void* const hBin)
{
    panner_data *pData = (panner_data*)(hBin);
    return pData->bFlipRoll == 1 ? -RAD2DEG(pData->roll) : RAD2DEG(pData->roll);
}

int panner_getFlipYaw(void* const hBin)
{
    panner_data *pData = (panner_data*)(hBin);
    return pData->bFlipYaw;
}

int panner_getFlipPitch(void* const hBin)
{
    panner_data *pData = (panner_data*)(hBin);
    return pData->bFlipPitch;
}

int panner_getFlipRoll(void* const hBin)
{
    panner_data *pData = (panner_data*)(hBin);
    return pData->bFlipRoll;
}

int panner_getProcessingDelay()
{
    return 12*HOP_SIZE;
}





    
    
