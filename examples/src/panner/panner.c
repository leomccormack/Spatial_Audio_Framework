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
 * @file panner.c
 * @brief A frequency-dependent 3D panner based on the Vector-base Amplitude
 *        Panning (VBAP) method [1], with an optional spread control [2].
 *
 * Depending on the listening room, it may be beneficial to employ amplitude-
 * normalised gains for low frequencies, and energy-normalised gains for high
 * frequencies. Therefore, this VBAP implementation also uses the method
 * described in [3], to do just that.
 *
 * @see [1] Pulkki, V. (1997). Virtual sound source positioning using vector
 *          base amplitude panning. Journal of the audio engineering society,
 *          45(6), 456-466.
 * @see [2] Pulkki, V. (1999). Uniform spreading of amplitude panned virtual
 *          sources. In Proceedings of the 1999 IEEE Workshop on Applications of
 *          Signal Processing to Audio and Acoustics. WASPAA'99 (Cat. No.
 *          99TH8452) (pp. 187-190). IEEE.
 * @see [3] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V.
 *          (2014). Gain normalisation in amplitude panning as a function of
 *          frequency and room reverberance. 55th International Conference of
 *          the AES. Helsinki, Finland.
 *
 * @author Leo McCormack
 * @date 25.09.2017
 * @license ISC
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

    /* default user parameters */
    panner_loadSourcePreset(SOURCE_CONFIG_PRESET_DEFAULT, pData->src_dirs_deg, &(pData->new_nSources), &(dummy)); /*check setStateInformation if you change default preset*/
    pData->nSources = pData->new_nSources;
    pData->DTT = 0.5f;
    pData->spread_deg = 0.0f;
    panner_loadLoudspeakerPreset(LOUDSPEAKER_ARRAY_PRESET_STEREO, pData->loudpkrs_dirs_deg, &(pData->new_nLoudpkrs), &(pData->output_nDims)); /*check setStateInformation if you change default preset*/
    pData->nLoudpkrs = pData->new_nLoudpkrs;
    pData->yaw = 0.0f;
    pData->pitch = 0.0f;
    pData->roll = 0.0f;
    pData->bFlipYaw = 0;
    pData->bFlipPitch = 0;
    pData->bFlipRoll = 0;
    
    /* time-frequency transform + buffers */
    pData->hSTFT = NULL;
    pData->inputFrameTD = (float**)malloc2d(MAX_NUM_INPUTS, PANNER_FRAME_SIZE, sizeof(float));
    pData->outputFrameTD = (float**)malloc2d(MAX_NUM_OUTPUTS, PANNER_FRAME_SIZE, sizeof(float));
    pData->inputframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_INPUTS, TIME_SLOTS, sizeof(float_complex));
    pData->outputframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_OUTPUTS, TIME_SLOTS, sizeof(float_complex));

    /* flags and gain table */
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
    pData->vbap_gtable = NULL;
    pData->recalc_M_rotFLAG = 1;
    pData->reInitGainTables = 1;
}

void panner_destroy
(
    void ** const phPan
)
{
    panner_data *pData = (panner_data*)(*phPan);

    if (pData != NULL) {
        /* not safe to free memory during intialisation/processing loop */
        while (pData->codecStatus == CODEC_STATUS_INITIALISING ||
               pData->procStatus == PROC_STATUS_ONGOING){
            SAF_SLEEP(10);
        }
        
        /* free afSTFT and buffers */
        if(pData->hSTFT !=NULL)
            afSTFT_destroy(&(pData->hSTFT));
        free(pData->inputFrameTD);
        free(pData->outputFrameTD);
        free(pData->inputframeTF);
        free(pData->outputframeTF);
        free(pData->vbap_gtable);
        free(pData->progressBarText);
        
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
    
    /* define frequency vector */
    pData->fs = sampleRate;
    afSTFT_getCentreFreqs(pData->hSTFT, (float)sampleRate, HYBRID_BANDS, pData->freqVector);
    
    /* calculate pValue per frequency */
    getPvalues(pData->DTT, pData->freqVector, HYBRID_BANDS, pData->pValue);

    /* reinitialise if needed */
    pData->recalc_M_rotFLAG = 1;
}

void panner_initCodec
(
    void* const hPan
)
{
    panner_data *pData = (panner_data*)(hPan);
    
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
    
    /* reinit TFT if needed */
    panner_initTFT(hPan);
    
    /* reinit gain tables */
    if(pData->reInitGainTables){
        panner_initGainTables(hPan);
        pData->reInitGainTables = 0;
    }
    
    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
    
}

void panner_process
(
    void        *  const hPan,
    const float *const * inputs,
    float       ** const outputs,
    int                  nInputs,
    int                  nOutputs,
    int                  nSamples
)
{
    panner_data *pData = (panner_data*)(hPan);
    int t, ch, ls, i, band, nSources, nLoudspeakers, N_azi, aziIndex, elevIndex, idx3d, idx2D;
    float aziRes, elevRes, pv_f, gains3D_sum_pvf, gains2D_sum_pvf, Rxyz[3][3], hypotxy;
    float src_dirs[MAX_NUM_INPUTS][2], pValue[HYBRID_BANDS], gains3D[MAX_NUM_OUTPUTS], gains2D[MAX_NUM_OUTPUTS];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);
    float_complex outputTemp[MAX_NUM_OUTPUTS][TIME_SLOTS];

    /* copy user parameters to local variables */
    memcpy(src_dirs, pData->src_dirs_deg, MAX_NUM_INPUTS*2*sizeof(float));
    memcpy(pValue, pData->pValue, HYBRID_BANDS*sizeof(float));
    nSources = pData->nSources;
    nLoudspeakers = pData->nLoudpkrs;

    /* apply panner */
    if ((nSamples == PANNER_FRAME_SIZE) && (pData->vbap_gtable != NULL) && (pData->codecStatus == CODEC_STATUS_INITIALISED) ) {
        pData->procStatus = PROC_STATUS_ONGOING;

        /* Load time-domain data */
        for(i=0; i < SAF_MIN(nSources,nInputs); i++)
            utility_svvcopy(inputs[i], PANNER_FRAME_SIZE, pData->inputFrameTD[i]);
        for(; i<MAX_NUM_INPUTS; i++)
            memset(pData->inputFrameTD[i], 0, PANNER_FRAME_SIZE * sizeof(float));

        /* Apply time-frequency transform (TFT) */
        afSTFT_forward_knownDimensions(pData->hSTFT, pData->inputFrameTD, PANNER_FRAME_SIZE, MAX_NUM_INPUTS, TIME_SLOTS, pData->inputframeTF);
        memset(FLATTEN3D(pData->outputframeTF), 0, HYBRID_BANDS*MAX_NUM_OUTPUTS*TIME_SLOTS * sizeof(float_complex));
        memset(outputTemp, 0, MAX_NUM_OUTPUTS*TIME_SLOTS * sizeof(float_complex));

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
                                gains3D_sum_pvf += powf(SAF_MAX(gains3D[ls], 0.0f), pv_f);
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
                    FLATTEN2D(pData->inputframeTF[band]), TIME_SLOTS, &cbeta,
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
                                gains2D_sum_pvf += powf(SAF_MAX(gains2D[ls], 0.0f), pv_f);
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
            cblas_sscal(/*re+im*/2*nLoudspeakers*TIME_SLOTS, 1.0f/sqrtf((float)nSources), (float*)FLATTEN2D(pData->outputframeTF[band]), 1);

        /* inverse-TFT and copy to output */
        afSTFT_backward_knownDimensions(pData->hSTFT, pData->outputframeTF, PANNER_FRAME_SIZE, MAX_NUM_OUTPUTS, TIME_SLOTS, pData->outputFrameTD);
        for (ch = 0; ch < SAF_MIN(nLoudspeakers, nOutputs); ch++)
            utility_svvcopy(pData->outputFrameTD[ch], PANNER_FRAME_SIZE, outputs[ch]);
        for (; ch < nOutputs; ch++)
            memset(outputs[ch], 0, PANNER_FRAME_SIZE*sizeof(float));

    }
    else
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, PANNER_FRAME_SIZE*sizeof(float));


    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}


/* Set Functions */

void panner_refreshSettings(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    pData->reInitGainTables = 1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
    panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
}

void panner_setSourceAzi_deg(void* const hPan, int index, float newAzi_deg)
{
    panner_data *pData = (panner_data*)(hPan);
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = SAF_MAX(newAzi_deg, -180.0f);
    newAzi_deg = SAF_MIN(newAzi_deg, 180.0f);
    if(pData->src_dirs_deg[index][0] != newAzi_deg){
        pData->src_dirs_deg[index][0] = newAzi_deg;
        pData->recalc_gainsFLAG[index] = 1;
        pData->recalc_M_rotFLAG = 1;
    }
}

void panner_setSourceElev_deg(void* const hPan, int index, float newElev_deg)
{
    panner_data *pData = (panner_data*)(hPan);
    newElev_deg = SAF_MAX(newElev_deg, -90.0f);
    newElev_deg = SAF_MIN(newElev_deg, 90.0f);
    if(pData->src_dirs_deg[index][1] != newElev_deg){
        pData->src_dirs_deg[index][1] = newElev_deg;
        pData->recalc_gainsFLAG[index] = 1;
        pData->recalc_M_rotFLAG = 1;
    }
}

void panner_setNumSources(void* const hPan, int new_nSources)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    /* determine if TFT must be reinitialised */
    new_nSources = new_nSources > MAX_NUM_INPUTS ? MAX_NUM_INPUTS : new_nSources;
    if(pData->nSources != new_nSources){
        pData->new_nSources = new_nSources;
        for(ch=pData->nSources; ch<pData->new_nSources; ch++)
            pData->recalc_gainsFLAG[ch] = 1;
        pData->recalc_M_rotFLAG = 1;
        panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
    }
}

void panner_setLoudspeakerAzi_deg(void* const hPan, int index, float newAzi_deg)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = SAF_MAX(newAzi_deg, -180.0f);
    newAzi_deg = SAF_MIN(newAzi_deg, 180.0f);
    if(pData->loudpkrs_dirs_deg[index][0] != newAzi_deg){
        pData->loudpkrs_dirs_deg[index][0] = newAzi_deg;
        pData->reInitGainTables=1;
        for(ch=0; ch<MAX_NUM_INPUTS; ch++)
            pData->recalc_gainsFLAG[ch] = 1;
        pData->recalc_M_rotFLAG = 1;
        panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
    }
}

void panner_setLoudspeakerElev_deg(void* const hPan, int index, float newElev_deg)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    newElev_deg = SAF_MAX(newElev_deg, -90.0f);
    newElev_deg = SAF_MIN(newElev_deg, 90.0f);
    if(pData->loudpkrs_dirs_deg[index][1] != newElev_deg){
        pData->loudpkrs_dirs_deg[index][1] = newElev_deg;
        pData->reInitGainTables=1;
        for(ch=0; ch<MAX_NUM_INPUTS; ch++)
            pData->recalc_gainsFLAG[ch] = 1;
        pData->recalc_M_rotFLAG = 1;
        panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
    }
}

void panner_setNumLoudspeakers(void* const hPan, int new_nLoudspeakers)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    
    new_nLoudspeakers  = new_nLoudspeakers > MAX_NUM_OUTPUTS ? MAX_NUM_OUTPUTS : new_nLoudspeakers;
    if(pData->new_nLoudpkrs != new_nLoudspeakers){
        pData->new_nLoudpkrs = new_nLoudspeakers;
        pData->reInitGainTables=1;
        for(ch=0; ch<MAX_NUM_INPUTS; ch++)
            pData->recalc_gainsFLAG[ch] = 1;
        pData->recalc_M_rotFLAG = 1;
        panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
    }
}

void panner_setOutputConfigPreset(void* const hPan, int newPresetID)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch, dummy;
    panner_loadLoudspeakerPreset(newPresetID, pData->loudpkrs_dirs_deg, &(pData->new_nLoudpkrs), &dummy);
    pData->reInitGainTables=1;
    for(ch=0; ch<MAX_NUM_INPUTS; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
    pData->recalc_M_rotFLAG = 1;
    panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
}

void panner_setInputConfigPreset(void* const hPan, int newPresetID)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch, dummy;
    panner_loadSourcePreset(newPresetID, pData->src_dirs_deg, &(pData->new_nSources), &dummy);
    for(ch=0; ch<pData->new_nSources; ch++)
        pData->recalc_gainsFLAG[ch] = 1;
    pData->recalc_M_rotFLAG = 1;
    panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
}

void panner_setDTT(void* const hPan, float newValue)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    if(pData->DTT != newValue){
        pData->DTT = newValue;
        getPvalues(pData->DTT, pData->freqVector, HYBRID_BANDS, pData->pValue);
        for(ch=0; ch<pData->new_nSources; ch++)
            pData->recalc_gainsFLAG[ch] = 1;
        pData->recalc_M_rotFLAG = 1;
        panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
    }
}

void panner_setSpread(void* const hPan, float newValue)
{
    panner_data *pData = (panner_data*)(hPan);
    int ch;
    if(pData->spread_deg!=newValue){
        pData->spread_deg = SAF_CLAMP(newValue, PANNER_SPREAD_MIN_VALUE, PANNER_SPREAD_MAX_VALUE);
        pData->reInitGainTables=1;
        for(ch=0; ch<MAX_NUM_INPUTS; ch++)
            pData->recalc_gainsFLAG[ch] = 1;
        pData->recalc_M_rotFLAG = 1;
        panner_setCodecStatus(hPan, CODEC_STATUS_NOT_INITIALISED);
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

int panner_getFrameSize(void)
{
    return PANNER_FRAME_SIZE;
}

CODEC_STATUS panner_getCodecStatus(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->codecStatus;
}

float panner_getProgressBar0_1(void* const hPan)
{
    panner_data *pData = (panner_data*)(hPan);
    return pData->progressBar0_1;
}

void panner_getProgressBarText(void* const hPan, char* text)
{
    panner_data *pData = (panner_data*)(hPan);
    memcpy(text, pData->progressBarText, PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

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
