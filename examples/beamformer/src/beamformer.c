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

/*
 * Filename: beamformer.c
 * ----------------------
 * Generates beamformers/virtual microphones in arbitrary directions. Several
 * different beam pattern types are included.
 *
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
    beamformer_data* pData = (beamformer_data*)malloc1d(sizeof(beamformer_data));
    *phBeam = (void*)pData;
    int i, j, ch, band;

    /* default user parameters */
    pData->beamOrder = 1;
    pData->nSH = pData->new_nSH = (pData->beamOrder + 1)*(pData->beamOrder + 1);
    for(i=0; i<MAX_NUM_BEAMS; i++){
        pData->beam_dirs_deg[i][0] = default_LScoords64_rad[i][0]*180.0f/M_PI;
        pData->beam_dirs_deg[i][1] = (default_LScoords64_rad[i][1] - M_PI/2.0f) < -M_PI/2.0f ?
        (M_PI/2.0f + default_LScoords64_rad[i][1]) :  (default_LScoords64_rad[i][1] - M_PI/2.0f);
        pData->beam_dirs_deg[i][1] *= 180.0f/M_PI;
    }
    pData->nBeams = pData->new_nBeams = 1;
    pData->beamType = BEAM_TYPE_HYPERCARDIOID;
    pData->chOrdering = CH_ACN;
    pData->norm = NORM_SN3D;
    
    /* afSTFT stuff */
//    pData->hSTFT = NULL;
//    pData->STFTInputFrameTF = malloc1d(MAX_NUM_SH_SIGNALS * sizeof(complexVector));
//    for(ch=0; ch< MAX_NUM_SH_SIGNALS; ch++) {
//        pData->STFTInputFrameTF[ch].re = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
//        pData->STFTInputFrameTF[ch].im = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
//    }
//    pData->tempHopFrameTD = (float**)malloc2d( MAX(MAX_NUM_SH_SIGNALS, MAX_NUM_LOUDSPEAKERS), HOP_SIZE, sizeof(float));
//    pData->STFTOutputFrameTF = malloc1d(MAX_NUM_LOUDSPEAKERS * sizeof(complexVector));
//    for(ch=0; ch< MAX_NUM_LOUDSPEAKERS; ch++) {
//        pData->STFTOutputFrameTF[ch].re = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
//        pData->STFTOutputFrameTF[ch].im = (float*)calloc1d(HYBRID_BANDS, sizeof(float));
//    }
    
    /* internal parameters */
    pData->new_nSH = pData->new_nSH = (pData->beamOrder+1)*(pData->beamOrder+1);
    
    /* flags */
    pData->reInitTFT = 1;
    for(ch=0; ch<MAX_NUM_BEAMS; ch++)
        pData->recalc_beamWeights[ch] = 1;
}

void beamformer_destroy
(
    void ** const phBeam
)
{
    beamformer_data *pData = (beamformer_data*)(*phBeam);
    int i, j, ch;
    
    if (pData != NULL) {
//        if(pData->hSTFT!=NULL)
//            afSTFTfree(pData->hSTFT);
//        if(pData->STFTInputFrameTF!=NULL){
//            for (ch = 0; ch< MAX_NUM_SH_SIGNALS; ch++) {
//                free(pData->STFTInputFrameTF[ch].re);
//                free(pData->STFTInputFrameTF[ch].im);
//            }
//        }
//        if((pData->STFTOutputFrameTF!=NULL)){
//            for (ch = 0; ch < MAX_NUM_LOUDSPEAKERS; ch++) {
//                free(pData->STFTOutputFrameTF[ch].re);
//                free(pData->STFTOutputFrameTF[ch].im);
//            }
//        }
//        free(pData->STFTInputFrameTF);
//        free(pData->STFTOutputFrameTF);
//        free2d((void**)pData->tempHopFrameTD, MAX(NUM_EARS, MAX_NUM_SH_SIGNALS));
//        free2d((void**)pData->tempHopFrameTD, MAX(MAX_NUM_LOUDSPEAKERS, MAX_NUM_SH_SIGNALS));

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
    int band, ch;
    
    /* define frequency vector */
    pData->fs = sampleRate;
//    for(band=0; band <HYBRID_BANDS; band++){
//        if(sampleRate == 44100)
//            pData->freqVector[band] =  (float)__afCenterFreq44100[band];
//        else /* Assume 48kHz */
//            pData->freqVector[band] =  (float)__afCenterFreq48e3[band];
//    }

    /* reinitialise if needed */
    beamformer_checkReInit(hBeam);
    
    /* defaults */
    memset(pData->beamWeights, 0, MAX_NUM_BEAMS*MAX_NUM_SH_SIGNALS*sizeof(float));
    //memset(pData->beamWeights_cmplx, 0, MAX_NUM_BEAMS*MAX_NUM_SH_SIGNALS*sizeof(float_complex));
    for(ch=0; ch<MAX_NUM_BEAMS; ch++)
        pData->recalc_beamWeights[ch] = 1;
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
    int n, t, ch, i, j, bi, band;
    int o[MAX_SH_ORDER+2];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    /* local copies of user parameters */
    int nBeams, beamOrder;
    NORM_TYPES norm;
    CH_ORDER chOrdering;
    
    /* reinitialise if needed */
#ifdef __APPLE__
    beamformer_checkReInit(hBeam);
#else
    if(pData->reInitTFT==1){
        pData->reInitTFT = 2;
        beamformer_initTFT(hBeam); /* always init before codec or hrtfs  */
        pData->reInitTFT = 0;
    }
#endif

    /* decode audio to loudspeakers or headphones */
    if( (nSamples == FRAME_SIZE) && (pData->reInitTFT==0) ) {
        /* copy user parameters to local variables */
        for(n=0; n<MAX_SH_ORDER+2; n++){  o[n] = n*n;  }
        beamOrder = pData->beamOrder;
        nBeams = pData->nBeams;
        norm = pData->norm;
        chOrdering = pData->chOrdering;
        
        /* Load time-domain data */
        switch(chOrdering){
            case CH_ACN:
                for(i=0; i < MIN(pData->nSH, nInputs); i++)
                    utility_svvcopy(inputs[i], FRAME_SIZE, pData->SHFrameTD[i]);
                for(; i<pData->nSH; i++)
                    memset(pData->SHFrameTD[i], 0, FRAME_SIZE * sizeof(float)); /* fill remaining channels with zeros */
                break;
            case CH_FUMA:   /* only for first-order, convert to ACN */
                if(nInputs>=4){
                    utility_svvcopy(inputs[0], FRAME_SIZE, pData->SHFrameTD[0]);
                    utility_svvcopy(inputs[1], FRAME_SIZE, pData->SHFrameTD[3]);
                    utility_svvcopy(inputs[2], FRAME_SIZE, pData->SHFrameTD[1]);
                    utility_svvcopy(inputs[3], FRAME_SIZE, pData->SHFrameTD[2]);
                    for(i=4; i<pData->nSH; i++)
                        memset(pData->SHFrameTD[i], 0, FRAME_SIZE * sizeof(float)); /* fill remaining channels with zeros */
                }
                else
                    for(i=0; i<pData->nSH; i++)
                        memset(pData->SHFrameTD[i], 0, FRAME_SIZE * sizeof(float));
                break;
        }
        
        /* account for input normalisation scheme */
        switch(norm){
            case NORM_N3D:  /* already in N3D, do nothing */
                break;
            case NORM_SN3D: /* convert to N3D */
                for (n = 0; n<beamOrder+1; n++)
                    for (ch = o[n]; ch<o[n+1]; ch++)
                        for(i = 0; i<FRAME_SIZE; i++)
                            pData->SHFrameTD[ch][i] *= sqrtf(2.0f*(float)n+1.0f);
                break;
            case NORM_FUMA: /* only for first-order, convert to N3D */
                for(i = 0; i<FRAME_SIZE; i++)
                    pData->SHFrameTD[0][i] *= sqrtf(2.0f);
                for (ch = 1; ch<4; ch++)
                    for(i = 0; i<FRAME_SIZE; i++)
                        pData->SHFrameTD[ch][i] *= sqrtf(3.0f);
                break;
        }
        
//        /* Apply time-frequency transform (TFT) */
//        for(t=0; t< TIME_SLOTS; t++) {
//            for(ch = 0; ch < pData->nSH; ch++)
//                utility_svvcopy(&(pData->SHFrameTD[ch][t*HOP_SIZE]), HOP_SIZE, pData->tempHopFrameTD[ch]);
//            afSTFTforward(pData->hSTFT, (float**)pData->tempHopFrameTD, (complexVector*)pData->STFTInputFrameTF);
//            for(band=0; band<HYBRID_BANDS; band++)
//                for(ch=0; ch < pData->nSH; ch++)
//                    pData->SHframeTF[band][ch][t] = cmplxf(pData->STFTInputFrameTF[ch].re[band], pData->STFTInputFrameTF[ch].im[band]);
//        }
        
        /* Main processing: */
        if(isPlaying){
            float* c_n;
            c_n = malloc1d((beamOrder+1)*sizeof(float));
            
            /* calculate beamforming coeffients */
            for(bi=0; bi<nBeams; bi++){
                if(pData->recalc_beamWeights[bi]){
                    memset(pData->beamWeights[bi], 0, MAX_NUM_SH_SIGNALS*sizeof(float));
                    //memset(pData->beamWeights_cmplx[bi], 0, MAX_NUM_SH_SIGNALS*sizeof(float_complex));
                    switch(pData->beamType){
                        case BEAM_TYPE_CARDIOID: beamWeightsCardioid2Spherical(beamOrder, c_n); break;
                        case BEAM_TYPE_HYPERCARDIOID: beamWeightsHypercardioid2Spherical(beamOrder, c_n); break;
                        case BEAM_TYPE_MAX_EV: beamWeightsMaxEV(beamOrder, c_n); break;
                    }
                    rotateAxisCoeffsReal(beamOrder, c_n, M_PI/2.0f - pData->beam_dirs_deg[bi][1]*M_PI/180.0f,
                                         pData->beam_dirs_deg[bi][0]*M_PI/180.0f, (float*)pData->beamWeights[bi]);
//                    for(i=0; i<pData->new_nSH; i++)
//                        pData->beamWeights_cmplx[bi][i] = cmplxf(pData->beamWeights[bi][i], 0.0f);
                    
                     pData->recalc_beamWeights[bi] = 0;
                }
            }
            free(c_n);
            
            /* apply beam weights (frequency-domain, for future adaptive beamformers) */
//            for(bi=0; bi<nBeams; bi++){
//                cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nBeams, TIME_SLOTS, pData->nSH, &calpha,
//                            pData->beamWeights_cmplx, MAX_NUM_SH_SIGNALS,
//                            pData->SHframeTF[band], TIME_SLOTS, &cbeta,
//                            pData->outputframeTF[band], TIME_SLOTS);
//            }
            
            /* apply beam weights */
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nBeams, FRAME_SIZE, pData->nSH, 1.0f,
                        (const float*)pData->beamWeights, MAX_NUM_SH_SIGNALS,
                        (const float*)pData->SHFrameTD, FRAME_SIZE, 0.0f,
                        (float*)pData->outputFrameTD, FRAME_SIZE);
 
            
        }
        else{
            //memset(pData->outputframeTF, 0, HYBRID_BANDS*MAX_NUM_BEAMS*TIME_SLOTS*sizeof(float_complex));
            memset(pData->outputFrameTD, 0, MAX_NUM_BEAMS*FRAME_SIZE*sizeof(float));
        }
        
//        /* inverse-TFT */
//        for(t = 0; t < TIME_SLOTS; t++) {
//            for(band = 0; band < HYBRID_BANDS; band++) {
//                for(ch = 0; ch < nLoudspeakers; ch++) {
//                    pData->STFTOutputFrameTF[ch].re[band] = crealf(pData->outputframeTF[band][ch][t]);
//                    pData->STFTOutputFrameTF[ch].im[band] = cimagf(pData->outputframeTF[band][ch][t]);
//                }
//            }
//
//            afSTFTinverse(pData->hSTFT, pData->STFTOutputFrameTF, pData->tempHopFrameTD);
//            for(ch = 0; ch < MIN(nBeams, nOutputs); ch++)
//                utility_svvcopy(pData->tempHopFrameTD[ch], HOP_SIZE, &(outputs[ch][t* HOP_SIZE]));
//            for (; ch < nOutputs; ch++)
//                memset(&(outputs[ch][t* HOP_SIZE]), 0, HOP_SIZE*sizeof(float));
//        }
        
        for(ch = 0; ch < MIN(nBeams, nOutputs); ch++)
            utility_svvcopy(pData->outputFrameTD[ch], FRAME_SIZE, outputs[ch]);
        for (; ch < nOutputs; ch++)
            memset(outputs[ch], 0, FRAME_SIZE*sizeof(float));
    }
    else
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch], 0, FRAME_SIZE*sizeof(float));
}


/* Set Functions */

void beamformer_refreshSettings(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    pData->reInitTFT = 1;
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
}

void beamformer_setBeamOrder(void  * const hBeam, int newValue)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int ch;
    pData->beamOrder = MIN(MAX(newValue,1), MAX_SH_ORDER);
    pData->new_nSH = (pData->beamOrder+1)*(pData->beamOrder+1);
    pData->reInitTFT = 1;
    for(ch=0; ch<MAX_NUM_BEAMS; ch++)
        pData->recalc_beamWeights[ch] = 1;
    /* FUMA only supports 1st order */
    if(pData->beamOrder!=BEAM_ORDER_FIRST && pData->chOrdering == CH_FUMA)
        pData->chOrdering = CH_ACN;
    if(pData->beamOrder!=BEAM_ORDER_FIRST && pData->norm == NORM_FUMA)
        pData->norm = NORM_SN3D;
}

void beamformer_setBeamAzi_deg(void* const hBeam, int index, float newAzi_deg)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = MAX(newAzi_deg, -180.0f);
    newAzi_deg = MIN(newAzi_deg, 180.0f);
    pData->beam_dirs_deg[index][0] = newAzi_deg;
    pData->recalc_beamWeights[index] = 1;
}

void beamformer_setBeamElev_deg(void* const hBeam, int index, float newElev_deg)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    newElev_deg = MAX(newElev_deg, -90.0f);
    newElev_deg = MIN(newElev_deg, 90.0f);
    pData->beam_dirs_deg[index][1] = newElev_deg;
    pData->recalc_beamWeights[index] = 1;
}

void beamformer_setNumBeams(void* const hBeam, int new_nBeams)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int ch;
    pData->new_nBeams = CLAMP(new_nBeams, 1, MAX_NUM_BEAMS);
    if(pData->nBeams != pData->new_nBeams){
        pData->reInitTFT = 1;
        for(ch=0; ch<MAX_NUM_BEAMS; ch++)
            pData->recalc_beamWeights[ch] = 1;
    }
}

void beamformer_setChOrder(void* const hBeam, int newOrder)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    if((CH_ORDER)newOrder != CH_FUMA || pData->beamOrder==BEAM_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->chOrdering = (CH_ORDER)newOrder;
}

void beamformer_setNormType(void* const hBeam, int newType)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    if((NORM_TYPES)newType != NORM_FUMA || pData->beamOrder==BEAM_ORDER_FIRST)/* FUMA only supports 1st order */
        pData->norm = (NORM_TYPES)newType;
}

void beamformer_setBeamType(void* const hBeam, int newID)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    int ch;
    pData->beamType = newID;
    for(ch=0; ch<MAX_NUM_BEAMS; ch++)
        pData->recalc_beamWeights[ch] = 1;
}

/* Get Functions */

int beamformer_getBeamOrder(void  * const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->beamOrder;
}

int beamformer_getNumberOfBands(void)
{
    return HYBRID_BANDS;
}

float beamformer_getBeamAzi_deg(void* const hBeam, int index)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->beam_dirs_deg[index][0];
}

float beamformer_getBeamElev_deg(void* const hBeam, int index)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->beam_dirs_deg[index][1];
}

int beamformer_getNumBeams(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->new_nBeams;
}

int beamformer_getMaxNumBeams()
{
    return MAX_NUM_BEAMS;
}

int  beamformer_getNSHrequired(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->nSH;
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

int beamformer_getBeamType(void* const hBeam)
{
    beamformer_data *pData = (beamformer_data*)(hBeam);
    return pData->beamType;
}



