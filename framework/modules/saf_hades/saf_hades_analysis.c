/*
 * This file is part of the saf_hades module.
 * Copyright (c) 2021 - Leo McCormack & Janani Fernandez
 *
 * The saf_hades module is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * The saf_hades module is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * See <http://www.gnu.org/licenses/> for a copy of the GNU General Public
 * License.
 */

/**
 * @file saf_hades_analysis.c
 * @ingroup HADES
 * @brief Source file for the HADES analysis (#SAF_HADES_MODULE)
 *
 * The framework for binaural rendering of Hearing-Assistive/Augmented-reality
 * Devices (HADES) is described further in [1].
 *
 * @see [1] paper submitted for review.
 *
 * @author Leo McCormack and Janani Fernandez
 * @date 01.02.2021
 * @license GNU GPLv2
 */

#include "saf_hades_analysis.h"
#include "saf_hades_internal.h"

#ifdef  SAF_ENABLE_HADES_MODULE

/* ========================================================================== */
/*                             HADES Analysis                               */
/* ========================================================================== */

void hades_analysis_create
(
    hades_analysis_handle* const phAna,
    float fs,
    HADES_FILTERBANKS fbOption,
    int hopsize,
    int blocksize,
    int hybridmode,
    float* h_array,
    float* grid_dirs_deg,
    int nGrid,
    int nMics,
    int h_len,
    HADES_DIFFUSENESS_ESTIMATORS diffOption,
    HADES_DOA_ESTIMATORS doaOption
)
{
    hades_analysis_data* a = (hades_analysis_data*)malloc1d(sizeof(hades_analysis_data));
    *phAna = (void*)a;
    int band, i, idx_max;
    float* w_tmp;
    float_complex *U, *E, *H_W;
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta = cmplxf(0.0f, 0.0f); /* blas */

    nMics = SAF_MIN(nMics, HADES_MAX_NMICS);
    assert(blocksize % hopsize == 0); /* Must be a multiple of hopsize */
    assert(blocksize<=HADES_MAX_BLOCKSIZE);

    /* User parameters */ 
    a->fs = fs;
    a->fbOpt = fbOption;
    a->hopsize = hopsize;
    a->blocksize = blocksize;
    a->hybridmode = hybridmode;
    a->h_array = malloc1d(nGrid*nMics*h_len*sizeof(float));
    memcpy(a->h_array, h_array, nGrid*nMics*h_len*sizeof(float));
    a->grid_dirs_deg = malloc1d(nGrid*2*sizeof(float));
    memcpy(a->grid_dirs_deg, grid_dirs_deg, nGrid*2*sizeof(float));
    a->nGrid = nGrid;
    a->nMics = nMics;
    a->h_len = h_len;
    a->diffOpt = diffOption;
    a->doaOpt = doaOption;
    a->covAvgCoeff = 1.0f - 1.0f/(4096.0f/a->blocksize);
    a->covAvgCoeff = SAF_CLAMP(a->covAvgCoeff, 0.0f, 0.99999f);

    /* Scale steering vectors so that the peak of loudest measurement is 1 */
    utility_simaxv(a->h_array, nGrid*nMics*h_len, &idx_max);
    cblas_sscal(nGrid*nMics*h_len, 1.0f/a->h_array[idx_max], a->h_array, 1);

    /* Initialise time-frequency transform  */
    a->timeSlots = a->blocksize/a->hopsize;
    switch(a->fbOpt){
        case HADES_USE_AFSTFT_LD: /* fall through */
        case HADES_USE_AFSTFT:
            afSTFT_create(&(a->hFB_enc), a->nMics, 0, a->hopsize, a->fbOpt==HADES_USE_AFSTFT_LD ? 1 : 0, a->hybridmode, AFSTFT_BANDS_CH_TIME);
            a->nBands = afSTFT_getNBands(a->hFB_enc);
            a->freqVector = malloc1d(a->nBands*sizeof(float));
            a->filterbankDelay = afSTFT_getProcDelay(a->hFB_enc);
            afSTFT_getCentreFreqs(a->hFB_enc, a->fs, a->nBands, a->freqVector);
            a->H_array = malloc1d(a->nBands*(a->nMics)*(a->nGrid)*sizeof(float_complex));
            a->H_array_w = malloc1d(a->nBands*(a->nMics)*(a->nGrid)*sizeof(float_complex));
            afSTFT_FIRtoFilterbankCoeffs(a->h_array, a->nGrid, a->nMics, a->h_len, a->hopsize, a->fbOpt==HADES_USE_AFSTFT_LD ? 1 : 0, a->hybridmode, a->H_array);
            break;
    }

    /* Initialise DoA estimator */
    utility_cseig_create(&(a->hEig), a->nMics);
    a->grid_dirs_xyz = malloc1d(a->nGrid*3*sizeof(float));
    unitSph2cart(a->grid_dirs_deg, a->nGrid, 1, a->grid_dirs_xyz);
    switch(a->doaOpt){
        case HADES_USE_MUSIC: hades_sdMUSIC_create(&(a->hDoA), a->nMics, a->grid_dirs_deg, a->nGrid); break;
    }

    /* Integration weights */
    a->W = calloc1d(a->nGrid*a->nGrid,sizeof(float_complex));
    if (cblas_sasum(a->nGrid, a->grid_dirs_deg+1, 2)/(float)a->nGrid<0.0001){
        for(i=0; i<a->nGrid; i++)
            a->W[i*(a->nGrid)+i] = cmplxf(1.0f, 0.0f);
    }
    else{
        w_tmp = malloc1d(a->nGrid*sizeof(float));
        getVoronoiWeights(a->grid_dirs_deg, a->nGrid, 0, w_tmp);
        for(i=0; i<a->nGrid; i++)
            a->W[i*(a->nGrid)+i] = cmplxf(w_tmp[i], 0.0f);
    }

    /* For spatial whitening of the spatial covariance matrix, such that it has an identity structure under diffuse-field conditions */
    a->T = (float_complex**)malloc2d(a->nBands, a->nMics*(a->nMics), sizeof(float_complex));
    a->DCM_array = malloc1d(a->nBands*(a->nMics)*(a->nMics)*sizeof(float_complex));
    U = malloc1d(a->nMics*(a->nMics)*sizeof(float_complex));
    E = malloc1d(a->nMics*(a->nMics)*sizeof(float_complex));
    H_W = malloc1d(a->nMics*(a->nGrid)*sizeof(float_complex));
    for(band=0; band<a->nBands; band++){
        /* Diffuse covariance matrix */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a->nMics, a->nGrid, a->nGrid, &calpha,
                    &(a->H_array[band*(a->nMics)*(a->nGrid)]), a->nGrid,
                    a->W, a->nGrid, &cbeta,
                    H_W, a->nGrid);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, a->nMics, a->nMics, a->nGrid, &calpha,
                    H_W, a->nGrid,
                    &(a->H_array[band*(a->nMics)*(a->nGrid)]), a->nGrid, &cbeta,
                    &(a->DCM_array[band*(a->nMics)*(a->nMics)]), a->nMics);
        cblas_sscal(/*re+im*/2*(a->nMics)*(a->nMics), 1.0f/(float)a->nGrid, (float*)&(a->DCM_array[band*(a->nMics)*(a->nMics)]), 1);

        /* Decomposition of the diffuse covariance matrix */
        utility_cseig(a->hEig, &(a->DCM_array[band*(a->nMics)*(a->nMics)]), a->nMics, 1, U, E, NULL);

        /* Compute spatial whitening matrix */
        for(i=0; i<a->nMics; i++)
            E[i*a->nMics+i] = cmplxf(sqrtf(1.0f/(crealf(E[i*a->nMics+i])+2.23e-10f)), 0.0f);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, a->nMics, a->nMics, a->nMics, &calpha,
                    E, a->nMics,
                    U, a->nMics, &cbeta,
                    a->T[band], a->nMics);

        /* Whiten the array steering vectors / anechoic relative transfer functions (RTFs) */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a->nMics, a->nGrid, a->nMics, &calpha,
                    a->T[band], a->nMics,
                    &(a->H_array[band*(a->nMics)*(a->nGrid)]), a->nGrid, &cbeta,
                    &(a->H_array_w[band*(a->nMics)*(a->nGrid)]), a->nGrid);
    }
    free(U);
    free(E);
    free(H_W);

    /* Run-time variables */
    a->inputBlock = (float**)malloc2d(a->nMics, a->blocksize, sizeof(float));
    a->Cx = malloc1d(a->nBands*sizeof(CxMic));
    a->V  = malloc1d((a->nMics)*(a->nMics)*sizeof(float_complex));
    a->Vn = malloc1d((a->nMics)*(a->nMics)*sizeof(float_complex));
    a->lambda = malloc1d((a->nMics)*sizeof(float));

    /* Flush run-time buffers with zeros */
    hades_analysis_reset((*phAna));
}

void hades_analysis_destroy
(
    hades_analysis_handle* const phAna
)
{
    hades_analysis_data *a = (hades_analysis_data*)(*phAna);

    if (a != NULL) {
        free(a->h_array);
        free(a->H_array);
        free(a->H_array_w);
        free(a->DCM_array);
        free(a->W);
        free(a->T);
        free(a->grid_dirs_xyz);
        free(a->grid_dirs_deg);

        /* Destroy time-frequency transform  */
        switch(a->fbOpt){
            case HADES_USE_AFSTFT_LD: /* fall through */
            case HADES_USE_AFSTFT:    afSTFT_destroy(&(a->hFB_enc)); break;
        }
        free(a->freqVector);

        /* Destroy DoA estimator */
        utility_cseig_destroy(&(a->hEig));
        switch(a->doaOpt){
            case HADES_USE_MUSIC:
                hades_sdMUSIC_destroy(&(a->hDoA));
                break;
        }

        /* Free run-time variables */
        free(a->inputBlock);
        free(a->Cx);
        free(a->V);
        free(a->Vn);
        free(a->lambda);

        free(a);
        a = NULL;
        (*phAna) = NULL;
    }
}

void hades_analysis_reset
(
    hades_analysis_handle const hAna
)
{
    hades_analysis_data *a;
    int band;
    if(hAna==NULL)
        return;
    a = (hades_analysis_data*)(hAna);

    for(band=0; band<a->nBands; band++)
        memset(a->Cx[band].Cx, 0, HADES_MAX_NMICS*HADES_MAX_NMICS*sizeof(float_complex));
}

void hades_analysis_apply
(
    hades_analysis_handle const hAna,
    float** input,
    int nChannels,
    int blocksize,
    void* const hPCon,
    void* const hSCon
)
{
    hades_analysis_data *a = (hades_analysis_data*)(hAna);
    hades_param_container_data *pcon = (hades_param_container_data*)(hPCon);
    hades_signal_container_data *scon = (hades_signal_container_data*)(hSCon);
    int i, j, k, ch, band, est_idx;
    float diffuseness;
    CxMic Cx_new, T_Cx, T_Cx_TH;
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta = cmplxf(0.0f, 0.0f); /* blas */

    assert(blocksize==a->blocksize);

    /* Load time-domain data */
    for(ch=0; ch<SAF_MIN(nChannels, a->nMics); ch++)
        cblas_scopy(blocksize, input[ch], 1, a->inputBlock[ch], 1);
    for(; ch<a->nMics; ch++)
        memset(a->inputBlock[ch], 0, blocksize*sizeof(float));

    /* Forward time-frequency transform */
    switch(a->fbOpt){
        case HADES_USE_AFSTFT_LD: /* fall through */
        case HADES_USE_AFSTFT:    afSTFT_forward_knownDimensions(a->hFB_enc, a->inputBlock, blocksize, a->nMics, a->timeSlots, scon->inTF); break;
    }

    /* Update covarience matrix per band */
    for(band=0; band<a->nBands; band++){
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, a->nMics, a->nMics, a->timeSlots, &calpha,
                    FLATTEN2D(scon->inTF[band]), a->timeSlots,
                    FLATTEN2D(scon->inTF[band]), a->timeSlots, &cbeta,
                    Cx_new.Cx, a->nMics);

        /* Make a copy for the signal container */
        cblas_ccopy(a->nMics*a->nMics, (float_complex*)Cx_new.Cx, 1, (float_complex*)scon->Cx[band].Cx, 1);

        /* Apply temporal averaging */
        cblas_sscal(/*re+im*/2*(a->nMics) * (a->nMics), SAF_CLAMP(a->covAvgCoeff, 0.0f, 0.999f), (float*)a->Cx[band].Cx, 1);
        cblas_saxpy(/*re+im*/2*(a->nMics) * (a->nMics), 1.0f-SAF_CLAMP(a->covAvgCoeff, 0.0f, 0.999f), (float*)Cx_new.Cx, 1, (float*)a->Cx[band].Cx, 1);
    }

    /* Spatial parameter estimation per band */
    for (band = 0; band < a->nBands; band++) {
        /* Apply diffuse whitening process */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a->nMics, a->nMics, a->nMics, &calpha,
                    a->T[band], a->nMics,
                    a->Cx[band].Cx, a->nMics, &cbeta,
                    T_Cx.Cx, a->nMics);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, a->nMics, a->nMics, a->nMics, &calpha,
                    T_Cx.Cx, a->nMics,
                    a->T[band], a->nMics, &cbeta,
                    T_Cx_TH.Cx, a->nMics);
        utility_cseig(a->hEig, T_Cx_TH.Cx, a->nMics, 1, a->V, NULL, a->lambda);

        /* Estimate diffuseness */
        diffuseness = 0.0f;
        switch(a->diffOpt){
            case HADES_USE_COMEDIE: diffuseness = hades_comedie(a->lambda, a->nMics); break;
        }

        /* Store diffuseness and source number estimates */
        pcon->diffuseness[band] = diffuseness;
        pcon->gains_dir[band] = pcon->gains_diff[band] = 1.0f; /* Default gains per band */

        /* Apply DoA estimator */
        switch(a->doaOpt){
            case HADES_USE_MUSIC:
                /* perform sphMUSIC on the noise subspace */
                for(i=0; i<a->nMics; i++)
                    for(j=0, k=1; j<a->nMics-1; j++, k++)
                        a->Vn[i*(a->nMics-1)+j] = a->V[i*(a->nMics)+k];
                hades_sdMUSIC_compute(a->hDoA, &(a->H_array_w[band*(a->nMics)*(a->nGrid)]), a->Vn, 1, NULL, &est_idx);
                break;
        }

        /* Store */
        pcon->doa_idx[band] = pcon->gains_idx[band] = est_idx; 
    }
}

const float* hades_analysis_getFrequencyVectorPtr
(
    hades_analysis_handle const hAna,
    int* nBands
)
{
    hades_analysis_data *a;
    if(hAna==NULL){
        if(nBands!=NULL)
            (*nBands) = 0;
        return NULL;
    }
    a = (hades_analysis_data*)(hAna);
    if(nBands!=NULL)
       (*nBands) = a->nBands;
    return (const float*)a->freqVector;
}

int hades_analysis_getNbands
(
    hades_analysis_handle const hAna
)
{
    return hAna == NULL ? 0 : ((hades_analysis_data*)(hAna))->nBands;
}

float* hades_analysis_getCovarianceAvagingCoeffPtr
(
    hades_analysis_handle const hAna
)
{
    hades_analysis_data *a;
    if(hAna==NULL)
        return NULL;
    a = (hades_analysis_data*)(hAna);
    return &(a->covAvgCoeff);
}

int hades_analysis_getProcDelay
(
    hades_analysis_handle const hAna
)
{
    return hAna == NULL ? 0 : ((hades_analysis_data*)(hAna))->filterbankDelay;
}


/* ========================================================================== */
/*                      Parameter and Signal Containers                       */
/* ========================================================================== */

void hades_param_container_create
(
    hades_param_container_handle* const phPCon,
    hades_analysis_handle const hAna
)
{
    hades_param_container_data* pcon = (hades_param_container_data*)malloc1d(sizeof(hades_param_container_data));
    *phPCon = (void*)pcon;
    hades_analysis_data *a = (hades_analysis_data*)(hAna);

    /* Copy data that is relevant to the container */
    pcon->nBands = a->nBands;

    /* Allocate parameter storage */
    pcon->diffuseness = malloc1d(pcon->nBands*sizeof(float));
    pcon->doa_idx = malloc1d(pcon->nBands*sizeof(int));
    pcon->gains_idx = malloc1d(pcon->nBands*sizeof(int));

    /* Allocate the optional parameter storage */
    pcon->gains_dir = calloc1d(pcon->nBands, sizeof(float));
    pcon->gains_diff = calloc1d(pcon->nBands, sizeof(float));
}

void hades_param_container_destroy
(
    hades_param_container_handle* const phPCon
)
{
    hades_param_container_data *pcon = (hades_param_container_data*)(*phPCon);

    if (pcon != NULL) {
        /* Free parameter storage */
        free(pcon->diffuseness);
        free(pcon->doa_idx);
        free(pcon->gains_idx);
        free(pcon->gains_dir);
        free(pcon->gains_diff);

        free(pcon);
        pcon = NULL;
        (*phPCon) = NULL;
    }
}

void hades_signal_container_create
(
    hades_signal_container_handle* const phSCon,
    hades_analysis_handle const hAna
)
{
    hades_signal_container_data* scon = (hades_signal_container_data*)malloc1d(sizeof(hades_signal_container_data));
    *phSCon = (void*)scon;
    hades_analysis_data *a = (hades_analysis_data*)(hAna);

    /* Copy data that is relevant to the container */
    scon->nMics = a->nMics;
    scon->nBands = a->nBands;
    scon->timeSlots = a->timeSlots;

    /* Copy of the NON-time-averaged covariance matrix per band */
    scon->Cx = malloc1d(a->nBands*sizeof(CxMic));

    /* Time-frequency frame */
    scon->inTF = (float_complex***)malloc3d(scon->nBands, scon->nMics, scon->timeSlots, sizeof(float_complex));
}

void hades_signal_container_destroy
(
    hades_signal_container_handle* const phSCon
)
{
    hades_signal_container_data *scon = (hades_signal_container_data*)(*phSCon);

    if (scon != NULL) {
        /* Free time-frequency frame */
        free(scon->Cx);
        free(scon->inTF);

        free(scon);
        scon = NULL;
        (*phSCon) = NULL;
    }
}

#endif /* SAF_ENABLE_HADES_MODULE */
