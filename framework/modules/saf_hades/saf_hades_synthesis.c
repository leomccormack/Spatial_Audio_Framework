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
 * @file saf_hades_synthesis.c
 * @ingroup HADES
 * @brief Source file for the HADES synthesis (#SAF_HADES_MODULE)
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

#include "saf_hades_synthesis.h"
#include "saf_hades_internal.h"

#ifdef  SAF_ENABLE_HADES_MODULE

/* ========================================================================== */
/*                            HADES Radial Editor                             */
/* ========================================================================== */

void hades_radial_editor_create
(
    hades_radial_editor_handle* const phREd,
    hades_analysis_handle const hAna
)
{
    hades_radial_editor_data* e = (hades_radial_editor_data*)malloc1d(sizeof(hades_radial_editor_data));
    *phREd = (hades_radial_editor_handle)e;
    hades_analysis_data *a = (hades_analysis_data*)(hAna);

    e->nBands = a->nBands;
    e->nGrid = a->nGrid;
    e->pGrid_dirs_deg = a->grid_dirs_deg;
    e->pGrid_dirs_xyz = a->grid_dirs_xyz;
}

void hades_radial_editor_destroy
(
    hades_radial_editor_handle* const phREd
)
{
    hades_radial_editor_data *e = (hades_radial_editor_data*)(*phREd);

    if (e != NULL) {
        free(e);
        e = NULL;
        (*phREd) = NULL;
    }
}

void hades_radial_editor_apply
(
    hades_radial_editor_handle const hREd,
    hades_param_container_handle  const hPCon,
    float dirGain_dB[360]
)
{
    hades_radial_editor_data *e = (hades_radial_editor_data*)(hREd);
    hades_param_container_data *pcon = (hades_param_container_data*)(hPCon);
    int band, edit_idx;
    float azi, gain_lin;

    for(band=0; band<e->nBands; band++){
        /* Determine edit index */
        azi = e->pGrid_dirs_deg[pcon->gains_idx[band]*2];
        azi = azi < 0.0f ? azi+360.0f : azi;             /* convert -180..180 if needed */
        edit_idx = SAF_CLAMP((int)(azi + 0.5f), 0, 359); /* round to nearest integer */

        /* Extra gain factor for the direct stream */
        gain_lin = powf(10.0f, SAF_CLAMP(dirGain_dB[edit_idx], -60.0f, 12.0f)/20.0f);
        pcon->gains_dir[band] *= gain_lin;
    }
}


/* ========================================================================== */
/*                             HADES Synthesis                              */
/* ========================================================================== */

void hades_synthesis_create
(
    hades_synthesis_handle* const phSyn,
    hades_analysis_handle const hAna,
    HADES_BEAMFORMER_TYPE beamOption,
    int enableCM,
    int refIndices[2],
    hades_binaural_config* binConfig,
    HADES_HRTF_INTERP_OPTIONS interpOption
)
{
    hades_synthesis_data* s = (hades_synthesis_data*)malloc1d(sizeof(hades_synthesis_data));
    *phSyn = (hades_synthesis_handle)s;
    hades_analysis_data *a = (hades_analysis_data*)(hAna);
    int band;
    float_complex* H_W;
    hades_binaural_config* bConfig;
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta = cmplxf(0.0f, 0.0f); /* blas */

    /* User configuration parameters */
    s->beamOption = beamOption;
    s->enableCM = enableCM;
    s->refIndices[0] = refIndices[0];
    s->refIndices[1] = refIndices[1];
    s->interpOption = interpOption;
 
    /* Default user parameters */
    s->eq = malloc1d(a->nBands * sizeof(float));
    s->streamBalance = malloc1d(a->nBands * sizeof(float));
    for(band = 0; band<a->nBands; band++){
        s->eq[band] = 1.0;             /* Flat EQ */
        s->streamBalance[band] = 1.0f; /* 50/50 direct/ambient balance (i.e., no biasing) */
    }
    s->synAvgCoeff = 1.0f - 1.0f/(4096.0f/a->blocksize); /* How much averaging of current mixing matrices with the previous mixing matrices */

    /* Things relevant to the synthesiser, which are copied from the analyser to keep things aligned */
    s->fbOpt = a->fbOpt;
    s->nBands = a->nBands;
    s->hopsize = a->hopsize;
    s->blocksize = a->blocksize;
    s->nGrid = a->nGrid;
    s->nMics = a->nMics;
    s->H_array = malloc1d(s->nBands*(s->nMics)*(s->nGrid)*sizeof(float_complex));
    memcpy(s->H_array, a->H_array, s->nBands*(s->nMics)*(s->nGrid)*sizeof(float_complex));
    s->DCM_array = malloc1d(s->nBands*(s->nMics)*(s->nMics)*sizeof(float_complex));
    memcpy(s->DCM_array, a->DCM_array, s->nBands*(s->nMics)*(s->nMics)*sizeof(float_complex));
    s->W = malloc1d(s->nGrid*(s->nGrid)*sizeof(float_complex));
    memcpy(s->W, a->W, s->nGrid*(s->nGrid)*sizeof(float_complex));
    s->grid_dirs_deg = malloc1d((s->nGrid)*2*sizeof(float));
    memcpy(s->grid_dirs_deg, a->grid_dirs_deg, (s->nGrid)*2*sizeof(float));
    s->grid_dirs_xyz = (float**)malloc2d((s->nGrid), 3, sizeof(float));
    memcpy(FLATTEN2D(s->grid_dirs_xyz), a->grid_dirs_xyz, (s->nGrid)*3*sizeof(float));
    s->timeSlots = a->timeSlots;
    s->freqVector = malloc1d(s->nBands*sizeof(float));
    memcpy(s->freqVector, a->freqVector, s->nBands*sizeof(float));

    /* Time-frequency transform */
    switch(s->fbOpt){
        case HADES_USE_AFSTFT_LD: afSTFT_create(&(s->hFB_dec), 0, NUM_EARS, s->hopsize, 1, a->hybridmode, AFSTFT_BANDS_CH_TIME); break;
        case HADES_USE_AFSTFT:    afSTFT_create(&(s->hFB_dec), 0, NUM_EARS, s->hopsize, 0, a->hybridmode, AFSTFT_BANDS_CH_TIME); break;
    }
 
    /* Copy binaural configuration */
    s->binConfig = malloc1d(sizeof(hades_binaural_config));
    bConfig = s->binConfig;
    bConfig->lHRIR = binConfig->lHRIR;
    bConfig->nHRIR = binConfig->nHRIR;
    bConfig->hrir_fs = binConfig->hrir_fs;
    bConfig->hrirs = malloc1d(bConfig->nHRIR * NUM_EARS * (bConfig->lHRIR) * sizeof(float));
    memcpy(bConfig->hrirs, binConfig->hrirs, bConfig->nHRIR * NUM_EARS * (bConfig->lHRIR) * sizeof(float));
    bConfig->hrir_dirs_deg = malloc1d(bConfig->nHRIR*2*sizeof(float));
    memcpy(bConfig->hrir_dirs_deg, binConfig->hrir_dirs_deg, bConfig->nHRIR*2*sizeof(float));

    /* Pre-process HRTFs, interpolate them for the scanning grid */
    s->H_bin = calloc1d(s->nBands*NUM_EARS*(s->nGrid),sizeof(float_complex));
    hades_getInterpolatedHRTFs(hAna, interpOption, bConfig, a->grid_dirs_deg, s->nGrid, s->H_bin);

    /* Diffuse rendering variables */
    s->DCM_bin_norm = malloc1d(a->nBands*NUM_EARS*NUM_EARS*sizeof(float_complex));
    H_W = malloc1d(NUM_EARS*(a->nGrid)*sizeof(float_complex));
    s->diffEQ = malloc1d(a->nBands*sizeof(float));
    for(band=0; band<s->nBands; band++){
        /* Binaural diffuse coherence matrix (not normalised yet!) */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, a->nGrid, a->nGrid, &calpha,
                    &(s->H_bin[band*NUM_EARS*(a->nGrid)]), a->nGrid,
                    a->W, a->nGrid, &cbeta,
                    H_W, a->nGrid);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, NUM_EARS, NUM_EARS, a->nGrid, &calpha,
                    H_W, a->nGrid,
                    &(s->H_bin[band*NUM_EARS*(a->nGrid)]), a->nGrid, &cbeta,
                    &(s->DCM_bin_norm[band*NUM_EARS*NUM_EARS]), NUM_EARS);
        cblas_sscal(/*re+im*/2*NUM_EARS*NUM_EARS, 1.0f/(float)a->nGrid, (float*)&(s->DCM_bin_norm[band*NUM_EARS*NUM_EARS]), 1);

        /* Compute EQ required to bring the overall diffuse-field magnitude response of the array to that of the HRTFs instead: */
        /* sqrt(trace(H_bin_dcm(:,:,band))/(trace(H_grid_dcm(refIndices,refIndices,band))+eps)) */
        s->diffEQ[band] = sqrtf((crealf(s->DCM_bin_norm[band*NUM_EARS*NUM_EARS]) + crealf(s->DCM_bin_norm[band*NUM_EARS*NUM_EARS + 3])) /
                                (crealf(s->DCM_array[band*(s->nMics)*(s->nMics) + s->refIndices[0]*(s->nMics) + s->refIndices[0]]) +
                                 crealf(s->DCM_array[band*(s->nMics)*(s->nMics) + s->refIndices[1]*(s->nMics) + s->refIndices[1]]) + 2.23e-10f));
        s->diffEQ[band] = SAF_MIN(s->diffEQ[band], 3.0f); /* Cap at a maximum of +9dB */

        /* Normalise the binaural diffuse coherence matrix */
        cblas_sscal(/*re+im*/2*NUM_EARS*NUM_EARS, 1.0f/(crealf(s->DCM_bin_norm[band*NUM_EARS*NUM_EARS]) + crealf(s->DCM_bin_norm[band*NUM_EARS*NUM_EARS + 3]) + 2.23e-10f),
                    (float*)&(s->DCM_bin_norm[band*NUM_EARS*NUM_EARS]), 1);
    }
    free(H_W);

    /* Run-time variables */
    utility_cpinv_create(&(s->hPinv), s->nMics, s->nMics);
    utility_cglslv_create(&(s->hLinSolve), s->nMics, s->nMics);
    cdf4sap_cmplx_create(&(s->hCDF), s->nMics, NUM_EARS);
    s->As   = malloc1d(s->nMics*sizeof(float_complex));
    s->As_l = malloc1d(s->nMics*sizeof(float_complex));
    s->As_r = malloc1d(s->nMics*sizeof(float_complex));
    s->Q_diff = malloc1d(NUM_EARS*(s->nMics)*sizeof(float_complex));
    s->Q_dir  = malloc1d(NUM_EARS*(s->nMics)*sizeof(float_complex));
    s->Q      = malloc1d(NUM_EARS*(s->nMics)*sizeof(float_complex));
    s->Cy = malloc1d(NUM_EARS*NUM_EARS*sizeof(float_complex));
    s->new_M = malloc1d(NUM_EARS*(s->nMics)*sizeof(float_complex));
    s->M  = (float_complex**)malloc2d(s->nBands, NUM_EARS*(s->nMics), sizeof(float_complex));

    /* Run-time audio buffers */
    s->outTF = (float_complex***)malloc3d(s->nBands, NUM_EARS, s->timeSlots, sizeof(float_complex));
    s->outTD = (float**)malloc2d(NUM_EARS, s->blocksize, sizeof(float));

    /* Flush run-time buffers with zeros */
    hades_synthesis_reset((*phSyn));
}

void hades_synthesis_destroy
(
    hades_synthesis_handle* const phSyn
)
{
    hades_synthesis_data *s = (hades_synthesis_data*)(*phSyn);

    if (s != NULL) {
        /* Free user parameters */
        free(s->eq);
        free(s->streamBalance);
        free(s->binConfig);

        /* Free things copied from analyser */
        free(s->H_array);
        free(s->DCM_array);
        free(s->W);
        free(s->grid_dirs_deg);
        free(s->grid_dirs_xyz);
        free(s->freqVector);

        /* Free time-frequency transform */
        switch(s->fbOpt){
            case HADES_USE_AFSTFT_LD: /* fall through */
            case HADES_USE_AFSTFT:    afSTFT_destroy(&(s->hFB_dec)); break;
        }

        /* HRTF and diffuse rendering variables */
        free(s->H_bin);
        free(s->DCM_bin_norm);
        free(s->diffEQ);

        /* Run-time variables */
        utility_cpinv_destroy(&(s->hPinv));
        utility_cglslv_destroy(&(s->hLinSolve));
        cdf4sap_cmplx_destroy(&(s->hCDF));
        free(s->As);
        free(s->As_l);
        free(s->As_r);
        free(s->Q_diff);
        free(s->Q_dir);
        free(s->Q);
        free(s->Cy);
        free(s->new_M);
        free(s->M);

        /* Run-time audio buffers */
        free(s->outTF);
        free(s->outTD);

        free(s);
        s = NULL;
        (*phSyn) = NULL;
    }
}

void hades_synthesis_reset
(
    hades_synthesis_handle const hSyn
)
{
    hades_synthesis_data *s;
    if(hSyn==NULL)
        return;
    s = (hades_synthesis_data*)(hSyn);

    /* Zero buffers, matrices etc. */
    switch(s->fbOpt){
        case HADES_USE_AFSTFT_LD: /* fall through */
        case HADES_USE_AFSTFT:    afSTFT_clearBuffers(s->hFB_dec); break;
    }
    memset(FLATTEN2D(s->M), 0, s->nBands*NUM_EARS*(s->nMics)*sizeof(float_complex));
}

void hades_synthesis_apply
(
    hades_synthesis_handle const hSyn,
    hades_param_container_handle  const hPCon,
    hades_signal_container_handle const hSCon,
    int nChannels,
    int blocksize,
    float** output
)
{
    hades_synthesis_data *s = (hades_synthesis_data*)(hSyn);
    hades_param_container_data *pcon = (hades_param_container_data*)(hPCon);
    hades_signal_container_data *scon = (hades_signal_container_data*)(hSCon);
    int i, j, ch, nMics, band, doa_idx, gain_idx;
    float a, b, diffuseness, synAvgCoeff, streamBalance, eq, gain_dir, gain_diff, trace_M, reg_M, sum_As, targetEnergy;
    float_complex g_l, g_r, h_dir[NUM_EARS], AsH_invCx_As;
    float_complex Cx[HADES_MAX_NMICS*HADES_MAX_NMICS], conj_As[HADES_MAX_NMICS], AsH_invCx[HADES_MAX_NMICS*HADES_MAX_NMICS];
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta = cmplxf(0.0f, 0.0f); /* blas */

    nMics = s->nMics;
    synAvgCoeff = SAF_CLAMP((s->synAvgCoeff), 0.0f, 0.99f);

    /* Loop over bands and compute the mixing matrices */
    for (band = 0; band < s->nBands; band++) {
        /* Pull estimated (and possibly modified) spatial parameters for this band */
        diffuseness = pcon->diffuseness[band];
        saf_assert(diffuseness>-0.0001f && diffuseness < 1.00001f, "Erroneous parameter analysis");
        doa_idx = pcon->doa_idx[band];
        gain_idx = pcon->gains_idx[band];
        gain_dir  = pcon->gains_dir[band];
        gain_diff = pcon->gains_diff[band];

        /* Optional biasing (e.g. to conduct de-reverberation or to emphasise reverberation) */
        streamBalance = SAF_CLAMP(s->streamBalance[band], 0.0f, 2.0f);
        eq = s->eq[band];
        if(streamBalance<1.0f){
            a = streamBalance;        /* pump more direct energy into output */
            b = 1.0f;                 /* pass ambient stream as normal */
        }
        else {
            a = 1.0f;                 /* pass source stream as normal */
            b = 2.0f - streamBalance; /* pump less ambient energy into output */
        }
        a *= gain_dir;
        b *= gain_diff;

        /* Source array steering vector for the estimated DoAs */
        for(i=0; i<nMics; i++)
            s->As[i] = s->H_array[band*nMics*(s->nGrid) + i*(s->nGrid) + doa_idx];

        /* Anechoic relative transfer functions (RTFs) */
        for(i=0; i<nMics; i++){
            s->As_l[i] = ccdivf(s->As[i], s->As[s->refIndices[0]]);
            s->As_r[i] = ccdivf(s->As[i], s->As[s->refIndices[1]]);
        }

        /* HRTF for this reproduction DoA */
        h_dir[0] = s->H_bin[band*NUM_EARS*(s->nGrid) + 0*(s->nGrid) + gain_idx];
        h_dir[1] = s->H_bin[band*NUM_EARS*(s->nGrid) + 1*(s->nGrid) + gain_idx];
        g_l = ccdivf(h_dir[0], s->As[s->refIndices[0]]); /* (Relative transfer functions) */
        g_r = ccdivf(h_dir[1], s->As[s->refIndices[1]]);
        if(cabsf(g_l)>4.0f || cabsf(g_r)>4.0f) /* if >12dB, then bypass: */
            g_l = g_r = cmplxf(1.0f, 0.0f);

        /* Diffuse mixing matrix (if the sound-field is analysed to be more diffuse, then we mix in more of just the reference sensors) */
        memset(s->Q_diff, 0, NUM_EARS*nMics*sizeof(float_complex));
        s->Q_diff[0*nMics+s->refIndices[0]] = cmplxf(s->diffEQ[band], 0.0f);
        s->Q_diff[1*nMics+s->refIndices[1]] = cmplxf(s->diffEQ[band], 0.0f);

        /* Source mixing matrix (beamforming towards the estimated DoAs) */
        switch(s->beamOption){
            case HADES_BEAMFORMER_NONE: /* No beamforming required */ break;
            case HADES_BEAMFORMER_FILTER_AND_SUM:
                /* Normalise the beamformers to unity gain in the look direction */
                utility_cpinv(s->hPinv, s->As_l, nMics, 1, s->Q_dir);
                utility_cpinv(s->hPinv, s->As_r, nMics, 1, s->Q_dir + nMics);

                /* Now bring their response from being w.r.t the array to being w.r.t the HRTF instead */
                cblas_cscal(nMics, &g_l, s->Q_dir, 1);
                cblas_cscal(nMics, &g_r, s->Q_dir + nMics, 1);
                break;

            case HADES_BEAMFORMER_BMVDR:
                /* prep */
                cblas_ccopy(nMics*nMics, scon->Cx[band].Cx, 1, Cx, 1);
                trace_M = 0.0f;
                for(i=0; i<nMics; i++)
                    trace_M += crealf(Cx[i*nMics+i]);
                sum_As = cblas_scasum(nMics, s->As, 1);

                /* Compute beamforming weights if checks pass */
                if( trace_M < 0.0001f || sum_As < 0.0001f)
                    memset(s->Q_dir, 0, NUM_EARS*nMics*sizeof(float_complex));
                else{
                    /* Regularise Cx */
                    reg_M = (trace_M/(float)nMics) * 10.0f + 0.0001f;
                    for(i=0; i<nMics; i++)
                        Cx[i*nMics+i] = craddf(Cx[i*nMics+i], reg_M);

                    /* Compute MVDR weights w.r.t the reference sensor at each ear, [As^H Cx^-1 As]^-1 As^H Cx^-1  */
                    for(j=0; j<NUM_EARS; j++){
                        /* Solve As^H Cx-1 */
                        utility_cvconj(j==0 ? s->As_l : s->As_r, nMics, conj_As);
                        utility_cglslv(s->hLinSolve, Cx, nMics, conj_As, 1, AsH_invCx);

                        /* Compute As^H Cx-1 As */
                        utility_cvvdot(AsH_invCx, j==0 ? s->As_l : s->As_r, nMics, NO_CONJ, &AsH_invCx_As);
                        AsH_invCx_As = craddf(AsH_invCx_As, 0.00001f);

                        /* The solution */
                        AsH_invCx_As = ccdivf(cmplxf(1.0f, 0.0f), AsH_invCx_As);
                        cblas_cscal(nMics, &AsH_invCx_As, AsH_invCx, 1);
                        cblas_ccopy(nMics, AsH_invCx, 1, s->Q_dir + j*nMics, 1);
                    }

                    /* Now bring their response from being w.r.t the array to instead being w.r.t the HRTF */
                    cblas_cscal(nMics, &g_l, s->Q_dir, 1);
                    cblas_cscal(nMics, &g_r, s->Q_dir + nMics, 1);
                }
                break;
        }

        /* Prototype mixing matrix */
        if(s->beamOption==HADES_BEAMFORMER_NONE){
            /* No beamforming (just pass through the reference signals) */
            memset(s->Q, 0, NUM_EARS*nMics*sizeof(float_complex));
            s->Q[0*nMics+s->refIndices[0]] = cmplxf(1.0f, 0.0f);
            s->Q[1*nMics+s->refIndices[1]] = cmplxf(1.0f, 0.0f);
            //s->Q[0*nMics+s->refIndices[0]] = cmplxf(s->diffEQ[band], 0.0f);
            //s->Q[1*nMics+s->refIndices[1]] = cmplxf(s->diffEQ[band], 0.0f);
        }
        else{
            /* Mix in the beamforming weights, conforming to the assumed direct-diffuse model */
            cblas_ccopy(NUM_EARS*nMics, s->Q_dir, 1, s->Q, 1);
            cblas_sscal(/*re+im*/2*NUM_EARS*nMics, eq*a*(1.0f-diffuseness), (float*)s->Q, 1);
            cblas_saxpy(/*re+im*/2*NUM_EARS*nMics, eq*b*diffuseness, (float*)s->Q_diff, 1, (float*)s->Q, 1);
        }

        /* Target output signal energy (used for the covariance matching) */
        targetEnergy = 0.0f;
        for(i=0; i<nMics; i++)
            targetEnergy += crealf(scon->Cx[band].Cx[i*nMics+i]);
        targetEnergy = eq*0.25f*targetEnergy * s->diffEQ[band];

        /* Final mixing matrix */
        if(s->enableCM && targetEnergy>0.0001f){
            /* "Direct" contributions to the target spatial covariance matrix */
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, NUM_EARS, NUM_EARS, 1, &calpha,
                        h_dir, 1,
                        h_dir, 1, &cbeta,
                        s->Cy, NUM_EARS);
            cblas_sscal(/*re+im*/2*NUM_EARS*NUM_EARS, eq*a*(1.0f-diffuseness)*targetEnergy, (float*)s->Cy, 1);

            /* "Diffuse" contributions to the target spatial covariance matrix */
            cblas_saxpy(/*re+im*/2*NUM_EARS*NUM_EARS, eq*b*diffuseness*targetEnergy, (float*)&(s->DCM_bin_norm[band*NUM_EARS*NUM_EARS]), 1, (float*)s->Cy, 1);

            /* Solve the covariance matching problem */
            formulate_M_and_Cr_cmplx(s->hCDF, (float_complex*)scon->Cx[band].Cx, s->Cy, s->Q, 1, 0.1f, s->new_M, NULL);
        }
        else
            cblas_ccopy(NUM_EARS*nMics, s->Q, 1, s->new_M, 1);

        /* Optional Equalisation */
        cblas_sscal(/*re+im*/2*NUM_EARS*nMics, eq, (float*)s->new_M, 1);

        /* Temporal averaging of mixing matrices */
        cblas_sscal(/*re+im*/2*NUM_EARS*nMics, synAvgCoeff, (float*)s->M[band], 1);
        cblas_saxpy(/*re+im*/2*NUM_EARS*nMics, 1.0f-synAvgCoeff, (float*)s->new_M, 1, (float*)s->M[band], 1);
    }

    /* Apply mixing matrices */
    for(band=0; band<s->nBands; band++){
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NUM_EARS, s->timeSlots, nMics, &calpha,
                    s->M[band], nMics,
                    FLATTEN2D(scon->inTF[band]), s->timeSlots, &cbeta,
                    FLATTEN2D(s->outTF[band]), s->timeSlots);
    }

    /* inverse time-frequency transform */
    switch(s->fbOpt){
        case HADES_USE_AFSTFT_LD: /* fall through */
        case HADES_USE_AFSTFT:    afSTFT_backward_knownDimensions(s->hFB_dec, s->outTF, blocksize, NUM_EARS, s->timeSlots, s->outTD);   break;
    }

    /* Copy to output */
    for(ch=0; ch<SAF_MIN(nChannels, NUM_EARS); ch++)
        memcpy(output[ch], s->outTD[ch], blocksize*sizeof(float));
    for(; ch<nChannels; ch++)
        memset(output[ch], 0, blocksize*sizeof(float));
}

float* hades_synthesis_getEqPtr
(
    hades_synthesis_handle const hSyn,
    int* nBands
)
{
    hades_synthesis_data *s;
    if(hSyn==NULL){
        if(nBands!=NULL)
            (*nBands) = 0;
        return NULL;
    }
    s = (hades_synthesis_data*)(hSyn);
    if(nBands!=NULL)
       (*nBands) = s->nBands;
    return s->eq;
}

float* hades_synthesis_getStreamBalancePtr
(
    hades_synthesis_handle const hSyn,
    int* nBands
)
{
    hades_synthesis_data *s;
    if(hSyn==NULL){
        if(nBands!=NULL)
            (*nBands) = 0;
        return NULL;
    }
    s = (hades_synthesis_data*)(hSyn);
    if(nBands!=NULL)
       (*nBands) = s->nBands;
    return s->streamBalance;
}

float* hades_synthesis_getSynthesisAveragingCoeffPtr
(
    hades_synthesis_handle const hSyn
)
{
    hades_synthesis_data *s;
    if(hSyn==NULL)
        return NULL;
    s = (hades_synthesis_data*)(hSyn);
    return &(s->synAvgCoeff);
}
int hades_synthesis_getProcDelay
(
    hades_synthesis_handle const hSyn
)
{
    if(hSyn==NULL)
        return 0;
    return 0; /* Accounted for in hades_analysis_getProcDelay() */
}

#endif /* SAF_ENABLE_HADES_MODULE */
