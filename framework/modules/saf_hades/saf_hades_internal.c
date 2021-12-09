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
 * @file saf_hades_internal.c
 * @ingroup HADES
 * @brief Internal source for the HADES module (#SAF_HADES_MODULE)
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

#include "saf_hades_internal.h"

#ifdef  SAF_ENABLE_HADES_MODULE

void hades_getInterpolatedHRTFs
(
    hades_analysis_handle const hAna,
    HADES_HRTF_INTERP_OPTIONS interpOption,
    hades_binaural_config* binConfig,
    float* target_dirs_deg,
    int nTargetDirs,
    float_complex* hrtf_interp
)
{
    hades_analysis_data *a = (hades_analysis_data*)(hAna);
    int band, i, j, ntable, ntri;
    int* idx;
    float* itds_s, *interpTable, *w;
    float_complex*** hrtf_fb;    /* nBands x NUM_EARS x N_dirs */

    /* Pass HRIRs through the filterbank */
    hrtf_fb = (float_complex***)malloc3d(a->nBands, NUM_EARS, binConfig->nHRIR, sizeof(float_complex));
    switch(a->fbOpt){
        case HADES_USE_AFSTFT_LD: HRIRs2HRTFs_afSTFT(binConfig->hrirs, binConfig->nHRIR, binConfig->lHRIR, a->hopsize, 1, a->hybridmode, FLATTEN3D(hrtf_fb)); break;
        case HADES_USE_AFSTFT:    HRIRs2HRTFs_afSTFT(binConfig->hrirs, binConfig->nHRIR, binConfig->lHRIR, a->hopsize, 0, a->hybridmode, FLATTEN3D(hrtf_fb)); break;
    }

    /* Integration weights */
    if (cblas_sasum(nTargetDirs, target_dirs_deg+1, 2)/(float)nTargetDirs<0.0001)
        w = NULL;
    else{
        w = malloc1d(nTargetDirs*sizeof(float));
        getVoronoiWeights(target_dirs_deg, nTargetDirs, 0, w);
    }

    /* estimate the ITDs for each HRIR */
    itds_s = malloc1d(binConfig->nHRIR*sizeof(float));
    estimateITDs(binConfig->hrirs, binConfig->nHRIR, binConfig->lHRIR, binConfig->hrir_fs, itds_s);

    /* Apply HRTF interpolation */
    switch(interpOption){
        case HADES_HRTF_INTERP_NEAREST:
            /* Quantise to nearest hrir direction */
            idx = malloc1d(nTargetDirs*sizeof(int));
            findClosestGridPoints(binConfig->hrir_dirs_deg, binConfig->nHRIR, target_dirs_deg, nTargetDirs, 1, idx, NULL, NULL);
            for(band=0; band<a->nBands; band++)
                for(i=0; i<NUM_EARS; i++)
                    for(j=0; j<nTargetDirs; j++)
                        hrtf_interp[band*NUM_EARS*nTargetDirs + i*nTargetDirs + j] = hrtf_fb[band][i][idx[j]];

            /* Diffuse-field EQ without phase-simplification */
            diffuseFieldEqualiseHRTFs(nTargetDirs, NULL, NULL, a->nBands, w, 1, 0, hrtf_interp);
            free(idx);
            break;

        case HADES_HRTF_INTERP_TRIANGULAR:
            /* Diffuse-field EQ with phase-simplification */
            diffuseFieldEqualiseHRTFs(binConfig->nHRIR, itds_s, a->freqVector, a->nBands, w, 1, 1, FLATTEN3D(hrtf_fb));

            /* Interpolation table */
            interpTable = NULL;
            generateVBAPgainTable3D_srcs(target_dirs_deg, nTargetDirs, binConfig->hrir_dirs_deg, binConfig->nHRIR, 0, 0, 0.0f, &interpTable, &ntable, &ntri);
            VBAPgainTable2InterpTable(interpTable, nTargetDirs, binConfig->nHRIR);

            /* Interpolate */
            interpHRTFs(FLATTEN3D(hrtf_fb), itds_s, a->freqVector, interpTable, binConfig->nHRIR, a->nBands, nTargetDirs, hrtf_interp);

            /* Clean-up */
            free(interpTable);
            break;
    }

    /* Clean-up */
    free(itds_s);
    free(hrtf_fb);
    free(w);
}

/** Internal data structure for sdMUSIC */
typedef struct _hades_sdMUSIC_data {
    int nMics, nDirs;
    float_complex* VnA;
    float* grid_dirs_xyz;
    float* abs_VnA;
    float* pSpec;
    float* pSpecInv;
    float* P_minus_peak;
    float* VM_mask;

}hades_sdMUSIC_data;

void hades_sdMUSIC_create
(
    void ** const phMUSIC,
    int nMics,
    float* grid_dirs_deg,
    int nDirs
)
{
    *phMUSIC = malloc1d(sizeof(hades_sdMUSIC_data));
    hades_sdMUSIC_data *h = (hades_sdMUSIC_data*)(*phMUSIC);

    h->nMics = nMics;
    h->nDirs = nDirs;

    /* store cartesian coords of scanning directions (for optional peak finding) */
    h->grid_dirs_xyz = malloc1d(h->nDirs * 3 * sizeof(float));
    unitSph2cart(grid_dirs_deg, h->nDirs, 1, h->grid_dirs_xyz);

    /* for run-time */
    h->VnA = malloc1d(h->nMics * (h->nDirs) * sizeof(float_complex));
    h->abs_VnA = malloc1d(h->nMics * (h->nDirs) * sizeof(float));
    h->pSpec = malloc1d(h->nDirs*sizeof(float));
    h->pSpecInv = malloc1d(h->nDirs*sizeof(float));
    h->P_minus_peak = malloc1d(h->nDirs*sizeof(float));
    h->VM_mask = malloc1d(h->nDirs*sizeof(float));
}

void hades_sdMUSIC_destroy
(
    void ** const phMUSIC
)
{
    hades_sdMUSIC_data *h = (hades_sdMUSIC_data*)(*phMUSIC);

    if (h != NULL) {
        free(h->grid_dirs_xyz);
        free(h->VnA);
        free(h->abs_VnA);
        free(h->pSpec);
        free(h->pSpecInv);
        free(h->P_minus_peak);
        free(h->VM_mask);
        free(h);
        h = NULL;
        *phMUSIC = NULL;
    }
}

void hades_sdMUSIC_compute
(
    void* const hMUSIC,
    float_complex* A_grid, /* nMics x nGrid */
    float_complex *Vn, /* nMics x (nMics - nSrcs) */
    int nSrcs,
    float* P_music,
    int* peak_inds
)
{
    hades_sdMUSIC_data *h = (hades_sdMUSIC_data*)(hMUSIC);
    int i, k, VnD2, peak_idx;
    float kappa, scale;
    float VM_mean[3];
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta = cmplxf(0.0f, 0.0f);

    VnD2 = h->nMics - nSrcs; /* noise subspace second dimension length */

    /* derive the pseudo-spectrum value for each grid direction */
    cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, h->nDirs, VnD2, h->nMics, &calpha,
                A_grid, h->nDirs,
                Vn, VnD2, &cbeta,
                h->VnA, VnD2);
    utility_cvabs(h->VnA, (h->nDirs)*VnD2, h->abs_VnA);
    for (i = 0; i < (h->nDirs); i++)
        h->pSpecInv[i] = cblas_sdot(VnD2, &(h->abs_VnA[i*VnD2]), 1, &(h->abs_VnA[i*VnD2]), 1);
    //h->pSpec[i] = 1.0f / (h->pSpecInv[i] + 2.23e-10f);
    utility_svrecip(h->pSpecInv, h->nDirs, h->pSpec);

    /* Output pseudo-spectrum */
    if(P_music!=NULL)
        cblas_scopy(h->nDirs, h->pSpec, 1, P_music, 1);

    /* Peak-finding */
    if(peak_inds!=NULL){
        kappa = 50.0f;
        scale = kappa/(2.0f*SAF_PI*expf(kappa)-expf(-kappa));
        cblas_scopy(h->nDirs, h->pSpec, 1, h->P_minus_peak, 1);

        /* Loop over the number of sources */
        for(k=0; k<nSrcs; k++){
            utility_simaxv(h->P_minus_peak, h->nDirs, &peak_idx);
            peak_inds[k] = peak_idx;
            if(k==nSrcs-1)
                break;
            VM_mean[0] = h->grid_dirs_xyz[peak_idx*3];
            VM_mean[1] = h->grid_dirs_xyz[peak_idx*3+1];
            VM_mean[2] = h->grid_dirs_xyz[peak_idx*3+2];

            /* Apply mask for next iteration */
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, h->nDirs, 1, 3, 1.0f,
                        h->grid_dirs_xyz, 3,
                        (const float*)VM_mean, 3, 0.0f,
                        h->VM_mask, 1);
            cblas_sscal(h->nDirs, kappa, h->VM_mask, 1);
            for(i=0; i<h->nDirs; i++)
                h->VM_mask[i] = expf(h->VM_mask[i]);             /* VM distribution */
            cblas_sscal(h->nDirs, scale, h->VM_mask, 1);
            for(i=0; i<h->nDirs; i++)
                h->VM_mask[i] = 1.0f/(0.00001f+(h->VM_mask[i])); /* inverse VM distribution */
            utility_svvmul(h->P_minus_peak, h->VM_mask, h->nDirs, h->P_minus_peak);
        }
    }
}

float hades_comedie
(
    float* lambda,
    int N
)
{
    int i;
    float Nord, sum, g_0, mean_ev, g, sumAbsDiff;

    Nord = sqrtf((float)N)-1.0f;
    sum = 0.0f;
    for(i=0; i<N; i++)
        sum += lambda[i];
    if(sum < 0.0001f) /* FLT_EPS*FLT_MIN -/+ range */
        return 1.0f;
    else{
        g_0 = 2.0f*(powf(Nord+1.0f,2.0f)-1.0f);
        mean_ev = (1.0f/(powf(Nord+1.0f,2.0f)))*sum;
        sumAbsDiff = 0.0f;
        for(i=0; i<N; i++)
            sumAbsDiff += fabsf(lambda[i]-mean_ev);
        g = (1.0f/mean_ev)*sumAbsDiff;
        /* due to numerical error small (10e-7) negative numbers were occuring
         * sometimes for the single plane-wave case; hence bounding it to >=0 */
        return SAF_MAX(1.0f-g/g_0, 0.0f);
    }
}

#endif /* SAF_ENABLE_HADES_MODULE */
