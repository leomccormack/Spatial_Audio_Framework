/*
 * Copyright 2021 Leo McCormack
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
 * @file: spreader.h
 * @brief spreader
 * @author Leo McCormack
 * @date 07.04.2021
 */

#include "spreader_internal.h"

void spreader_create
(
    void ** const phSpr
)
{
    spreader_data* pData = (spreader_data*)malloc1d(sizeof(spreader_data));
    *phSpr = (void*)pData;
    int band, t;

    printf(SAF_VERSION_LICENSE_STRING);

    /* user parameters */
    pData->nSources = 1;
    pData->procMode = SPREADER_MODE_NAIVE;
    pData->useDefaultHRIRsFLAG = 1;
    pData->covAvgCoeff = 0.85f;
    memset(pData->src_spread, 0, SPREADER_MAX_NUM_SOURCES*sizeof(float));
    memset(pData->src_dirs_deg, 0, SPREADER_MAX_NUM_SOURCES*2*sizeof(float));

    /* time-frequency transform + buffers */
    pData->hSTFT = NULL;
    pData->inputFrameTD = (float**)malloc2d(MAX_NUM_INPUTS, FRAME_SIZE, sizeof(float));
    pData->outframeTD = (float**)malloc2d(MAX_NUM_OUTPUTS, FRAME_SIZE, sizeof(float));
    pData->inputframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_INPUTS, TIME_SLOTS, sizeof(float_complex));
    pData->protoframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_OUTPUTS, TIME_SLOTS, sizeof(float_complex));
    pData->decorframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_OUTPUTS, TIME_SLOTS, sizeof(float_complex));
    pData->spreadframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_OUTPUTS, TIME_SLOTS, sizeof(float_complex));
    pData->outputframeTF = (float_complex***)malloc3d(HYBRID_BANDS, MAX_NUM_OUTPUTS, TIME_SLOTS, sizeof(float_complex));
    
    /* Internal */
    pData->Q = pData->nGrid = pData->h_len = 0;
    pData->h_fs = 0.0f;
    pData->h_grid = NULL;
    pData->H_grid = NULL;
    for(band=0; band<HYBRID_BANDS; band++)
        pData->HHH[band] = NULL;
    pData->grid_dirs_deg = NULL;
    pData->grid_dirs_xyz = NULL;
    pData->weights = NULL;
    pData->hDecor = NULL;
    pData->Cy = NULL;
    pData->Cproto = NULL;
    pData->new_M = NULL;
    pData->prev_M = NULL;
    pData->new_Mr = NULL;
    pData->prev_Mr = NULL;
    pData->interp_M = NULL;
    pData->interp_Mr = NULL;
    pData->interp_Mr_cmplx = NULL;
    pData->inFrame_t = NULL;
    for(t=0; t<TIME_SLOTS; t++){
        pData->interpolatorFadeIn[t] = ((float)t+1.0f)/(float)TIME_SLOTS;
        pData->interpolatorFadeOut[t] = 1.0f - ((float)t+1.0f)/(float)TIME_SLOTS;
    }

    /* Optimal mixing */
    pData->hCdf = NULL;
    pData->hCdf_res = NULL;
    pData->Qmix = NULL;
    pData->Qmix_cmplx = NULL;
    pData->Cr = NULL;
    pData->Cr_cmplx = NULL;

    /* flags/status */
    pData->new_nSources = pData->nSources;
    pData->progressBar0_1 = 0.0f;
    pData->progressBarText = malloc1d(PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
    strcpy(pData->progressBarText,"");
    pData->codecStatus = CODEC_STATUS_NOT_INITIALISED;
    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

void spreader_destroy
(
    void ** const phSpr
)
{
    spreader_data *pData = (spreader_data*)(*phSpr);
    int band;

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
        free(pData->outframeTD);
        free(pData->inputframeTF);
        free(pData->decorframeTF);
        free(pData->spreadframeTF);
        free(pData->outputframeTF);

        /* internal */
        free(pData->h_grid);
        free(pData->H_grid);
        for(band=0; band<HYBRID_BANDS; band++)
            free(pData->HHH[band]);
        free(pData->grid_dirs_deg);
        free(pData->grid_dirs_xyz);
        free(pData->weights);
        free(pData->hDecor);
        free(pData->angles);
        free(pData->Cy);
        free(pData->Cproto);
        free(pData->new_M);
        free(pData->prev_M);
        free(pData->new_Mr);
        free(pData->prev_Mr);
        free(pData->interp_M);
        free(pData->interp_Mr);
        free(pData->interp_Mr_cmplx);
        free(pData->inFrame_t);

        /* Optimal mixing */
        cdf4sap_cmplx_destroy(&(pData->hCdf));
        cdf4sap_destroy(&(pData->hCdf_res));
        free(pData->Qmix);
        free(pData->Qmix_cmplx);
        free(pData->Cr);
        free(pData->Cr_cmplx);

        free(pData->progressBarText);
         
        free(pData);
        pData = NULL;
    }
}

void spreader_init
(
    void * const hSpr,
    int          sampleRate
)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    
    /* define frequency vector */
    pData->fs = sampleRate;
    afSTFT_getCentreFreqs(pData->hSTFT, (float)sampleRate, HYBRID_BANDS, pData->freqVector);
}

void spreader_initCodec
(
    void* const hSpr
)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    int q, band, ng, nSources;
    float_complex scaleC;
    SAF_SOFA_ERROR_CODES error;
    float_complex H_tmp[MAX_NUM_CHANNELS];
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    if (pData->codecStatus != CODEC_STATUS_NOT_INITIALISED)
        return; /* re-init not required, or already happening */
    while (pData->procStatus == PROC_STATUS_ONGOING){
        /* re-init required, but we need to wait for the current processing loop to end */
        pData->codecStatus = CODEC_STATUS_INITIALISING; /* indicate that we want to init */
        SAF_SLEEP(10);
    }

    nSources = pData->new_nSources;
    
    /* for progress bar */
    pData->codecStatus = CODEC_STATUS_INITIALISING;
    strcpy(pData->progressBarText,"Initialising");
    pData->progressBar0_1 = 0.0f;

    /* Load measurements (e.g. HRIRs, Microphone array IRs etc.) */
    if(pData->useDefaultHRIRsFLAG){
        /* Load default HRIR data */
        pData->Q = NUM_EARS;
        pData->nGrid = __default_N_hrir_dirs;
        pData->h_len = __default_hrir_len;
        pData->h_fs = (float)__default_hrir_fs;
        pData->h_grid = realloc1d(pData->h_grid, pData->nGrid * (pData->Q) * (pData->h_len) * sizeof(float));
        memcpy(pData->h_grid, (float*)__default_hrirs, pData->nGrid * (pData->Q) * (pData->h_len) * sizeof(float));
        pData->grid_dirs_deg = realloc1d(pData->grid_dirs_deg, pData->nGrid * 2 * sizeof(float));
        memcpy(pData->grid_dirs_deg, (float*)__default_hrir_dirs_deg, pData->nGrid * 2 * sizeof(float));
    }
    else{
        saf_sofa_container sofa;
        error = saf_sofa_open(&sofa, pData->sofa_filepath);
        if(error!=SAF_SOFA_OK){
            /* if failed, then load default data instead */
            pData->useDefaultHRIRsFLAG = 1;
            spreader_initCodec(hSpr);
        }
        pData->h_fs = (int)sofa.DataSamplingRate;
        pData->h_len = sofa.DataLengthIR;
        pData->nGrid = sofa.nSources;
        pData->h_grid = realloc1d(pData->h_grid, pData->nGrid*(pData->Q)*(pData->h_len)*sizeof(float));
        memcpy(pData->h_grid, sofa.DataIR, pData->nGrid*(pData->Q)*(pData->h_len)*sizeof(float));
        pData->grid_dirs_deg = realloc1d(pData->grid_dirs_deg, pData->nGrid*2*sizeof(float));
        cblas_scopy(pData->nGrid, sofa.SourcePosition, 3, pData->grid_dirs_deg, 2); /* azi */
        cblas_scopy(pData->nGrid, &sofa.SourcePosition[1], 3, &pData->grid_dirs_deg[1], 2); /* elev */
        saf_sofa_close(&sofa);
    }

    /* Initialise time-frequency transform and decorrelators */
    afSTFT_destroy(&(pData->hSTFT));
    afSTFT_create(&(pData->hSTFT), nSources, pData->Q, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    int orders[4] = {20, 15, 6, 6}; /* 20th order up to 700Hz, 15th->2.4kHz, 6th->4kHz, 3rd->12kHz, NONE(only delays)->Nyquist */
    //float freqCutoffs[4] = {600.0f, 2.6e3f, 4.5e3f, 12e3f};
    float freqCutoffs[4] = {900.0f, 6.8e3f, 12e3f, 24e3f};
    const int maxDelay = 12;
    latticeDecorrelator_destroy(&(pData->hDecor));
    latticeDecorrelator_create(&(pData->hDecor), pData->fs, HOP_SIZE, pData->freqVector, HYBRID_BANDS, pData->Q, orders, freqCutoffs, 4, maxDelay, 0, 0.75f);

    /* Convert to filterbank coefficients and pre-compute outer products */
    pData->H_grid = realloc1d(pData->H_grid, HYBRID_BANDS*(pData->Q)*pData->nGrid*sizeof(float_complex));
    afSTFT_FIRtoFilterbankCoeffs(pData->h_grid, pData->nGrid, pData->Q, pData->h_len, HOP_SIZE, 0, 1, pData->H_grid);
    pData->weights = realloc1d(pData->weights, pData->nGrid*sizeof(float));
    getVoronoiWeights(pData->grid_dirs_deg, pData->nGrid, 0, pData->weights);
    cblas_sscal(pData->nGrid, 1.0f/FOURPI, pData->weights, 1);
    for(band=0; band<HYBRID_BANDS; band++){
        pData->HHH[band] = (float_complex**)realloc2d((void**)pData->HHH[band], pData->nGrid, pData->Q * (pData->Q), sizeof(float_complex));
        for(ng=0; ng<pData->nGrid; ng++){
            for(q=0; q<pData->Q; q++)
                H_tmp[q] = pData->H_grid[band*(pData->Q)*pData->nGrid + q*pData->nGrid + ng];
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, pData->Q, pData->Q, 1, &calpha,
                        H_tmp, 1,
                        H_tmp, 1, &cbeta,
                        pData->HHH[band][ng], pData->Q);
            scaleC = cmplxf(pData->weights[ng], 0.0f);
            cblas_cscal(pData->Q * (pData->Q), &scaleC, pData->HHH[band][ng], 1);
        }
    }
    pData->angles = realloc1d(pData->angles, pData->nGrid*sizeof(float));

    /* OM structures */
    cdf4sap_cmplx_destroy(&(pData->hCdf));
    cdf4sap_cmplx_create(&(pData->hCdf), pData->Q, pData->Q);
    cdf4sap_destroy(&(pData->hCdf_res));
    cdf4sap_create(&(pData->hCdf_res), pData->Q, pData->Q);
    pData->Qmix = realloc1d(pData->Qmix, pData->Q*(pData->Q)*sizeof(float));
    memset(pData->Qmix, 0, pData->Q*(pData->Q)*sizeof(float));
    pData->Qmix_cmplx = realloc1d(pData->Qmix_cmplx, pData->Q*(pData->Q)*sizeof(float_complex));
    memset(pData->Qmix_cmplx, 0, pData->Q*(pData->Q)*sizeof(float_complex));
    for(q=0; q<pData->Q; q++){
        pData->Qmix[q*(pData->Q)+q] = 1.0f;
        pData->Qmix_cmplx[q*(pData->Q)+q] = cmplxf(1.0f, 0.0f);
    }
    pData->Cr = realloc1d(pData->Cr, pData->Q*(pData->Q)*sizeof(float));
    pData->Cr_cmplx = realloc1d(pData->Cr_cmplx, pData->Q*(pData->Q)*sizeof(float_complex));

    /* mixing matrices and buffers */
    pData->Cy = (float_complex**)realloc2d((void**)pData->Cy, HYBRID_BANDS, (pData->Q)*(pData->Q), sizeof(float_complex));
    memset(FLATTEN2D(pData->Cy), 0, HYBRID_BANDS * (pData->Q)*(pData->Q) * sizeof(float_complex));
    pData->Cproto = (float_complex**)realloc2d((void**)pData->Cproto, HYBRID_BANDS, (pData->Q)*(pData->Q), sizeof(float_complex));
    memset(FLATTEN2D(pData->Cproto), 0, HYBRID_BANDS * (pData->Q)*(pData->Q) * sizeof(float_complex));
    pData->new_M = (float_complex**)realloc2d((void**)pData->new_M, HYBRID_BANDS, (pData->Q)*(pData->Q), sizeof(float_complex));
    pData->prev_M = (float_complex**)realloc2d((void**)pData->prev_M, HYBRID_BANDS, (pData->Q)*(pData->Q), sizeof(float_complex));
    memset(FLATTEN2D(pData->prev_M), 0, HYBRID_BANDS * (pData->Q)*(pData->Q) * sizeof(float_complex));
    pData->new_Mr = (float**)realloc2d((void**)pData->new_Mr, HYBRID_BANDS, (pData->Q)*(pData->Q), sizeof(float));
    pData->prev_Mr = (float**)realloc2d((void**)pData->prev_Mr, HYBRID_BANDS, (pData->Q)*(pData->Q), sizeof(float));
    memset(FLATTEN2D(pData->prev_Mr), 0, HYBRID_BANDS * (pData->Q)*(pData->Q) * sizeof(float));
    pData->interp_M = realloc1d(pData->interp_M, (pData->Q)*(pData->Q) * sizeof(float_complex));
    pData->interp_Mr = realloc1d(pData->interp_Mr, (pData->Q)*(pData->Q) * sizeof(float));
    pData->interp_Mr_cmplx = realloc1d(pData->interp_Mr_cmplx, (pData->Q)*(pData->Q) * sizeof(float_complex));
    memset(pData->interp_Mr_cmplx, 0, (pData->Q)*(pData->Q) * sizeof(float_complex));
    pData->inFrame_t = realloc1d(pData->inFrame_t, (pData->Q) * sizeof(float_complex));

    /* New config */
    pData->nSources = pData->new_nSources;

    /* done! */
    strcpy(pData->progressBarText,"Done!");
    pData->progressBar0_1 = 1.0f;
    pData->codecStatus = CODEC_STATUS_INITIALISED;
}

void spreader_process
(
    void  *  const hSpr,
    float ** const inputs,
    float ** const outputs,
    int            nInputs,
    int            nOutputs,
    int            nSamples
)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    int q, src, ng, ch, i, j, band, t, nSources, Q, centre_ind;
    float trace;
    float src_dir_xyz[3], CprotoDiag[MAX_NUM_OUTPUTS*MAX_NUM_OUTPUTS];
    float_complex scaleC, tmp;
    float_complex tmpFrame[MAX_NUM_CHANNELS][TIME_SLOTS], H_tmp[MAX_NUM_CHANNELS], Cy_new[MAX_NUM_CHANNELS*MAX_NUM_CHANNELS];
    float_complex E_dir[MAX_NUM_CHANNELS*MAX_NUM_CHANNELS], V[MAX_NUM_OUTPUTS*MAX_NUM_OUTPUTS], D[MAX_NUM_OUTPUTS*MAX_NUM_OUTPUTS];
    float_complex Cproto[MAX_NUM_OUTPUTS*MAX_NUM_OUTPUTS];
    SPREADER_PROC_MODES procMode;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    /* copy user parameters to local variables */
    procMode = pData->procMode;
    nSources = pData->nSources;
    Q = pData->Q;

    /* apply binaural panner */
    if ((nSamples == FRAME_SIZE) && (pData->codecStatus==CODEC_STATUS_INITIALISED) ){
        pData->procStatus = PROC_STATUS_ONGOING;

        /* Load time-domain data */
        for(i=0; i < MIN(nSources,nInputs); i++)
            utility_svvcopy(inputs[i], FRAME_SIZE, pData->inputFrameTD[i]);
        for(; i<nSources; i++)
            memset(pData->inputFrameTD[i], 0, FRAME_SIZE * sizeof(float));

        /* Apply time-frequency transform (TFT) */
        afSTFT_forward(pData->hSTFT, pData->inputFrameTD, FRAME_SIZE, pData->inputframeTF);

        /* Zero output buffer */
        for(band=0; band<HYBRID_BANDS; band++)
            memset(FLATTEN2D(pData->outputframeTF[band]), 0, Q*TIME_SLOTS*sizeof(float_complex));

        /* Loop over sources */
        assert(nSources==1); // need to add support for multi-source
        for(src=0; src<nSources; src++){
            /* Find the "spread" indices */
            unitSph2cart(pData->src_dirs_deg[src], 1, 1, src_dir_xyz);
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, pData->nGrid, 3, 1.0f,
                        src_dir_xyz, 3,
                        pData->grid_dirs_xyz, 3, 0.0f,
                        pData->angles, 1);
            for(i=0; i<pData->nGrid; i++)
                pData->angles[i] = acosf(MIN(pData->angles[i], 0.999999f));
            utility_siminv(pData->angles, pData->nGrid, &centre_ind);

            /* Define Prototype signals */
            switch(procMode){
                case SPREADER_MODE_NAIVE:/* fall through */
                case SPREADER_MODE_OM:
                    for(band=0; band<HYBRID_BANDS; band++){
                        memset(FLATTEN2D(pData->protoframeTF[band]), 0, Q*TIME_SLOTS*sizeof(float_complex));
                        for(ng=0,j=0; ng<pData->nGrid; ng++){
                            if(pData->angles[ng] <= pData->src_spread[ng]/2.0f*SAF_PI/180.0f){
                                for(q=0; q<Q; q++)
                                    H_tmp[q] = pData->H_grid[band*Q*pData->nGrid + q*pData->nGrid + ng];
                                cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Q, TIME_SLOTS, 1, &calpha,
                                            H_tmp, 1,
                                            pData->inputframeTF[band][src], TIME_SLOTS, &cbeta,
                                            tmpFrame, TIME_SLOTS);
                                cblas_caxpy(Q*TIME_SLOTS, &calpha, tmpFrame, 1, FLATTEN2D(pData->protoframeTF[band]), 1);
                                j++;
                            }
                        }
                        if(j==0){  /* If no directions found in the spread area, then just include the nearest one */
                            for(q=0; q<Q; q++)
                                H_tmp[q] = pData->H_grid[band*Q*pData->nGrid + q*pData->nGrid + centre_ind];
                            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Q, TIME_SLOTS, 1, &calpha,
                                        H_tmp, 1,
                                        pData->inputframeTF[band][src], TIME_SLOTS, &cbeta,
                                        tmpFrame, TIME_SLOTS);
                            cblas_caxpy(Q*TIME_SLOTS, &calpha, tmpFrame, 1, FLATTEN2D(pData->protoframeTF[band]), 1);
                        }
                    }
                    break;
                case SPREADER_MODE_EVD:
                    for(band=0; band<HYBRID_BANDS; band++)
                        for(q=0; q<Q; q++)
                            memcpy(pData->protoframeTF[band][q], pData->inputframeTF[band][src], TIME_SLOTS*sizeof(float_complex));
                    break;
            }

            /* Main processing */
            if(procMode==SPREADER_MODE_NAIVE) {
                /* If naive mode, then we're already done... */
                for(band=0; band<HYBRID_BANDS; band++)
                    memcpy(FLATTEN2D(pData->spreadframeTF[band]), FLATTEN2D(pData->protoframeTF[band]), Q*TIME_SLOTS*sizeof(float_complex));
            }
            else{
                /* Apply decorrelation of prototype signals */
                latticeDecorrelator_apply(pData->hDecor, pData->protoframeTF, TIME_SLOTS, pData->decorframeTF);

                /* Define target covariance matrices */
                for(band=0; band<HYBRID_BANDS; band++){
                    /* Sum the H_array outer product matrices for the whole spreading area */
                    memset(Cy_new, 0, Q*Q*sizeof(float_complex));
                    for(ng=0, j=0; ng<pData->nGrid; ng++){
                        if(pData->angles[ng] <= pData->src_spread[src]/2.0f*SAF_PI/180.0f){
                            cblas_caxpy(Q*Q, &calpha, pData->HHH[band][ng], 1, Cy_new, 1);
                            j++;
                        }
                    }
                    if(j==0) /* If no directions found in the spread area, then just include the nearest one */
                        cblas_caxpy(Q*Q, &calpha, pData->HHH[band][centre_ind], 1, Cy_new, 1);

                    /* The OM-solution requires some target energies too */
                    if (procMode == SPREADER_MODE_OM){
                        /* Normalise Cy */
                        trace = 0.0f;
                        for(q=0; q<Q; q++)
                            trace += crealf(Cy_new[q*Q+q]);
                        scaleC = cmplxf(1.0f/(trace+2.23e-9f), 0.0f);
                        cblas_cscal(Q*Q, &scaleC, Cy_new, 1);

                        /* Compute signals for the centre of the spread */
                        for(q=0; q<Q; q++)
                            H_tmp[q] = pData->H_grid[band*Q*pData->nGrid + q*pData->nGrid + centre_ind];
                        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Q, TIME_SLOTS, 1, &calpha,
                                    H_tmp, 1,
                                    pData->inputframeTF[band][src], TIME_SLOTS, &cbeta,
                                    tmpFrame, TIME_SLOTS);

                        /* Introduce their channel energies into the target */
                        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Q, Q, TIME_SLOTS, &calpha,
                                    tmpFrame, TIME_SLOTS,
                                    tmpFrame, TIME_SLOTS, &cbeta,
                                    E_dir, Q);
                        trace = 0.0f;
                        for(q=0; q<Q; q++)
                            trace += crealf(E_dir[q*Q+q]);
                        scaleC = cmplxf(trace, 0.0f);
                        cblas_cscal(Q*Q, &scaleC, Cy_new, 1);
                    }

                    /* Average over time */
                    scaleC = cmplxf(pData->covAvgCoeff, 0.0f);
                    cblas_cscal(Q*Q, &scaleC, pData->Cy[band], 1);
                    scaleC = cmplxf(1.0f-pData->covAvgCoeff, 0.0f);
                    cblas_caxpy(Q*Q, &scaleC, Cy_new, 1, pData->Cy[band], 1);
                }

                /* Formulate mixing matrices */
                switch(procMode){
                    case SPREADER_MODE_NAIVE: assert(0); break;
                    case SPREADER_MODE_EVD:
                        for(band=0; band<HYBRID_BANDS; band++){
                            utility_cseig(pData->Cy[band], Q, 1, V, D, NULL);
                            for(i=0; i<Q; i++)
                                for(j=0; j<Q; j++)
                                    D[i*Q+j] = i==j ? csqrtf(D[i*Q+j]) : cmplxf(0.0f, 0.0f);
                            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Q, Q, Q, &calpha,
                                        V, Q,
                                        D, Q, &cbeta,
                                        pData->new_M[band], Q);
                        }
                        break;
                    case SPREADER_MODE_OM:
                        for(band=0; band<HYBRID_BANDS; band++){
                            /* Prototype covariance matrix averaged over time */
                            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, Q, Q, TIME_SLOTS, &calpha,
                                        FLATTEN2D(pData->protoframeTF[band]), TIME_SLOTS,
                                        FLATTEN2D(pData->protoframeTF[band]), TIME_SLOTS, &cbeta,
                                        Cproto, Q);
                            scaleC = cmplxf(pData->covAvgCoeff, 0.0f);
                            cblas_cscal(Q*Q, &scaleC, pData->Cproto[band], 1);
                            scaleC = cmplxf(1.0f-pData->covAvgCoeff, 0.0f);
                            cblas_caxpy(Q*Q, &scaleC, Cproto, 1, pData->Cproto[band], 1);

                            /* Diagonalise and diagonal loading */
                            cblas_ccopy(Q*Q, pData->Cproto[band], 1, Cproto, 1);
                            for(i=0; i<Q; i++){
                                for(j=0; j<Q; j++){
                                    if(i==j)
                                        Cproto[i*Q+i] = craddf(Cproto[i*Q+i], 0.000001f);
                                    CprotoDiag[i*Q+j] = i==j ? crealf(Cproto[i*Q+i]) : 0.0f;
                                }
                            }
                            formulate_M_and_Cr_cmplx(pData->hCdf, Cproto, pData->Cy[band], pData->Qmix_cmplx, 0, 0.2f, pData->new_M[band], pData->Cr_cmplx);
                            for(i=0; i<Q*Q; i++)
                                pData->Cr[i] = crealf(pData->Cr_cmplx[i]);
                            formulate_M_and_Cr(pData->hCdf_res, CprotoDiag, pData->Cr, pData->Qmix, 0, 0.2f, pData->new_Mr[band], NULL);
                        }
                        break;
                }

                /* Apply mixing matrices */
                for(band=0; band<HYBRID_BANDS; band++){
                    for(t=0; t<TIME_SLOTS; t++){
                        scaleC = cmplxf(pData->interpolatorFadeIn[t], 0.0f);
                        utility_cvsmul(pData->new_M[band], &scaleC, Q*Q, pData->interp_M);
                        scaleC = cmplxf(pData->interpolatorFadeOut[t], 0.0f);
                        cblas_caxpy(Q*Q, &scaleC, pData->prev_M[band], 1, pData->interp_M, 1);
                        for(j=0; j<Q; j++)
                            pData->inFrame_t[j] = procMode == SPREADER_MODE_EVD ? pData->decorframeTF[band][j][t] : pData->protoframeTF[band][j][t];
                        for(i=0; i<Q; i++)
                            utility_cvvdot((float_complex*)(&(pData->interp_M[i*Q])), pData->inFrame_t, Q, NO_CONJ, &(pData->spreadframeTF[band][i][t]));
                    }
                    /* Also mix in the residual part */
                    if(procMode == SPREADER_MODE_OM){
                        for(t=0; t<TIME_SLOTS; t++){
                            utility_svsmul(pData->new_Mr[band], &(pData->interpolatorFadeIn[t]), Q*Q, pData->interp_Mr);
                            cblas_saxpy(Q*Q, pData->interpolatorFadeOut[t], pData->prev_Mr[band], 1, pData->interp_Mr, 1);
                            for(j=0; j<Q; j++)
                                pData->inFrame_t[j] = pData->decorframeTF[band][j][t];
                            cblas_scopy(Q*Q, pData->interp_Mr, 1, (float*)pData->interp_Mr_cmplx, 2);
                            for(i=0; i<Q; i++){
                                utility_cvvdot((float_complex*)(&(pData->interp_Mr_cmplx[i*Q])), pData->inFrame_t, Q, NO_CONJ, &tmp);
                                pData->spreadframeTF[band][i][t] = ccaddf(pData->spreadframeTF[band][i][t], tmp);
                            }
                        }
                    }
                }
            }

            /* Add the spread frame to the output frame, then move onto the next source... */
            for(band=0; band<HYBRID_BANDS; band++)
               cblas_caxpy(Q*TIME_SLOTS, &calpha, FLATTEN2D(pData->spreadframeTF[band]), 1, FLATTEN2D(pData->outputframeTF[band]), 1);

            /* For next frame */
            cblas_ccopy(HYBRID_BANDS*Q*Q, FLATTEN2D(pData->new_M), 1, FLATTEN2D(pData->prev_M), 1);
            cblas_scopy(HYBRID_BANDS*Q*Q, FLATTEN2D(pData->new_Mr), 1, FLATTEN2D(pData->prev_Mr), 1);
        }

        /* inverse-TFT */
        afSTFT_backward(pData->hSTFT, pData->outputframeTF, FRAME_SIZE, pData->outframeTD);

        /* Copy to output buffer */
        for (ch = 0; ch < MIN(Q, nOutputs); ch++)
            utility_svvcopy(pData->outframeTD[ch], FRAME_SIZE, outputs[ch]);
        for (; ch < nOutputs; ch++)
            memset(outputs[ch], 0, FRAME_SIZE*sizeof(float));
    }
    else{
        for (ch=0; ch < nOutputs; ch++)
            memset(outputs[ch],0, FRAME_SIZE*sizeof(float));
    }

    pData->procStatus = PROC_STATUS_NOT_ONGOING;
}

/* Set Functions */

void spreader_refreshSettings(void* const hSpr)
{
    spreader_setCodecStatus(hSpr, CODEC_STATUS_NOT_INITIALISED);
}

void spreader_setSourceAzi_deg(void* const hSpr, int index, float newAzi_deg)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    if(newAzi_deg>180.0f)
        newAzi_deg = -360.0f + newAzi_deg;
    newAzi_deg = MAX(newAzi_deg, -180.0f);
    newAzi_deg = MIN(newAzi_deg, 180.0f);
    if(pData->src_dirs_deg[index][0]!=newAzi_deg){
        pData->src_dirs_deg[index][0] = newAzi_deg;
    }
}

void spreader_setSourceElev_deg(void* const hSpr, int index, float newElev_deg)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    newElev_deg = MAX(newElev_deg, -90.0f);
    newElev_deg = MIN(newElev_deg, 90.0f);
    if(pData->src_dirs_deg[index][1] != newElev_deg){
        pData->src_dirs_deg[index][1] = newElev_deg;
    }
}

void spreader_setNumSources(void* const hSpr, int new_nSources)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    pData->new_nSources = CLAMP(new_nSources, 1, SPREADER_MAX_NUM_SOURCES);
    spreader_setCodecStatus(hSpr, CODEC_STATUS_NOT_INITIALISED);
}

void spreader_setUseDefaultHRIRsflag(void* const hSpr, int newState)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    if((!pData->useDefaultHRIRsFLAG) && (newState)){
        pData->useDefaultHRIRsFLAG = newState;
        spreader_refreshSettings(hSpr);  // re-init and re-calc
    }
}

void spreader_setSofaFilePath(void* const hSpr, const char* path)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    
    pData->sofa_filepath = malloc1d(strlen(path) + 1);
    strcpy(pData->sofa_filepath, path);
    pData->useDefaultHRIRsFLAG = 0;
    spreader_refreshSettings(hSpr);  // re-init and re-calc
}


/* Get Functions */

int spreader_getFrameSize(void)
{
    return FRAME_SIZE;
}

CODEC_STATUS spreader_getCodecStatus(void* const hSpr)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    return pData->codecStatus;
}

float spreader_getProgressBar0_1(void* const hSpr)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    return pData->progressBar0_1;
}

void spreader_getProgressBarText(void* const hSpr, char* text)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    memcpy(text, pData->progressBarText, PROGRESSBARTEXT_CHAR_LENGTH*sizeof(char));
}

float spreader_getSourceAzi_deg(void* const hSpr, int index)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    return pData->src_dirs_deg[index][0];
}

float spreader_getSourceElev_deg(void* const hSpr, int index)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    return pData->src_dirs_deg[index][1];
}

int spreader_getNumSources(void* const hSpr)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    return pData->new_nSources;
}

int spreader_getMaxNumSources()
{
    return MAX_NUM_INPUTS;
}

int spreader_getNumOutputs(void* const hSpr)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    return pData->Q;
}

int spreader_getNDirs(void* const hSpr)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    return pData->nGrid;
}

float spreader_getIRAzi_deg(void* const hSpr, int index)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    if(pData->grid_dirs_deg!=NULL)
        return pData->grid_dirs_deg[index*2+0];
    else
        return 0.0f;
}

float spreader_getIRElev_deg(void* const hSpr, int index)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    if(pData->grid_dirs_deg!=NULL)
        return pData->grid_dirs_deg[index*2+1];
    else
        return 0.0f;
}

int spreader_getIRlength(void* const hSpr)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    return pData->h_len;
}

int spreader_getIRsamplerate(void* const hSpr)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    return pData->h_fs;
}

int spreader_getUseDefaultHRIRsflag(void* const hSpr)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    return pData->useDefaultHRIRsFLAG;
}

char* spreader_getSofaFilePath(void* const hSpr)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    if(pData->sofa_filepath!=NULL)
        return pData->sofa_filepath;
    else
        return "no_file";
}


int spreader_getDAWsamplerate(void* const hSpr)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    return pData->fs;
}

int spreader_getProcessingDelay()
{
    return 12*HOP_SIZE;
}
