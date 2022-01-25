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
 * @file array2sh_internal.c
 * @brief Spatially encodes spherical microphone array signals into spherical
 *        harmonic signals (aka: Ambisonic signals) utilising theoretical
 *        encoding filters.
 *
 * The algorithms within array2sh were pieced together and developed in
 * collaboration with Symeon Delikaris-Manias and Angelo Farina.
 * A detailed explanation of the algorithms within array2sh can be found in [1].
 * Also included, is a diffuse-field equalisation option for frequencies past
 * aliasing, developed in collaboration with Archontis Politis, 8.02.2019
 *
 * @note Since the algorithms are based on theory, only array designs where
 *       there are analytical solutions available are supported. i.e. only
 *       spherical or cylindrical arrays, which have phase-matched sensors.
 *       For more information, the reader is referred to [2,3].
 * @test test__saf_example_array2sh()
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., "Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis," in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 * @see [2] Williams EG. Fourier acoustics: sound radiation and nearfield
 *          acoustical holography. Elsevier; 1999 Jun 10.
 * @see [3] Rafaely B. Fundamentals of spherical array processing. Berlin:
 *          Springer; 2015 Feb 18.
 *
 * @author Leo McCormack
 * @date 13.09.2017
 * @license ISC
 */

#include "array2sh_internal.h"

/**
 * Takes the bNs computed up to N+1, and replicates them to be of length
 * (N+1)^2 (replicating the 1st order bNs 3 times, 2nd -> 5 times etc.)
 */
static void array2sh_replicate_order
(
    void* const hA2sh,
    int order
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    int band, n, i;
    int o[MAX_SH_ORDER+2];
    
    for(n=0; n<order+2; n++)
        o[n] = n*n;
    for(band=0; band<HYBRID_BANDS; band++)
        for(n=0; n < order+1; n++)
            for(i=o[n]; i < o[n+1]; i++)
                pData->bN_inv_R[band][i] = pData->bN_inv[band][n];
}

void array2sh_initTFT
(
    void* const hA2sh
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    int new_nSH, nSH;
    
    new_nSH = (pData->new_order+1)*(pData->new_order+1);
    nSH = (pData->order+1)*(pData->order+1);
    if(pData->hSTFT==NULL)
        afSTFT_create(&(pData->hSTFT), arraySpecs->newQ, new_nSH, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    else if(arraySpecs->newQ != arraySpecs->Q || nSH != new_nSH){
        afSTFT_channelChange(pData->hSTFT, arraySpecs->newQ, new_nSH);
        afSTFT_clearBuffers(pData->hSTFT); 
        pData->reinitSHTmatrixFLAG = 1; /* filters will need to be updated too */
    }
    arraySpecs->Q = arraySpecs->newQ;
}

void array2sh_calculate_sht_matrix
(
    void* const hA2sh
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    int i, j, band, n, order, nSH;
    double alpha, beta, g_lim, regPar;
    double kr[HYBRID_BANDS], kR[HYBRID_BANDS];
    float* Y_mic, *pinv_Y_mic;
    float_complex* pinv_Y_mic_cmplx, *diag_bN_inv_R;
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta  = cmplxf(0.0f, 0.0f);
    
    /* prep */
    order = pData->new_order;
    nSH = (order+1)*(order+1);
    arraySpecs->R = SAF_MIN(arraySpecs->R, arraySpecs->r);
    for(band=0; band<HYBRID_BANDS; band++){
        kr[band] = 2.0*SAF_PId*(pData->freqVector[band])*(arraySpecs->r)/pData->c;
        kR[band] = 2.0*SAF_PId*(pData->freqVector[band])*(arraySpecs->R)/pData->c;
    }
    
    /* Spherical harmponic weights for each sensor direction */
    Y_mic = malloc1d(nSH*(arraySpecs->Q)*sizeof(float));
    getRSH(order, (float*)arraySpecs->sensorCoords_deg, arraySpecs->Q, Y_mic); /* nSH x Q */
    pinv_Y_mic = malloc1d( arraySpecs->Q * nSH *sizeof(float));
    utility_spinv(NULL, Y_mic, nSH, arraySpecs->Q, pinv_Y_mic);
    pinv_Y_mic_cmplx =  malloc1d((arraySpecs->Q) * nSH *sizeof(float_complex));
    for(i=0; i<(arraySpecs->Q)*nSH; i++)
        pinv_Y_mic_cmplx[i] = cmplxf(pinv_Y_mic[i], 0.0f);
    
    /* ------------------------------------------------------------------------------ */
    /* Encoding filters based on the regularised inversion of the modal coefficients: */
    /* ------------------------------------------------------------------------------ */
    if ( (pData->filterType==FILTER_SOFT_LIM) || (pData->filterType==FILTER_TIKHONOV) ){
        /* Compute modal responses */
        free(pData->bN);
        pData->bN = malloc1d((HYBRID_BANDS)*(order+1)*sizeof(double_complex));
        switch(arraySpecs->arrayType){
            case ARRAY_CYLINDRICAL:
                switch (arraySpecs->weightType){
                    case WEIGHT_RIGID_OMNI:   cylModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_RIGID, pData->bN); break;
                    case WEIGHT_RIGID_CARD:   saf_print_error("weightType is not supported"); break;
                    case WEIGHT_RIGID_DIPOLE: saf_print_error("weightType is not supported"); break;
                    case WEIGHT_OPEN_OMNI:    cylModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN, pData->bN);  break;
                    case WEIGHT_OPEN_CARD:    saf_print_error("weightType is not supported"); break;
                    case WEIGHT_OPEN_DIPOLE:  saf_print_error("weightType is not supported"); break;
                }
                break;
            case ARRAY_SPHERICAL:
                switch (arraySpecs->weightType){
                    case WEIGHT_OPEN_OMNI:   sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN, 1.0, pData->bN); break;
                    case WEIGHT_OPEN_CARD:   sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.5, pData->bN); break;
                    case WEIGHT_OPEN_DIPOLE: sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.0, pData->bN); break;
                    case WEIGHT_RIGID_OMNI:
                    case WEIGHT_RIGID_CARD:
                    case WEIGHT_RIGID_DIPOLE:
                        /* if sensors are flushed with the rigid baffle: */
                        if(arraySpecs->R == arraySpecs->r )
                            sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_RIGID, 1.0, pData->bN);

                        /* if sensors protrude from the rigid baffle: */
                        else{
                            if (arraySpecs->weightType == WEIGHT_RIGID_OMNI)
                                sphScattererModalCoeffs(order, kr, kR, HYBRID_BANDS, pData->bN);
                            else if (arraySpecs->weightType == WEIGHT_RIGID_CARD)
                                sphScattererDirModalCoeffs(order, kr, kR, HYBRID_BANDS, 0.5, pData->bN);
                            else if (arraySpecs->weightType == WEIGHT_RIGID_DIPOLE)
                                sphScattererDirModalCoeffs(order, kr, kR, HYBRID_BANDS, 0.0, pData->bN);
                        }
                        break;
                }
                break;
        }
        
        for(band=0; band<HYBRID_BANDS; band++)
            for(n=0; n < order+1; n++)
                pData->bN[band*(order+1)+n] = ccdiv(pData->bN[band*(order+1)+n], cmplx(4.0*SAF_PId, 0.0)); /* 4pi term */

        /* direct inverse */
        regPar = pData->regPar;
        for(band=0; band<HYBRID_BANDS; band++)
            for(n=0; n < order+1; n++)
                pData->bN_modal[band][n] = ccdiv(cmplx(1.0,0.0), (pData->bN[band*(order+1)+n]));
        
        /* regularised inverse */
        if (pData->filterType == FILTER_SOFT_LIM){
            /* Bernschutz, B., Porschmann, C., Spors, S., Weinzierl, S., Versterkung, B., 2011. Soft-limiting der
             modalen amplitudenverst?rkung bei sph?rischen mikrofonarrays im plane wave decomposition verfahren.
             Proceedings of the 37. Deutsche Jahrestagung fur Akustik (DAGA 2011) */
            g_lim = sqrt(arraySpecs->Q)*pow(10.0,(regPar/20.0));
            for(band=0; band<HYBRID_BANDS; band++)
                for(n=0; n < order+1; n++)
                    pData->bN_inv[band][n] = crmul(pData->bN_modal[band][n], (2.0*g_lim*cabs(pData->bN[band*(order+1)+n]) / SAF_PId)
                                                     * atan(SAF_PId / (2.0*g_lim*cabs(pData->bN[band*(order+1)+n]))) );
        }
        else if(pData->filterType == FILTER_TIKHONOV){
            /* Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording with higher order ambisonics-objective
             measurements and validation of spherical microphone. In Audio Engineering Society Convention 120. */
            alpha = sqrt(arraySpecs->Q)*pow(10.0,(regPar/20.0));
            for(band=0; band<HYBRID_BANDS; band++){
                for(n=0; n < order+1; n++){
                    beta = sqrt((1.0-sqrt(1.0-1.0/ pow(alpha,2.0)))/(1.0+sqrt(1.0-1.0/pow(alpha,2.0))));
                    pData->bN_inv[band][n] = ccdiv(conj(pData->bN[band*(order+1)+n]), cmplx((pow(cabs(pData->bN[band*(order+1)+n]), 2.0) + pow(beta, 2.0)),0.0));
                }
            }
        }
        
        /* diag(filters) * Y */
        array2sh_replicate_order(hA2sh, order); /* replicate orders */
        
        diag_bN_inv_R = calloc1d(nSH*nSH, sizeof(float_complex));
        for(band=0; band<HYBRID_BANDS; band++){
            for(i=0; i<nSH; i++)
                diag_bN_inv_R[i*nSH+i] = cmplxf((float)creal(pData->bN_inv_R[band][i]), (float)cimag(pData->bN_inv_R[band][i]));
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, (arraySpecs->Q), nSH, &calpha,
                        diag_bN_inv_R, nSH,
                        pinv_Y_mic_cmplx, nSH, &cbeta,
                        pData->W[band], MAX_NUM_SENSORS);
        }
        free(diag_bN_inv_R);
    }
    
    /* ------------------------------------------------------------- */
    /* Encoding filters based on a linear-phase filter-bank approach */
    /* ------------------------------------------------------------- */
    else if ( (pData->filterType==FILTER_Z_STYLE) || (pData->filterType==FILTER_Z_STYLE_MAXRE) ) {
        /* Zotter, F. A Linear-Phase Filter-Bank Approach to Process Rigid Spherical Microphone Array Recordings. */
        double normH;
        float f_lim[MAX_SH_ORDER+1];
        double H[HYBRID_BANDS][MAX_SH_ORDER+1];
        double_complex Hs[HYBRID_BANDS][MAX_SH_ORDER+1];
        
        /* find suitable cut-off frequencies */
        switch (arraySpecs->weightType){
            case WEIGHT_OPEN_OMNI:   sphArrayNoiseThreshold(order, arraySpecs->Q, arraySpecs->r, pData->c, ARRAY_CONSTRUCTION_OPEN, 1.0, pData->regPar, f_lim); break;
            case WEIGHT_OPEN_CARD:   sphArrayNoiseThreshold(order, arraySpecs->Q, arraySpecs->r, pData->c, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.5, pData->regPar, f_lim); break;
            case WEIGHT_OPEN_DIPOLE: sphArrayNoiseThreshold(order, arraySpecs->Q, arraySpecs->r, pData->c, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.0, pData->regPar, f_lim); break;
            case WEIGHT_RIGID_OMNI:
            case WEIGHT_RIGID_CARD:
            case WEIGHT_RIGID_DIPOLE:
                /* Currently no support for estimating the noise cut-off frequencies for rigid scatterers. */
                sphArrayNoiseThreshold(order, arraySpecs->Q, arraySpecs->r, pData->c, ARRAY_CONSTRUCTION_RIGID, 1.0, pData->regPar, f_lim); break;
        }
        
        /* design prototype filterbank */
        for(band=0; band<HYBRID_BANDS; band++){
            normH = 0.0;
            for (n=0; n<order+1; n++){
                if (n==0)
                    H[band][n] = 1.0/(1.0+ pow((double)(pData->freqVector[band]/f_lim[n]),2.0));
                else if (n==order)
                    H[band][n] = pow((double)(pData->freqVector[band]/f_lim[n-1]), (double)order+1.0 )  /
                                 (1.0 + pow((double)(pData->freqVector[band]/f_lim[n-1]), (double)order+1.0));
                else
                    H[band][n] = pow((double)(pData->freqVector[band]/f_lim[n-1]), (double)n+1.0 )  /
                                 (1.0 + pow((double)(pData->freqVector[band]/f_lim[n-1]), (double)n+1.0)) *
                                 (1.0 / (1.0 + pow((double)(pData->freqVector[band]/f_lim[n]), (double)n+2.0)));
                normH += H[band][n];
            }
            /* normalise */
            for (n=0; n<order+1; n++)
                H[band][n] = H[band][n]/normH;
        }
                
        /* compute inverse radial response */ 
        free(pData->bN);
        pData->bN = malloc1d((HYBRID_BANDS)*(order+1)*sizeof(double_complex));
        switch(arraySpecs->arrayType){
            case ARRAY_CYLINDRICAL:
                switch (arraySpecs->weightType){
                    case WEIGHT_RIGID_OMNI:   cylModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_RIGID, pData->bN); break;
                    case WEIGHT_RIGID_CARD:   /* not supported */ break;
                    case WEIGHT_RIGID_DIPOLE: /* not supported */ break;
                    case WEIGHT_OPEN_OMNI:    cylModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN, pData->bN);  break;
                    case WEIGHT_OPEN_CARD:    /* not supported */ break;
                    case WEIGHT_OPEN_DIPOLE:  /* not supported */ break;
                }
                break;
            case ARRAY_SPHERICAL:
                switch (arraySpecs->weightType){
                    case WEIGHT_OPEN_OMNI:   sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN, 1.0, pData->bN); break;
                    case WEIGHT_OPEN_CARD:   sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.5, pData->bN); break;
                    case WEIGHT_OPEN_DIPOLE: sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.0, pData->bN); break;
                    case WEIGHT_RIGID_OMNI:
                    case WEIGHT_RIGID_CARD:
                    case WEIGHT_RIGID_DIPOLE:
                        /* if sensors are flushed with the rigid baffle: */
                        if(arraySpecs->R == arraySpecs->r )
                            sphModalCoeffs(order, kr, HYBRID_BANDS, ARRAY_CONSTRUCTION_RIGID, 1.0, pData->bN);
                        
                        /* if sensors protrude from the rigid baffle: */
                        else{
                            if (arraySpecs->weightType == WEIGHT_RIGID_OMNI)
                                sphScattererModalCoeffs(order, kr, kR, HYBRID_BANDS, pData->bN);
                            else if (arraySpecs->weightType == WEIGHT_RIGID_CARD)
                                sphScattererDirModalCoeffs(order, kr, kR, HYBRID_BANDS, 0.5, pData->bN);
                            else if (arraySpecs->weightType == WEIGHT_RIGID_DIPOLE)
                                sphScattererDirModalCoeffs(order, kr, kR, HYBRID_BANDS, 0.0, pData->bN);
                        }
                        break;
                }
                break;
        }
        
        /* direct inverse (only required for GUI) */
        for(band=0; band<HYBRID_BANDS; band++)
            for(n=0; n < order+1; n++)
                pData->bN_modal[band][n] = ccdiv(cmplx(4.0*SAF_PId, 0.0), pData->bN[band*(order+1)+n]);

        /* phase shift */
        for(band=0; band<HYBRID_BANDS; band++)
            for (n=0; n<order+1; n++)
                Hs[band][n] = ccmul(cexp(cmplx(0.0, kr[band])), ccdiv(cmplx(4.0*SAF_PId, 0.0), pData->bN[band*(order+1)+n]));
        
        /* apply max-re order weighting and diffuse equalisation (not the same as "array2sh_apply_diff_EQ") */
        float* wn;
        double W[MAX_SH_ORDER+1][MAX_SH_ORDER+1];
        double EN, scale;
        int nSH_n;
        memset(W, 0, (MAX_SH_ORDER+1)*(MAX_SH_ORDER+1)*sizeof(double));
        for (n=0; n<order+1; n++){
            nSH_n = (n+1)*(n+1);
            wn = calloc1d(nSH_n*nSH_n, sizeof(float));
            if(pData->filterType==FILTER_Z_STYLE)
                for (i=0; i<n+1; i++)
                    wn[(i*i)*nSH_n+(i*i)] = 1.0f;
            else if(pData->filterType==FILTER_Z_STYLE_MAXRE)
                getMaxREweights(n, 1, wn);
            scale = 0.0;
            for (i=0; i<n+1; i++)
                scale += (double)(2*i+1)*pow((double)wn[(i*i)*nSH_n + (i*i)], 2.0);
            for (i=0; i<n+1; i++)
                W[i][n] = (double)wn[(i*i)*nSH_n + (i*i)]/ sqrt(scale);
            free(wn);
        }
        EN=W[0][n-1];
        for (n=0; n<order+1; n++)
            for (i=0; i<order+1; i++)
                W[i][n] /= EN;
        
        /* apply bandpass filterbank to the inverse array response to regularise it */
        double HW[HYBRID_BANDS];
        double H_np[HYBRID_BANDS][MAX_SH_ORDER+1];
        double W_np[MAX_SH_ORDER+1];
        for (n=0; n<order+1; n++){
            for(band=0; band< HYBRID_BANDS; band++)
                for (i=n, j=0; i<order+1; i++, j++)
                    H_np[band][j] = H[band][i];
            for (i=n, j=0; i<order+1; i++, j++)
                W_np[j] = W[n][i];
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, HYBRID_BANDS, 1, order+1-n, 1.0,
                        (const double*)H_np, MAX_SH_ORDER+1,
                        (const double*)W_np, MAX_SH_ORDER+1, 0.0,
                        (double*)HW, 1);
            for(band=0; band<HYBRID_BANDS; band++)
                pData->bN_inv[band][n] = crmul(Hs[band][n], HW[band]);
        }
        
        /* diag(filters) * Y */
        array2sh_replicate_order(hA2sh, order); /* replicate orders */
        diag_bN_inv_R = calloc1d(nSH*nSH, sizeof(float_complex));
        for(band=0; band<HYBRID_BANDS; band++){
            for(i=0; i<nSH; i++)
                diag_bN_inv_R[i*nSH+i] = cmplxf((float)creal(pData->bN_inv_R[band][i]), (float)cimag(pData->bN_inv_R[band][i])); /* double->single */
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, (arraySpecs->Q), nSH, &calpha,
                        diag_bN_inv_R, nSH,
                        pinv_Y_mic_cmplx, nSH, &cbeta,
                        pData->W[band], MAX_NUM_SENSORS);
        }
        free(diag_bN_inv_R);
        
    }
     
    pData->order = order;
    
    if(pData->enableDiffEQpastAliasing)
        array2sh_apply_diff_EQ(hA2sh);
    
    free(Y_mic);
    free(pinv_Y_mic);
    free(pinv_Y_mic_cmplx);
}

/* Based on a MatLab script by Archontis Politis, 2019 */
void array2sh_apply_diff_EQ(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    int i, j, band, array_order, idxf_alias, nSH;
    float f_max, kR_max, f_alias, f_f_alias;
    double_complex* dM_diffcoh_s;
    const double_complex calpha = cmplx(1.0, 0.0); const double_complex cbeta  = cmplx(0.0, 0.0);
    double kr[HYBRID_BANDS];
    double_complex L_diff_fal[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];
    double_complex L_diff[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];
    double_complex E_diff[MAX_NUM_SH_SIGNALS][MAX_NUM_SENSORS];
    double_complex W_diffEQ[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];
    double_complex W_tmp[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];
    double* dM_diffcoh; 
    
    if(arraySpecs->arrayType==ARRAY_CYLINDRICAL)
        return; /* unsupported */
    
    /* prep */
    nSH = (pData->order+1)*(pData->order+1);
    dM_diffcoh = malloc1d((arraySpecs->Q)*(arraySpecs->Q)* (HYBRID_BANDS) * sizeof(double_complex));
    dM_diffcoh_s = malloc1d((arraySpecs->Q)*(arraySpecs->Q) * sizeof(double_complex));
    f_max = 20e3f;
    kR_max = 2.0f*SAF_PI*f_max*(arraySpecs->r)/pData->c;
    array_order = SAF_MIN((int)(ceilf(2.0f*kR_max)+0.01f), 28); /* Cap at around 28, as Bessels at 30+ can be numerically unstable */
    for(band=0; band<HYBRID_BANDS; band++)
        kr[band] = 2.0*SAF_PId*(pData->freqVector[band])*(arraySpecs->r)/pData->c;
    
    /* Get theoretical diffuse coherence matrix */
    switch(arraySpecs->arrayType){
        case ARRAY_CYLINDRICAL:
            return; /* Unsupported */
            break;
        case ARRAY_SPHERICAL:
            switch (arraySpecs->weightType){
                case WEIGHT_RIGID_OMNI: /* Does not handle the case where kr != kR ! */
                    sphDiffCohMtxTheory(array_order, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, ARRAY_CONSTRUCTION_RIGID, 1.0, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
                case WEIGHT_RIGID_CARD:
                    sphDiffCohMtxTheory(array_order, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL, 0.5, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
                case WEIGHT_RIGID_DIPOLE:
                    sphDiffCohMtxTheory(array_order, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL, 0.0, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
                case WEIGHT_OPEN_OMNI:
                    sphDiffCohMtxTheory(array_order, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, ARRAY_CONSTRUCTION_OPEN, 1.0, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
                case WEIGHT_OPEN_CARD:
                    sphDiffCohMtxTheory(array_order, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.5, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
                case WEIGHT_OPEN_DIPOLE:
                    sphDiffCohMtxTheory(array_order, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.0, kr, HYBRID_BANDS, dM_diffcoh);
                    break;
            }
            break;
    }
    
    /* determine band index for the spatial aliasing limit */
    f_alias = sphArrayAliasLim(arraySpecs->r, pData->c, pData->order);
    idxf_alias = 1;
    f_f_alias = 1e13f;
    for(band=0; band<HYBRID_BANDS; band++){
        if( fabsf(pData->freqVector[band]-f_alias) < f_f_alias){
            f_f_alias = fabsf(pData->freqVector[band]-f_alias);
            idxf_alias = band;
        }
    }
    
    /* baseline */
    for(i=0; i<arraySpecs->Q; i++)
        for(j=0; j<arraySpecs->Q; j++)
            dM_diffcoh_s[i*(arraySpecs->Q)+j] = cmplx(dM_diffcoh[i*(arraySpecs->Q)* (HYBRID_BANDS) + j*(HYBRID_BANDS) + (idxf_alias)], 0.0);
    for(i=0; i<nSH; i++)
        for(j=0; j<arraySpecs->Q; j++)
            W_tmp[i][j]= cmplx((double)crealf(pData->W[idxf_alias][i][j]), (double)cimagf(pData->W[idxf_alias][i][j]));
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, (arraySpecs->Q), (arraySpecs->Q), &calpha,
                W_tmp, MAX_NUM_SENSORS,
                dM_diffcoh_s, (arraySpecs->Q), &cbeta,
                E_diff, MAX_NUM_SENSORS);
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, nSH, (arraySpecs->Q), &calpha,
                E_diff, MAX_NUM_SENSORS,
                W_tmp, MAX_NUM_SENSORS, &cbeta,
                L_diff_fal, MAX_NUM_SH_SIGNALS);
    for(i=0; i<nSH; i++)
        L_diff_fal[i][i] = crmul(L_diff_fal[i][i], 1.0/(4.0*SAF_PId)); /* only care about the diagonal entries */
    
    /* diffuse-field equalise bands above aliasing. */
    for(band = SAF_MAX(idxf_alias,0)+1; band<HYBRID_BANDS; band++){
        for(i=0; i<arraySpecs->Q; i++)
            for(j=0; j<arraySpecs->Q; j++)
                dM_diffcoh_s[i*(arraySpecs->Q)+j] = cmplx(dM_diffcoh[i*(arraySpecs->Q)* (HYBRID_BANDS) + j*(HYBRID_BANDS) + (band)], 0.0);
        for(i=0; i<nSH; i++)
            for(j=0; j<arraySpecs->Q; j++)
                W_tmp[i][j]= cmplx((double)crealf(pData->W[band][i][j]), (double)cimagf(pData->W[band][i][j]));
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, (arraySpecs->Q), (arraySpecs->Q), &calpha,
                    W_tmp, MAX_NUM_SENSORS,
                    dM_diffcoh_s, (arraySpecs->Q), &cbeta,
                    E_diff, MAX_NUM_SENSORS);
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nSH, nSH, (arraySpecs->Q), &calpha,
                    E_diff, MAX_NUM_SENSORS,
                    W_tmp, MAX_NUM_SENSORS, &cbeta,
                    L_diff, MAX_NUM_SH_SIGNALS);
        for(i=0; i<nSH; i++)
            for(j=0; j<nSH; j++)
                L_diff[i][j] = i==j? csqrt(cradd(ccdiv(L_diff_fal[i][j], crmul(L_diff[i][j], 1.0/(4.0*SAF_PId))), 2.23e-10)): cmplx(0.0,0.0);
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, (arraySpecs->Q), nSH, &calpha,
                    L_diff, MAX_NUM_SH_SIGNALS,
                    W_tmp, MAX_NUM_SENSORS, &cbeta,
                    W_diffEQ, MAX_NUM_SENSORS);
        for(i=0; i<nSH; i++)
            for(j=0; j<arraySpecs->Q; j++)
                pData->W[band][i][j] = cmplxf((float)creal(W_diffEQ[i][j]), (float)cimag(W_diffEQ[i][j]));
    }
    
    pData->evalStatus = EVAL_STATUS_NOT_EVALUATED;
    
    free(dM_diffcoh);
    free(dM_diffcoh_s);
}

void array2sh_calculate_mag_curves(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    int band, n;
    
    for(band = 0; band <HYBRID_BANDS; band++){
        for(n = 0; n <pData->order+1; n++){
            pData->bN_inv_dB[band][n] = 20.0f * (float)log10(cabs(pData->bN_inv[band][n]));
            pData->bN_modal_dB[band][n] = 20.0f * (float)log10(cabs(pData->bN_modal[band][n]));
        }
    }
}

void array2sh_evaluateSHTfilters(void* hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    array2sh_arrayPars* arraySpecs = (array2sh_arrayPars*)(pData->arraySpecs);
    int band, i, j, simOrder, order, nSH;
    double kr[HYBRID_BANDS];
    double kR[HYBRID_BANDS];
    float* Y_grid_real;
    float_complex* Y_grid, *H_array, *Wshort;
     
    saf_assert(pData->W != NULL, "The initCodec function must have been called prior to calling array2sh_evaluateSHTfilters()");
    
    strcpy(pData->progressBarText,"Simulating microphone array");
    pData->progressBar0_1 = 0.35f;
    
    /* simulate the current array by firing 812 plane-waves around the surface of a theoretical version of the array
     * and ascertaining the transfer function for each */
    simOrder = (int)(2.0f*SAF_PI*MAX_EVAL_FREQ_HZ*(arraySpecs->r)/pData->c)+1;
    for(band=0; band<HYBRID_BANDS; band++){
        kr[band] = 2.0*SAF_PId*(pData->freqVector[band])*(arraySpecs->r)/pData->c;
        kR[band] = 2.0*SAF_PId*(pData->freqVector[band])*(arraySpecs->R)/pData->c;
    }
    H_array = malloc1d((HYBRID_BANDS) * (arraySpecs->Q) * 812*sizeof(float_complex));
    switch(arraySpecs->arrayType){
        case ARRAY_SPHERICAL:
            switch(arraySpecs->weightType){
                default:
                case WEIGHT_RIGID_OMNI:
                    simulateSphArray(simOrder, kr, kR, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID, 1.0, H_array);
                    break;
                case WEIGHT_RIGID_CARD:
                    simulateSphArray(simOrder, kr, kR, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL, 0.5, H_array);
                    break;
                case WEIGHT_RIGID_DIPOLE:
                    simulateSphArray(simOrder, kr, kR, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL, 0.0, H_array);
                    break;
                case WEIGHT_OPEN_OMNI:
                    simulateSphArray(simOrder, kr, NULL, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN, 1.0, H_array);
                    break;
                case WEIGHT_OPEN_CARD:
                    simulateSphArray(simOrder, kr, NULL, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.5, H_array);
                    break;
                case WEIGHT_OPEN_DIPOLE:
                    simulateSphArray(simOrder, kr, NULL, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q,
                                     (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, 0.0, H_array);
                    break;
            }
            break;
            
        case ARRAY_CYLINDRICAL:
            switch(arraySpecs->weightType){
                default:
                case WEIGHT_RIGID_OMNI:
                case WEIGHT_RIGID_CARD:
                case WEIGHT_RIGID_DIPOLE:
                    simulateCylArray(simOrder, kr, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID, H_array);
                    break;
                case WEIGHT_OPEN_DIPOLE:
                case WEIGHT_OPEN_CARD:
                case WEIGHT_OPEN_OMNI:
                    simulateCylArray(simOrder, kr, HYBRID_BANDS, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN, H_array);
                    break;
            }
            break;
    }
    
    strcpy(pData->progressBarText,"Evaluating encoding performance");
    pData->progressBar0_1 = 0.8f;
    
    /* generate ideal (real) spherical harmonics to compare with */
    order = pData->order;
    nSH = (order+1)*(order+1);
    Y_grid_real = malloc1d(nSH*812*sizeof(float));
    getRSH(order, (float*)__geosphere_ico_9_0_dirs_deg, 812, Y_grid_real);
    Y_grid = malloc1d(nSH*812*sizeof(float_complex));
    for(i=0; i<nSH*812; i++)
        Y_grid[i] = cmplxf(Y_grid_real[i], 0.0f); /* "evaluateSHTfilters" function requires complex data type */
    
    /* compare the spherical harmonics obtained from encoding matrix 'W' with the ideal patterns */
    Wshort = malloc1d(HYBRID_BANDS*nSH*(arraySpecs->Q)*sizeof(float_complex));
    for(band=0; band<HYBRID_BANDS; band++)
        for(i=0; i<nSH; i++)
            for(j=0; j<(arraySpecs->Q); j++)
                Wshort[band*nSH*(arraySpecs->Q) + i*(arraySpecs->Q) + j] = pData->W[band][i][j];
    evaluateSHTfilters(order, Wshort, arraySpecs->Q, HYBRID_BANDS, H_array, 812, Y_grid, pData->cSH, pData->lSH);

    free(Y_grid_real);
    free(Y_grid);
    free(H_array);
    free(Wshort);
}

void array2sh_createArray(void ** const hPars)
{
    array2sh_arrayPars* pars = (array2sh_arrayPars*)malloc1d(sizeof(array2sh_arrayPars));
    *hPars = (void*)pars;
}

void array2sh_destroyArray(void ** const hPars)
{
    array2sh_arrayPars *pars = (array2sh_arrayPars*)(*hPars);
    if(pars!=NULL) {
        free(pars);
        pars=NULL;
    }
}
 
void array2sh_initArray
(
    void* const hPars,
    ARRAY2SH_MICROPHONE_ARRAY_PRESETS preset,
    int* arrayOrder,
    int firstInitFlag
)
{
    array2sh_arrayPars *pars = (array2sh_arrayPars*)(hPars);
    int ch, i, Q;
    
    switch(preset){
        default:
        case MICROPHONE_ARRAY_PRESET_DEFAULT:
            (*arrayOrder) = 1;
            Q = 4; /* number of mics */
            pars->r = 0.042f; /* array radius */
            pars->R = 0.042f; /* radius of the sensors (incase they protrude from the surface of the array), (only for rigid arrays) */
            pars->arrayType = ARRAY_SPHERICAL;    /* spherical or cylindrical */
            pars->weightType = WEIGHT_RIGID_OMNI; /* open or rigid, and directivity of the sensors (only for open arrays) */
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __default_SENSORcoords64_rad[ch][i]; /* spherical coordinates of the sensors, in radians */
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_AALTO_HYDROPHONE:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.173f;
            pars->R = 0.173f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_OMNI;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Aalto_Hydrophone_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_SENNHEISER_AMBEO:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.014f;
            pars->R = 0.014f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Sennheiser_Ambeo_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_CORE_SOUND_TETRAMIC:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.02f;
            pars->R = 0.02f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Core_Sound_TetraMic_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_ZOOM_H3VR_PRESET:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.012f;
            pars->R = 0.012f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Zoom_H3VR_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_SOUND_FIELD_SPS200:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.02f;
            pars->R = 0.02f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Sound_field_SPS200_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_ZYLIA_1D:
            (*arrayOrder) = 3;
            Q = 19;
            pars->r = 0.049f;
            pars->R = 0.049f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID_OMNI;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Zylia1D_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_EIGENMIKE32:
            (*arrayOrder) = 4;
            Q = 32;
            pars->r = 0.042f;
            pars->R = 0.042f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID_OMNI;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Eigenmike32_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
        case MICROPHONE_ARRAY_PRESET_DTU_MIC:
            (*arrayOrder) = 6;
            Q = 52;
            pars->r = 0.05f;
            pars->R = 0.05f; 
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID_OMNI;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __DTU_mic_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
                }
            }
            break;
    }
    
    /* Fill remaining slots with default coords */
    for(; ch<MAX_NUM_SENSORS_IN_PRESET; ch++){
        for(i=0; i<2; i++){
            pars->sensorCoords_rad[ch][i] = __default_SENSORcoords64_rad[ch][i];
            pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/SAF_PI);
        }
    }
    
    /* For dynamically changing the number of TFT channels */
    if(firstInitFlag==1){
        pars->Q = Q;
        pars->newQ = pars->Q;
    }
    else
        pars->newQ = Q;
}

