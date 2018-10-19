    /*
 Copyright 2017-2018 Leo McCormack
 
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
 *     array2sh_internal.c
 * Description:
 *     Spatially encodes spherical or cylindrical sensor array signals into spherical harmonic
 *     signals utilising theoretical encoding filters.
 *     The algorithms within array2sh were pieced together and developed in collaboration
 *     with Symeon Delikaris-Manias.
 *     A more detailed explanation of the algorithms in array2sh can be found in:
 *     McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki, V.,
 *     “Real-time conversion of sensor array signals into spherical harmonic signals with
 *     applications to spatially localised sub-band sound-field analysis,” in Audio
 *     Engineering Society Convention 144, Audio Engineering Society, 2018.
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_sh
 * Author, date created:
 *     Leo McCormack, 13.09.2017
 */

#include "array2sh_internal.h"

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
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    
    if(pData->hSTFT==NULL)
        afSTFTinit(&(pData->hSTFT), HOP_SIZE, arraySpecs->newQ, pData->new_nSH, 0, 1);
    else
        afSTFTchannelChange(pData->hSTFT, arraySpecs->newQ, pData->new_nSH);

    arraySpecs->Q = arraySpecs->newQ;
    pData->nSH = pData->new_nSH;
    pData->reinitSHTmatrixFLAG = 1; /* filters need to be updated too */
}

void array2sh_calculate_sht_matrix
(
    void* const hA2sh
)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    int i, band, n, order, nSH;
    double f_n, alpha, beta, g_lim, regPar;
    double kr[HYBRID_BANDS-1], kR[HYBRID_BANDS-1];
    float* Y_mic, *pinv_Y_mic;
    float_complex* pinv_Y_mic_cmplx, *diag_bN_inv_R;
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta  = cmplxf(0.0f, 0.0f);
    
    order = pData->new_order;
    nSH = (order+1)*(order+1);
    arraySpecs->R = MAX(arraySpecs->R, arraySpecs->r);
    
    /* Compute modal coefficients */
    for(band=0; band<HYBRID_BANDS-1; band++)
        kr[band] = 2.0*M_PI*(pData->freqVector[band+1/* ignore DC */])*(arraySpecs->r)/pData->c;
    free(pData->bN);
    pData->bN = malloc((HYBRID_BANDS-1)*(order+1)*sizeof(double_complex));
    switch(arraySpecs->arrayType){
        case ARRAY_CYLINDRICAL:
            switch (arraySpecs->weightType){
                case WEIGHT_RIGID:       cylModalCoeffs(order, kr, HYBRID_BANDS-1, ARRAY_CONSTRUCTION_RIGID, pData->bN); break;
                case WEIGHT_OPEN_OMNI:   cylModalCoeffs(order, kr, HYBRID_BANDS-1, ARRAY_CONSTRUCTION_OPEN, pData->bN);  break;
                case WEIGHT_OPEN_CARD:   /* not supported */ break;
                case WEIGHT_OPEN_DIPOLE: /* not supported */ break;
            }
            break;
        case ARRAY_SPHERICAL:
            switch (arraySpecs->weightType){
                case WEIGHT_OPEN_OMNI:   sphModalCoeffs(order, kr, HYBRID_BANDS-1, ARRAY_CONSTRUCTION_OPEN, 0.0, pData->bN); break;
                case WEIGHT_OPEN_CARD:   sphModalCoeffs(order, kr, HYBRID_BANDS-1, ARRAY_CONSTRUCTION_OPEN, 0.5, pData->bN); break;
                case WEIGHT_OPEN_DIPOLE: sphModalCoeffs(order, kr, HYBRID_BANDS-1, ARRAY_CONSTRUCTION_OPEN, 1.0, pData->bN); break;
                case WEIGHT_RIGID:
                    if(arraySpecs->R == arraySpecs->r )
                        sphModalCoeffs(order, kr, HYBRID_BANDS-1, ARRAY_CONSTRUCTION_RIGID, 1.0, pData->bN);
                    else{
                        for(band=0; band<HYBRID_BANDS-1; band++)
                            kR[band] = 2.0*M_PI*(double)pData->freqVector[band+1] * (double)arraySpecs->R / (double)pData->c;
                        sphScattererModalCoeffs(order, kr, kR, HYBRID_BANDS-1, pData->bN);
                    }
                    break;
            }
            break;
    }
    for(band=0; band<HYBRID_BANDS-1; band++)
        for(n=0; n < order+1; n++)
            pData->bN[band*(order+1)+n] = ccdiv(pData->bN[band*(order+1)+n], cmplx(4.0*M_PI, 0.0f)); /* 4pi term */
    
    /* regularised inversion */
    regPar = pData->regPar;
    for(band=0; band<HYBRID_BANDS-1; band++)
        for(n=0; n < order+1; n++)
            pData->bN_modal[band+1][n] = ccdiv(cmplx(1.0,0.0), (pData->bN[band*(order+1)+n]));
    for(n=0; n < order+1; n++)
        pData->bN_modal[0][n] = cmplx(0.0f, 0.0f); /* zero DC */
    switch(pData->regType){
        case REG_DAS:
            for(band=0; band<HYBRID_BANDS-1; band++){
                f_n = 0.0;
                for(n=0; n < order+1; n++)
                    f_n += (2.0*(double)n+1.0) * pow(cabs(pData->bN[band*(order+1)+n]), 2.0);
                beta = (1.0/ pow((double)(MAX_SH_ORDER)+1.0,2.0)) * f_n;
                for(n=0; n < MAX_SH_ORDER+1; n++)
                    pData->bN_inv[band+1][n] = crmul(pData->bN_modal[band+1][n], pow(cabs(pData->bN[band*(order+1)+n]), 2.0) / beta);
            }
            break;
            
        case REG_SOFT_LIM:
            g_lim = sqrt(arraySpecs->Q)*pow(10.0,(regPar/20.0));
            for(band=0; band<HYBRID_BANDS-1; band++)
                for(n=0; n < order+1; n++)
                    pData->bN_inv[band+1][n] = crmul(pData->bN_modal[band+1][n], (2.0*g_lim*cabs(pData->bN[band*(order+1)+n]) / M_PI)
                                                   * atan(M_PI / (2.0*g_lim*cabs(pData->bN[band*(order+1)+n]))) );
            break;
            
        case REG_TIKHONOV:
            alpha = sqrt(arraySpecs->Q)*pow(10.0,(regPar/20.0));
            for(band=0; band<HYBRID_BANDS-1; band++){
                for(n=0; n < order+1; n++){
                    beta = sqrt((1.0-sqrt(1.0-1.0/ pow(alpha,2.0)))/(1.0+sqrt(1.0-1.0/pow(alpha,2.0)))); /* Moreau & Daniel */
                    pData->bN_inv[band+1][n] = ccdiv(conj(pData->bN[band*(order+1)+n]), cmplx((pow(cabs(pData->bN[band*(order+1)+n]), 2.0) + pow(beta, 2.0)),0.0));
                }
            }
            break;
    }
    for(n=0; n < order+1; n++)
        pData->bN_inv[0][n] = cmplx(0.0f, 0.0f); /* remove NaN at DC */
    array2sh_replicate_order(hA2sh, order); /* replicate orders */
    
    /* Generate encoding matrix per band */
    Y_mic = NULL;
    getRSH(order, (float*)arraySpecs->sensorCoords_deg, arraySpecs->Q, &Y_mic); /* nSH x Q */
#if 1
    pinv_Y_mic = malloc( arraySpecs->Q * nSH *sizeof(float));
    utility_spinv(Y_mic, nSH, arraySpecs->Q, pinv_Y_mic);
    pinv_Y_mic_cmplx =  malloc((arraySpecs->Q) * nSH *sizeof(float_complex));
    for(i=0; i<(arraySpecs->Q)*nSH; i++)
        pinv_Y_mic_cmplx[i] = cmplxf(pinv_Y_mic[i], 0.0f);
#else
    int j;
    float* YYT;
    YYT = malloc(nSH*nSH*sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, nSH, (arraySpecs->Q), 1.0f,
                Y_mic, (arraySpecs->Q),
                Y_mic, (arraySpecs->Q), 0.0f,
                YYT, nSH);
    pinv_Y_mic =  malloc((arraySpecs->Q) * nSH *sizeof(float_complex));
    utility_sslslv(YYT, nSH, Y_mic, (arraySpecs->Q), pinv_Y_mic);
    pinv_Y_mic_cmplx =  malloc((arraySpecs->Q) * nSH *sizeof(float_complex));
    for(i=0; i<nSH; i++)
        for(j=0; j<(arraySpecs->Q); j++)
            pinv_Y_mic_cmplx[j*nSH+i] = cmplxf(pinv_Y_mic[i*(arraySpecs->Q)+j], 0.0f);
    free(YYT);
#endif 
    diag_bN_inv_R = calloc(nSH*nSH, sizeof(float_complex));
    for(band=0; band<HYBRID_BANDS; band++){
        for(i=0; i<nSH; i++)
            diag_bN_inv_R[i*nSH+i] = cmplxf((float)creal(pData->bN_inv_R[band][i]), (float)cimag(pData->bN_inv_R[band][i]));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, (arraySpecs->Q), nSH, &calpha,
                    diag_bN_inv_R, nSH,
                    pinv_Y_mic_cmplx, nSH, &cbeta,
                    pData->W[band], MAX_NUM_SENSORS);
    }
    
    pData->order = order; 
    
    free(Y_mic);
    free(pinv_Y_mic);
    free(pinv_Y_mic_cmplx);
    free(diag_bN_inv_R);
}


void array2sh_calculate_mag_curves(void* const hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    int band, n;
    
    for(band = 0; band <HYBRID_BANDS-1; band++){
        for(n = 0; n <pData->order+1; n++){
            pData->bN_inv_dB[band][n] = 20.0f * (float)log10(cabs(pData->bN_inv[band+1][n])); /* Ignore DC */
            pData->bN_modal_dB[band][n] = 20.0f * (float)log10(cabs(pData->bN_modal[band+1][n]));
        }
    }
}

void array2sh_evaluateSHTfilters(void* hA2sh)
{
    array2sh_data *pData = (array2sh_data*)(hA2sh);
    arrayPars* arraySpecs = (arrayPars*)(pData->arraySpecs);
    int band, i, j, simOrder, order, nSH;
    double kr[HYBRID_BANDS-1];
    float* Y_grid_real;
    float_complex* Y_grid, *H_array, *Wshort;
    
	pData->evalReady = 0;
    assert(pData->W != NULL);
    
    /* simulate the current array by firing 812 plane-waves around the surface of a theoretical sphere
     * and ascertaining the transfer function for each */
    simOrder = (int)(2.0f*M_PI*MAX_EVAL_FREQ_HZ*(arraySpecs->r)/pData->c)+1;
    for(i=0; i<HYBRID_BANDS-1; i++)
        kr[i] = 2.0*M_PI*(pData->freqVector[i+1/* ignore DC */])*(arraySpecs->R)/pData->c;
    H_array = malloc((HYBRID_BANDS-1) * (arraySpecs->Q) * 812*sizeof(float_complex));
    switch(arraySpecs->arrayType){
        case ARRAY_SPHERICAL:
            switch(arraySpecs->weightType){
                default:
                case WEIGHT_RIGID:
                    simulateSphArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID, 0.0, H_array);
                    break;
                case WEIGHT_OPEN_OMNI:
                    simulateSphArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN, 0.0, H_array);
                    break;
                case WEIGHT_OPEN_DIPOLE:
                    simulateSphArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_DIRECTIONAL, 1.0, H_array);
                    break;
                case WEIGHT_OPEN_CARD:
                    simulateSphArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_DIRECTIONAL, 0.5, H_array);
                    break;
            }
            break;
            
        case ARRAY_CYLINDRICAL:
            switch(arraySpecs->weightType){
                default:
                case WEIGHT_RIGID:
                    simulateCylArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_RIGID, H_array);
                    break;
                case WEIGHT_OPEN_DIPOLE:
                case WEIGHT_OPEN_CARD:
                case WEIGHT_OPEN_OMNI:
                    simulateCylArray(simOrder, kr, HYBRID_BANDS-1, (float*)arraySpecs->sensorCoords_rad, arraySpecs->Q, (float*)__geosphere_ico_9_0_dirs_deg, 812, ARRAY_CONSTRUCTION_OPEN, H_array);
                    break;
            }
            break;
    }
    
    /* generate ideal (real) spherical harmonics to compare with */
    order = pData->order;
    nSH = (order+1)*(order+1);
    Y_grid_real = NULL;
    getRSH(order, (float*)__geosphere_ico_9_0_dirs_deg, 812, &Y_grid_real);
    Y_grid = malloc(nSH*812*sizeof(float_complex));
    for(i=0; i<nSH*812; i++)
        Y_grid[i] = cmplxf(Y_grid_real[i], 0.0f); /* "evaluateSHTfilters" function requires complex data type */
    
    /* compare the spherical harmonics obtained from encoding matrix 'W' with the ideal patterns */
    Wshort = malloc(HYBRID_BANDS*nSH*(arraySpecs->Q)*sizeof(float_complex));
    for(band=0; band<HYBRID_BANDS-1; band++)
        for(i=0; i<nSH; i++)
            for(j=0; j<(arraySpecs->Q); j++)
                Wshort[band*nSH*(arraySpecs->Q) + i*(arraySpecs->Q) + j] = pData->W[band+1/* skip DC */][i][j];
    evaluateSHTfilters(order, Wshort, arraySpecs->Q, HYBRID_BANDS-1, H_array, 812, Y_grid, pData->cSH, pData->lSH);

	pData->evalReady = 1;

    free(Y_grid_real);
    free(Y_grid);
    free(H_array);
    free(Wshort);
}

void array2sh_createArray(void ** const hPars)
{
    arrayPars* pars = (arrayPars*)malloc(sizeof(arrayPars));
    if (pars == NULL) { return;/*error*/ }
    *hPars = (void*)pars;
}

void array2sh_destroyArray(void ** const hPars)
{
    arrayPars *pars = (arrayPars*)(*hPars);
    
    if(pars!=NULL) {
        free(pars);
        pars=NULL;
    }
}
 
void array2sh_initArray(void* const hPars, PRESETS preset, int* arrayOrder, int firstInitFlag)
{
    arrayPars *pars = (arrayPars*)(hPars);
    int ch, i, Q;
    
    switch(preset){
        default:
        case PRESET_DEFAULT:
            (*arrayOrder) = 1;
            Q = 4; /* number of mics */
            pars->r = 0.042f; /* array radius */
            pars->R = 0.042f; /* radius of the sensors (incase they protrude from the surface of the array), (only for rigid arrays) */
            pars->arrayType = ARRAY_SPHERICAL; /* spherical or cylindrical */
            pars->weightType = WEIGHT_RIGID; /* open or rigid, and directivity of the sensors (only for open arrays) */
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __default_coords_rad[ch][i]; /* spherical coordinates of the sensors, in radians */
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#ifdef ENABLE_AALTO_HYDROPHONE_PRESET
        case PRESET_AALTO_HYDROPHONE:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.1555f;
            pars->R = 0.1555f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_OMNI;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Aalto_Hydrophone_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_SENNHEISER_AMBEO_PRESET
        case PRESET_SENNHEISER_AMBEO:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.02f;
            pars->R = 0.02f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Sennheiser_Ambeo_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_CORE_SOUND_TETRAMIC_PRESET
        case PRESET_CORE_SOUND_TETRAMIC:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.02f;
            pars->R = 0.02f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Core_Sound_TetraMic_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_SOUND_FIELD_SPS200_PRESET
        case PRESET_SOUND_FIELD_SPS200:
            (*arrayOrder) = 1;
            Q = 4;
            pars->r = 0.02f;
            pars->R = 0.02f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_OPEN_CARD;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Sound_field_SPS200_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_ZYLIA_1D_PRESET
        case PRESET_ZYLIA_1D:
            (*arrayOrder) = 3;
            Q = 19;
            pars->r = 0.049f;
            pars->R = 0.049f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Zylia1D_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_EIGENMIKE32_PRESET
        case PRESET_EIGENMIKE32:
            (*arrayOrder) = 4;
            Q = 32;
            pars->r = 0.042f;
            pars->R = 0.042f;
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __Eigenmike32_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif
#ifdef ENABLE_DTU_MIC_PRESET
        case PRESET_DTU_MIC:
            (*arrayOrder) = 6;
            Q = 52;
            pars->r = 0.05f;
            pars->R = 0.05f; 
            pars->arrayType = ARRAY_SPHERICAL;
            pars->weightType = WEIGHT_RIGID;
            for(ch=0; ch<Q; ch++){
                for(i=0; i<2; i++){
                    pars->sensorCoords_rad[ch][i] = __DTU_mic_coords_rad[ch][i];
                    pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
                }
            }
            break;
#endif 
    }
    
    /* Fill remaining slots with default coords */
    for(; ch<MAX_NUM_SENSORS; ch++){
        for(i=0; i<2; i++){
            pars->sensorCoords_rad[ch][i] = __default_SENSORcoords64_rad[ch][i];
            pars->sensorCoords_deg[ch][i] = pars->sensorCoords_rad[ch][i] * (180.0f/M_PI);
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




