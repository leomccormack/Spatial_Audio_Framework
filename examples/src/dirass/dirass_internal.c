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
 * @file dirass_internal.c
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

void dirass_setCodecStatus(void* const hDir, CODEC_STATUS newStatus)
{
    dirass_data *pData = (dirass_data*)(hDir);
    if(newStatus==CODEC_STATUS_NOT_INITIALISED){
        /* Pause until current initialisation is complete */
        while(pData->codecStatus == CODEC_STATUS_INITIALISING)
            SAF_SLEEP(10);
    }
    pData->codecStatus = newStatus;
}

void dirass_initAna(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    dirass_codecPars* pars = pData->pars;
    int i, j, N_azi, N_ele, nSH_order, order, nSH_sec, order_sec, order_up, nSH_up, geosphere_ico_freq, td_degree;
    float hfov, vfov, fi, aspectRatio;
    float *grid_x_axis, *grid_y_axis, *c_n;
    float_complex* A_xyz;
    
    order = pData->new_inputOrder;
    order_up = pData->new_upscaleOrder;
    nSH_order = (order+1)*(order+1);
    nSH_up = (order_up+1)*(order_up+1);

    strcpy(pData->progressBarText,"Preparing scanning grid");
    pData->progressBar0_1 = 0.4f;
    
    /* Determine scanning grid */
    switch(pData->gridOption){
        case T_DESIGN_3:           /* 6 points */
            td_degree = 3;
            pars->grid_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[td_degree-1];
            pars->grid_nDirs = __Tdesign_nPoints_per_degree[td_degree-1];
            break;
        case T_DESIGN_4:           /* 12 points */
            td_degree = 4;
            pars->grid_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[td_degree-1];
            pars->grid_nDirs = __Tdesign_nPoints_per_degree[td_degree-1];
            break;
        case T_DESIGN_6:           /* 24 points */
            td_degree = 6;
            pars->grid_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[td_degree-1];
            pars->grid_nDirs = __Tdesign_nPoints_per_degree[td_degree-1];
            break;
        case T_DESIGN_9:           /* 48 points */
            td_degree = 9;
            pars->grid_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[td_degree-1];
            pars->grid_nDirs = __Tdesign_nPoints_per_degree[td_degree-1];
            break;
        case T_DESIGN_13:          /* 94 points */
            td_degree = 13;
            pars->grid_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[td_degree-1];
            pars->grid_nDirs = __Tdesign_nPoints_per_degree[td_degree-1];
            break;
        case T_DESIGN_18:          /* 180 points */
            td_degree = 18;
            pars->grid_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[td_degree-1];
            pars->grid_nDirs = __Tdesign_nPoints_per_degree[td_degree-1];
            break;
        case GRID_GEOSPHERE_6:     /* 362 points */
            geosphere_ico_freq = 6;
            pars->grid_dirs_deg = (float*)__HANDLES_geosphere_ico_dirs_deg[geosphere_ico_freq];
            pars->grid_nDirs = __geosphere_ico_nPoints[geosphere_ico_freq];
            break;
        case T_DESIGN_30:          /* 480 points */
            pars->grid_dirs_deg = (float*)__Tdesign_degree_30_dirs_deg;
            pars->grid_nDirs = 480;
            break;
        case GRID_GEOSPHERE_8:     /* 642 points */
            geosphere_ico_freq = 8;
            pars->grid_dirs_deg = (float*)__HANDLES_geosphere_ico_dirs_deg[geosphere_ico_freq];
            pars->grid_nDirs = __geosphere_ico_nPoints[geosphere_ico_freq];
            break;
        case GRID_GEOSPHERE_9:     /* 812 points */
            geosphere_ico_freq = 9;
            pars->grid_dirs_deg = (float*)__HANDLES_geosphere_ico_dirs_deg[geosphere_ico_freq];
            pars->grid_nDirs = __geosphere_ico_nPoints[geosphere_ico_freq];
            break;
        case GRID_GEOSPHERE_10:    /* 1002 points */
            geosphere_ico_freq = 10;
            pars->grid_dirs_deg = (float*)__HANDLES_geosphere_ico_dirs_deg[geosphere_ico_freq];
            pars->grid_nDirs = __geosphere_ico_nPoints[geosphere_ico_freq];
            break;
        case GRID_GEOSPHERE_12:    /* 1442 points */
            geosphere_ico_freq = 12;
            pars->grid_dirs_deg = (float*)__HANDLES_geosphere_ico_dirs_deg[geosphere_ico_freq];
            pars->grid_nDirs = __geosphere_ico_nPoints[geosphere_ico_freq];
            break;
    }
    
    /* generate interpolation table for current display config */
    switch(pData->HFOVoption){
        default:
        case HFOV_360: hfov = 360.0f; break;
        case HFOV_180: hfov = 180.0f; break;
        case HFOV_90:  hfov = 90.0f;  break;
        case HFOV_60:  hfov = 60.0f;  break;
    }
    switch(pData->aspectRatioOption){
        default:
        case ASPECT_RATIO_2_1:  aspectRatio = 2.0f; break;
        case ASPECT_RATIO_16_9: aspectRatio = 16.0f/9.0f; break;
        case ASPECT_RATIO_4_3:  aspectRatio = 4.0f/3.0f; break;
    }
    N_azi = pData->dispWidth;
    N_ele = (int)((float)pData->dispWidth/aspectRatio + 0.5f);
    grid_x_axis = malloc1d(N_azi * sizeof(float));
    grid_y_axis = malloc1d(N_ele * sizeof(float));
    vfov = hfov/aspectRatio;
    for(fi = -hfov/2.0f, i = 0; i<N_azi; fi+=hfov/N_azi, i++)
        grid_x_axis[i] = fi;
    for(fi = -vfov/2.0f,  i = 0; i<N_ele; fi+=vfov/N_ele, i++)
        grid_y_axis[i] = fi;
    pars->interp_dirs_deg = realloc1d(pars->interp_dirs_deg, N_azi*N_ele*2*sizeof(float));
    pars->interp_dirs_rad = realloc1d(pars->interp_dirs_rad, N_azi*N_ele*2*sizeof(float));
    for(i = 0; i<N_ele; i++){
        for(j=0; j<N_azi; j++){
            pars->interp_dirs_deg[(i*N_azi + j)*2]   = grid_x_axis[j];
            pars->interp_dirs_deg[(i*N_azi + j)*2+1] = grid_y_axis[i];
            pars->interp_dirs_rad[(i*N_azi + j)*2] = grid_x_axis[j] * SAF_PI/180.0f;
            pars->interp_dirs_rad[(i*N_azi + j)*2+1] = grid_y_axis[i] * SAF_PI/180.0f;
        }
    }
    free(pars->interp_table);
    generateVBAPgainTable3D_srcs(pars->interp_dirs_deg, N_azi*N_ele, pars->grid_dirs_deg, pars->grid_nDirs, 0, 0, 0.0f, &(pars->interp_table), &(pars->interp_nDirs), &(pars->interp_nTri));
    VBAPgainTable2InterpTable(pars->interp_table, pars->interp_nDirs, pars->grid_nDirs);
    
    strcpy(pData->progressBarText,"Computing Sector coefficients");
    pData->progressBar0_1 = 0.85f;
    
    /* get beamforming matrices for sector velocity and sector patterns */
    order_sec = order-1;
    nSH_sec = (order_sec+1)*(order_sec+1);
    A_xyz = malloc1d(nSH_order*nSH_sec*3*sizeof(float_complex));
    c_n = malloc1d((order_sec+1)*sizeof(float));
    computeVelCoeffsMtx(order_sec, A_xyz);
    switch(pData->beamType){
        case STATIC_BEAM_TYPE_CARDIOID: beamWeightsCardioid2Spherical(order_sec, c_n); break;
        case STATIC_BEAM_TYPE_HYPERCARDIOID: beamWeightsHypercardioid2Spherical(order_sec, c_n); break;
        case STATIC_BEAM_TYPE_MAX_EV: beamWeightsMaxEV(order_sec, c_n); break;
    }
    pars->Cxyz = realloc1d(pars->Cxyz, pars->grid_nDirs * nSH_order * 3 * sizeof(float));
    pars->Cw = realloc1d(pars->Cw, pars->grid_nDirs * nSH_sec * sizeof(float));
    for(i=0; i<pars->grid_nDirs; i++){
        beamWeightsVelocityPatternsReal(order_sec, c_n, pars->grid_dirs_deg[i*2]*SAF_PI/180.0f,
                                        pars->grid_dirs_deg[i*2+1]*SAF_PI/180.0f, A_xyz, &(pars->Cxyz[i*nSH_order*3]));
        rotateAxisCoeffsReal(order_sec, c_n, SAF_PI/2.0f - pars->grid_dirs_deg[i*2+1]*SAF_PI/180.0f,
                             pars->grid_dirs_deg[i*2]*SAF_PI/180.0f, &(pars->Cw[i*nSH_sec]));
    }
    free(A_xyz);
    free(c_n);

    /* get regular beamforming weights */
    c_n = malloc1d((order+1)*sizeof(float));
    switch(pData->beamType){
        case STATIC_BEAM_TYPE_CARDIOID: beamWeightsCardioid2Spherical(order, c_n); break;
        case STATIC_BEAM_TYPE_HYPERCARDIOID: beamWeightsHypercardioid2Spherical(order, c_n); break;
        case STATIC_BEAM_TYPE_MAX_EV: beamWeightsMaxEV(order, c_n); break;
    }
    pars->w = realloc1d(pars->w, pars->grid_nDirs * nSH_order * sizeof(float));
    for(i=0; i<pars->grid_nDirs; i++){
        rotateAxisCoeffsReal(order, c_n, SAF_PI/2.0f - pars->grid_dirs_deg[i*2+1]*SAF_PI/180.0f,
                             pars->grid_dirs_deg[i*2]*SAF_PI/180.0f, &(pars->w[i*nSH_order]));
    }
    free(c_n);
 
    /* beamforming weights for upscaled */
    c_n = malloc1d((order_up+1)*sizeof(float)); 
    switch(pData->beamType){
        case STATIC_BEAM_TYPE_CARDIOID: beamWeightsCardioid2Spherical(order_up, c_n); break;
        case STATIC_BEAM_TYPE_HYPERCARDIOID: beamWeightsHypercardioid2Spherical(order_up, c_n); break;
        case STATIC_BEAM_TYPE_MAX_EV: beamWeightsMaxEV(order_up, c_n); break;
    } 
    pars->Uw = realloc1d(pars->Uw, pars->grid_nDirs * nSH_up * sizeof(float));
    for(i=0; i<pars->grid_nDirs; i++){
        rotateAxisCoeffsReal(order_up, c_n, SAF_PI/2.0f - pars->grid_dirs_deg[i*2+1]*SAF_PI/180.0f,
                             pars->grid_dirs_deg[i*2]*SAF_PI/180.0f, &(pars->Uw[i*nSH_up]));
    }
    free(c_n);
 
    /* reallocate memory */
    pars->Y_up = realloc1d(pars->Y_up, nSH_up * (pars->grid_nDirs)*sizeof(float));
    pars->est_dirs = realloc1d(pars->est_dirs, pars->grid_nDirs * 2 * sizeof(float));
    pars->ss = realloc1d(pars->ss, pars->grid_nDirs * DIRASS_FRAME_SIZE * sizeof(float));
    pars->ssxyz = realloc1d(pars->ssxyz, 3 * DIRASS_FRAME_SIZE * sizeof(float));
    pData->pmap = realloc1d(pData->pmap, pars->grid_nDirs*sizeof(float));
    pars->est_dirs_idx = realloc1d(pars->est_dirs_idx, pars->grid_nDirs*sizeof(int));
    pars->prev_intensity = realloc1d(pars->prev_intensity, pars->grid_nDirs*3*sizeof(float));
    pars->prev_energy = realloc1d(pars->prev_energy, pars->grid_nDirs*sizeof(float));
    memset(pars->prev_intensity, 0, pars->grid_nDirs*3*sizeof(float));
    memset(pars->prev_energy, 0, pars->grid_nDirs*sizeof(float)); 
    for(i=0; i<NUM_DISP_SLOTS; i++){
        pData->pmap_grid[i] = realloc1d(pData->pmap_grid[i], pars->interp_nDirs*sizeof(float));
        memset(pData->pmap_grid[i], 0, pars->interp_nDirs*sizeof(float));
    }
    
    pData->inputOrder = order;
    pData->upscaleOrder = order_up;
    
    free(grid_x_axis);
    free(grid_y_axis);
}
