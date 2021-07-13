/*
 * Copyright 2016-2018 Leo McCormack
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
 * @file powermap_internal.c
 * @brief A sound-field visualiser, which utilises spherical harmonic signals as
 *        input; note this code is a remnant from the work conducted in [1]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S. and Pulkki, V., 2017. Parametric
 *          acoustic camera for real-time sound capture, analysis and tracking.
 *          In Proceedings of the 20th International Conference on Digital Audio
 *          Effects (DAFx-17) (pp. 412-419)
 *
 * @author Leo McCormack
 * @date 26.04.2016
 * @license ISC
 */

#include "powermap.h"
#include "powermap_internal.h"

void powermap_setCodecStatus(void* const hPm, CODEC_STATUS newStatus)
{
    powermap_data *pData = (powermap_data*)(hPm);
    if(newStatus==CODEC_STATUS_NOT_INITIALISED){
        /* Pause until current initialisation is complete */
        while(pData->codecStatus == CODEC_STATUS_INITIALISING)
            SAF_SLEEP(10);
    }
    pData->codecStatus = newStatus;
}

void powermap_initAna(void* const hPm)
{
    powermap_data *pData = (powermap_data*)(hPm);
    powermap_codecPars* pars = pData->pars;
    int i, j, n, N_azi, N_ele, nSH_order, order;
    float scaleY, hfov, vfov, fi, aspectRatio;
    float* Y_grid_N, *grid_x_axis, *grid_y_axis;
    
    order = pData->new_masterOrder;
    
    /* Store Y_grid per order */
    int geosphere_ico_freq = 9;
    pars->grid_dirs_deg = (float*)__HANDLES_geosphere_ico_dirs_deg[geosphere_ico_freq];
    pars->grid_nDirs = __geosphere_ico_nPoints[geosphere_ico_freq];
    Y_grid_N = malloc1d(((order+1)*(order+1))*(pars->grid_nDirs)*sizeof(float));
    getRSH(order, pars->grid_dirs_deg, pars->grid_nDirs, Y_grid_N);
    for(n=1; n<=order; n++){
        nSH_order = (n+1)*(n+1);
        scaleY = 1.0f/(float)nSH_order;
        free(pars->Y_grid[n-1]);
        free(pars->Y_grid_cmplx[n-1]);
        pars->Y_grid[n-1] = malloc1d(nSH_order * (pars->grid_nDirs)*sizeof(float));
        pars->Y_grid_cmplx[n-1] = malloc1d(nSH_order * (pars->grid_nDirs)*sizeof(float_complex));
        memcpy(pars->Y_grid[n-1], Y_grid_N, nSH_order * (pars->grid_nDirs)*sizeof(float));
        utility_svsmul(pars->Y_grid[n-1], &scaleY, nSH_order * (pars->grid_nDirs), NULL);
        for(i=0; i<nSH_order; i++)
            for(j=0; j<pars->grid_nDirs; j++)
                pars->Y_grid_cmplx[n-1][i*(pars->grid_nDirs)+j] = cmplxf(pars->Y_grid[n-1][i*(pars->grid_nDirs)+j], 0.0f);
    }

    /* generate interpolation table for current display settings */
    switch(pData->HFOVoption){
        default:
        case HFOV_360: hfov = 360.0f; break;
    }
    switch(pData->aspectRatioOption){
        default:
        case ASPECT_RATIO_2_1: aspectRatio = 2.0f; break;
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
    free(pars->interp_dirs_deg);
    pars->interp_dirs_deg = malloc1d(N_azi*N_ele*2*sizeof(float));
    for(i = 0; i<N_ele; i++){
        for(j=0; j<N_azi; j++){
            pars->interp_dirs_deg[(i*N_azi + j)*2]   = grid_x_axis[j];
            pars->interp_dirs_deg[(i*N_azi + j)*2+1] = grid_y_axis[i];
        }
    }
    free(pars->interp_table);
    generateVBAPgainTable3D_srcs(pars->interp_dirs_deg, N_azi*N_ele, pars->grid_dirs_deg, pars->grid_nDirs, 0, 0, 0.0f, &(pars->interp_table), &(pars->interp_nDirs), &(pars->interp_nTri));
    VBAPgainTable2InterpTable(pars->interp_table, pars->interp_nDirs, pars->grid_nDirs);
    
    /* reallocate memory for storing the powermaps */
    free(pData->pmap);
    pData->pmap = malloc1d(pars->grid_nDirs*sizeof(float));
    free(pData->prev_pmap);
    pData->prev_pmap = calloc1d(pars->grid_nDirs, sizeof(float));
    for(i=0; i<NUM_DISP_SLOTS; i++){
        free(pData->pmap_grid[i]);
        pData->pmap_grid[i] = calloc1d(pars->interp_nDirs,sizeof(float));
    }
    
    pData->masterOrder = order;
    
    free(Y_grid_N);
    free(grid_x_axis);
    free(grid_y_axis);
}

void powermap_initTFT
(
    void* const hPm
)
{
    powermap_data *pData = (powermap_data*)(hPm);
    int nSH, new_nSH;
    
    nSH = (pData->masterOrder+1)*(pData->masterOrder+1);
    new_nSH = (pData->new_masterOrder+1)*(pData->new_masterOrder+1);
    if(pData->hSTFT==NULL)
        afSTFT_create(&(pData->hSTFT), new_nSH, 0, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    else if(nSH!=new_nSH){
        afSTFT_channelChange(pData->hSTFT, new_nSH, 0);
        afSTFT_clearBuffers(pData->hSTFT);
        memset(pData->Cx, 0 , MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS*HYBRID_BANDS*sizeof(float_complex));
    }
}
