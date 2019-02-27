/*
 Copyright 2019 Leo McCormack
 
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
 *     dirass_internal.c
 * Description:
 *     A sound-field visualiser based on the directional re-assignment of beamformer energy,
 *     utilising the DoA estimates extracted from spatially-localised active-intensity
 *     (SLAI) vectors; which correspond to the scanning grid directions.
 *     For more information on the method, refer to:
 *         McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of angular
 *         spectra based on a directional re-assignment approach for ambisonic sound-field
 *         visualisation". IEEE International Conference on Acoustics, Speech and Signal
 *         Processing (ICASSP).
 *
 * Dependencies:
 *     saf_utilities, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 21.02.2019
 */

#include "dirass.h"
#include "dirass_internal.h" 

void dirass_initAna(void* const hDir)
{
    dirass_data *pData = (dirass_data*)(hDir);
    codecPars* pars = pData->pars;
    int i, j, N_azi, N_ele, nSH_order, order, nSH_sec, order_sec;
    float scaleY, hfov, vfov, fi, aspectRatio;
    float *grid_x_axis, *grid_y_axis, *c_n;
    float_complex* A_xyz;
    
    order = pData->new_inputOrder;
    nSH_order = (order+1)*(order+1);
    
    /* Calculate spherical harmonic weights for each grid direction (Y_grid) */
    free(pars->Y_grid);
    int geosphere_ico_freq = 9;
    pars->grid_dirs_deg = (float*)__HANDLES_geosphere_ico_dirs_deg[geosphere_ico_freq];
    pars->grid_nDirs = __geosphere_ico_nPoints[geosphere_ico_freq];
    getRSH(order, pars->grid_dirs_deg, pars->grid_nDirs, &(pars->Y_grid));
    scaleY = 1.0f/(float)nSH_order;
    //utility_svsmul(pars->Y_grid, &scaleY, nSH_order * (pars->grid_nDirs), NULL);
    
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
    grid_x_axis = malloc(N_azi * sizeof(float));
    grid_y_axis = malloc(N_ele * sizeof(float));
    vfov = hfov/aspectRatio;
    for(fi = -hfov/2.0f, i = 0; i<N_azi; fi+=hfov/N_azi, i++)
        grid_x_axis[i] = fi;
    for(fi = -vfov/2.0f,  i = 0; i<N_ele; fi+=vfov/N_ele, i++)
        grid_y_axis[i] = fi;
    free(pars->interp_dirs_deg);
    pars->interp_dirs_deg = malloc(N_azi*N_ele*2*sizeof(float));
    for(i = 0; i<N_ele; i++){
        for(j=0; j<N_azi; j++){
            pars->interp_dirs_deg[(i*N_azi + j)*2]   = grid_x_axis[j];
            pars->interp_dirs_deg[(i*N_azi + j)*2+1] = grid_y_axis[i];
        }
    }
    free(pars->interp_table);
    generateVBAPgainTable3D_srcs(pars->interp_dirs_deg, N_azi*N_ele, pars->grid_dirs_deg, pars->grid_nDirs, 0, 0, 0.0f, &(pars->interp_table), &(pars->interp_nDirs), &(pars->interp_nTri));
    VBAPgainTable2InterpTable(pars->interp_table, pars->interp_nDirs, pars->grid_nDirs);
    
    /* get sector matrices (hyper-cardioids) */
    order_sec = order-1;
    nSH_sec = (order_sec+1)*(order_sec+1);
    A_xyz = malloc(nSH_order*nSH_sec*3*sizeof(float_complex));
    c_n = malloc((order_sec+1)*sizeof(float));
    computeVelCoeffsMtx(order_sec, A_xyz);
    beamWeightsHypercardioid2Spherical(order_sec, c_n);
    free(pars->Cxyz);
    free(pars->Cw);
    pars->Cxyz = malloc(pars->grid_nDirs * nSH_order * 3 * sizeof(float));
    pars->Cw = malloc(pars->grid_nDirs * nSH_sec * sizeof(float));
    for(i=0; i<pars->grid_nDirs; i++){
        beamWeightsVelocityPatternsReal(order_sec, c_n, pars->grid_dirs_deg[i*2]*M_PI/180.0f,
                                        pars->grid_dirs_deg[i*2+1]*M_PI/180.0f, A_xyz, &(pars->Cxyz[i*nSH_order*3]));
        rotateAxisCoeffsReal(order_sec, c_n, M_PI/2.0f - pars->grid_dirs_deg[i*2+1]*M_PI/180.0f,
                             pars->grid_dirs_deg[i*2]*M_PI/180.0f, &(pars->Cw[i*nSH_sec]));
    }
    free(A_xyz);
    free(c_n);
    
//    float c_nTEST[3];
//    float_complex A_xyzTEST[16][9][3];
//    float CxyzTEST[812][16][3];
//    float CwTEST[812][9];
//    memcpy(c_nTEST, c_n, 3*sizeof(float));
//    memcpy(A_xyzTEST, A_xyz, 16*9*3*sizeof(float_complex));
//    memcpy(CxyzTEST, Cxyz, 812*16*3*sizeof(float));
//    memcpy(CwTEST, Cw, 812*9*sizeof(float));
//
    
//    /* reallocate memory for storing the dirasss */
//    free(pData->pmap);
//    pData->pmap = malloc(pars->grid_nDirs*sizeof(float));
//    free(pData->prev_pmap);
//    pData->prev_pmap = calloc(pars->grid_nDirs, sizeof(float));
//    for(i=0; i<NUM_DISP_SLOTS; i++){
//        free(pData->pmap_grid[i]);
//        pData->pmap_grid[i] = calloc(pars->interp_nDirs,sizeof(float));
//    }
//    
//    pData->inputOrder = order;
    
    free(grid_x_axis);
    free(grid_y_axis);
}
 








