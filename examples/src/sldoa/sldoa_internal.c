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
 * @file sldoa_internal.c
 * @brief A spatially-localised active-intensity (SLAI) based direction-of-
 *        arrival estimator (SLDoA)
 *
 * VBAP gain patterns are imposed on the spherical harmonic signals, such that
 * the DoA can be estimated in a spatially-constrained region; thus mitigating
 * the effect of interferes and reflections arriving from other directions.
 * The DoA is estimated per sector for each frequency band.
 *
 * The algorithms within sldoa were developed in collaboration with Symeon
 * Delikaris-Manias and Angelo Farina, and are explained in more detail in [1,2]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 * @see [2] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., "Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis," in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 *
 * @author Leo McCormack
 * @date 18.10.2017
 * @license ISC
 */

#include "sldoa.h"
#include "sldoa_internal.h"
#include "sldoa_database.h"

void sldoa_setCodecStatus(void* const hSld, CODEC_STATUS newStatus)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    if(newStatus==CODEC_STATUS_NOT_INITIALISED){
        /* Pause until current initialisation is complete */
        while(pData->codecStatus == CODEC_STATUS_INITIALISING)
            SAF_SLEEP(10);
    }
    pData->codecStatus = newStatus;
}

void sldoa_initAna(void* const hSld)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    int i, n, j, k, order, nSectors, nSH, grid_N_vbap_gtable, grid_nGroups, maxOrder;
    float* sec_dirs_deg, *grid_vbap_gtable, *w_SG, *pinv_Y, *grid_vbap_gtable_T;
    float secPatterns[4][NUM_GRID_DIRS];

    maxOrder = pData->new_masterOrder;
    
    grid_vbap_gtable_T = malloc1d(ORDER2NUMSECTORS(maxOrder) * NUM_GRID_DIRS * sizeof(float));
    
    for(i=0, order=2; order<=maxOrder; i++,order++){
        nSectors = ORDER2NUMSECTORS(order);
        nSH = (order+1)*(order+1);
        
        /* define sector directions */
        sec_dirs_deg = malloc1d(nSectors*2*sizeof(float));
        memcpy(sec_dirs_deg, __HANDLES_SphCovering_dirs_deg[nSectors-1], nSectors*2*sizeof(float));
        
        /* generate VBAP gain table */
        generateVBAPgainTable3D_srcs((float*)pData->grid_dirs_deg, NUM_GRID_DIRS, sec_dirs_deg, nSectors, 0, 0, 0.0f,
                                     &(grid_vbap_gtable), &(grid_N_vbap_gtable), &(grid_nGroups));
        
        /* convert to amplitude preserving gains */
        VBAPgainTable2InterpTable(grid_vbap_gtable, NUM_GRID_DIRS, nSectors);
        
        /* transpose */
        for(n=0; n<nSectors; n++)
            for(j=0; j<NUM_GRID_DIRS; j++)
                grid_vbap_gtable_T[n*NUM_GRID_DIRS+j] = grid_vbap_gtable[j*nSectors+n];
        
        /* generate sector coefficients */
        if(pData->secCoeffs[i]!=NULL)
            free(pData->secCoeffs[i]);
        pData->secCoeffs[i] = malloc1d(4 * (nSH*nSectors) * sizeof(float_complex));
        w_SG = malloc1d(4 * (nSH) * sizeof(float));
        pinv_Y = malloc1d(NUM_GRID_DIRS*nSH*sizeof(float));
        for(n=0; n<nSectors; n++){ 
            utility_svvmul(&(grid_vbap_gtable_T[n*NUM_GRID_DIRS]), pData->grid_Y[0], NUM_GRID_DIRS, secPatterns[0]);
            for(j=0; j<3; j++)
                utility_svvmul(&(grid_vbap_gtable_T[n*NUM_GRID_DIRS]), pData->grid_Y_dipoles_norm[j], NUM_GRID_DIRS, secPatterns[j+1]);
            utility_spinv(NULL, &(pData->grid_Y[0][0]), nSH, NUM_GRID_DIRS, pinv_Y);
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, nSH, NUM_GRID_DIRS, 1.0f,
                        &(secPatterns[0][0]), NUM_GRID_DIRS,
                        pinv_Y, nSH, 0.0f,
                        w_SG, nSH);
            
            /* stack the sector coefficients */
            for(j=0; j<4; j++)
                for(k=0; k<nSH; k++)
                    pData->secCoeffs[i][j*(nSectors*nSH)+n*nSH+k] = cmplxf(w_SG[j*nSH+k], 0.0f);
        }
        free(w_SG);
        free(pinv_Y);
        free(grid_vbap_gtable);
        free(sec_dirs_deg);
    }
    
    free(grid_vbap_gtable_T);
    
    pData->masterOrder = maxOrder;
}

void sldoa_initTFT
(
    void* const hSld
)
{
    sldoa_data *pData = (sldoa_data*)(hSld);
    int nSH, new_nSH;
    
    nSH = (pData->masterOrder+1)*(pData->masterOrder+1);
    new_nSH = (pData->new_masterOrder+1)*(pData->new_masterOrder+1);
    if(pData->hSTFT==NULL)
        afSTFT_create(&(pData->hSTFT), new_nSH, 0, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    else if(nSH!=new_nSH){
        afSTFT_channelChange(pData->hSTFT, new_nSH, 0);
        afSTFT_clearBuffers(pData->hSTFT);
    }
}

void sldoa_estimateDoA
(
    float_complex** SHframeTF,
    int anaOrder,
    float_complex* secCoeffs,
    float doa[MAX_NUM_SECTORS][TIME_SLOTS][2],
    float energy[MAX_NUM_SECTORS][TIME_SLOTS]
)
{
    int n, ch, i, j, nSectors, analysisOrder, nSH;
    float_complex secSig[4][TIME_SLOTS];
    float_complex* sec_c;
    float secEnergy[TIME_SLOTS], secIntensity[3][TIME_SLOTS], secAzi[TIME_SLOTS], secElev[TIME_SLOTS];
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta = cmplxf(0.0f, 0.0f);
    
    /* prep */
    memset(doa,0,MAX_NUM_SECTORS*TIME_SLOTS*2*sizeof(float));
    memset(energy,0,MAX_NUM_SECTORS*TIME_SLOTS*sizeof(float));
    analysisOrder = SAF_MAX(SAF_MIN(MAX_SH_ORDER, anaOrder),1);
    nSectors = ORDER2NUMSECTORS(analysisOrder);
    nSH = (analysisOrder+1)*(analysisOrder+1);
    sec_c = malloc1d(4*nSH*sizeof(float_complex));
    
    /* calculate energy and DoA for each sector */
    for( n=0; n<nSectors; n++){
        if(anaOrder==1 || secCoeffs == NULL) /* standard first order active-intensity based DoA estimation */
            for (i=0; i<4; i++)
                memcpy(secSig[i], SHframeTF[i], TIME_SLOTS * sizeof(float_complex));
        else{ /* spatially localised active-intensity based DoA estimation */
            for (i=0; i<4; i++)
                for (j=0; j<nSH; j++)
                    sec_c[i*nSH+j] = secCoeffs[i*(nSectors*nSH)+n*nSH+j];
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, TIME_SLOTS, nSH, &calpha,
                        sec_c, nSH,
                        FLATTEN2D(SHframeTF), TIME_SLOTS, &cbeta,
                        secSig, TIME_SLOTS);
        }
        
        /* convert N3D to SN3D */
        for (ch = 1; ch<4; ch++)
            for(i = 0; i<TIME_SLOTS; i++)
                secSig[ch][i] = crmulf(secSig[ch][i], 1.0f/sqrtf(3.0f));
        
        /* calculate sector energy and intensity vector */
        memset(secEnergy, 0, TIME_SLOTS*sizeof(float));
        for (i=0; i<4; i++)
            for (j=0; j<TIME_SLOTS; j++)
                secEnergy[j] += 0.5f*powf(cabsf(secSig[i][j]), 2.0f);
        for (i=0; i<3; i++)
            for (j=0; j<TIME_SLOTS; j++)
                secIntensity[i][j] = crealf(ccmulf(conjf(secSig[0][j]), secSig[1+i][j]));
        
        /* extract DoA */
        for (j=0; j<TIME_SLOTS; j++){
            secAzi[j]  = atan2f( secIntensity[0][j], secIntensity[2][j] );
            secElev[j] = atan2f( secIntensity[1][j], sqrtf( powf(secIntensity[2][j], 2.0f) + powf(secIntensity[0][j], 2.0f)) );
        }
        
        /* store energy and DoA estimate */
        for (j=0; j<TIME_SLOTS; j++){
            doa[n][j][0] = secAzi[j];
            doa[n][j][1] = secElev[j];
            energy[n][j] = secEnergy[j]*1e6f;
        }
    }
    
    free(sec_c);
}

