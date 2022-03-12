/*
 * This file is part of the saf_hades module unit tests.
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
 * @file test__hades_module.c
 * @brief Unit tests for the SAF HADES module
 * @author Leo McCormack and Janani Fernandez
 * @date 09.04.2021
 * @license GNU GPLv2
 */

#include "saf_test.h"

#ifdef SAF_ENABLE_HADES_MODULE

void test__hades(void){
    hades_analysis_handle hAna = NULL;          /* Analysis handle */
    hades_synthesis_handle hSyn = NULL;         /* Synthesis handle */
    hades_param_container_handle hPCon = NULL;  /* Parameter container handle */
    hades_signal_container_handle hSCon = NULL; /* Signal container handle */
    hades_binaural_config binConfig;
    SAF_SOFA_ERROR_CODES error;
    saf_sofa_container sofa;
    int i, ch, nDirs, nMics;
    int refIndices[2];
    float* grid_dirs_deg;

    /* Config */
    const int fs = 48000;
    const int sigLen = fs*2;
    const int hopsize = 64;
    const int blocksize = 256;
    const int hybridmode = 0;

    /* Analysis */
    error = saf_sofa_open(&sofa, "/Users/mccorml1/Documents/git/matlab/h_array/h_array_horiz1deg_357.sofa", SAF_SOFA_READER_OPTION_DEFAULT);
    if(error!=SAF_SOFA_OK)
        return; /* SOFA File does not exist, so skip this unit test. */
    nDirs = sofa.nSources;
    nMics = sofa.nReceivers;
    grid_dirs_deg = malloc1d(nDirs*2*sizeof(float));
    cblas_scopy(nDirs, sofa.SourcePosition, 3, grid_dirs_deg, 2); /* azi */
    cblas_scopy(nDirs, &sofa.SourcePosition[1], 3, &grid_dirs_deg[1], 2); /* elev */
    hades_analysis_create(&hAna, (float)fs, HADES_USE_AFSTFT_LD, hopsize, blocksize, hybridmode, /* for time-frequency transform */
                          sofa.DataIR, grid_dirs_deg, nDirs, nMics, sofa.DataLengthIR,    /* for the array measurements */
                          HADES_USE_COMEDIE, HADES_USE_MUSIC);                            /* for parameter analysis */
    saf_sofa_close(&sofa);

    /* Parameter/signal containers */
    hades_param_container_create(&hPCon, hAna);
    hades_signal_container_create(&hSCon, hAna);

    /* Synthesis */
    binConfig.hrir_fs = __default_hrir_fs;
    binConfig.lHRIR = __default_hrir_len;
    binConfig.nHRIR = __default_N_hrir_dirs;
    binConfig.hrirs = (float*)__default_hrirs;
    binConfig.hrir_dirs_deg = (float*)__default_hrir_dirs_deg;
    refIndices[0] = 1;
    refIndices[1] = 5;
    hades_synthesis_create(&hSyn, hAna, HADES_BEAMFORMER_BMVDR, SAF_TRUE, refIndices, &binConfig, HADES_HRTF_INTERP_NEAREST);

    /* Define input audio */
    float** inSigMIC;
    inSigMIC = (float**)malloc2d(nMics, sigLen, sizeof(float));
    rand_m1_1(FLATTEN2D(inSigMIC), nMics*sigLen);

    /* Main loop */
    float **inSigMIC_block, **outSigBIN_block, **outSigBIN;
    inSigMIC_block = (float**)malloc2d(nMics, blocksize, sizeof(float));
    outSigBIN_block = (float**)malloc2d(NUM_EARS, blocksize, sizeof(float));
    outSigBIN = (float**)malloc2d(NUM_EARS, sigLen, sizeof(float));
    for(i=0; i < (int)((float)sigLen/(float)blocksize); i++){
        /* Copy input to buffer */
        for(ch=0; ch<nMics; ch++)
            memcpy(inSigMIC_block[ch], &inSigMIC[ch][i*blocksize], blocksize*sizeof(float));

        /* Analysis/synthesis */
        hades_analysis_apply(hAna, inSigMIC_block, nMics, blocksize, hPCon, hSCon);
        hades_synthesis_apply(hSyn, hPCon, hSCon, NUM_EARS, blocksize, outSigBIN_block);

        /* Copy buffer to output */
        for(ch=0; ch<NUM_EARS; ch++)
            memcpy(&outSigBIN[ch][i*blocksize], outSigBIN_block[ch], blocksize*sizeof(float));
    }

    /* Clean-up */
    hades_analysis_destroy(&hAna);
    hades_param_container_destroy(&hPCon);
    hades_signal_container_destroy(&hSCon);
    hades_synthesis_destroy(&hSyn);
    free(grid_dirs_deg);
    free(inSigMIC);
    free(inSigMIC_block);
    free(outSigBIN_block);
    free(outSigBIN);
}


#endif /* SAF_ENABLE_HADES_MODULE */
