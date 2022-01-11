/*
 * This file is part of the saf_tracker module unit tests.
 * Copyright (c) 2020-2021 - Leo McCormack
 *
 * The saf_tracker module is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * The saf_tracker module is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * See <http://www.gnu.org/licenses/> for a copy of the GNU General Public
 * License.
 */

/**
 * @file test__tracker_module.c
 * @brief Unit tests for the SAF tracker module
 * @author Leo McCormack
 * @date 27.04.2020
 * @license GNU GPLv2
 */

#include "saf_test.h"

#ifdef SAF_ENABLE_TRACKER_MODULE

void test__tracker3d(void){
    int hop, i, j, k, nSH, nGrid, rand_idx, dropouts;
    int inds[2];
    int* target_IDs;
    void* hT3d, *hMUSIC;
    float measNoiseSD_deg, noiseSpecDen_deg, scale, rand01;
    float est_dirs_deg[2][2], est_dirs_xyz[2][3];
    float* grid_dirs_deg;
    float** insigs, **inputSH, **inputSH_noise, **inputSH_hop, **Y, **Cx, **V, **Vn;
    float_complex** Vn_cmplx;
    int nTargets;
    float *target_dirs_xyz, *target_var_xyz;

    /* Test configuration */
    const float acceptedTolerance = 0.005f;
    const int order = 2;
    const float fs = 48e3;
    const int hopsize = 128;
    const int sigLen = (int)fs*5;
    const int nSources = 2; /* cannot be changed, hard-coded for 2 */
    const float src_dirs_deg[2][2] = { {-35.0f, 30.0f}, {120.0f, 0.0f} };

    /* Configure the tracker */
    tracker3d_config tpars;
    /* Number of Monte Carlo samples/particles. The more complex the
     * distribution is, the more particles required (but also, the more
     * computationally expensive the tracker becomes). */
    tpars.Np = 20;
    tpars.ARE_UNIT_VECTORS = 1;
    tpars.maxNactiveTargets = 2; /* about 2 higher than expected is good */
    /* Likelihood of an estimate being noise/clutter */
    tpars.noiseLikelihood = 0.2f; /* between [0..1] */
    /* Measurement noise - e.g. to assume that estimates within the range +/-20
     * degrees belong to the same target, set SDmnoise_deg = 20 */
    measNoiseSD_deg = 20.0f;
    tpars.measNoiseSD = 1.0f-cosf(measNoiseSD_deg*SAF_PI/180.0f); /* Measurement noise standard deviation */
    /* Noise spectral density - not fully understood. But it influences the
     * smoothness of the target tracks */
    noiseSpecDen_deg = 1.0f;
    tpars.noiseSpecDen = 1.0f-cosf(noiseSpecDen_deg*SAF_PI/180.0f);  /* Noise spectral density */
    /* FLAG - whether to allow for multiple target deaths in the same tracker
     * prediction step */
    tpars.ALLOW_MULTI_DEATH = 1;
    /* Probability of birth and death */
    tpars.init_birth = 0.5f; /* value between [0 1] - Prior probability of birth */
    tpars.alpha_death = 20.0f; /* always >= 1; 1 is good */      /* 20-> means death is very unlikely... */
    tpars.beta_death = 1.0f; /* always >= 1; 1 is good */
    /* Elapsed time (in seconds) between observations */
    tpars.dt = 1.0f/(fs/(float)hopsize); /* Hop length of frames */
    /* Whether or not to allow multiple active sources for each update */
    /* Real-time tracking is based on the particle with highest weight. A
     * one-pole averaging filter is used to smooth the weights over time. */
    tpars.W_avg_coeff = 0.5f;
    /* Force kill targets that are close to another target. In these cases, the
     * target that has been 'alive' for the least amount of time, is killed */
    tpars.FORCE_KILL_TARGETS = 1;
    tpars.forceKillDistance = 0.2f;
    /* Mean position priors x,y,z (assuming directly in-front) */
    tpars.M0[0] = 1.0f; tpars.M0[1] = 0.0f; tpars.M0[2] = 0.0f;
    /* Mean Velocity priors x,y,z (assuming stationary) */
    tpars.M0[3] = 0.0f; tpars.M0[4] = 0.0f; tpars.M0[5] = 0.0f;
    /* Target velocity - e.g. to assume that a target can move 20 degrees in
    * two seconds along the horizontal, set V_azi = 20/2 */
    const float Vazi_deg = 3.0f;  /* Velocity of target on azimuthal plane */
    const float Vele_deg = 3.0f;  /* Velocity of target on median plane */
    memset(tpars.P0, 0, 6*6*sizeof(float));
    /* Variance PRIORs of estimates along the x,y,z axes, respectively. Assuming
     * coordinates will lay on the unit sphere +/- x,y,z, so a range of 2, and
     * hence a variance of 2^2: */
    tpars.P0[0][0] = 4.0f; tpars.P0[1][1] = 4.0f; tpars.P0[2][2] = 4.0f;
    /* Velocity PRIORs of estimates the x,y,z axes */
    tpars.P0[3][3] = 1.0f-cosf(Vazi_deg*SAF_PI/180.0f); /* x */
    tpars.P0[4][4] = tpars.P0[3][3];                    /* y */
    tpars.P0[5][5] = 1.0f-cosf(Vele_deg*SAF_PI/180.0f); /* z */
    /* PRIOR probabilities of noise. (Assuming the noise is uniformly
     * distributed in the entire spatial grid). */
    tpars.cd = 1.0f/(4.0f*SAF_PI);

    /* Create tracker */
    tracker3d_create(&hT3d, tpars);

    /* Create spherical harmonic input signals */
    insigs = (float**)malloc2d(nSources, sigLen, sizeof(float));
    rand_m1_1(FLATTEN2D(insigs), nSources*sigLen);
    nSH = ORDER2NSH(order);
    Y = (float**)malloc2d(nSH, nSources, sizeof(float));
    getRSH(order, (float*)src_dirs_deg, nSources, FLATTEN2D(Y));
    scale = 1.0f/(float)nSources;
    utility_svsmul(FLATTEN2D(Y), &scale, nSH*nSources, NULL);
    inputSH = (float**)malloc2d(nSH, sigLen, sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nSH, sigLen, nSources, 1.0f,
                FLATTEN2D(Y), nSources,
                FLATTEN2D(insigs), sigLen, 0.0f,
                FLATTEN2D(inputSH), sigLen);

    /* Add some noise */
    inputSH_noise = (float**)malloc2d(nSH, sigLen, sizeof(float));
    rand_m1_1(FLATTEN2D(inputSH_noise), nSH*sigLen);
    scale = 0.05f;
    utility_svsmul(FLATTEN2D(inputSH_noise), &scale, nSH*sigLen, NULL);
    cblas_saxpy(nSH*sigLen, 1.0f, FLATTEN2D(inputSH_noise), 1, FLATTEN2D(inputSH), 1);

    /* Create DoA estimator */
    nGrid = 240; /* number of points in a t-design of degree 21 */
    grid_dirs_deg = (float*)__Tdesign_degree_21_dirs_deg;
    sphMUSIC_create(&hMUSIC, order, grid_dirs_deg, nGrid);

    /* Memory allocations */
    inputSH_hop = (float**)malloc2d(nSH, hopsize, sizeof(float));
    Cx = (float**)malloc2d(nSH, nSH, sizeof(float));
    V = (float**)malloc2d(nSH, nSH, sizeof(float));
    Vn = (float**)malloc2d(nSH, (nSH-nSources), sizeof(float)); /* noise subspace */
    Vn_cmplx = (float_complex**)malloc2d(nSH, (nSH-nSources), sizeof(float_complex));
    target_dirs_xyz = NULL;
    target_var_xyz = NULL;
    target_IDs = NULL;

    /* Loop over hops */
    dropouts = 0;
    for(hop=0; hop<(int)((float)sigLen/(float)hopsize); hop++){
        /* Grab current hop */
        for(i=0; i<nSH; i++)
            memcpy(inputSH_hop[i], &inputSH[i][hop*hopsize], hopsize*sizeof(float));

        /* Eigenvalue decomposition and truncation of eigen vectors to obtain
         * noise subspace (based on source number) */
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nSH, nSH, hopsize, 1.0f,
                    FLATTEN2D(inputSH_hop), hopsize,
                    FLATTEN2D(inputSH_hop), hopsize, 0.0f,
                    FLATTEN2D(Cx), nSH);
        utility_sseig(NULL, FLATTEN2D(Cx), nSH, 1, FLATTEN2D(V), NULL, NULL);
        for(i=0; i<nSH; i++)
            for(j=0, k=nSources; j<nSH-nSources; j++, k++)
                Vn[i][j] = V[i][k];
        for(i=0; i<nSH; i++)
            for(j=0; j<nSH-nSources; j++)
                Vn_cmplx[i][j] = cmplxf(Vn[i][j], 0.0f);

        /* DoA estimation */
        sphMUSIC_compute(hMUSIC, FLATTEN2D(Vn_cmplx), nSources, NULL, (int*)inds);
        est_dirs_deg[0][0] = grid_dirs_deg[inds[0]*2+0];
        est_dirs_deg[0][1] = grid_dirs_deg[inds[0]*2+1];
        est_dirs_deg[1][0] = grid_dirs_deg[inds[1]*2+0];
        est_dirs_deg[1][1] = grid_dirs_deg[inds[1]*2+1];
        unitSph2cart((float*)est_dirs_deg, nSources, 1, (float*)est_dirs_xyz);

        /* Pick an estimate at random */
        rand_0_1(&rand01, 1);
        rand_idx = (int)(rand01*(float)nSources);

        /* Feed tracker */
        tracker3d_step(hT3d, (float*)&est_dirs_xyz[rand_idx], 1, &target_dirs_xyz, &target_var_xyz, &target_IDs, &nTargets);

        /* Give the tracker a couple of steps, and then assert that it is keeping track of these two targets */
        if(hop>10){
            if(nTargets==nSources){
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[0][0] - target_dirs_xyz[0*3+0]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[0][0] - target_dirs_xyz[1*3+0]) <= acceptedTolerance);
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[0][1] - target_dirs_xyz[0*3+1]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[0][1] - target_dirs_xyz[1*3+1]) <= acceptedTolerance);
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[0][2] - target_dirs_xyz[0*3+2]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[0][2] - target_dirs_xyz[1*3+2]) <= acceptedTolerance);
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[1][0] - target_dirs_xyz[0*3+0]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[1][0] - target_dirs_xyz[1*3+0]) <= acceptedTolerance);
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[1][1] - target_dirs_xyz[0*3+1]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[1][1] - target_dirs_xyz[1*3+1]) <= acceptedTolerance);
                TEST_ASSERT_TRUE( fabsf(est_dirs_xyz[1][2] - target_dirs_xyz[0*3+2]) <= acceptedTolerance ||
                                  fabsf(est_dirs_xyz[1][2] - target_dirs_xyz[1*3+2]) <= acceptedTolerance);
            }
            else
                dropouts++; /* Should be very unlikely, (as the probably of death set to be so low), but it can still happen... */
        }
    }
    TEST_ASSERT_TRUE(dropouts<12);

    /* Clean-up */
    tracker3d_destroy(&hT3d);
    sphMUSIC_destroy(&hMUSIC);
    free(target_dirs_xyz);
    free(target_var_xyz);
    free(target_IDs);
    free(insigs);
    free(inputSH);
    free(inputSH_noise);
    free(inputSH_hop);
    free(Y);
    free(Cx);
    free(V);
    free(Vn);
    free(Vn_cmplx);
}

#endif /* SAF_ENABLE_TRACKER_MODULE */
