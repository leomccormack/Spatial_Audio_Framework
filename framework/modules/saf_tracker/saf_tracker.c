/*
 * This file is part of the saf_tracker module.
 * Copyright (c) 2020 - Leo McCormack
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
 * @file saf_tracker.c
 * @brief Particle filtering based tracker
 *
 * Based on the RBMCDA Matlab toolbox (GPLv2 license) by Simo Särkkä and Jouni
 * Hartikainen (Copyright (C) 2003-2008):
 *     https://users.aalto.fi/~ssarkka/#softaudio
 *
 * And also inspired by the work of Sharath Adavanne, Archontis Politis, Joonas
 * Nikunen, and Tuomas Virtanen (GPLv2 license):
 *     https://github.com/sharathadavanne/multiple-target-tracking
 *     
 * @author Leo McCormack
 * @date 12.08.2020
 */

#include "saf_tracker_internal.h"
#include "saf_tracker.h"

void tracker3d_create
(
    void ** const phT3d,
    tracker3d_config tpars 
)
{
    tracker3d_data* pData = (tracker3d_data*)malloc1d(sizeof(tracker3d_data));
    *phT3d = (void*)pData;
    int i;
    float sd_xyz, q_xyz, W0;
    float Qc[6][6];

    /* Store user configuration */
    pData->tpars = tpars;

    /* Measurement noise PRIORs along the x,y,z axis, respectively  */
    sd_xyz = 1.0f-cosf(tpars.measNoiseSD_deg*SAF_PI/180.0f);
    memset(pData->R, 0, 3*3*sizeof(float));
    pData->R[0][0] = powf(sd_xyz,2.0f);
    pData->R[1][1] = powf(sd_xyz,2.0f);
    pData->R[2][2] = powf(sd_xyz,2.0f);

    /* Noise spectral density along x, y, z axis qx,y,z in combination with
     * sd_xyz decides how smooth the target tracks are. */
    q_xyz = 1.0f-cosf(tpars.noiseSpecDen_deg*SAF_PI/180.0f);

    /* Dynamic and measurement models */
    const float F[6][6] =
     {  {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}  };
    memset(Qc, 0, 6*6*sizeof(float));
    Qc[3][3] = q_xyz;
    Qc[4][4] = q_xyz;
    Qc[5][5] = q_xyz;
    lti_disc((float*)F, 6, 6, NULL, (float*)Qc, tpars.dt, (float*)pData->A, (float*)pData->Q);
    memset(pData->H, 0, 3*6*sizeof(float));
    pData->H[0][0] = 1.0f;
    pData->H[1][1] = 1.0f;
    pData->H[2][2] = 1.0f;

    /* Create particles */
    pData->SS = malloc1d(tpars.Np * sizeof(voidPtr));
    pData->W0 = 1.0f/(float)tpars.Np;
    for(i=0; i<tpars.Np; i++)
        tracker3d_particleCreate(&(pData->SS[i]), pData->W0, tpars.dt);

    /* Possible combinations of active sources */
    if(tpars.MULTI_ACTIVE){
        assert(0); // not implemented yet!
    }

    /* Starting values */
    for(i=0; i<TRACKER3D_MAX_NUM_EVENTS; i++){
        pData->evta[i] = NULL;
        pData->str[i] = NULL;
    }

    pData->incrementTime = 0;
}

void tracker3d_destroy
(
    void ** const phT3d
)
{
    tracker3d_data *pData = (tracker3d_data*)(*phT3d);

    if (pData != NULL) {

        free(pData);
        pData = NULL;
    }
}
  
void tracker3d_step
(
    void* const hT3d,
    float* newObs_xyz,
    int nObs,
    float** target_xyz,
    int* nTargets
)
{
    tracker3d_data *pData = (tracker3d_data*)(hT3d);
    int ob;

    pData->incrementTime = 0;
    pData->incrementTime++;

#ifdef TRACKER_VERBOSE
    printf("%s\n", "Prediction step");
#endif
    tracker3d_predict(hT3d, pData->incrementTime);
    pData->incrementTime = 0;

#ifdef TRACKER_VERBOSE
    printf("%s\n", "Update step");
#endif
    tracker3d_update(hT3d, newObs_xyz[0], pData->incrementTime);

    for(ob=0; ob<nObs; ob++){

    }


}

