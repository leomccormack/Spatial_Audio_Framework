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
 * Based on the RBMCDA [1] Matlab toolbox (GPLv2 license) by Simo Särkkä and
 * Jouni Hartikainen (Copyright (C) 2003-2008):
 *     https://users.aalto.fi/~ssarkka/#softaudio
 *
 * And also inspired by the work of Sharath Adavanne, Archontis Politis, Joonas
 * Nikunen, and Tuomas Virtanen (GPLv2 license):
 *     https://github.com/sharathadavanne/multiple-target-tracking
 *
 * @see [1] Särkkä, S., Vehtari, A. and Lampinen, J., 2004, June. Rao-
 *          Blackwellized Monte Carlo data association for multiple target
 *          tracking. In Proceedings of the seventh international conference on
 *          information fusion (Vol. 1, pp. 583-590). I.
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
    float sd_xyz, q_xyz;
    float Qc[6][6];

    assert(tpars.Np<=TRACKER3D_MAX_NUM_PARTICLES);

    /* Store user configuration */
    pData->tpars = tpars;

    /* Measurement noise PRIORs along the x,y,z axis, respectively  */
    sd_xyz = tpars.measNoiseSD;
    memset(pData->R, 0, 3*3*sizeof(float));
    pData->R[0][0] = powf(sd_xyz,2.0f);
    pData->R[1][1] = powf(sd_xyz,2.0f);
    pData->R[2][2] = powf(sd_xyz,2.0f);

    /* Noise spectral density along x, y, z axis qx,y,z which, (in combination
     * with sd_xyz), dictates how smooth the target tracks are. */
    q_xyz = tpars.noiseSpecDen;

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
    pData->SS_resamp = malloc1d(tpars.Np * sizeof(voidPtr));
    pData->W0 = 1.0f/(float)tpars.Np;
    for(i=0; i<tpars.Np; i++){
        tracker3d_particleCreate(&(pData->SS[i]), pData->W0, tpars.dt);
        tracker3d_particleCreate(&(pData->SS_resamp[i]), pData->W0, tpars.dt);
    }
    
    /* Event starting values */
    for(i=0; i<TRACKER3D_MAX_NUM_EVENTS; i++){
        pData->evta[i] = NULL;
        tracker3d_particleCreate(&(pData->str[i]), pData->W0, tpars.dt);
    }
    pData->incrementTime = 0;
}

void tracker3d_destroy
(
    void ** const phT3d
)
{
    tracker3d_data *pData = (tracker3d_data*)(*phT3d);
    int i;

    if (pData != NULL) {

        for(i=0; i<pData->tpars.Np; i++){
            tracker3d_particleDestroy(&pData->SS[i]);
            tracker3d_particleDestroy(&pData->SS_resamp[i]);
        }
        free(pData->SS);
        free(pData->SS_resamp);

        for(i=0; i<TRACKER3D_MAX_NUM_EVENTS; i++){
            free(pData->evta[i]);
            tracker3d_particleDestroy(&pData->str[i]);
        }

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
    int** target_IDs,
    int* nTargets
)
{
    tracker3d_data *pData = (tracker3d_data*)(hT3d);
    int i, kt, ob, maxIdx, nt;
    float Neff, maxVal, tmpW;
    int s[TRACKER3D_MAX_NUM_PARTICLES];
    MCS_data* S_max;
#ifdef TRACKER_VERBOSE
    char c_str[256], tmp[256];
#endif

    pData->incrementTime++; 

    /* Loop over measurements */
    for(ob=0; ob<nObs; ob++){
#ifdef TRACKER_VERY_VERBOSE
        printf("%s\n", "Update step");
#endif
        tracker3d_update(hT3d, &newObs_xyz[ob*3], pData->incrementTime);

#ifdef TRACKER_VERY_VERBOSE
		printf("%s\n", "Prediction step");
#endif
		for (kt = 0; kt < pData->incrementTime; kt++)
			tracker3d_predict(hT3d, 1);

        pData->incrementTime = 0;

        /* Resample if needed */
        Neff = eff_particles(pData->SS, pData->tpars.Np);
        if (Neff < (float)pData->tpars.Np/4.0f){
#ifdef TRACKER_VERBOSE
            printf("%s\n", "Resampling");
#endif
            resampstr(pData->SS, pData->tpars.Np, s);
            for(i=0; i<pData->tpars.Np; i++)
                tracker3d_particleCopy(pData->SS[s[i]], pData->SS_resamp[i]);
            for(i=0; i<pData->tpars.Np; i++){
                tracker3d_particleCopy(pData->SS_resamp[i], pData->SS[i]);
                ((MCS_data*)pData->SS[i])->W = ((MCS_data*)pData->SS[i])->W0;
            }
        }

        /* Apply (optional) temporal smoothing of the particle importance weights */
        if(pData->tpars.W_avg_coeff>0.0f){
            for(i=0; i<pData->tpars.Np; i++){
                ((MCS_data*)pData->SS[i])->W = ((MCS_data*)pData->SS[i])->W * (1.0f-pData->tpars.W_avg_coeff) +
                       ((MCS_data*)pData->SS[i])->W_prev * pData->tpars.W_avg_coeff;
                ((MCS_data*)pData->SS[i])->W_prev = ((MCS_data*)pData->SS[i])->W;
               // ((MCS_data*)pData->SS[i])->W = tmpW;
            }
        }
    }

    /* Find most significant particle.. */
    maxVal = FLT_MIN;
    maxIdx = -1;
    for(i=0; i<pData->tpars.Np; i++){
        if(maxVal<((MCS_data*)pData->SS[i])->W){
            maxVal = ((MCS_data*)pData->SS[i])->W;
            maxIdx = i;
        }
    }
    assert(maxIdx!=-1);
    S_max = (MCS_data*)pData->SS[maxIdx];
 
    /* Output */
    if(S_max->nTargets==0){
        (*target_xyz) = NULL;
        (*target_IDs) = NULL;
        (*nTargets) = 0;
#ifdef TRACKER_VERBOSE
        printf("%s\n", "No targets");
#endif
    }
    else{
        (*target_xyz) = realloc1d((*target_xyz), S_max->nTargets*3*sizeof(float));
        (*target_IDs) = realloc1d((*target_IDs), S_max->nTargets*sizeof(int));
        (*nTargets) = S_max->nTargets;

#ifdef TRACKER_VERBOSE
        memset(c_str, 0, 256*sizeof(char)); 
#endif
        for(nt=0; nt<S_max->nTargets; nt++){
#ifdef TRACKER_VERBOSE
            sprintf(tmp, "ID_%d: [%.5f,%.5f,%.5f] ", S_max->targetIDs[nt], S_max->M[nt].m0, S_max->M[nt].m1, S_max->M[nt].m2);
            strcat(c_str, tmp);
#endif 
            (*target_IDs)[nt] = S_max->targetIDs[nt];
            memcpy(&(*target_xyz)[nt*3], S_max->M[nt].M, 3*sizeof(float));
        }
#ifdef TRACKER_VERBOSE
        printf("%s\n", c_str);
#endif
    }
}

