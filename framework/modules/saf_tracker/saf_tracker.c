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
 * @ingroup Tracker
 * @brief Particle filtering based 3D multi-target tracker (#SAF_TRACKER_MODULE)
 *
 * Based on the RBMCDA [1] MATLAB toolbox (GPLv2 license) by Simo Sa"rkka" and
 * Jouni Hartikainen (Copyright (C) 2003-2008):
 *     https://users.aalto.fi/~ssarkka/#softaudio
 *
 * More information regarding this specific implementation can be found in [2]
 *
 * @see [1] Sa"rkka", S., Vehtari, A. and Lampinen, J., 2004, June. Rao-
 *          Blackwellized Monte Carlo data association for multiple target
 *          tracking. In Proceedings of the seventh international conference on
 *          information fusion (Vol. 1, pp. 583-590). I.
 * @see [2] McCormack, L., Politis, A. Sa"rkka", S., and Pulkki, V., 2021.
 *          Real-Time Tracking of Multiple Acoustical Sources Utilising
 *          Rao-Blackwellised Particle Filtering. In 29th European Signal
 *          Processing Conference (EUSIPCO), (pp. 206-210).
 *
 * @author Leo McCormack
 * @date 12.08.2020
 * @license GNU GPLv2
 */

#include "saf_tracker_internal.h"
#include "saf_tracker.h"

#ifdef  SAF_ENABLE_TRACKER_MODULE

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

    /* Store user configuration */
    pData->tpars = tpars;

    /* Parameter checking */
    pData->tpars.Np = SAF_CLAMP(pData->tpars.Np, 1, TRACKER3D_MAX_NUM_PARTICLES);
    saf_assert(pData->tpars.ARE_UNIT_VECTORS == 0 || pData->tpars.ARE_UNIT_VECTORS == 1, "ARE_UNIT_VECTORS is a bool");
    pData->tpars.init_birth = SAF_CLAMP(pData->tpars.init_birth, 0.0f, 0.99f);
    pData->tpars.alpha_death = SAF_CLAMP(pData->tpars.alpha_death, 1.0f, 20.0f);
    pData->tpars.beta_death = SAF_CLAMP(pData->tpars.beta_death, 1.0f, 20.0f);
    pData->tpars.dt = SAF_MAX(pData->tpars.dt, 0.0001f);
    pData->tpars.cd = SAF_MAX(pData->tpars.cd, 0.0001f);
    pData->tpars.W_avg_coeff = SAF_CLAMP(pData->tpars.W_avg_coeff, 0.0f, 0.99f);
    pData->tpars.noiseSpecDen = SAF_MAX(pData->tpars.noiseSpecDen, 0.0001f);
    pData->tpars.noiseLikelihood = SAF_CLAMP(pData->tpars.noiseLikelihood, 0.0f, 0.99f);
    pData->tpars.measNoiseSD = SAF_MAX(pData->tpars.measNoiseSD, 0.001f);

    /* Measurement noise PRIORs along the x,y,z axis, respectively  */
    sd_xyz = pData->tpars.measNoiseSD;
    memset(pData->R, 0, 3*3*sizeof(float));
    pData->R[0][0] = powf(sd_xyz,2.0f);
    pData->R[1][1] = powf(sd_xyz,2.0f);
    pData->R[2][2] = powf(sd_xyz,2.0f);

    /* Noise spectral density along x, y, z axis qx,y,z which, (in combination
     * with sd_xyz), dictates how smooth the target tracks are. */
    q_xyz = pData->tpars.noiseSpecDen;

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
    lti_disc((float*)F, 6, 6, NULL, (float*)Qc, pData->tpars.dt, (float*)pData->A, (float*)pData->Q);
    memset(pData->H, 0, 3*6*sizeof(float));
    pData->H[0][0] = 1.0f;
    pData->H[1][1] = 1.0f;
    pData->H[2][2] = 1.0f;

    /* Create kalman filter helper struct */
    kf_update6_create(&(pData->hKF6));

    /* Create particles */
    pData->SS = malloc1d(pData->tpars.Np * sizeof(voidPtr));
    pData->SS_resamp = malloc1d(pData->tpars.Np * sizeof(voidPtr));
    pData->W0 = 1.0f/(float)pData->tpars.Np;
    for(i=0; i<pData->tpars.Np; i++){
        tracker3d_particleCreate(&(pData->SS[i]), pData->W0, pData->tpars.dt);
        tracker3d_particleCreate(&(pData->SS_resamp[i]), pData->W0, pData->tpars.dt);
    }
    
    /* Event starting values */
    for(i=0; i<TRACKER3D_MAX_NUM_EVENTS; i++){
        pData->evta[i] = -1;
        tracker3d_particleCreate(&(pData->str[i]), pData->W0, pData->tpars.dt);
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

        kf_update6_destroy(&(pData->hKF6));

        for(i=0; i<pData->tpars.Np; i++){
            tracker3d_particleDestroy(&pData->SS[i]);
            tracker3d_particleDestroy(&pData->SS_resamp[i]);
        }
        free(pData->SS);
        free(pData->SS_resamp);

        for(i=0; i<TRACKER3D_MAX_NUM_EVENTS; i++)
            tracker3d_particleDestroy(&pData->str[i]);

        free(pData);
        pData = NULL;
    }
}

void tracker3d_reset
(
    void* const hT3d
)
{
    tracker3d_data *pData = (tracker3d_data*)(hT3d);
    int i;
 
    pData->incrementTime = 0;
    for(i=0; i<pData->tpars.Np; i++)
        tracker3d_particleReset(pData->SS[i]);
}
  
void tracker3d_step
(
    void* const hT3d,
    float* newObs_xyz,
    int nObs,
    float** target_pos_xyz,
    float** target_var_xyz,
    int** target_IDs,
    int* nTargets
)
{
    tracker3d_data *pData = (tracker3d_data*)(hT3d);
    int i, kt, ob, maxIdx, nt;
    float Neff;
    int s[TRACKER3D_MAX_NUM_PARTICLES];
    MCS_data* S_max;
#if 0
    int nt2;
    float w_sum;
    MCS_data* S_tmp;
#endif
#ifdef TRACKER_VERBOSE
    char c_str[256], tmp[256];
    memset(c_str, 0, 256*sizeof(char));
#endif

    pData->incrementTime++;
 
    /* Loop over measurements */
    if(nObs>0 && newObs_xyz!=NULL){
        for(ob=0; ob<nObs; ob++){
            /* Predict and update steps */
            for (kt = 0; kt < pData->incrementTime; kt++)
                tracker3d_predict(hT3d, 1);
            tracker3d_update(hT3d, &newObs_xyz[ob*3], pData->incrementTime);
            pData->incrementTime = 0;

            /* Resample if needed */
            Neff = eff_particles(pData->SS, pData->tpars.Np);
            if (Neff < (float)pData->tpars.Np/4.0f){
#ifdef TRACKER_VERBOSE
                printf("%s\n", "Resampling");
#endif
                maxIdx = tracker3d_getMaxParticleIdx(hT3d);
                for(i=0; i<pData->tpars.Np; i++)
                    s[i] = maxIdx;
                //resampstr(pData->SS, pData->tpars.Np, s);

                for(i=0; i<pData->tpars.Np; i++)
                    tracker3d_particleCopy(pData->SS[s[i]], pData->SS_resamp[i]);
                for(i=0; i<pData->tpars.Np; i++){
                    tracker3d_particleCopy(pData->SS_resamp[i], pData->SS[i]);
                    ((MCS_data*)pData->SS[i])->W = ((MCS_data*)pData->SS[i])->W0;
                }
            }

            /* Apply (optional) temporal smoothing of the particle importance weights */
            if(pData->tpars.W_avg_coeff>0.0001f){
                for(i=0; i<pData->tpars.Np; i++){
                    ((MCS_data*)pData->SS[i])->W = ((MCS_data*)pData->SS[i])->W * (1.0f-pData->tpars.W_avg_coeff) +
                           ((MCS_data*)pData->SS[i])->W_prev * pData->tpars.W_avg_coeff;
                    ((MCS_data*)pData->SS[i])->W_prev = ((MCS_data*)pData->SS[i])->W;
                }
            }
        }
    } 

    /* Find most dominant particle.. */
    maxIdx = tracker3d_getMaxParticleIdx(hT3d);
    S_max = (MCS_data*)pData->SS[maxIdx];
 
    /* Output */
    if(S_max->nTargets==0){
        free((*target_pos_xyz));
        free((*target_var_xyz));
        free((*target_IDs));
        (*target_pos_xyz) = NULL;
        (*target_var_xyz) = NULL;
        (*target_IDs) = NULL;
        (*nTargets) = 0;
#ifdef TRACKER_VERBOSE
        printf("%s\n", "No targets");
#endif
    }
    else{
        (*target_pos_xyz) = realloc1d((*target_pos_xyz), S_max->nTargets*3*sizeof(float));
        (*target_var_xyz) = realloc1d((*target_var_xyz), S_max->nTargets*3*sizeof(float));
        (*target_IDs) = realloc1d((*target_IDs), S_max->nTargets*sizeof(int));
        (*nTargets) = S_max->nTargets;

        /* Loop over targets */
        for(nt=0; nt<S_max->nTargets; nt++){
#ifdef TRACKER_VERBOSE
            sprintf(tmp, "ID_%d: [%.5f,%.5f,%.5f] ", S_max->targetIDs[nt], S_max->M[nt].m0, S_max->M[nt].m1, S_max->M[nt].m2);
            strcat(c_str, tmp);
#endif
            /* Target IDs are based on those defined by the most dominant particle */
            (*target_IDs)[nt] = S_max->targetIDs[nt];
            (*target_pos_xyz)[nt*3]   = S_max->M[nt].m0;
            (*target_pos_xyz)[nt*3+1] = S_max->M[nt].m1;
            (*target_pos_xyz)[nt*3+2] = S_max->M[nt].m2;
            (*target_var_xyz)[nt*3]   = S_max->P[nt].p00;
            (*target_var_xyz)[nt*3+1] = S_max->P[nt].p11;
            (*target_var_xyz)[nt*3+2] = S_max->P[nt].p22;

# if 0
            /* Apply the corresponding importance weight */
            w_sum = S_max->W;
            (*target_pos_xyz)[nt*3]   *= S_max->W;
            (*target_pos_xyz)[nt*3+1] *= S_max->W;
            (*target_pos_xyz)[nt*3+2] *= S_max->W;
            (*target_var_xyz)[nt*3]   *= S_max->W;
            (*target_var_xyz)[nt*3+1] *= S_max->W;
            (*target_var_xyz)[nt*3+2] *= S_max->W;

            /* Loop over all of the other particles - importance sampling */
            for(int p = 0; p<pData->tpars.Np; p++){
                if(p!=maxIdx){
                    S_tmp = (MCS_data*)pData->SS[p];
                    for(nt2=0; nt2<S_tmp->nTargets; nt2++){
                        if((*target_IDs)[nt] == S_tmp->targetIDs[nt2]){
                            w_sum += S_tmp->W;
                            (*target_pos_xyz)[nt*3]   += (S_tmp->M[nt2].m0 * S_tmp->W);
                            (*target_pos_xyz)[nt*3+1] += (S_tmp->M[nt2].m1 * S_tmp->W);
                            (*target_pos_xyz)[nt*3+2] += (S_tmp->M[nt2].m2 * S_tmp->W);
                            (*target_var_xyz)[nt*3]   += (S_tmp->P[nt2].p00 * S_tmp->W);
                            (*target_var_xyz)[nt*3+1] += (S_tmp->P[nt2].p11 * S_tmp->W);
                            (*target_var_xyz)[nt*3+2] += (S_tmp->P[nt2].p22 * S_tmp->W);
                        }
                    }
                }
            }

            /* Renormalise based on the combined importance weights */
            (*target_pos_xyz)[nt*3]   /= w_sum;
            (*target_pos_xyz)[nt*3+1] /= w_sum;
            (*target_pos_xyz)[nt*3+2] /= w_sum;
            (*target_var_xyz)[nt*3]   /= w_sum;
            (*target_var_xyz)[nt*3+1] /= w_sum;
            (*target_var_xyz)[nt*3+2] /= w_sum;
#endif
        }
#ifdef TRACKER_VERBOSE
        printf("%s\n", c_str);
#endif
    } 
}

#endif /* SAF_ENABLE_TRACKER_MODULE */
