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
 * @file saf_tracker_internal.c
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
 * @date 04.07.2020
 */

#include "saf_tracker.h"
#include "saf_tracker_internal.h"

/* https://stackoverflow.com/questions/24460776/fast-incomplete-gamma-function */
float poisson_z(float x, float mu){
    static const float twoThirds = 2.0f/3.0f;
    float w = sqrt((x+0.5f)/mu) - 1.0f;
    float coeff = w>=0.0f ? 0.085f : 0.15f;
    return (x-mu+twoThirds)/sqrtf(mu*(1.0f+w*(0.68f+w*coeff)));
}

float gamma_cdf(float x, float gam, float beta)
{
    /* Convert to standard form */
    x = (x)/beta;

    /* Compute the probability using the imcomplete gamma function */
    return poisson_z(x,gam);
}

void tracker3d_predict
(
    voidPtr* SS,
    tracker3d_config* tpars,
    int Tinc
)
{
    int i, j, nDead;
    int* dead;
    float dt0, dt1, p_death;
    MCS_data* S;

    nDead = 0;
    dead = NULL;

    for (i=0; i<tpars->Np; i++){
       // dead=[];
        S = (MCS_data*)SS[i];

        /* Loop over targets */
        for (j=0; j<S->nTargets; j++){

            /* No target has died yet or multiple targets are allowed to die in
             * one prediction step */
            if (nDead==0 || tpars->ALLOW_MULTI_DEATH){
                /* Probability of death */
                dt0 = (float)S->Tcount[j] * S->dt;
                dt1 = dt0 + S->dt * (float)Tinc;
                if (dt0 == 0)
                    p_death = gamma_cdf(dt1, tpars->alpha_death, tpars->beta_death);
                else
                    p_death = 1 - (1-gamma_cdf(dt1, tpars->alpha_death,tpars->beta_death)) /
                        (1-gamma_cdf(dt0,tpars->alpha_death,tpars->beta_death));


                /* Force probability of death to 1, if this target is too close
                   another target that has been alive longer. */
//                if tpars.FORCE_KILL_TARGETS
//                    for k=1:S{i}.nTargets
//                        if k~=j
//                            %dist = sqrt(sum((S{i}.M{j}(1:3,1)-S{i}.M{k}(1:3,1)).^2));
//                            %angle_diff = acos( dot(S{i}.M{j}(1:3,1), S{i}.M{k}(1:3,1))./(norm(S{i}.M{j}(1:3,1))*norm(S{i}.M{k}(1:3,1) )));
//                            angle_diff = atan2(norm(cross(S{i}.M{j}(1:3,1),S{i}.M{k}(1:3,1))),dot(S{i}.M{j}(1:3,1),S{i}.M{k}(1:3,1)));
//                            if angle_diff < tpars.forceKillAngle_rad && S{i}.Tcount(j) <= S{i}.Tcount(k)
//                                p_death = 1;
//                            end
//                        end
//                    end
//                end

//                if (rand < p_death)
//                    /* Target dies */
//                    dead = [dead j];
//                end
            }

//            if tpars.ALLOW_MULTI_DEATH
//                if all(j ~= dead)
//                    % Kalman Filter prediction for the target if alive
//                    [S{i}.M{j},S{i}.P{j}] = kf_predict(S{i}.M{j},S{i}.P{j},tpars.A,tpars.Q,B,u);
//                end
//            else
//                if  (isempty(dead)) || (j ~= dead)
//                    % Kalman Filter prediction for the target if alive
//                    [S{i}.M{j},S{i}.P{j}] = kf_predict(S{i}.M{j},S{i}.P{j},tpars.A,tpars.Q,B,u);
//                end
//            end

        }
//        % Remove the dead target
//        if tpars.ALLOW_MULTI_DEATH
//            if ~isempty(dead), S{i}.D = S{i}.targetIDs(dead,1); else, S{i}.D = []; end
//            for j=1:length(dead)
//                ind = find((1:S{i}.nTargets)~=dead(j)); % Remove index
//                ev_strs{i} = sprintf('Target %d died ',dead(j));
//                S{i}.nTargets = S{i}.nTargets-1;
//                S{i}.M = S{i}.M(ind);
//                S{i}.P = S{i}.P(ind);
//                S{i}.Tcount = S{i}.Tcount(ind);
//                S{i}.targetIDs = S{i}.targetIDs(ind,1);
//                S{i}.B = 0;
//            end
//        else
//            if ~isempty(dead)
//                ind = find((1:S{i}.nTargets)~=dead); % Remove index
//                ev_strs{i} = sprintf('Target %d died ',dead);
//                S{i}.nTargets = S{i}.nTargets-1;
//                S{i}.D = S{i}.targetIDs(dead,1);
//                S{i}.M = S{i}.M(ind);
//                S{i}.P = S{i}.P(ind);
//                S{i}.Tcount = S{i}.Tcount(ind);
//                S{i}.targetIDs = S{i}.targetIDs(ind,1);
//                S{i}.B = 0;
//            else
//                S{i}.D = 0;
//            end
//        end

    }
}


void tracker3d_particleCreate
(
     void** phPart,
     float W0,
     float M0[6],
     float P0[6][6],
     float dt
)
{
    *phPart = malloc1d(sizeof(MCS_data));
    MCS_data *p = (MCS_data*)(*phPart);

    p->W = W0;
    p->W_prev = W0;
    p->W0 = W0;
    p->nTargets = 0;
    p->nActive = 0;
    p->B = 0;
    p->D = 0;
    p->dt = dt;
    p->M = NULL;
    p->P = NULL;
    memcpy(p->M0.M, M0, 6*sizeof(float));
    memcpy(p->P0.P, P0, 6*6*sizeof(float));
    p->targetIDs = NULL;
    p->activeIDs = NULL;
    p->Tcount = NULL;
}

void tracker3d_particleDestroy
(
     void** phPart
)
{
    MCS_data *p = (MCS_data*)(*phPart); 

    if(p!=NULL){
        free(p->M);
        free(p->P);
        free(p->targetIDs);
        free(p->activeIDs);
        free(p->Tcount);
        p=NULL;
        *phPart = NULL;
    }
}

/* Copyright (C) 2002, 2003 Simo Särkkä (GPLv2) */
void lti_disc
(
    float* F,
    int len_N,
    int len_L,
    float* opt_L,
    float* opt_Qc,
    float dt,
    float* A,
    float* Q
)
{
    int i, j;
    float* L, *Qc;
    float* Fdt;
    float** L_Qc, **L_Qc_LT, **Phi, **ZE, **B, **AB, **AB1_T, **AB2_T, **Q_T;

    /* Defaults: */
    if(opt_L==NULL){ /* Identity */
        L = calloc1d(len_N*len_L, sizeof(float));
        for(i=0; i<MIN(len_N, len_L); i++)
            L[i*len_L+i] = 1.0f;
    }
    else
        L = opt_L;
    if(opt_Qc==NULL) /* zeros */
        Qc = calloc1d(len_L*len_L, sizeof(float));
    else
        Qc = opt_Qc;
    Fdt = malloc1d(len_N*len_N*sizeof(float));

    /* Closed form integration of transition matrix */
    for(i=0; i<len_N*len_N; i++)
        Fdt[i] = F[i]*dt;
    gexpm(Fdt, len_N, 0, A);

    /* Closed form integration of covariance by matrix fraction decomposition */
    L_Qc = (float**)malloc2d(len_N, len_L, sizeof(float));
    L_Qc_LT = (float**)malloc2d(len_N, len_N, sizeof(float));
    Phi = (float**)calloc2d(len_N*2, len_N*2, sizeof(float)); // ca
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, len_N, len_L, len_L, 1.0f,
                L, len_L,
                Qc, len_L, 0.0f,
                FLATTEN2D(L_Qc), len_L);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, len_N, len_N, len_L, 1.0f,
                FLATTEN2D(L_Qc), len_L,
                L, len_L, 0.0f,
                FLATTEN2D(L_Qc_LT), len_N);
    for(i=0; i<len_N; i++){
        for(j=0; j<len_N; j++){
            Phi[i][j] = F[i*len_N+j];
            Phi[i][j+len_N] = L_Qc_LT[i][j];
            Phi[i+len_N][j+len_N] = -F[j*len_N+i];
        }
    }
    utility_svsmul(FLATTEN2D(Phi), &dt, (len_N*2)*(len_N*2), NULL);
    ZE = (float**)calloc2d(len_N*2, len_N, sizeof(float));
    for(i=0; i<len_N; i++)
        ZE[i+len_N][i] = 1.0f;
    B = (float**)malloc2d(len_N*2, len_N*2, sizeof(float));
    AB = (float**)malloc2d(len_N*2, len_N, sizeof(float));
    gexpm(FLATTEN2D(Phi), len_N*2, 0, FLATTEN2D(B)); 
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, len_N*2, len_N, len_N*2, 1.0f,
                FLATTEN2D(B), len_N*2,
                FLATTEN2D(ZE), len_N, 0.0f,
                FLATTEN2D(AB), len_N);
    AB1_T = (float**)malloc2d(len_N, len_N, sizeof(float));
    AB2_T = (float**)malloc2d(len_N, len_N, sizeof(float));
    Q_T = (float**)malloc2d(len_N, len_N, sizeof(float));
    for(i=0; i<len_N; i++){
        for(j=0; j<len_N; j++){
            AB1_T[j][i] = AB[i][j];
            AB2_T[j][i] = AB[i+len_N][j];
        }
    }
    utility_sglslv(FLATTEN2D(AB2_T), len_N, FLATTEN2D(AB1_T), len_N, FLATTEN2D(Q_T));

    /* transpose back */
    for(i=0; i<len_N; i++)
        for(j=0; j<len_N; j++)
            Q[i*len_N+j] = Q_T[j][i];
 
    /* clean-up */
    if(opt_L==NULL)
        free(L);
    if(opt_Qc==NULL)
        free(Qc);
    free(Fdt);
    free(L_Qc);
    free(L_Qc_LT);
    free(Phi);
    free(ZE);
    free(B);
    free(AB);
    free(AB1_T);
    free(AB2_T);
    free(Q_T);
}
