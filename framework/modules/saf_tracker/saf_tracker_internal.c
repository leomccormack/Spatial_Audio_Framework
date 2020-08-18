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
 

//%LTI_DISC  Discretize LTI ODE with Gaussian Noise
//%
//% Syntax:
//%   [A,Q] = lti_disc(F,L,Qc,dt)
//%
//% In:
//%   F  - NxN Feedback matrix
//%   L  - NxL Noise effect matrix        (optional, default identity)
//%   Qc - LxL Diagonal Spectral Density  (optional, default zeros)
//%   dt - Time Step                      (optional, default 1)
//%
//% Out:
//%   A - Transition matrix
//%   Q - Discrete Process Covariance
//%
//% Description:
//%   Discretize LTI ODE with Gaussian Noise. The original
//%   ODE model is in form
//%
//%     dx/dt = F x + L w,  w ~ N(0,Qc)
//%
//%   Result of discretization is the model
//%
//%     x[k] = A x[k-1] + q, q ~ N(0,Q)
//%
//%   Which can be used for integrating the model
//%   exactly over time steps, which are multiples
//%   of dt.

/* Copyright (C) 2002, 2003 Simo Särkkä (GPLv2) */
void lti_disc
(
    float* F,      // NxN
    int len_N,
    int len_L,
    float* opt_L,  // NxL
    float* opt_Qc, // LxL
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
