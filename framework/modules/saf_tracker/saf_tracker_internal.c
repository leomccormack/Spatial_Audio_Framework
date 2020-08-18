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
    int i;
    float* L, *Qc;

    /* Defaults: */
    if(opt_L==NULL){ /* Identity */
        L = calloc1d(len_N*len_N, sizeof(float));
        for(i=0; i<len_N; i++)
            L[i] = 1.0f;
    }
    else
        L = opt_L;
    if(opt_Qc==NULL) /* zeros */
        Qc = calloc1d(len_N*len_N, sizeof(float));
    else
        Qc = opt_Qc;

    

//  %
//  % Closed form integration of transition matrix
//  %
//  A = expm(F*dt);
//
//  %
//  % Closed form integration of covariance
//  % by matrix fraction decomposition
//  %
//  n   = size(F,1);
//  Phi = [F L*Qc*L'; zeros(n,n) -F'];
//  AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];
//  Q   = AB(1:n,:)/AB((n+1):(2*n),:);

}
