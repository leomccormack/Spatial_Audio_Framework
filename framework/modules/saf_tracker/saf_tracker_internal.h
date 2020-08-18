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
 * @file saf_tracker_internal.h
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

#ifndef __SAF_TRACKER_INTERNAL_H_INCLUDED__
#define __SAF_TRACKER_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "saf_tracker.h"
#include "saf.h"
#include "saf_externals.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */ 
    
/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */



/* ========================================================================== */
/*                            Internal Structures                             */
/* ========================================================================== */

/** Void pointer (just to improve code readability when working with arrays of
 * handles) */
typedef void* voidPtr;

typedef struct _M6 {
    union {
        struct { float m0, m1, m2, m3, m4, m5; };
        float M[6];
    };
} M6;

typedef struct _P66 {
    union {
        struct {
            float p00, p01, p02, p03, p04, p05,
                  p10, p11, p12, p13, p14, p15,
                  p20, p21, p22, p23, p24, p25,
                  p30, p31, p32, p33, p34, p35,
                  p40, p41, p42, p43, p44, p45,
                  p50, p51, p52, p53, p54, p55;
        };
        float P[6][6];
    };
} P66;

/** Monte Carlo Sample (particle) structure */
typedef struct _MCS {
    float W;        /**< Importance weight */
    float W_prev;   /**< Previous important weight */
    float W0;       /**< PRIOR importance weight */
    M6 M0;          /**< 0,1,2: Position of sound source PRIORs (x,y,z), 3,4,5:
                     *   Mean velocity PRIORs (x,y,z) */
    P66 P0;         /**< Diagonal matrix, 0,1,2: Variance PRIORs of estimates
                     *   along the x,y,z axes; 3,4,5 Velocity PRIORs of
                     *   estimates along the x,y,z axes */
    int nTargets;   /**< Number of targets being tracked */
    int nActive;    /**< Number of active targets */
    int B;          /**< Birth indicator (non-zero indicates new birth with
                     *   index of the value of 'B') */
    int D;          /**< Death indicator (non-zero indicates new death with
                     *   index of the value of 'D') */
    float dt;       /**< Elapsed time inbetween each observation/measurment */
    M6* M;          /**< Current target means; nTargets x ([6]) */
    P66* P;         /**< Current target variances; nTargets x ([6][6]) */
    int* targetIDs; /**< Unique ID assigned to each target; nTargets x 1 */
    int* activeIDs; /**< IDs of targets currently active; nActive x 1 */
    int* Tcount;    /**< Time elapsed since birth of target: Tcount * dt;
                     *   nTargets x 1 */

} MCS_data;

/** Main structure for tracker3d */
typedef struct _tracker3d
{
    /** User parameters struct */
    tracker3d_config tpars;

    /* Internal */
    voidPtr* S;     /**< The particles; tpars.Np x 1 */
    float R[3][3];  /**< Diagonal matrix, measurement noise PRIORs along the
                     *   x,y,z axes */
    float A[6][6];  /**< Transition matrix */
    float Q[6][6];  /**< Discrete Process Covariance */
    float H[3][6];  /**< Measurement matrix */

    
} tracker3d_data;
     

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */


void tracker3d_particleCreate(void** phPart,
                              float W0,
                              float M0[6],
                              float P0[6][6],
                              float dt);

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
 );


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_TRACKER_INTERNAL_H_INCLUDED__ */
