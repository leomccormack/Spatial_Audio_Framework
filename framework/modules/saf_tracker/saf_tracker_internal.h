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
 * @date 12.08.2020
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


#define TRACKER_VERBOSE

/* ========================================================================== */
/*                            Internal Structures                             */
/* ========================================================================== */

/** Void pointer (just to improve code readability when working with arrays of
 * handles) */
typedef void* voidPtr;

/** Mean union struct for 3-D */
typedef struct _M6 {
    union {
        struct { float m0, m1, m2, m3, m4, m5; };
        float M[6];
    };
} M6;

/** Variance union struct for 3-D */
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

/** Monte-Carlo Sample (particle) structure */
typedef struct _MCS {
    float W;         /**< Importance weight */
    float W_prev;    /**< Previous important weight */
    float W0;        /**< PRIOR importance weight */
    M6 M0;           /**< 0,1,2: Position of sound source PRIORs (x,y,z), 3,4,5:
                      *   Mean velocity PRIORs (x,y,z) */
    P66 P0;          /**< Diagonal matrix, 0,1,2: Variance PRIORs of estimates
                      *   along the x,y,z axes; 3,4,5 Velocity PRIORs of
                      *   estimates along the x,y,z axes */
    int nTargets;    /**< Number of targets being tracked */
    int nActive;     /**< Number of active targets */
    float dt;        /**< Elapsed time inbetween each observation/measurment */
    M6* M;           /**< Current target means; nTargets x ([6]) */
    P66* P;          /**< Current target variances; nTargets x ([6][6]) */
    int* targetIDs;  /**< Unique ID assigned to each target; nTargets x 1 */
    int* activeIDs;  /**< IDs of targets currently active; nActive x 1 */
    int* Tcount;     /**< Time elapsed since birth of target: Tcount * dt;
                      *   nTargets x 1 */
#ifdef TRACKER_VERBOSE
    char evstr[256]; /** Event string */
#endif

} MCS_data;

/** Main structure for tracker3d */
typedef struct _tracker3d
{
    /** User parameters struct */
    tracker3d_config tpars;

    /* Internal */
    voidPtr* SS;    /**< The particles; tpars.Np x 1 */
    float R[3][3];  /**< Diagonal matrix, measurement noise PRIORs along the
                     *   x,y,z axes */
    float A[6][6];  /**< Transition matrix */
    float Q[6][6];  /**< Discrete Process Covariance */
    float H[3][6];  /**< Measurement matrix */
    int incrementTime;

    
} tracker3d_data;
     

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 *
 */
void tracker3d_particleCreate(void** phPart,
                              float W0,
                              float M0[6],
                              float P0[6][6],
                              float dt);

/**
 * Prediction step
 */
void tracker3d_predict(void* const hT3d, 
                       int Tinc);


/* ========================================================================== */
/*                              RBMCDA Functions                              */
/* ========================================================================== */

/**
 * Perform Kalman Filter prediction step
 *
 * The model is:
 *    x[k] = A*x[k-1] + B*u[k-1] + q,  q ~ N(0,Q).
 * The predicted state is distributed as follows:
 *    p(x[k] | x[k-1]) = N(x[k] | A*x[k-1], Q[k-1])
 *
 * The predicted mean x-[k] and covariance P-[k] are calculated with the
 * following equations:
 *    m-[k] = A*x[k-1]
 *    P-[k] = A*P[k-1]*A' + Q.
 *
 * @note This has been hard-coded for N=6 and without 'B' and 'u', but a general
 *       version should be quite straight-forward (just not needed yet)...
 *
 * @param[in,out] M Mean state estimate of previous step
 * @param[in,out] P State covariance of previous step
 * @param[in]     A Transition matrix of discrete model
 * @param[in]     Q Process noise of discrete model
 *
 * Original Copyright (C) 2002-2006 Simo Särkkä, 2007 Jouni Hartikainen (GPLv2)
 */
void kf_predict6(float M[6],
                 float P[6][6],
                 float A[6][6],
                 float Q[6][6]);

/**
 * Cumulative density function of a Gamma distribution
 *
 * @param[in] x    Locations where to evaluate the CDF
 * @param[in] gam  Parameter of the distribution
 * @param[in] beta Parameter of the distribution
 * @param[in] mu   Mean of the distribution
 * @returns Cumulative density at the given locations
 *
 * Original Copyright (C) 2003 Simo Särkkä, 2008 Jouni Hartikainen (GPLv2)
 */
float gamma_cdf(float x,
                float gam,
                float beta,
                float mu);

/**
 * LTI_DISC  Discretize LTI ODE with Gaussian Noise
 *
 * Discretize LTI ODE with Gaussian Noise. The original ODE model is in form:
 *    dx/dt = F x + L w,  w ~ N(0,Qc)
 *
 * Result of discretization is the model:
 *    x[k] = A x[k-1] + q, q ~ N(0,Q)
 *
 * Which can be used for integrating the model exactly over time steps, which
 * are multiples of dt.
 *
 * @param[in]  F      Square feedback matrix; FLAT: len_N x len_N
 * @param[in]  len_N  Size of square matrix 'F'
 * @param[in]  len_Q  Size of square matrix 'opt_Qc'
 * @param[in]  opt_L  Noise effect matrix (optional, set to NULL for default
 *                    values); FLAT: len_N x len_Q
 * @param[in]  opt_Qc Diagonal Spectral Density (optional, set to NULL for
 *                    default values); FLAT: len_Q x len_Q
 * @param[in]  dt     Time Step
 * @param[out] A      Transition matrix; FLAT: len_N x len_N
 * @param[out] Q      Discrete Process Covariance; FLAT: len_N x len_N
 *
 * Original Copyright (C) 2002, 2003 Simo Särkkä (GPLv2)
 */
void lti_disc(/* Input Arguments */
              float* F,
              int len_N,
              int len_Q,
              float* opt_L,
              float* opt_Qc,
              float dt,
              /* Output Arguments */
              float* A,
              float* Q);


/**
 * Multivariate Gaussian PDF
 *
 * Calculate values of PDF (Probability Density Function) of multivariate
 * Gaussian distribution:
 *    N(X | M, S)
 *
 * @note This has been hard-coded for N=3, but a general version should be quite
 *       straight-forward (just not needed yet)...
 *
 * @param[in] X Values
 * @param[in] M Mean of distibution or values as 3x3 matrix.
 * @param[in] S 3x3 covariance matrix
 * @returns Probability of X.
 *
 * Copyright (C) 2002 Simo Särkkä (GPLv2)
 */
float gauss_pdf3(float X[3], float M[3], float S[3][3])


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_TRACKER_INTERNAL_H_INCLUDED__ */
