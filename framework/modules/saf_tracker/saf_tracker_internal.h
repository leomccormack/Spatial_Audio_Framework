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

#ifndef __SAF_TRACKER_INTERNAL_H_INCLUDED__
#define __SAF_TRACKER_INTERNAL_H_INCLUDED__

#include "saf_tracker.h"
#include "../../modules/saf_utilities/saf_utilities.h"
#include "saf_externals.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef  SAF_ENABLE_TRACKER_MODULE

/** Spits out tracker status info to the terminal */
//#define TRACKER_VERBOSE
#ifdef TRACKER_VERBOSE
/** Spits out even more tracker information */
//# define TRACKER_VERY_VERBOSE
#endif

/** Maximum number of targets that can be tracked */
#define TRACKER3D_MAX_NUM_TARGETS ( 24 )
/** Maximum number of possible events during update */
#define TRACKER3D_MAX_NUM_EVENTS ( 24 )
/** Maximum number of particles */
#define TRACKER3D_MAX_NUM_PARTICLES ( 100 )

/* ========================================================================== */
/*                            Internal Structures                             */
/* ========================================================================== */

/** Void pointer (just to improve code readability when working with arrays of
 *  handles) */
typedef void* voidPtr;

/** Union struct for 3-D mean values */
typedef struct _M6 {
    union {
        struct { float m0, m1, m2, m3, m4, m5; };
        float M[6];
    };
} M6;

/** Union struct for 3-D variance values */
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
    float W_prev;    /**< Previous importance weight */
    float W0;        /**< PRIOR importance weight */ 
    int nTargets;    /**< Number of targets being tracked */
    float dt;        /**< Elapsed time between each observation/measurement */
    M6* M;           /**< Current target means; nTargets x ([6]) */
    P66* P;          /**< Current target variances; nTargets x ([6][6]) */
    int* targetIDs;  /**< Unique ID assigned to each target; nTargets x 1 */
    int* Tcount;     /**< Time elapsed since birth of target (Tcount * dt);
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
    void* hKF6;         /**< kf_update6 handle */
    voidPtr* SS;        /**< The particles; tpars.Np x 1 */
    voidPtr* SS_resamp; /**< Resampled particles; tpars.Np x 1 */
    float R[3][3];      /**< Diagonal matrix, measurement noise PRIORs along the
                         *   x,y,z axes */
    float A[6][6];      /**< Transition matrix */
    float Q[6][6];      /**< Discrete Process Covariance */
    float H[3][6];      /**< Measurement matrix */
    int incrementTime;  /**< Number steps of "tpars.dt" to increment time by */
    float W0;           /**< PRIOR importance weight */

    /* Events */
#ifdef TRACKER_VERBOSE
    char evt[TRACKER3D_MAX_NUM_EVENTS][256]; /**< Event descriptions */
#endif
    int evta[TRACKER3D_MAX_NUM_EVENTS];    /**< Event targets */
    float evp[TRACKER3D_MAX_NUM_EVENTS];   /**< Event priors */
    float evl[TRACKER3D_MAX_NUM_EVENTS];   /**< Event likelhoods*/
    float imp[TRACKER3D_MAX_NUM_EVENTS];   /**< Event distributions */
    voidPtr str[TRACKER3D_MAX_NUM_EVENTS]; /**< Structure after each event */

} tracker3d_data;
     

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Creates an instance of a particle / Monte-Carlo Sample
 *
 * @param[in] phPart (&) address of particle structure
 * @param[in] W0     Importance weight PRIOR
 * @param[in] dt     Time step
 */
void tracker3d_particleCreate(void** phPart,
                              float W0, 
                              float dt);

/**
 * Resets a particle structure to defaults
 *
 * @param[in] hPart Particle structure
 */
void tracker3d_particleReset(void* hPart);

/**
 * Copies particle structure "hPart1" into structure "hPart2"
 *
 * @param[in] hPart1 Particle structure 1
 * @param[in] hPart2 Particle structure 2
 */
void tracker3d_particleCopy(void* hPart1,
                            void* hPart2);

/**
 * Destroys an instance of a particle / Monte-Carlo Sample
 *
 * @param[in] phPart (&) address of particle structure
 */
void tracker3d_particleDestroy(void** phPart);

/**
 * Prediction step
 *
 * @param[in] hT3d tracker3d handle
 * @param[in] Tinc Number of time steps to increment by
 */
void tracker3d_predict(void* const hT3d, 
                       int Tinc);

/**
 * Prediction update
 *
 * @param[in] hT3d tracker3d handle
 * @param[in] Y    New observation/measurement; 3 x 1
 * @param[in] Tinc Number of time steps to increment by
 */
void tracker3d_update(void* const hT3d,
                      float* Y,
                      int Tinc);

/**
 * Returns the index of the most important particle
 *
 * @param[in] hT3d tracker3d handle
 * @returns index of the most important particle
 */
int tracker3d_getMaxParticleIdx(void* const hT3d);


/* ========================================================================== */
/*                              RBMCDA Functions                              */
/* ========================================================================== */

/**
 * Stratified resampling - returns a new set of indices according to the
 * probabilities P
 *
 * Sorted re-sampling is slower but has slightly smaller variance. Stratified
 * resampling is unbiased, almost as fast as deterministic resampling, and has
 * only slightly larger variance.
 *
 * In stratified resampling indices are sampled using random numbers [1]
 *    u_j~U[(j-1)/n,j/n],
 * where n is length of P. Compare this to simple random resampling where
 *    u_j~U[0,1].
 *
 * @warning This function assumes that the weights have been normalised!
 *
 * @param[in]  SS Array of particle structures; NP x 1
 * @param[in]  NP Number of particle structures
 * @param[out] s  Resampled indices; NP x 1
 *
 * @see [1] Kitagawa, G., Monte Carlo Filter and Smoother for Non-Gaussian
 *          Nonlinear State Space Models, Journal of Computational and Graphical
 *          Statistics, 5(1):1-25, 1996.
 *
 * Original Copyright (c) 2003-2004 Aki Vehtari (GPLv2)
 */
void resampstr(voidPtr* SS,
               int NP,
               int* s);

/**
 * Estimate the number of effective particles
 *
 * @warning This function assumes that the weights have been normalised!
 *
 * @param[in] SS Array of particle structures; NP x 1
 * @param[in] NP Number of particle structures
 * @returns Number of effective particles
 *
 * Original Copyright (C) 2003 Simo Särkkä, 2008 Jouni Hartikainen (GPLv2)
 */
float eff_particles(voidPtr* SS,
                    int NP);

/**
 * Normalises the weights of the given particles
 *
 * @param[in,out] SS Array of particle structures; NP x 1
 * @param[in]     NP Number of particle structures
 *
 * Original Copyright (C) 2008 Jouni Hartikainen (GPLv2)
 */
void normalise_weights(voidPtr* SS,
                       int NP);

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

/** Creates helper structure for kf_update6() */
void kf_update6_create(void ** const phUp6);

/** Destroys helper structure for kf_update6() */
void kf_update6_destroy(void ** const phUp6);

/**
 * Kalman Filter update step
 *
 * Kalman Filter model is:
 *    x[k] = A*x[k-1] + B*u[k-1] + q,  q ~ N(0,Q)
 *    y[k] = H*x[k]   + r,             r ~ N(0,R)
 *
 * Prediction step of Kalman filter computes predicted mean m-[k] and covariance
 * P-[k] of state:
 *    p(x[k] | y[1:k-1]) = N(x[k] | m-[k], P-[k])
 *
 * See for instance kf_predict6() how m-[k] and P-[k] are calculated.
 *
 * Update step computes the posterior mean m[k] and covariance P[k] of state
 * given new measurement:
 *    p(x[k] | y[1:k]) = N(x[k] | m[k], P[k])
 *
 * Innovation distribution is defined as:
 *    p(y[k] | y[1:k-1]) = N(y[k] | IM[k], IS[k])
 *
 * Updated mean x[k] and covarience P[k] are given by the following equations
 * (not the only possible ones):
 *    v[k] = y[k] - H[k]*m-[k]
 *    S[k] = H[k]*P-[k]*H[k]' + R[k]
 *    K[k] = P-[k]*H[k]'*[S[k]]^(-1)
 *    m[k] = m-[k] + K[k]*v[k]
 *    P[k] = P-[k] - K[k]*S[k]*K[k]'
 *
 * @note This has been hard-coded for N=6 and without 'K', 'IM' and 'IS', but a
 *       general version should be quite straight-forward (just not needed yet).
 *
 * @param[in]  hUp6  Handle for helper structure
 * @param[in]  X     Nx1 mean state estimate after prediction step
 * @param[in]  P     NxN state covariance after prediction step
 * @param[in]  y     Dx1 measurement vector.
 * @param[in]  H     Measurement matrix.
 * @param[in]  R     Measurement noise covariance.
 * @param[out] X_out Updated state mean
 * @param[out] P_out Updated state covariance
 * @param[out] LH    (&) Predictive probability (likelihood) of measurement.
 *
 * Original Copyright (C) 2002, 2003 Simo Särkkä, 2007 Jouni Hartikainen (GPLv2)
 */
void kf_update6(void * const hUp6,
                float X[6],
                float P[6][6],
                float y[3],
                float H[3][6],
                float R[3][3],
                float X_out[6],
                float P_out[6][6],
                float* LH);

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
 * @param[in] hUp6  Handle for helper structure
 * @param[in] X     Values
 * @param[in] M     Mean of distibution or values as 3x3 matrix.
 * @param[in] S     3x3 covariance matrix
 * @returns Probability of X.
 *
 * Original Copyright (C) 2002 Simo Särkkä (GPLv2)
 */
float gauss_pdf3(/* Input Arguments */
                 void * const hUp6,
                 float X[3],
                 float M[3],
                 float S[3][3]);

/**
 * Draws samples from a given one dimensional discrete distribution
 *
 * @param[in] P     Discrete distribution; len_P x 1
 * @param[in] len_P length of P
 *
 * Original Copyright (C) 2002 Simo Särkkä, 2008 Jouni Hartikainen (GPLv2)
 */
int categ_rnd(float* P,
              int len_P);

#endif /* SAF_ENABLE_TRACKER_MODULE */


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_TRACKER_INTERNAL_H_INCLUDED__ */
