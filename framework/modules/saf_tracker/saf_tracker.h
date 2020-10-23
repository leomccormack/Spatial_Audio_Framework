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
 *@addtogroup Tracker
 *@{
 * @file saf_tracker.h
 * @brief Particle filtering based tracker (#SAF_TRACKER_MODULE)
 *
 * Based on the RBMCDA [1] Matlab toolbox (GPLv2 license) by Simo Sa"rkka" and
 * Jouni Hartikainen (Copyright (C) 2003-2008):
 *     https://users.aalto.fi/~ssarkka/#softaudio
 *
 * And also inspired by the work of Sharath Adavanne, Archontis Politis, Joonas
 * Nikunen, and Tuomas Virtanen (GPLv2 license):
 *     https://github.com/sharathadavanne/multiple-target-tracking
 *
 * @see [1] Sa"rkka", S., Vehtari, A. and Lampinen, J., 2004, June. Rao-
 *          Blackwellized Monte Carlo data association for multiple target
 *          tracking. In Proceedings of the seventh international conference on
 *          information fusion (Vol. 1, pp. 583-590). I.
 *
 * @author Leo McCormack
 * @date 12.08.2020
 */

#ifndef __SAF_TRACKER_H_INCLUDED__
#define __SAF_TRACKER_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef  SAF_ENABLE_TRACKER_MODULE

/* ========================================================================== */
/*                               Public Structs                               */
/* ========================================================================== */

/** User parameters for tracker3d */
typedef struct _tracker3d_config{
    int Np;                  /**< Number of Monte Carlo samples/particles. The
                              *   more complex the distribution is, the more
                              *   particles required (but also the more
                              *   computationally expensive the tracker
                              *   becomes). */
    int ARE_UNIT_VECTORS;    /**< 1: if the Cartesian coordinates are given
                              *   as unit vectors, 0: if not */
    int maxNactiveTargets;   /**< Maximum number of simultaneous targets
                              *   permitted */
    float noiseLikelihood;   /**< Likelihood of an estimate being noise/clutter
                              *   between [0..1] */
    float measNoiseSD;       /**< Measurement noise standard deviation.
                              *   Estimates within this standard deviation
                              *   range belong to the same target */
    float noiseSpecDen;      /**< Noise spectral density; influences the
                              *   smoothness of the traget tracks */
    int ALLOW_MULTI_DEATH;   /**< FLAG whether to allow for multiple target
                              *   deaths in the same tracker prediction step.
                              *   1: enabled, 0: disabled */
    float init_birth;        /**< PRIOR probability of birth [0 1] */
    float alpha_death;       /**< Coefficient influencing the likelihood that
                              *   a target will die; always >= 1 */
    float beta_death;        /**< Coefficient influencing the likelihood that
                              *   a target will die; always >= 1 */
    float dt;                /**< Elapsed time (in seconds) between
                              *   observations/measurements */
    float W_avg_coeff;       /**< Real-time tracking is based on the particle
                              *   with highest weight. A one-pole averaging
                              *   filter is used to smooth these weights over
                              *   time [0..0.999] */
    int FORCE_KILL_TARGETS;  /**< FLAG force kill targets that are too close to
                              *   one another. In these cases, the target which
                              *   has been 'alive' for the least amount of
                              *   of time is killed; 1: enabled, 0: disabled */
    float forceKillDistance; /**< Euclidian distance at which to kill targets
                              *   that come too close to other (older) targets
                              *   (<=). */
    float M0[6];             /**< [0,1,2] Position of sound source PRIORs
                              *   (x,y,z), [3,4,5] Mean velocity PRIORs (x,y,z)
                              *   Note: If there is no reasonable prior
                              *   position, then you can set higher variances.
                              */
    float P0[6][6];          /**< Diagonal matrix, [0,1,2] Variance PRIORs of
                              *   estimates along the x,y,z axes; [3,4,5]
                              *   Velocity PRIORs of estimates along the x,y,z
                              *   axes */
    float cd;                /**< PRIOR probability of noise. */
    
}tracker3d_config;


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the mighty tracker3d
 *
 * @test test__tracker3d()
 *
 * @param[in] phT3d (&) address of tracker3d handle
 * @param[in] tpars Tracker configuration/user parameter struct
 */
void tracker3d_create(void** const phT3d,
                      tracker3d_config tpars);

/**
 * Destroys an instance of the mighty tracker3d
 *
 * @param[in] phT3d (&) address of tracker3d handle
 */
void tracker3d_destroy(void** const phT3d);

/**
 * Resets an instance of the mighty tracker3d
 *
 * @param[in] hT3d tracker3d handle
 */
void tracker3d_reset(void* const hT3d);
    
/**
 * Tracker time step to update & predict current target locations and to parse
 * new measurements/observations
 *
 * @note It is assumed that this will be called every step in time (tpars.dt).
 *       If there are no new observations/measurements then still call this
 *       function, but set newObs_xyz=NULL and/or, nObs=0.
 *
 * @param[in]  hT3d       tracker3d handle
 * @param[in]  newObs_xyz New observations/measurements; nObs x 3
 * @param[in]  nObs       Number of new observations/measurements
 * @param[out] target_xyz (&) Current target locations; *nTargets x 3
 * @param[out] target_IDs (&) (Unique) target IDs; *nTargets x 1
 * @param[out] nTargets   (&) Current number of targets being tracked
 */
void tracker3d_step(void* const hT3d,
                    float* newObs_xyz,
                    int nObs,
                    float** target_xyz,
                    int** target_IDs,
                    int* nTargets);

#endif /* SAF_ENABLE_TRACKER_MODULE */


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_TRACKER_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup Tracker */
