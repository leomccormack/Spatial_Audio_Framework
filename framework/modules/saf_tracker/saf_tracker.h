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
 * @file saf_tracker.h
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

#ifndef __SAF_TRACKER_H_INCLUDED__
#define __SAF_TRACKER_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**< User parameters struct */
typedef struct _tracker3d_config{
    int Np;                   /**< Number of Monte Carlo samples/particles. */
    int maxNactiveTargets;    /**< Max numer of simultaneous targets */
    float Vazi_deg;           /**< Velocity of target on azimuthal plane */
    float Vele_deg;           /**< Velocity of target on median plane */
    float noiseLikelihood;    /**< Likelihood of an estimate being noise/clutter
                               *   between [0..1] */
    float measNoiseSD_deg;    /**< Measurement noise standard deviation
                               *   estimates within the range +/-20 degrees
                               *   belong to the same target */
    float noiseSpecDen_deg;   /**< Noise spectral density; influences the
                               *   smoothness of the traget tracks */
    int ALLOW_MULTI_DEATH;    /**< FLAG whether to allow for multiple target
                               *   deaths in the same tracker prediction step */
    float init_birth;         /**< Prior probability of birth [0 1] */
    float alpha_death;        /**< Prior probability of death; always >= 1 */
    float beta_death;         /**< Prior probability of death; always >= 1 */
    float dt;                 /**< Hop length of frames in seconds */
    int MULTI_ACTIVE;         /**< FLAG whether or not to allow multiple active
                               *   sources for each update */
    float W_avg_coeff;        /**< Real-time tracking is based on the particle
                               *   with highest weight. A one-pole averaging
                               *   filter is used to smooth these weights over
                               *   time */
    int FORCE_KILL_TARGETS;   /**< FLAG force kill targets that are close to
                               *   another target. In these cases, the target
                               *   that has been 'alive' for the least amount
                               *   of time, is killed */
    float forceKillAngle_deg; /**< Angle at which to kill targets, in degrees */
    
}tracker3d_config;


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the mighty tracker3dlib
 *
 * @param[in] phT3d (&) address of tracker3dlib handle
 */
void tracker3dlib_create(void** const phT3d,
                         tracker3d_config tpars);

/**
 * Destroys an instance of the mighty tracker3dlib
 *
 * @param[in] phT3d (&) address of tracker3dlib handle
 */
void tracker3dlib_destroy(void** const phT3d);

/**
 * Initialises an instance of tracker3dlib with default settings
 *
 * @param[in] hT3d       tracker3dlib handle
 * @param[in] samplerate Host samplerate.
 */
void tracker3d_init(void* const hT3d,
                    int samplerate);
    
/** 
 */
void tracker3dlib_predict(void* const hT3d,
                          float** const inputs,
                          float** const outputs,
                          int nInputs,
                          int nOutputs,
                          int nSamples);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_TRACKER_H_INCLUDED__ */
