/*
 * Copyright 2017-2018 Leo McCormack
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

/**
 * @file rotator_internal.h
 * @brief A basic spherical harmonic/ Ambisonic signals rotator, based on the
 *        recursive approach detailed in [1]
 *
 * @test test__saf_example_rotator()
 *
 * @see [1] Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real
 *          Spherical Harmonics. Direct Determination by Recursion Page:
 *          Additions and Corrections. Journal of Physical Chemistry A, 102(45),
 *          9099?9100.
 *
 * @author Leo McCormack
 * @date 02.11.2017
 * @license ISC
 */

#ifndef __ROTATOR_INTERNAL_H_INCLUDED__
#define __ROTATOR_INTERNAL_H_INCLUDED__

#include "rotator.h"       /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#if !defined(ROTATOR_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define ROTATOR_FRAME_SIZE ( FRAME_SIZE ) /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define ROTATOR_FRAME_SIZE ( 64 )         /**< Framesize, in time-domain samples */
# endif
#endif

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/** Available Ambisonic channel ordering conventions */
typedef enum {
    M_ROT_READY = 1,           /**< M_rot is ready */
    M_ROT_RECOMPUTE_EULER,     /**< Use Euler angles to recompute M_rot */
    M_ROT_RECOMPUTE_QUATERNION /**< Use Quaternions to recompute M_rot */
} M_ROT_STATUS;

/** Main struct for the rotator */
typedef struct _rotator
{
    /* Internal buffers */
    float inputFrameTD[MAX_NUM_SH_SIGNALS][ROTATOR_FRAME_SIZE];         /**< Input frame of signals */
    float prev_inputFrameTD[MAX_NUM_SH_SIGNALS][ROTATOR_FRAME_SIZE];    /**< Previous frame of signals */
    float tempFrame[MAX_NUM_SH_SIGNALS][ROTATOR_FRAME_SIZE];            /**< Temporary frame */
    float tempFrame_fadeOut[MAX_NUM_SH_SIGNALS][ROTATOR_FRAME_SIZE];    /**< Temporary frame with linear interpolation (fade-out) applied */
    float outputFrameTD[MAX_NUM_SH_SIGNALS][ROTATOR_FRAME_SIZE];        /**< Output frame of SH signals */
    float outputFrameTD_fadeIn[MAX_NUM_SH_SIGNALS][ROTATOR_FRAME_SIZE]; /**< Output frame of SH signals with linear interpolation (fade-in) applied */

    /* Internal variables */
    float interpolator_fadeIn[ROTATOR_FRAME_SIZE];       /**< Linear Interpolator (fade-in) */
    float interpolator_fadeOut[ROTATOR_FRAME_SIZE];      /**< Linear Interpolator (fade-out) */
    float M_rot[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];      /**< Current SH rotation matrix [1] */
    float prev_M_rot[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS]; /**< Previous SH rotation matrix [1] */
    M_ROT_STATUS M_rot_status;      /**< see #M_ROT_STATUS */
    int fs;                         /**< Host sampling rate, in Hz */

    /* user parameters */
    quaternion_data Q;              /**< Quaternion used for rotation */
    int bFlipQuaternion;            /**< 1: invert quaternion, 0: no inversion */
    float yaw;                      /**< yaw (Euler) rotation angle, in degrees */
    float roll;                     /**< roll (Euler) rotation angle, in degrees */
    float pitch;                    /**< pitch (Euler) rotation angle, in degrees */
    int bFlipYaw;                   /**< flag to flip the sign of the yaw rotation angle */
    int bFlipPitch;                 /**< flag to flip the sign of the pitch rotation angle */
    int bFlipRoll;                  /**< flag to flip the sign of the roll rotation angle */
    int useRollPitchYawFlag;        /**< rotation order flag, 1: r-p-y, 0: y-p-r */
    CH_ORDER chOrdering;            /**< Ambisonic channel order convention (see #CH_ORDER) */
    NORM_TYPES norm;                /**< Ambisonic normalisation convention (see #NORM_TYPES) */
    SH_ORDERS inputOrder;           /**< current input/output SH order */ 
    
} rotator_data;
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __ROTATOR_INTERNAL_H_INCLUDED__ */
