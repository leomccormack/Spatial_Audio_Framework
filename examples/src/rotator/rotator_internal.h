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
 */

#ifndef __ROTATOR_INTERNAL_H_INCLUDED__
#define __ROTATOR_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "rotator.h"
#include "saf.h"
#include "saf_externals.h" /* to also include saf dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#ifndef FRAME_SIZE
# define FRAME_SIZE ( 64 ) 
#endif
#ifndef DEG2RAD
# define DEG2RAD(x) (x * SAF_PI / 180.0f)
#endif
#ifndef RAD2DEG
# define RAD2DEG(x) (x * 180.0f / SAF_PI)
#endif

/** Available Ambisonic channel ordering conventions */
typedef enum {
    M_ROT_READY = 1,           /**< M_rot is ready */
    M_ROT_RECOMPUTE_EULER,     /**< Use Euler angles to recompute M_rot */
    M_ROT_RECOMPUTE_QUATERNION /**< Use Quaternions to recompute M_rot */

} M_ROT_STATUS;
    
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main struct for the rotator
 */
typedef struct _rotator
{
    /* internal */
    float inputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float prev_inputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float tempFrame[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float outputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float interpolator[FRAME_SIZE];
    float M_rot[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];
    float prev_M_rot[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];
    M_ROT_STATUS M_rot_status;
    int fs;

    /* user parameters */
    quaternion_data Q;                    /**< Quaternion used for rotation */
    int bFlipQuaternion;                  /**< 1: invert quaternion, 0: no inversion */
    float yaw, roll, pitch;               /**< rotation angles in degrees */
    int bFlipYaw, bFlipPitch, bFlipRoll;  /**< flags to flip the sign of the Euler rotation angles */
    CH_ORDER chOrdering;                  /**< only ACN is supported */
    NORM_TYPES norm;                      /**< N3D or SN3D */
    SH_ORDERS inputOrder;                 /**< current input/output order int order;*/
    int useRollPitchYawFlag;              /**< rotation order flag, 1: r-p-y, 0: y-p-r */
    
} rotator_data;
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __ROTATOR_INTERNAL_H_INCLUDED__ */
