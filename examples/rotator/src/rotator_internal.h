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
 * @brief  A simple spherical harmonic domain rotator, based on the recursive
 *         approach detailed in [1].
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

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define FRAME_SIZE ( 64 )
#define MAX_SH_ORDER ( ROTATOR_MAX_SH_ORDER )
#define MAX_NUM_SH_SIGNALS ( (MAX_SH_ORDER + 1)*(MAX_SH_ORDER + 1)  )    /* (L+1)^2 */
#ifndef DEG2RAD
  #define DEG2RAD(x) (x * SAF_PI / 180.0f)
#endif
#ifndef RAD2DEG
  #define RAD2DEG(x) (x * 180.0f / SAF_PI)
#endif
    
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main struct for the rotator
 */
typedef struct _rotator
{
    /* FIFO buffers */
    int FIFO_idx;
    float inFIFO[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float outFIFO[MAX_NUM_SH_SIGNALS][FRAME_SIZE];

    /* internal */
    float inputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float prev_inputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float tempFrame[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float outputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float interpolator[FRAME_SIZE];
    float M_rot[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];
    float prev_M_rot[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];
    int recalc_M_rotFLAG;

    /* user parameters */
    float yaw, roll, pitch;               /**< rotation angles in degrees */
    int bFlipYaw, bFlipPitch, bFlipRoll;  /**< flag to flip the sign of the individual rotation angles */
    ROTATOR_CH_ORDER chOrdering;          /**< only ACN is supported */
    ROTATOR_NORM_TYPES norm;              /**< N3D or SN3D */
    ROTATOR_INPUT_ORDERS inputOrder;      /**< current input/output order int order;*/
    int useRollPitchYawFlag;              /**< rotation order flag, 1: r-p-y, 0: y-p-r */
    
} rotator_data;
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __ROTATOR_INTERNAL_H_INCLUDED__ */
