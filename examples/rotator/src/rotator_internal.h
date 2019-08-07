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

/*
 * Filename: rotator_internal.h
 * ----------------------------
 * A simple spherical harmonic domain rotator, based on the recursive approach
 * detailed in [1].
 *
 * Dependencies:
 *     saf_utilities, saf_sh
 * Author, date created:
 *     Leo McCormack, 02.11.2017
 *
 * [1] Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical
 *     Harmonics. Direct Determination by Recursion Page: Additions and
 *     Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100.
 */

#ifndef __ROTATOR_INTERNAL_H_INCLUDED__
#define __ROTATOR_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "rotator.h"
#define SAF_ENABLE_SH /* for spherical harmonic domain rotation matrices */
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define MAX_SH_ORDER ( ROTATOR_MAX_SH_ORDER )
#define MAX_NUM_SH_SIGNALS ( (MAX_SH_ORDER + 1)*(MAX_SH_ORDER + 1)  )    /* (L+1)^2 */
#ifndef DEG2RAD
  #define DEG2RAD(x) (x * PI / 180.0f)
#endif
#ifndef RAD2DEG
  #define RAD2DEG(x) (x * 180.0f / PI)
#endif
    
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

typedef struct _rotator
{
    float inputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float prev_inputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float tempFrame[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float outputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    
    /* internal */
    float interpolator[FRAME_SIZE];
    float M_rot[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];
    float prev_M_rot[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];
    int recalc_M_rotFLAG;

    /* user parameters */
    float yaw, roll, pitch;      /* rotation angles in degrees */
    int bFlipYaw, bFlipPitch, bFlipRoll;  /* flag to flip the sign of the individual rotation angles */
    CH_ORDER chOrdering;         /* only ACN is supported */
    NORM_TYPES norm;             /* N3D or SN3D */
    INPUT_ORDERS inputOrder;     /* current input/output order int order;*/
    int useRollPitchYawFlag;     /* rotation order flag, 1: r-p-y, 0: y-p-r */
    
} rotator_data;
    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __ROTATOR_INTERNAL_H_INCLUDED__ */
