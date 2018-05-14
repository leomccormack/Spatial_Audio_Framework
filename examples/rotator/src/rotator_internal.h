/*
 Copyright 2017-2018 Leo McCormack
 
 Permission to use, copy, modify, and/or distribute this software for any purpose with or
 without fee is hereby granted, provided that the above copyright notice and this permission
 notice appear in all copies.
 
 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 OR PERFORMANCE OF THIS SOFTWARE.
*/
/*
 * Filename:
 *     rotator_internal.h  
 * Description:
 *     A simple spherical harmonic domain rotator.
 * Dependencies:
 *     saf_utilities, saf_sh
 * Author, date created:
 *     Leo McCormack, 02.11.2017
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
#endif

#define NUM_SH_SIGNALS ( (SH_ORDER + 1)*(SH_ORDER + 1)  )    /* (L+1)^2 */
    
#ifndef DEG2RAD
  #define DEG2RAD(x) (x * PI / 180.0f)
#endif
#ifndef RAD2DEG
  #define RAD2DEG(x) (x * 180.0f / PI)
#endif
    
typedef enum _CH_ORDER{
    CH_ACN = 1
}CH_ORDER;

typedef enum _NORM_TYPES{
    NORM_N3D = 1,
    NORM_SN3D
}NORM_TYPES;

typedef struct _rotator
{
    float inputFrameTD[NUM_SH_SIGNALS][FRAME_SIZE];
    float outputFrameTD[NUM_SH_SIGNALS][FRAME_SIZE];

    /* user parameters */
    float yaw, roll, pitch;
    int bFlipYaw, bFlipPitch, bFlipRoll;
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    
} rotator_data;
    
#ifdef __cplusplus
}
#endif


#endif /* __ROTATOR_INTERNAL_H_INCLUDED__ */




















