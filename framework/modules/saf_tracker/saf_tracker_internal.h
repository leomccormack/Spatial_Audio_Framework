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

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */ 
    
/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */



/* ========================================================================== */
/*                            Internal Structures                             */
/* ========================================================================== */

/**
 * Main structure for tracker3d
 */
typedef struct _tracker3dlib
{
    tracker3d_config tpars;
    
} tracker3dlib_data;
     

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */



#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_TRACKER_INTERNAL_H_INCLUDED__ */
