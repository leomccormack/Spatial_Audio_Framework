/*
 * Copyright 2020 Leo McCormack
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
 * @file ambi_roomsim_internal.h
 * @brief A simple shoebox room Ambisonic encoder.
 *
 * @author Leo McCormack
 * @date 10.08.2020
 */

#ifndef __AMBI_ROOMSIM_INTERNAL_H_INCLUDED__
#define __AMBI_ROOMSIM_INTERNAL_H_INCLUDED__

#include "ambi_roomsim.h"  /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#ifndef FRAME_SIZE
# define FRAME_SIZE ( 128 )
#endif

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for ambi_roomsim. Contains variables for audio buffers, internal
 * variables, user parameters
 */
typedef struct _ambi_roomsim
{
    /* Internals */
    float inputFrameTD[MAX_NUM_INPUTS][FRAME_SIZE];
    float prev_inputFrameTD[MAX_NUM_INPUTS][FRAME_SIZE];
    float tempFrame[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float outputFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float fs;

    /* Internal */
    void* hIms;               /**< Image source implementation handle */
    int sourceIDs[ROOM_SIM_MAX_NUM_SOURCES];     /**< Unique IDs per source in the simulation */
    int receiverIDs[ROOM_SIM_MAX_NUM_RECEIVERS]; /**< Unique IDs per receiver in the simulation */
    float** src_sigs;         /**< Source signal buffers; ROOM_SIM_MAX_NUM_SOURCES x FRAME_SIZE */
    float*** rec_sh_outsigs;  /**< Receiver signal buffers; ROOM_SIM_MAX_NUM_RECEIVERS x MAX_NUM_SH_SIGNALS x FRAME_SIZE */
    int reinit_room;          /**< Flag, 1: re-init required, 0: not required*/
    int new_sh_order;         /**< New receiver SH order */
    int new_nSources;         /**< New number of sources */
    int new_nReceivers;       /**< New number of receivers */
    
    /* user parameters */
    int sh_order;             /**< Current SH order of receivers */
    int enableReflections;    /**< 0: disabled, 1: enabled */
    int refl_order;           /**< Current maximum image source reflection order */
    int nSources;             /**< Current number of sources */
    int nReceivers;           /**< Current number of receivers */
    float room_dims[3];       /**< Room dimensions along the x,y,z axes in meters */
    float abs_wall[6];        /**< Absorption coefficients per wall, in the order in which the axis intersect walls: +x -x +y -y +z -z */
    float src_pos[ROOM_SIM_MAX_NUM_SOURCES][3];   /**< Current source Cartesian coordinates, meters */
    float rec_pos[ROOM_SIM_MAX_NUM_RECEIVERS][3]; /**< Current receiver Cartesian coordinates, meters */
    CH_ORDER chOrdering;      /**< Ambisonic channel order convention (see #CH_ORDER) */
    NORM_TYPES norm;          /**< Ambisonic normalisation convention (see #NORM_TYPES) */
    
} ambi_roomsim_data;
    

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __AMBI_ROOMSIM_INTERNAL_H_INCLUDED__ */
