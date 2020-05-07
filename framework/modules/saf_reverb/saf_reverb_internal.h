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
 * @file saf_reverb_internal.h
 * @brief Internal part of the reverb processing module (saf_reverb)
 *
 * ...
 *
 * @author Leo McCormack
 * @date 06.05.2020
 */

#ifndef __REVERB_INTERNAL_H_INCLUDED__
#define __REVERB_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h> 
#include <string.h>
#include <assert.h>
#include "saf_reverb.h"
#include "../saf_utilities/saf_utilities.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef NUM_EARS
# define NUM_EARS 2
#endif

typedef void* voidPtr;

typedef struct _position_xyz {
    union {
        struct { float x, y, z; };
        float v[3];
    };
} position_xyz;

typedef struct _ims_core_workspace
{
    /* Locals */
    int room[3];
    float d_max;
    position_xyz src, rec;

    /* Intermediates */
    float Nx, Ny, Nz;
    int lengthVec;
    float* II, *JJ, *KK;
    float* s_x, *s_y, *s_z, *s_d, *s_t, *s_att;


}ims_core_workspace;


typedef struct _ims_scene_data
{
    /* Copy of input arguments */
    int room_dimensions[3];
    float c_ms;
    int nBands;
    float** abs_wall;

    /* Source and receiver positions */
    position_xyz* src_xyz;
    position_xyz* rec_xyz;
    long* src_IDs;
    long* rec_IDs;
    long nSources;
    long nRecievers;

    /* Internal */
    float* band_centerfreqs;
    float max_time_s;
    voidPtr** ims_core_work;

} ims_scene_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

void ims_shoebox_coreWorkspaceCreate(void** hWork);

/*
Calculates an echogram of a rectangular space using ISM.
%   IMS_CORE calculates an echogram of a rectangular space using
%   the Image-Source Method, for a given source and receiver. Input
%   argument room should be a structure with the following fields:
%   room-->length, width, height, absorption. room.absoprtion is a 2x3 matrix
%   with absorption coefficients (broadband) for each of the walls on the
%   respective planes [x+ y+ z+; x- y- z-].
%
%   source and receiver are structures holding the coordinates of the
%   source/receiver as: source.coord = [Sx Sy Sz]. There are plans to
%   include directivity coefficients in the source structure.
%
%   Coordinates of source/receiver are specified from the left ground corner
%   of the room:
%                ^x
%              __|__    _
%             |  |  |   |
%             |  |  |   |
%          y<----.  |   | l
%             |     |   |
%             |     |   |
%             o_____|   -
%
%             |-----|
%                w
%
*/
void ims_shoebox_core(void* hWork,
                      int room[3],
                      position_xyz src,
                      position_xyz rec,
                      float maxTime,
                      float c_ms);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __REVERB_INTERNAL_H_INCLUDED__ */
