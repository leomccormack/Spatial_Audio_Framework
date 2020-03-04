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
 * @file saf_vbap_internal.h
 * @brief Internal part of the "saf_vbap" module
 *
 * VBAP functions largely derived from the MATLAB library by Archontis Politis,
 * found in [1].
 *
 * @see [1] https://github.com/polarch/Vector-Base-Amplitude-Panning
 *
 * @author Leo McCormack
 * @date 02.10.2017
 */

#ifndef __SAF_VBAP_INTERNAL_H_INCLUDED__
#define __SAF_VBAP_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "saf_vbap.h"
#include "../saf_utilities/saf_utilities.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/** In degrees, if no ls_dirs have elevation +/- this value. Dummies are placed
 * at +/- 90 elevation.  */
#define ADD_DUMMY_LIMIT ( 60.0f )
/** if omitLargeTriangles==1, triangles with an aperture larger than this are
 * discarded */
#define APERTURE_LIMIT_DEG ( 180.0f )

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Computes the 3D convex-hull of a spherical grid of loudspeaker directions
 *
 * @param[in]  ls_dirs_deg        Loudspeaker directions in DEGREES; FLAT: L x 2
 * @param[in]  L                  Number of loudspeakers
 * @param[in]  omitLargeTriangles '0' normal triangulation, '1' remove large
 *                                triangles
 * @param[out] out_vertices       (&) loudspeaker directions in cartesian
 *                                coordinates; FLAT: L x 3
 * @param[out] numOutVertices     (&) number of loudspeakers
 * @param[out] out_faces          (&) true loudspeaker triangle indices;
 *                                FLAT: numOutFaces x 3
 * @param[out] numOutFaces        (&) number of true loudspeaker triangles
 */
void findLsTriplets(/* Input Arguments */
                    float* ls_dirs_deg,
                    int L,
                    int omitLargeTriangles,
                    /* Output Arguments */
                    float** out_vertices,
                    int* numOutVertices,
                    int** out_faces,
                    int* numOutFaces);

/**
 * Inverts a loudspeaker matrix
 *
 * @param[in]  U_spkr       Loudspeaker directions as cartesian coordinates
 *                          (unit length); FLAT: L x 3
 * @param[in]  ls_groups    True loudspeaker triangle indices; FLAT: N_group x 3
 * @param[in]  N_group      Number of true loudspeaker triangles
 * @param[out] layoutInvMtx (&) inverted 3x3 loudspeaker matrices per group;
 *                          FLAT: N_group x 9
 */
void invertLsMtx3D(/* Input Arguments */
                   float* U_spkr,
                   int* ls_groups,
                   int N_group,
                   /* Output Arguments */
                   float** layoutInvMtx);

/**
 * Computes a set of points that surround the source direction with a specific
 * degree of spread.
 *
 * @param[in]  src_azi_rad  Source azimuth, in RADIANS
 * @param[in]  src_elev_rad Source elevation, in RADIANS
 * @param[in]  spread       Spread in DEGREES
 * @param[in]  num_src      Number of auxiliary sources to use for spreading
 * @param[in]  num_rings_3d Number of concentric rings of num_src each to
 *                          generate inside the spreading surface
 * @param[out] U_spread     Spread directions Cartesian coordinates;
 *                          FLAT: (num_src*num_rings_3d+1) x 3
 */
void getSpreadSrcDirs3D(/* Input Arguments */
                        float src_azi_rad,
                        float src_elev_rad,
                        float spread,
                        int num_src,
                        int num_rings_3d,
                        /* Output Arguments */
                        float* U_spread);

/**
 * Calculates 3D VBAP gains for pre-calculated loudspeaker triangles and
 * predefined source directions
 *
 * @param[in]  src_dirs     Source directions; FLAT: src_num x 2
 * @param[in]  src_num      Number of sources
 * @param[in]  ls_num       Number of loudspeakers
 * @param[in]  ls_groups    True loudspeaker triangle indices; FLAT: nFaces x 3
 * @param[in]  nFaces       Number of true loudspeaker triangles
 * @param[in]  spread       Spreading in degrees, 0: VBAP, >0: MDAP
 * @param[in]  layoutInvMtx Inverted 3x3 loudspeaker matrix flattened;
 *                          FLAT: nFaces x 9
 * @param[out] GainMtx      (&) Loudspeaker VBAP gain table;
 *                          FLAT: src_num x ls_num
 */
void vbap3D(/* Input Arguments */
            float* src_dirs,
            int src_num,
            int ls_num,
            int* ls_groups,
            int nFaces,
            float spread,
            float* layoutInvMtx,
            /* Output Arguments */
            float** GainMtx);

/**
 * Calculates loudspeaker pairs for a circular grid of loudspeaker directions
 *
 * @param[in]  ls_dirs_deg Loudspeaker/source directions; FLAT: L x 1
 * @param[in]  L           Number of loudspeakers
 * @param[out] out_pairs   (&) loudspeaker pair indices; FLAT: numOutPairs x 2
 * @param[out] numOutPairs (&) number of loudspeaker pairs
 */
void findLsPairs(/* Input Arguments */
                 float* ls_dirs_deg,
                 int L,
                 /* Output Arguments */
                 int** out_pairs,
                 int* numOutPairs);

/**
 * Inverts the loudspeaker matrix
 *
 * @param[in]  U_spkr       Loudspeaker directions in cartesian (xy)
 *                          coordinates; FLAT: L x 2
 * @param[in]  ls_pairs     Loudspeaker pair indices; FLAT: N_pairs x 3
 * @param[in]  N_pairs      Number of loudspeaker pairs
 * @param[out] layoutInvMtx (&) inverted 2x2 loudspeaker matrix flattened;
 *                          FLAT: N_group x 4
 */
void invertLsMtx2D(/* Input Arguments */
                   float* U_spkr,
                   int* ls_pairs,
                   int N_pairs,
                   /* Output Arguments */
                   float** layoutInvMtx);

/**
 * Calculates 2D VBAP gains for pre-calculated loudspeaker pairs and predefined
 * source positions
 *
 * @param[in]  src_dirs     Source directions in DEGREES; FLAT: src_num x 1
 * @param[in]  src_num      Number of sources
 * @param[in]  ls_num       Number of loudspeakers
 * @param[in]  ls_pairs     Loudspeaker pair indices; FLAT: N_pairs x 2
 * @param[in]  N_pairs      Number of loudspeaker pairs
 * @param[in]  layoutInvMtx Inverted 2x2 loudspeaker matrix flattened;
 *                          FLAT: N_pairs x 4
 * @param[out] GainMtx      (&) Loudspeaker VBAP gain table;
 *                          FLAT: src_num x ls_num
 */
void vbap2D(/* Input Arguments */
            float* src_dirs,
            int src_num,
            int ls_num,
            int* ls_pairs,
            int N_pairs,
            float* layoutInvMtx,
            /* Output Arguments */
            float** GainMtx);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_VBAP_INTERNAL_H_INCLUDED__ */
