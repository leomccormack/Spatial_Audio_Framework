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
 * Filename: saf_vbap_internal.h
 * -----------------------------
 * VBAP functions largely derived from the MATLAB library by Archontis Politis,
 * found here: https://github.com/polarch/Vector-Base-Amplitude-Panning
 *
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 02.10.2017
 */
 
#if defined(SAF_ENABLE_VBAP) && !defined(__SAF_VBAP_INTERNAL_H_INCLUDED__)
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
 
/* in degrees, if no ls_dirs have elevation +/- this value. Dummies are placed
 * at +/- 90 elevation.  */
#define ADD_DUMMY_LIMIT ( 60.0f )
/* if omitLargeTriangles==1, triangles with an aperture larger than this are
 * discarded */
#define APERTURE_LIMIT_DEG ( 180.0f )
    
    
/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/*
 * Function: findLsTriplets
 * ------------------------
 * Computes the 3D convex-hull of a spherical grid of loudspeaker directions
 *
 * Input Arguments:
 *     ls_dirs_deg        - Loudspeaker directions in DEGREES; FLAT: L x 2
 *     L                  - number of loudspeakers
 *     omitLargeTriangles - 0: normal triangulation, 1: remove large triangles
 * Output Arguments:
 *     out_vertices       - & loudspeaker directions in cartesian coordinates;
 *                          FLAT: L x 3
 *     numOutVertices     - & number of loudspeakers
 *     out_faces          - & true loudspeaker triangle indices;
 *                          FLAT: numOutFaces x 3
 *     numOutFaces        - & number of true loudspeaker triangles
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
    
/*
 * Function: invertLsMtx3D
 * ------------------------
 * Inverts the loudspeaker matrix
 *
 * Input Arguments:
 *     U_spkr       - loudspeaker directions as cartesian coordinates (unit
 *                    length); FLAT: L x 3
 *     ls_groups    - true loudspeaker triangle indices; FLAT: N_group x 3
 *     N_group      - number of true loudspeaker triangles
 * Output Arguments:
 *     layoutInvMtx - & inverted 3x3 loudspeaker matrices per group;
 *                    FLAT: N_group x 9
 */
void invertLsMtx3D(/* Input Arguments */
                   float* U_spkr,
                   int* ls_groups,
                   int N_group,
                   /* Output Arguments */
                   float** layoutInvMtx);
    
/*
 * Function: getSpreadSrcDirs3D
 * ----------------------------
 * Returns a set of points that surround the source direction with a specific
 * degree of spread.
 *
 * Input Arguments:
 *     src_azi_rad  - source azimuth, in RADIANS
 *     src_elev_rad - source elevation, in RADIANS
 *     spread       - spread in DEGREES
 *     num_src      - number of auxiliary sources to use for spreading
 *     num_rings_3d - number of concentric rings of num_src each to generate
 *                    inside the spreading surface
 * Output Arguments:
 *     U_spread     - spread directions Cartesian coordinates;
 *                    FLAT: (num_src*num_rings_3d+1) x 3
 */
void getSpreadSrcDirs3D(/* Input Arguments */
                        float src_azi_rad,
                        float src_elev_rad,
                        float spread,
                        int num_src,
                        int num_rings_3d,
                        /* Output Arguments */
                        float* U_spread);
    
/*
 * Function: vbap3D
 * ----------------
 * Calculates 3D VBAP gains for pre-calculated loudspeaker triangles and
 * predefined source directions
 *
 * Input Arguments:
 *     src_dirs     - source directions; FLAT: src_num x 2
 *     src_num      - number of sources
 *     ls_num       - number of loudspeakers
 *     ls_groups    - true loudspeaker triangle indices; FLAT: nFaces x 3
 *     nFaces       - number of true loudspeaker triangles
 *     spread       - spreading in degrees, 0: VBAP, >0: MDAP
 *     layoutInvMtx - inverted 3x3 loudspeaker matrix flattened;
 *                    FLAT: nFaces x 9
 * Output Arguments:
 *     GainMtx      - & Loudspeaker VBAP gain table; FLAT: src_num x ls_num
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
    
/*
 * Function: findLsPairs
 * ---------------------
 * Calculates loudspeaker pairs for a circular grid of loudspeaker directions
 *
 * Input Arguments:
 *     ls_dirs_deg - loudspeaker/source directions; FLAT: L x 1
 *     L           - number of loudspeakers
 * Output Arguments:
 *     out_pairs   - & loudspeaker pair indices; FLAT: numOutPairs x 2
 *     numOutPairs - & number of loudspeaker pairs
 */
void findLsPairs(/* Input Arguments */
                 float* ls_dirs_deg,
                 int L,
                 /* Output Arguments */
                 int** out_pairs,
                 int* numOutPairs);
    
/*
 * Function: invertLsMtx2D
 * ---------------------
 * Inverts the loudspeaker matrix
 *
 * Input Arguments:
 *     U_spkr       - loudspeaker directions in cartesian (xy) coordinates;
 *                    FLAT: L x 2
 *     ls_pairs     - loudspeaker pair indices; FLAT: N_pairs x 3
 *     N_pairs      - number of loudspeaker pairs
 * Output Arguments:
 *     layoutInvMtx - & inverted 2x2 loudspeaker matrix flattened;
 *                    FLAT: N_group x 4
 */
void invertLsMtx2D(/* Input Arguments */
                   float* U_spkr,
                   int* ls_pairs,
                   int N_pairs,
                   /* Output Arguments */
                   float** layoutInvMtx);
    
/*
 * Function: vbap2D
 * ---------------------
 * Calculates 2D VBAP gains for pre-calculated loudspeaker pairs and predefined
 * source positions
 *
 * Input Arguments:
 *     src_dirs     - source directions in DEGREES; FLAT: src_num x 1
 *     src_num      - number of sources
 *     ls_num       - number of loudspeakers
 *     ls_pairs     - loudspeaker pair indices; FLAT: N_pairs x 2
 *     N_pairs      - number of loudspeaker pairs
 *     layoutInvMtx - inverted 2x2 loudspeaker matrix flattened;
 *                    FLAT: N_pairs x 4
 * Output Arguments:
 *     GainMtx      - & Loudspeaker VBAP gain table; FLAT: src_num x ls_num
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

#endif /* defined(SAF_ENABLE_VBAP) && !defined(__SAF_VBAP_INTERNAL_H_INCLUDED__) */
