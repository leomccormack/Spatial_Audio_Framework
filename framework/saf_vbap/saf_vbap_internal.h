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
 *     saf_vbap_internal.h
 * Description:
 *    vbap functions largely derived from the MATLAB library by Archontis Politis,
 *    found here: https://github.com/polarch/Vector-Base-Amplitude-Panning
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 02.10.2017
 */
 
#ifndef __SAF_VBAP_INTERNAL_H_INCLUDED__
#define __SAF_VBAP_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "saf_vbap.h"
#include "saf_utilities.h"

#ifdef __cplusplus
extern "C" {
#endif
 
#define ADD_DUMMY_LIMIT ( 60.0f )             /* in degrees, if no ls_dirs have elevation +/- this value. Dummies are placed at +/- 90 elevation.  */
#define MAX_NUM_FACES ( 5000 )                /* avoids infinite loops in the 3d convexhull main loop */
#define APERTURE_LIMIT_DEG ( 180.0f )         /* if omitLargeTriangles==1, triangles with an aperture larger than this are discarded */
#ifndef M_PI
  #define M_PI ( 3.14159265359f )
#endif
    
/* Calculates the 3D convex-hull of a spherical grid of loudspeaker directions */
void findLsTriplets(float* ls_dirs_deg,       /* loudspeaker/source directions; FLAT: L x 2 */
                    int L,                    /* number of loudspeakers */
                    int omitLargeTriangles,   /* 0: normal triangultion, 1: remove large triangles */
                    float** out_vertices,     /* & loudspeaker directions in cartesian coordinates; FLAT: L x 3 */
                    int* numOutVertices,      /* & number of loudspeakers */
                    int** out_faces,          /* & true loudspeaker triangle indices; FLAT: numOutFaces x 3 */
                    int* numOutFaces);        /* & number of true loudspeaker triangles */
    
/* Inverts the loudspeaker matrix */
void invertLsMtx3D(float* U_spkr,             /* loudspeaker directions in cartesian coordinates; FLAT: L x 3 */ 
                   int* ls_groups,            /* true loudspeaker triangle indices; FLAT: N_group x 3 */
                   int N_group,               /* number of true loudspeaker triangles */
                   float** layoutInvMtx);     /* & inverted 3x3 loudspeaker matrices per group; FLAT: N_group x 9 */
    
/* Calculates 3D VBAP gains for pre-calculated loudspeaker triangles and predefined source positions */
void vbap3D(float* src_dirs,                  /* source directions; FLAT: src_num x 2 */
            int src_num,                      /* number of sources */
            int ls_num,                       /* number of loudspeakers */
            int* ls_groups,                   /* true loudspeaker triangle indices; FLAT: nFaces x 3 */
            int nFaces,                       /* number of true loudspeaker triangles */
            float* layoutInvMtx,              /* inverted 3x3 loudspeaker matrix flattened; FLAT: nFaces x 9 */
            float** GainMtx);                 /* & Loudspeaker VBAP gain table; FLAT: src_num x ls_num */
    
/* Calculates loudspeaker pairs for a circular grid of loudspeaker directions  */
void findLsPairs(float* ls_dirs_deg,          /* loudspeaker/source directions; FLAT: L x 1 */
                 int L,                       /* number of loudspeakers */
                 int** out_pairs,             /* & loudspeaker pair indices; FLAT: numOutPairs x 2 */
                 int* numOutPairs);           /* & number of loudspeaker pairs */
    
/* Inverts the loudspeaker matrix */
void invertLsMtx2D(float* U_spkr,             /* loudspeaker directions in cartesian (xy) coordinates; FLAT: L x 2 */ 
                   int* ls_pairs,             /* loudspeaker pair indices; FLAT: N_pairs x 3 */
                   int N_pairs,               /* number of loudspeaker pairs */
                   float** layoutInvMtx);     /* & inverted 2x2 loudspeaker matrix flattened; FLAT: N_group x 4 */
    
/* Calculates 2D VBAP gains for pre-calculated loudspeaker pairs and predefined source positions */
void vbap2D(float* src_dirs,                  /* source directions; FLAT: src_num x 1 */
            int src_num,                      /* number of sources */
            int ls_num,                       /* number of loudspeakers */
            int* ls_pairs,                    /* loudspeaker pair indices; FLAT: N_pairs x 2 */
            int N_pairs,                      /* number of loudspeaker pairs */
            float* layoutInvMtx,              /* inverted 2x2 loudspeaker matrix flattened; FLAT: N_pairs x 4 */
            float** GainMtx);                 /* & Loudspeaker VBAP gain table; FLAT: src_num x ls_num */
    
#ifdef __cplusplus
}
#endif

#endif /* __SAF_VBAP_INTERNAL_H_INCLUDED__ */




















