/*
 * Copyright 2016-2018 Leo McCormack
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
 *@addtogroup Utilities
 *@{
 * @file saf_utility_loudspeaker_presets.h
 * @brief A collection of loudspeaker array directions and (nearly) uniform
 *        spherical grids
 *
 * @author Leo McCormack
 * @date 11.07.2016
 * @license ISC
 */

#ifndef __SAF_LOUDSPEAKER_PRESETS_H_INCLUDED__
#define __SAF_LOUDSPEAKER_PRESETS_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
#include <stdio.h>  
    
/* ========================================================================== */
/*                      Loudspeaker Array Configurations                      */
/* ========================================================================== */
/*
 * NOTE: All loudspeaker array sensor directions are given in degrees, and in
 * the [azimuth, elevation] convention. [0 0] is looking directly in-front
 * (positive x-axis), with positive elevations looking upwards (positive z-axis)
 * and positive azimuth angles looking leftwards (positive y-axis).
 * In other words: the convention used by SAF follows the "right-hand-rule".
 */

/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a mono setup */
extern const float __mono_dirs_deg[1][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a stereo setup */
extern const float __stereo_dirs_deg[2][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a 5.x setup */
extern const float __5pX_dirs_deg[5][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a 7.x setup */
extern const float __7pX_dirs_deg[7][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a 8.x setup */
extern const float __8pX_dirs_deg[8][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a 9.x setup */
extern const float __9pX_dirs_deg[9][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a 10.x setup */
extern const float __10pX_dirs_deg[10][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a 11.x setup */
extern const float __11pX_dirs_deg[11][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a 7.4.x setup */
extern const float __11pX_7_4_dirs_deg[11][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a 13.x setup */
extern const float __13pX_dirs_deg[13][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a 22.x setup */
extern const float __22pX_dirs_deg[22][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for a 9+10+3.2 setup
 * BS 2051 recommedation: https://www.itu.int/rec/R-REC-BS.2051/en */
extern const float __9_10_3p2_dirs_deg[24][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for the multi-channel
 * anechoic chamber (MCC), at Aalto University. */
extern const float __Aalto_MCC_dirs_deg[45][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for the multi-channel
 * anechoic chamber (MCC) sub-set, at Aalto University. */
extern const float __Aalto_MCCsubset_dirs_deg[37][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for the audio-visual
 * listening room (Apaja), at Aalto University. */
extern const float __Aalto_Apaja_dirs_deg[29][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for the ITU standard
 * listening room (LR), at Aalto University. */
extern const float __Aalto_LR_dirs_deg[13][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for the Audio Visual
 * Immersion Lab (AVIL), at the Technical University of Denmark (DTU) */
extern const float __DTU_AVIL_dirs_deg[64][2];
/**
 * Loudspeaker directions [azimuth, Elevation] in degrees, for the 22.x setup,
 * at Zylia Labs */
extern const float __Zylia_Lab_dirs_deg[22][2];
/**
 * Default Loudspeaker directions [azimuth, Elevation] in RADIANS! */
extern const float __default_LScoords64_rad[64][2];
    
    
/* ========================================================================== */
/*                                  T-designs                                 */
/* ========================================================================== */

/** Number of directions in a minimum Tdesign of degree: 1 */
extern const int   __Tdesign_degree_1_nPoints;
/** Number of directions in a minimum Tdesign of degree: 2 */
extern const int   __Tdesign_degree_2_nPoints;
/** Number of directions in a minimum Tdesign of degree: 3 */
extern const int   __Tdesign_degree_3_nPoints;
/** Number of directions in a minimum Tdesign of degree: 4 */
extern const int   __Tdesign_degree_4_nPoints;
/** Number of directions in a minimum Tdesign of degree: 5 */
extern const int   __Tdesign_degree_5_nPoints;
/** Number of directions in a minimum Tdesign of degree: 6 */
extern const int   __Tdesign_degree_6_nPoints;
/** Number of directions in a minimum Tdesign of degree: 7 */
extern const int   __Tdesign_degree_7_nPoints;
/** Number of directions in a minimum Tdesign of degree: 8 */
extern const int   __Tdesign_degree_8_nPoints;
/** Number of directions in a minimum Tdesign of degree: 9 */
extern const int   __Tdesign_degree_9_nPoints;
/** Number of directions in a minimum Tdesign of degree: 10 */
extern const int   __Tdesign_degree_10_nPoints;
/** Number of directions in a minimum Tdesign of degree: 11 */
extern const int   __Tdesign_degree_11_nPoints;
/** Number of directions in a minimum Tdesign of degree: 12 */
extern const int   __Tdesign_degree_12_nPoints;
/** Number of directions in a minimum Tdesign of degree: 13 */
extern const int   __Tdesign_degree_13_nPoints;
/** Number of directions in a minimum Tdesign of degree: 14 */
extern const int   __Tdesign_degree_14_nPoints;
/** Number of directions in a minimum Tdesign of degree: 15 */
extern const int   __Tdesign_degree_15_nPoints;
/** Number of directions in a minimum Tdesign of degree: 16 */
extern const int   __Tdesign_degree_16_nPoints;
/** Number of directions in a minimum Tdesign of degree: 17 */
extern const int   __Tdesign_degree_17_nPoints;
/** Number of directions in a minimum Tdesign of degree: 18 */
extern const int   __Tdesign_degree_18_nPoints;
/** Number of directions in a minimum Tdesign of degree: 19 */
extern const int   __Tdesign_degree_19_nPoints;
/** Number of directions in a minimum Tdesign of degree: 20 */
extern const int   __Tdesign_degree_20_nPoints;
/** Number of directions in a minimum Tdesign of degree: 21 */
extern const int   __Tdesign_degree_21_nPoints;
/** Number of directions in a minimum Tdesign of degree: 30 */
extern const int   __Tdesign_degree_30_nPoints;
/** Number of directions in a minimum Tdesign of degree: 40 */
extern const int   __Tdesign_degree_40_nPoints;
/** Number of directions in a minimum Tdesign of degree: 50 */
extern const int   __Tdesign_degree_50_nPoints;
/** Number of directions in a minimum Tdesign of degree: 60 */
extern const int   __Tdesign_degree_60_nPoints;
/** Number of directions in a minimum Tdesign of degree: 70 */
extern const int   __Tdesign_degree_70_nPoints;
/** Number of directions in a minimum Tdesign of degree: 80 */
extern const int   __Tdesign_degree_80_nPoints;
/** Number of directions in a minimum Tdesign of degree: 90 */
extern const int   __Tdesign_degree_90_nPoints;
/** Number of directions in a minimum Tdesign of degree: 100 */
extern const int   __Tdesign_degree_100_nPoints;
/** Number of directions in a minimum Tdesign of degree: 124 */
extern const int   __Tdesign_degree_124_nPoints;
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 1
 * @see McLaren's Improved Snub Cube and Other New Spherical Designs in Three
 *      Dimensions", R. H. Hardin and N. J. A. Sloane, Discrete and
 *      Computational Geometry, 15 (1996), pp. 429-441. */
extern const float __Tdesign_degree_1_dirs_deg[2][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 2 */
extern const float __Tdesign_degree_2_dirs_deg[4][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 3 */
extern const float __Tdesign_degree_3_dirs_deg[6][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 4 */
extern const float __Tdesign_degree_4_dirs_deg[12][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 5 */
extern const float __Tdesign_degree_5_dirs_deg[12][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 6 */
extern const float __Tdesign_degree_6_dirs_deg[24][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 7 */
extern const float __Tdesign_degree_7_dirs_deg[24][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 8 */
extern const float __Tdesign_degree_8_dirs_deg[36][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 9 */
extern const float __Tdesign_degree_9_dirs_deg[48][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 10 */
extern const float __Tdesign_degree_10_dirs_deg[60][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 11 */
extern const float __Tdesign_degree_11_dirs_deg[70][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 12 */
extern const float __Tdesign_degree_12_dirs_deg[84][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 13 */
extern const float __Tdesign_degree_13_dirs_deg[94][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 14 */
extern const float __Tdesign_degree_14_dirs_deg[108][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 15 */
extern const float __Tdesign_degree_15_dirs_deg[120][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 16 */
extern const float __Tdesign_degree_16_dirs_deg[144][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 17 */
extern const float __Tdesign_degree_17_dirs_deg[156][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 18 */
extern const float __Tdesign_degree_18_dirs_deg[180][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 19 */
extern const float __Tdesign_degree_19_dirs_deg[204][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 20 */
extern const float __Tdesign_degree_20_dirs_deg[216][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 21 */
extern const float __Tdesign_degree_21_dirs_deg[240][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 30
 * @see Gra"f, M., & Potts, D. (2011). On the computation of spherical designs
 *      by a new optimization approach based on fast spherical Fourier
 *      transforms. Numerische Mathematik, 119(4), 699-724.
 * Taken from:
 * github.com/chris-hld/spaudiopy/blob/master/spaudiopy/n_designs_1_124.mat
 */
extern const float __Tdesign_degree_30_dirs_deg[480][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 40 */
extern const float __Tdesign_degree_40_dirs_deg[840][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 50 */
extern const float __Tdesign_degree_50_dirs_deg[1296][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 60 */
extern const float __Tdesign_degree_60_dirs_deg[1860][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 70 */
extern const float __Tdesign_degree_70_dirs_deg[2520][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 80 */
extern const float __Tdesign_degree_80_dirs_deg[3276][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 90 */
extern const float __Tdesign_degree_90_dirs_deg[4140][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 100
 */
extern const float __Tdesign_degree_100_dirs_deg[5100][2];
/**
 * Directions [azimuth, Elevation] in degrees, for minimum Tdesign degree: 124
 */
extern const float __Tdesign_degree_124_dirs_deg[7812][2];
/**
 * minimum T-design HANDLES (up to degree 21 only).
 *
 * Access as, e.g.
 * \code{.c}
 *   const int tdesign_degree = 7;
 *   float* tdesign_dirs_deg = __HANDLES_Tdesign_dirs_deg [tdesign_degree-1];
 *   int num_Tdesign_dirs = __Tdesign_nPoints_per_degree  [tdesign_degree-1];
 * \endcode
 */
extern const float* __HANDLES_Tdesign_dirs_deg[21];
/**
 * Number of points in each t-design (up to degree 21 only).
 *
 * Access as, e.g.
 * \code{.c}
 *   const int tdesign_degree = 7;
 *   float* tdesign_dirs_deg = __HANDLES_Tdesign_dirs_deg [tdesign_degree-1];
 *   int num_Tdesign_dirs = __Tdesign_nPoints_per_degree  [tdesign_degree-1];
 * \endcode
 */
extern const int __Tdesign_nPoints_per_degree[21];


/* ========================================================================== */
/*                              Sphere Coverings                              */
/* ========================================================================== */

/** Directions [azimuth, Elevation] in degrees, for sphere covering: 4 dirs
 * @see Belger, M. (1989). JH Conway, NJA Sloane. Sphere packings, lattices and
 *      groups. Springer Verlag New York--Berlin--Heidelberg--London--Paris
 *      Tokyo 1988, 663 pages, 112 illustrations, DM 178.00, ISBN 0--387--96617
 *      --X. Crystal Research and Technology, 24(1), 90--90. */
extern const float __SphCovering_4_dirs_deg[4][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 5 dirs */
extern const float __SphCovering_5_dirs_deg[5][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 6 dirs */
extern const float __SphCovering_6_dirs_deg[6][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 7 dirs */
extern const float __SphCovering_7_dirs_deg[7][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 8 dirs */
extern const float __SphCovering_8_dirs_deg[8][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 9 dirs */
extern const float __SphCovering_9_dirs_deg[9][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 10 dirs */
extern const float __SphCovering_10_dirs_deg[10][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 11 dirs */
extern const float __SphCovering_11_dirs_deg[11][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 12 dirs */
extern const float __SphCovering_12_dirs_deg[12][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 13 dirs */
extern const float __SphCovering_13_dirs_deg[13][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 14 dirs */
extern const float __SphCovering_14_dirs_deg[14][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 15 dirs */
extern const float __SphCovering_15_dirs_deg[15][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 16 dirs */
extern const float __SphCovering_16_dirs_deg[16][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 17 dirs */
extern const float __SphCovering_17_dirs_deg[17][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 18 dirs */
extern const float __SphCovering_18_dirs_deg[18][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 19 dirs */
extern const float __SphCovering_19_dirs_deg[19][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 20 dirs */
extern const float __SphCovering_20_dirs_deg[20][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 21 dirs */
extern const float __SphCovering_21_dirs_deg[21][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 22 dirs */
extern const float __SphCovering_22_dirs_deg[22][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 23 dirs */
extern const float __SphCovering_23_dirs_deg[23][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 24 dirs */
extern const float __SphCovering_24_dirs_deg[24][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 25 dirs */
extern const float __SphCovering_25_dirs_deg[25][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 26 dirs */
extern const float __SphCovering_26_dirs_deg[26][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 27 dirs */
extern const float __SphCovering_27_dirs_deg[27][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 28 dirs */
extern const float __SphCovering_28_dirs_deg[28][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 29 dirs */
extern const float __SphCovering_29_dirs_deg[29][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 30 dirs */
extern const float __SphCovering_30_dirs_deg[30][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 31 dirs */
extern const float __SphCovering_31_dirs_deg[31][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 32 dirs */
extern const float __SphCovering_32_dirs_deg[32][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 33 dirs */
extern const float __SphCovering_33_dirs_deg[33][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 34 dirs */
extern const float __SphCovering_34_dirs_deg[34][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 35 dirs */
extern const float __SphCovering_35_dirs_deg[35][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 36 dirs */
extern const float __SphCovering_36_dirs_deg[36][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 37 dirs */
extern const float __SphCovering_37_dirs_deg[37][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 38 dirs */
extern const float __SphCovering_38_dirs_deg[38][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 39 dirs */
extern const float __SphCovering_39_dirs_deg[39][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 40 dirs */
extern const float __SphCovering_40_dirs_deg[40][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 41 dirs */
extern const float __SphCovering_41_dirs_deg[41][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 42 dirs */
extern const float __SphCovering_42_dirs_deg[42][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 43 dirs */
extern const float __SphCovering_43_dirs_deg[43][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 44 dirs */
extern const float __SphCovering_44_dirs_deg[44][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 45 dirs */
extern const float __SphCovering_45_dirs_deg[45][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 46 dirs */
extern const float __SphCovering_46_dirs_deg[46][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 47 dirs */
extern const float __SphCovering_47_dirs_deg[47][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 48 dirs */
extern const float __SphCovering_48_dirs_deg[48][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 49 dirs */
extern const float __SphCovering_49_dirs_deg[49][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 50 dirs */
extern const float __SphCovering_50_dirs_deg[50][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 51 dirs */
extern const float __SphCovering_51_dirs_deg[51][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 52 dirs */
extern const float __SphCovering_52_dirs_deg[52][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 53 dirs */
extern const float __SphCovering_53_dirs_deg[53][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 54 dirs */
extern const float __SphCovering_54_dirs_deg[54][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 55 dirs */
extern const float __SphCovering_55_dirs_deg[55][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 56 dirs */
extern const float __SphCovering_56_dirs_deg[56][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 57 dirs */
extern const float __SphCovering_57_dirs_deg[57][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 58 dirs */
extern const float __SphCovering_58_dirs_deg[58][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 59 dirs */
extern const float __SphCovering_59_dirs_deg[59][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 60 dirs */
extern const float __SphCovering_60_dirs_deg[60][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 61 dirs */
extern const float __SphCovering_61_dirs_deg[61][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 62 dirs */
extern const float __SphCovering_62_dirs_deg[62][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 63 dirs */
extern const float __SphCovering_63_dirs_deg[63][2];
/**
 * Directions [azimuth, Elevation] in degrees, for sphere covering: 64 dirs */
extern const float __SphCovering_64_dirs_deg[64][2];
/**
 * Sphere covering handles ( between 4..64 points only)
 *
 * Access as, e.g.
 * \code{.c}
 *   const int numPoints = 44; // between 4..64 points
 *   float* sphCov_dirs_deg = __HANDLES_SphCovering_dirs_deg [numPoints-1];
 * \endcode
 */
extern const float* __HANDLES_SphCovering_dirs_deg[64];


/* ========================================================================== */
/*                                 Geosphere                                  */
/* ========================================================================== */

/** Number of directions in a ico geosphere of degree: 0 */
extern const int   __geosphere_ico_0_0_nPoints;
/** Number of directions in a ico geosphere of degree: 1 */
extern const int   __geosphere_ico_1_0_nPoints;
/** Number of directions in a ico geosphere of degree: 2 */
extern const int   __geosphere_ico_2_0_nPoints;
/** Number of directions in a ico geosphere of degree: 3 */
extern const int   __geosphere_ico_3_0_nPoints;
/** Number of directions in a ico geosphere of degree: 4 */
extern const int   __geosphere_ico_4_0_nPoints;
/** Number of directions in a ico geosphere of degree: 5 */
extern const int   __geosphere_ico_5_0_nPoints;
/** Number of directions in a ico geosphere of degree: 6 */
extern const int   __geosphere_ico_6_0_nPoints;
/** Number of directions in a ico geosphere of degree: 7 */
extern const int   __geosphere_ico_7_0_nPoints;
/** Number of directions in a ico geosphere of degree: 8 */
extern const int   __geosphere_ico_8_0_nPoints;
/** Number of directions in a ico geosphere of degree: 9 */
extern const int   __geosphere_ico_9_0_nPoints;
/** Number of directions in a ico geosphere of degree: 10 */
extern const int   __geosphere_ico_10_0_nPoints;
/** Number of directions in a ico geosphere of degree: 11 */
extern const int   __geosphere_ico_11_0_nPoints;
/** Number of directions in a ico geosphere of degree: 12 */
extern const int   __geosphere_ico_12_0_nPoints;
/** Number of directions in a ico geosphere of degree: 13 */
extern const int   __geosphere_ico_13_0_nPoints;
/** Number of directions in a ico geosphere of degree: 14 */
extern const int   __geosphere_ico_14_0_nPoints;
/** Number of directions in a ico geosphere of degree: 15 */
extern const int   __geosphere_ico_15_0_nPoints;
/** Number of directions in a ico geosphere of degree: 16 */
extern const int   __geosphere_ico_16_0_nPoints;
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 0
 * @see 3LD - Library for Loudspeaker Layout Design; release 2, 2006/03/15
 *      Copyright (c) 2006 Florian Hollerweger (floholl_AT_sbox.tugraz.at) and
 *      (c) 2002 Darren Weber. */
extern const float __geosphere_ico_0_0_dirs_deg[12][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 1 */
extern const float __geosphere_ico_1_0_dirs_deg[32][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 2 */
extern const float __geosphere_ico_2_0_dirs_deg[42][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 3 */
extern const float __geosphere_ico_3_0_dirs_deg[92][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 4 */
extern const float __geosphere_ico_4_0_dirs_deg[162][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 5 */
extern const float __geosphere_ico_5_0_dirs_deg[252][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 6 */
extern const float __geosphere_ico_6_0_dirs_deg[362][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 7 */
extern const float __geosphere_ico_7_0_dirs_deg[492][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 8 */
extern const float __geosphere_ico_8_0_dirs_deg[642][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 9 */
extern const float __geosphere_ico_9_0_dirs_deg[812][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 10 */
extern const float __geosphere_ico_10_0_dirs_deg[1002][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 11 */
extern const float __geosphere_ico_11_0_dirs_deg[1212][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 12 */
extern const float __geosphere_ico_12_0_dirs_deg[1442][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 13 */
extern const float __geosphere_ico_13_0_dirs_deg[1692][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 14 */
extern const float __geosphere_ico_14_0_dirs_deg[1962][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 15 */
extern const float __geosphere_ico_15_0_dirs_deg[2252][2];
/**
 * Directions [azimuth, Elevation] in degrees, for ico geosphere, degree: 16 */
extern const float __geosphere_ico_16_0_dirs_deg[2562][2];

/** Number of directions in a oct geosphere of degree: 0 */
extern const int   __geosphere_oct_0_0_nPoints;
/** Number of directions in a oct geosphere of degree: 1 */
extern const int   __geosphere_oct_1_0_nPoints;
/** Number of directions in a oct geosphere of degree: 2 */
extern const int   __geosphere_oct_2_0_nPoints;
/** Number of directions in a oct geosphere of degree: 3 */
extern const int   __geosphere_oct_3_0_nPoints;
/** Number of directions in a oct geosphere of degree: 4 */
extern const int   __geosphere_oct_4_0_nPoints;
/** Number of directions in a oct geosphere of degree: 5 */
extern const int   __geosphere_oct_5_0_nPoints;
/** Number of directions in a oct geosphere of degree: 6 */
extern const int   __geosphere_oct_6_0_nPoints;
/** Number of directions in a oct geosphere of degree: 7 */
extern const int   __geosphere_oct_7_0_nPoints;
/** Number of directions in a oct geosphere of degree: 8 */
extern const int   __geosphere_oct_8_0_nPoints;
/** Number of directions in a oct geosphere of degree: 9 */
extern const int   __geosphere_oct_9_0_nPoints;
/** Number of directions in a oct geosphere of degree: 10 */
extern const int   __geosphere_oct_10_0_nPoints;
/** Number of directions in a oct geosphere of degree: 11 */
extern const int   __geosphere_oct_11_0_nPoints;
/** Number of directions in a oct geosphere of degree: 12 */
extern const int   __geosphere_oct_12_0_nPoints;
/** Number of directions in a oct geosphere of degree: 13 */
extern const int   __geosphere_oct_13_0_nPoints;
/** Number of directions in a oct geosphere of degree: 14 */
extern const int   __geosphere_oct_14_0_nPoints;
/** Number of directions in a oct geosphere of degree: 15 */
extern const int   __geosphere_oct_15_0_nPoints;
/** Number of directions in a oct geosphere of degree: 16 */
extern const int   __geosphere_oct_16_0_nPoints;
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 0
 * @see 3LD - Library for Loudspeaker Layout Design; release 2, 2006/03/15
 *      Copyright (c) 2006 Florian Hollerweger (floholl_AT_sbox.tugraz.at) and
 *      (c) 2002 Darren Weber. */

extern const float __geosphere_oct_0_0_dirs_deg[6][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 1 */
extern const float __geosphere_oct_1_0_dirs_deg[14][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 2 */
extern const float __geosphere_oct_2_0_dirs_deg[18][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 3 */
extern const float __geosphere_oct_3_0_dirs_deg[38][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 4 */
extern const float __geosphere_oct_4_0_dirs_deg[66][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 5 */
extern const float __geosphere_oct_5_0_dirs_deg[102][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 6 */
extern const float __geosphere_oct_6_0_dirs_deg[146][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 7 */
extern const float __geosphere_oct_7_0_dirs_deg[198][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 8 */
extern const float __geosphere_oct_8_0_dirs_deg[258][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 9 */
extern const float __geosphere_oct_9_0_dirs_deg[326][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 10 */
extern const float __geosphere_oct_10_0_dirs_deg[402][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 11 */
extern const float __geosphere_oct_11_0_dirs_deg[486][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 12 */
extern const float __geosphere_oct_12_0_dirs_deg[578][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 13 */
extern const float __geosphere_oct_13_0_dirs_deg[678][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 14 */
extern const float __geosphere_oct_14_0_dirs_deg[786][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 15 */
extern const float __geosphere_oct_15_0_dirs_deg[902][2];
/**
 * Directions [azimuth, Elevation] in degrees, for oct geosphere, degree: 16 */
extern const float __geosphere_oct_16_0_dirs_deg[1026][2];
/**
 * 3LD geosphere HANDLES (freq = [0 0], [1 0],..., [16 0])
 *
 * Access as, e.g.
 * \code{.c}
 *   const int degree = 10; // between 0..16
 *   float* geo_dirs_degree_10 = __HANDLES_geosphere_ico_dirs_deg[degree];
 *   int num_dirs_degree_10 = __geosphere_ico_nPoints[degree];
 * \endcode
 */
extern const float* __HANDLES_geosphere_ico_dirs_deg[17];
/**
 * 3LD geosphere HANDLES (freq = [0 0], [1 0],..., [16 0])
 *
 * Access as, e.g.
 * \code{.c}
 *   const int degree = 10; // between 0..16
 *   float* geo_dirs_degree_10 = __HANDLES_geosphere_oct_dirs_deg[degree];
 *   int num_dirs_degree_10 = __geosphere_oct_nPoints[degree];
 * \endcode
 */
extern const float* __HANDLES_geosphere_oct_dirs_deg[17];
/**
 * 3LD geosphere number of points (freq = [0 0], [1 0],..., [16 0])
 *
 * Access as, e.g.
 * \code{.c}
 *   const int degree = 10; // between 0..16
 *   float* geo_dirs_degree_10 = __HANDLES_geosphere_ico_dirs_deg[degree];
 *   int num_dirs_degree_10 = __geosphere_ico_nPoints[degree];
 * \endcode
 */
extern const int __geosphere_ico_nPoints[17];
/**
 * 3LD geosphere number of points (freq = [0 0], [1 0],..., [16 0])
 *
 * Access as, e.g.
 * \code{.c}
 *   const int degree = 10; // between 0..16
 *   float* geo_dirs_degree_10 = __HANDLES_geosphere_oct_dirs_deg[degree];
 *   int num_dirs_degree_10 = __geosphere_oct_nPoints[degree];
 * \endcode
 */
extern const int __geosphere_oct_nPoints[17];
    

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_LOUDSPEAKER_PRESETS_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup Utilities */
