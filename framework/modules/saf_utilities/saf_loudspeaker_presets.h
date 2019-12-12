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

/*
 * Filename: saf_loudspeaker_presets.h
 * -----------------------------------
 * Comprises the directions of loudspeaker arrays and (nearly) uniform spherical
 * grids.
 *
 * Dependencies:
 *     none
 * Author, date created:
 *     Leo McCormack, 11.07.2016
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
    
extern const float __mono_dirs_deg[1][2];
extern const float __stereo_dirs_deg[2][2];
extern const float __5pX_dirs_deg[5][2];
extern const float __7pX_dirs_deg[7][2];
extern const float __8pX_dirs_deg[8][2];
extern const float __9pX_dirs_deg[9][2];
extern const float __10pX_dirs_deg[10][2];
extern const float __11pX_dirs_deg[11][2];
extern const float __11pX_7_4_dirs_deg[11][2];
extern const float __13pX_dirs_deg[13][2];
extern const float __22pX_dirs_deg[22][2];
extern const float __Aalto_MCC_dirs_deg[44][2];
extern const float __Aalto_Apaja_dirs_deg[29][2];
extern const float __Aalto_LR_dirs_deg[13][2];
extern const float __DTU_AVIL_dirs_deg[64][2];
extern const float __Zylia_Lab_dirs_deg[22][2];
extern const float default_LScoords64_rad[64][2];
    
    
/* ========================================================================== */
/*                        (Nearly) Uniform Arrangements                       */
/* ========================================================================== */
    
/* Minimum T-designs for given degrees
 *     McLaren's Improved Snub Cube and Other New Spherical Designs in Three
 *     Dimensions", R. H. Hardin and N. J. A. Sloane, Discrete and Computational
 *     Geometry, 15 (1996), pp. 429-441. */
extern const float __Tdesign_degree_1_dirs_deg[2][2];
extern const float __Tdesign_degree_2_dirs_deg[4][2];
extern const float __Tdesign_degree_3_dirs_deg[6][2];
extern const float __Tdesign_degree_4_dirs_deg[12][2];
extern const float __Tdesign_degree_5_dirs_deg[12][2];
extern const float __Tdesign_degree_6_dirs_deg[24][2];
extern const float __Tdesign_degree_7_dirs_deg[24][2];
extern const float __Tdesign_degree_8_dirs_deg[36][2];
extern const float __Tdesign_degree_9_dirs_deg[48][2];
extern const float __Tdesign_degree_10_dirs_deg[60][2];
extern const float __Tdesign_degree_11_dirs_deg[70][2];
extern const float __Tdesign_degree_12_dirs_deg[84][2];
extern const float __Tdesign_degree_13_dirs_deg[94][2];
extern const float __Tdesign_degree_14_dirs_deg[108][2];
extern const float __Tdesign_degree_15_dirs_deg[120][2];
extern const float __Tdesign_degree_16_dirs_deg[144][2];
extern const float __Tdesign_degree_17_dirs_deg[156][2];
extern const float __Tdesign_degree_18_dirs_deg[180][2];
extern const float __Tdesign_degree_19_dirs_deg[204][2];
extern const float __Tdesign_degree_20_dirs_deg[216][2];
extern const float __Tdesign_degree_21_dirs_deg[240][2];
/*     Gräf, M., & Potts, D. (2011). On the computation of spherical designs by
 *     a new optimization approach based on fast spherical Fourier transforms.
 *     Numerische Mathematik, 119(4), 699-724. */
extern const float __Tdesign_degree_30_dirs_deg[480][2];
extern const float __Tdesign_degree_40_dirs_deg[840][2];
extern const float __Tdesign_degree_50_dirs_deg[1296][2];
extern const float __Tdesign_degree_100_dirs_deg[5100][2];

/* minimum T-design HANDLES (up to degree 21 only) */
extern const float* __HANDLES_Tdesign_dirs_deg[21];
extern const int __Tdesign_nPoints_per_degree[21];
    
/* Sphere Covering
 *     Belger, M. (1989). JH Conway, NJA Sloane. Sphere packings, lattices and
 *     groups. Springer Verlag New York‐Berlin‐Heidelberg‐London‐Paris Tokyo
 *     1988, 663 pages, 112 illustrations, DM 178.00, ISBN 0‐387‐96617‐X.
 *     Crystal Research and Technology, 24(1), 90-90. */
extern const float __SphCovering_4_dirs_deg[4][2];
extern const float __SphCovering_5_dirs_deg[5][2];
extern const float __SphCovering_6_dirs_deg[6][2];
extern const float __SphCovering_7_dirs_deg[7][2];
extern const float __SphCovering_8_dirs_deg[8][2];
extern const float __SphCovering_9_dirs_deg[9][2];
extern const float __SphCovering_10_dirs_deg[10][2];
extern const float __SphCovering_11_dirs_deg[11][2];
extern const float __SphCovering_12_dirs_deg[12][2];
extern const float __SphCovering_13_dirs_deg[13][2];
extern const float __SphCovering_14_dirs_deg[14][2];
extern const float __SphCovering_15_dirs_deg[15][2];
extern const float __SphCovering_16_dirs_deg[16][2];
extern const float __SphCovering_17_dirs_deg[17][2];
extern const float __SphCovering_18_dirs_deg[18][2];
extern const float __SphCovering_19_dirs_deg[19][2];
extern const float __SphCovering_20_dirs_deg[20][2];
extern const float __SphCovering_21_dirs_deg[21][2];
extern const float __SphCovering_22_dirs_deg[22][2];
extern const float __SphCovering_23_dirs_deg[23][2];
extern const float __SphCovering_24_dirs_deg[24][2];
extern const float __SphCovering_25_dirs_deg[25][2];
extern const float __SphCovering_26_dirs_deg[26][2];
extern const float __SphCovering_27_dirs_deg[27][2];
extern const float __SphCovering_28_dirs_deg[28][2];
extern const float __SphCovering_29_dirs_deg[29][2];
extern const float __SphCovering_30_dirs_deg[30][2];
extern const float __SphCovering_31_dirs_deg[31][2];
extern const float __SphCovering_32_dirs_deg[32][2];
extern const float __SphCovering_33_dirs_deg[33][2];
extern const float __SphCovering_34_dirs_deg[34][2];
extern const float __SphCovering_35_dirs_deg[35][2];
extern const float __SphCovering_36_dirs_deg[36][2];
extern const float __SphCovering_37_dirs_deg[37][2];
extern const float __SphCovering_38_dirs_deg[38][2];
extern const float __SphCovering_39_dirs_deg[39][2];
extern const float __SphCovering_40_dirs_deg[40][2];
extern const float __SphCovering_41_dirs_deg[41][2];
extern const float __SphCovering_42_dirs_deg[42][2];
extern const float __SphCovering_43_dirs_deg[43][2];
extern const float __SphCovering_44_dirs_deg[44][2];
extern const float __SphCovering_45_dirs_deg[45][2];
extern const float __SphCovering_46_dirs_deg[46][2];
extern const float __SphCovering_47_dirs_deg[47][2];
extern const float __SphCovering_48_dirs_deg[48][2];
extern const float __SphCovering_49_dirs_deg[49][2];
extern const float __SphCovering_50_dirs_deg[50][2];
extern const float __SphCovering_51_dirs_deg[51][2];
extern const float __SphCovering_52_dirs_deg[52][2];
extern const float __SphCovering_53_dirs_deg[53][2];
extern const float __SphCovering_54_dirs_deg[54][2];
extern const float __SphCovering_55_dirs_deg[55][2];
extern const float __SphCovering_56_dirs_deg[56][2];
extern const float __SphCovering_57_dirs_deg[57][2];
extern const float __SphCovering_58_dirs_deg[58][2];
extern const float __SphCovering_59_dirs_deg[59][2];
extern const float __SphCovering_60_dirs_deg[60][2];
extern const float __SphCovering_61_dirs_deg[61][2];
extern const float __SphCovering_62_dirs_deg[62][2];
extern const float __SphCovering_63_dirs_deg[63][2];
extern const float __SphCovering_64_dirs_deg[64][2];
    
/* sphere covering handles ( between 4..64 points only) */
extern const float* __HANDLES_SphCovering_dirs_deg[64];
    
/* 3LD geosphere
 *     3LD - Library for Loudspeaker Layout Design; release 2, 2006/03/15
 *     Copyright (c) 2006 Florian Hollerweger (floholl_AT_sbox.tugraz.at) and
 *     (c) 2002 Darren Weber. */
extern const float __geosphere_ico_0_0_dirs_deg[12][2];
extern const float __geosphere_ico_1_0_dirs_deg[32][2];
extern const float __geosphere_ico_2_0_dirs_deg[42][2];
extern const float __geosphere_ico_3_0_dirs_deg[92][2];
extern const float __geosphere_ico_4_0_dirs_deg[162][2];
extern const float __geosphere_ico_5_0_dirs_deg[252][2];
extern const float __geosphere_ico_6_0_dirs_deg[362][2];
extern const float __geosphere_ico_7_0_dirs_deg[492][2];
extern const float __geosphere_ico_8_0_dirs_deg[642][2];
extern const float __geosphere_ico_9_0_dirs_deg[812][2];
extern const float __geosphere_ico_10_0_dirs_deg[1002][2];
extern const float __geosphere_ico_11_0_dirs_deg[1212][2];
extern const float __geosphere_ico_12_0_dirs_deg[1442][2];
extern const float __geosphere_ico_13_0_dirs_deg[1692][2];
extern const float __geosphere_ico_14_0_dirs_deg[1962][2];
extern const float __geosphere_ico_15_0_dirs_deg[2252][2];
extern const float __geosphere_ico_16_0_dirs_deg[2562][2];
extern const float __geosphere_oct_0_0_dirs_deg[6][2];
extern const float __geosphere_oct_1_0_dirs_deg[14][2];
extern const float __geosphere_oct_2_0_dirs_deg[18][2];
extern const float __geosphere_oct_3_0_dirs_deg[38][2];
extern const float __geosphere_oct_4_0_dirs_deg[66][2];
extern const float __geosphere_oct_5_0_dirs_deg[102][2];
extern const float __geosphere_oct_6_0_dirs_deg[146][2];
extern const float __geosphere_oct_7_0_dirs_deg[198][2];
extern const float __geosphere_oct_8_0_dirs_deg[258][2];
extern const float __geosphere_oct_9_0_dirs_deg[326][2];
extern const float __geosphere_oct_10_0_dirs_deg[402][2];
extern const float __geosphere_oct_11_0_dirs_deg[486][2];
extern const float __geosphere_oct_12_0_dirs_deg[578][2];
extern const float __geosphere_oct_13_0_dirs_deg[678][2];
extern const float __geosphere_oct_14_0_dirs_deg[786][2];
extern const float __geosphere_oct_15_0_dirs_deg[902][2];
extern const float __geosphere_oct_16_0_dirs_deg[1026][2];
    
/* 3LD geosphere HANDLES (freq = [0 0], [1 0],..., [16 0]) */
/*     3LD - Library for Loudspeaker Layout Design; release 2, 2006/03/15
 *     Copyright (c) 2006 Florian Hollerweger (floholl_AT_sbox.tugraz.at) and
 *     (c) 2002 Darren Weber. */
extern const float* __HANDLES_geosphere_ico_dirs_deg[17];
extern const float* __HANDLES_geosphere_oct_dirs_deg[17];
extern const int __geosphere_ico_nPoints[17];
extern const int __geosphere_oct_nPoints[17];
    

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_LOUDSPEAKER_PRESETS_H_INCLUDED__ */

