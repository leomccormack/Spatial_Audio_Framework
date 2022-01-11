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
 *@addtogroup VBAP
 *@{
 * @file saf_vbap.h
 * @brief Main header for the VBAP/MDAP module (#SAF_VBAP_MODULE)
 *
 * VBAP functions largely derived from the MATLAB library found in [1].
 *
 * @see [1] https://github.com/polarch/Vector-Base-Amplitude-Panning
 *          Copyright (c) 2015, Archontis Politis, BSD-3-Clause License
 *
 * @author Leo McCormack
 * @date 02.10.2017
 * @license ISC
 */

#ifndef __SAF_VBAP_H_INCLUDED__
#define __SAF_VBAP_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

/**
 * Generates a 3-D VBAP [1] gain table based on specified source and loudspeaker
 * directions, with optional spreading [2]
 *
 * @note gtable is returned as NULL if the triangulation fails. The VBAP gains
 *       are also ENERGY normalised; i.e. sum(gains^2) = 1
 *
 * @param[in]  src_dirs_deg       Source directions in degrees; FLAT: S x 2
 * @param[in]  S                  Number of Sources
 * @param[in]  ls_dirs_deg        Loudspeaker directions in degrees; FLAT: L x 2
 * @param[in]  L                  Number of loudspeakers
 * @param[in]  omitLargeTriangles '0' normal triangulation, '1' remove large
 *                                triangles
 * @param[in]  enableDummies      '0' disabled, '1' enabled, and dummies are
 *                                placed at +/-90 elevation if required
 * @param[in]  spread             Spreading factor in degrees, 0: VBAP, >0: MDAP
 * @param[out] gtable             (&) The 3D VBAP gain table energy normalised;
 *                                FLAT: N_gtable x L
 * @param[out] N_gtable           (&) number of points in the gain table
 * @param[out] nTriangles         (&) number of loudspeaker triangles
 *
 * @see [1] Pulkki, V. (1997). Virtual sound source positioning using vector
 *          base amplitude panning. Journal of the audio engineering society,
 *          45(6), 456-466.
 * @see [2] Pulkki, V. (1999). Uniform spreading of amplitude panned virtual
 *          sources. In Proceedings of the 1999 IEEE Workshop on Applications of
 *          Signal Processing to Audio and Acoustics. WASPAA'99 (Cat. No.
 *          99TH8452) (pp. 187-190). IEEE.
 */
void generateVBAPgainTable3D_srcs(/* Input arguments */
                                  float* src_dirs_deg,
                                  int S,
                                  float* ls_dirs_deg,
                                  int L,
                                  int omitLargeTriangles,
                                  int enableDummies,
                                  float spread,
                                  /* Output arguments */
                                  float** gtable,
                                  int* N_gtable,
                                  int* nTriangles);

/**
 * Generates a 3-D VBAP gain table based on specified loudspeaker directions,
 * with optional spreading [2]
 *
 * This function generates the VBAP gains for a grid: -180:az_res_deg:180
 * azimuths and -90:el_res_deg:90 elevations, which should be accessed as:
 * \code{.c}
 *   N_azi = (int)(360.0f / aziRes + 0.5f) + 1;
 *   aziIndex = (int)(matlab_fmodf(AZI + 180.0f, 360.0f)/az_res_deg + 0.5f);
 *   elevIndex = (int)((ELEV + 90.0f) / el_res_deg + 0.5f);
 *   idx3d = elevIndex * N_azi + aziIndex;
 *   for (ls = 0; ls < L; ls++)
 *       gains3D[ls] =  gtable[idx3d*L+ls];}
 * \endcode
 *
 * where 'gains3D' are the loudspeaker gains to pan the source to [AZI ELEV],
 * using the nearest grid point
 *
 * @note 'gtable' is returned as NULL if the triangulation fails. The VBAP
 *       gains are also ENERGY normalised; i.e. sum(gains^2) = 1
 *
 * @param[in]  ls_dirs_deg        Loudspeaker directions in degrees; FLAT: L x 2
 * @param[in]  L                  Number of loudspeakers
 * @param[in]  az_res_deg         Azimuthal resolution in degrees
 * @param[in]  el_res_deg         Elevation resolution in degrees
 * @param[in]  omitLargeTriangles '0' normal triangulation, '1' remove large
 *                                triangles
 * @param[in]  enableDummies      '0' disabled, '1' enabled. Dummies are placed
 *                                at +/-90 elevation if required
 * @param[in]  spread             Spreading factor in degrees, 0: VBAP, >0: MDAP
 * @param[out] gtable             (&) The 3D VBAP gain table ENERGY normalised;
 *                                FLAT: N_gtable x L
 * @param[out] N_gtable           (&) number of points in the gain table
 * @param[out] nTriangles         (&) number of loudspeaker triangles
 *
 * @see [1] Pulkki, V. (1997). Virtual sound source positioning using vector
 *          base amplitude panning. Journal of the audio engineering society,
 *          45(6), 456-466.
 * @see [2] Pulkki, V. (1999). Uniform spreading of amplitude panned virtual
 *          sources. In Proceedings of the 1999 IEEE Workshop on Applications of
 *          Signal Processing to Audio and Acoustics. WASPAA'99 (Cat. No.
 *          99TH8452) (pp. 187-190). IEEE.
 */
void generateVBAPgainTable3D(/* Input arguments */
                             float* ls_dirs_deg,
                             int L,
                             int az_res_deg,
                             int el_res_deg,
                             int omitLargeTriangles,
                             int enableDummies,
                             float spread,
                             /* Output arguments */
                             float** gtable,
                             int* N_gtable,
                             int* nTriangles);

/**
 * Compresses a VBAP gain table to use less memory and CPU (by removing the
 * elements which are just zero)
 *
 * Handy for large grid sizes for interpolation purposes. Therefore, the gains
 * are also re-normalised to have the AMPLITUDE-preserving property.
 * If vbap_gtable is generated by generateVBAPgainTable3D(), then the compressed
 * tables should be accessed as:
 * \code{.c}
 *   N_azi = (int)(360.0f / aziRes + 0.5f) + 1;
 *   aziIndex = (int)(matlab_fmodf(AZI + 180.0f, 360.0f)/az_res_deg + 0.5f);
 *   elevIndex = (int)((ELEV + 90.0f) / el_res_deg + 0.5f);
 *   idx3d = elevIndex * N_azi + aziIndex;
 *   for (i = 0; i < 3; i++){
 *       gains[i] =  vbap_gtableComp[idx3d*3+i];
 *       idx[i]   =  vbap_gtableIdx[idx3d*3+i];
 *   }
 * \endcode
 *
 * where 'gains' are then the gains for loudspeakers('idx') to pan the source to
 * [AZI ELEV], using the nearest grid point
 *
 * @note The VBAP gains are AMPLITUDE normalised; i.e. sum(gains) = 1
 *
 * @param[in]  vbap_gtable     The 3D VBAP gain table; FLAT: nTable x nDirs
 * @param[in]  nTable          number of points in the gain table
 * @param[in]  nDirs           number of loudspeakers
 * @param[out] vbap_gtableComp The compressed 3D VBAP gain table amplitude-
 *                             normalised; FLAT: nTable x 3
 * @param[out] vbap_gtableIdx  The indices for the compressed 3D VBAP gain
 *                             table; FLAT: nTable x 3
 */
void compressVBAPgainTable3D(/* Input arguments */
                             float* vbap_gtable,
                             int nTable,
                             int nDirs,
                             /* Output arguments */
                             float* vbap_gtableComp,
                             int* vbap_gtableIdx);

/**
 * Renormalises a VBAP gain table (in-place) so it may be utilised for
 * interpolation of data (for example, powermaps or HRTFs)
 *
 * @note The VBAP gains are AMPLITUDE normalised; i.e. sum(gains) = 1.
 *
 * @param[in,out] vbap_gtable The 3D VBAP gain table; FLAT: nTable x nDirs
 * @param[in]     nTable      Number of points in the gain table
 * @param[in]     nDirs       Number of loudspeaker directions
 */
void VBAPgainTable2InterpTable(float* vbap_gtable,
                               int nTable,
                               int nDirs);

/**
 * Generates a 2-D VBAP gain table based on specified source and loudspeaker
 * directions
 *
 * @note source and loudspeaker directions are required to be inter-leaved with
 *       zeros, i.e. [src_az1, 0; src_az2, 0; src_az3, 0;]. The VBAP gains are
 *       also ENERGY normalised; i.e. sum(gains^2) = 1
 *
 * @param[in]  src_dirs_deg Source directions in degrees (elev assumed to be 0
 *                          for all); FLAT: S x 2
 * @param[in]  S            Number of Sources
 * @param[in]  ls_dirs_deg  Loudspeaker directions in degrees (elev assumed to
 *                          be 0 for all); FLAT: L x 2
 * @param[in]  L            Number of loudspeakers
 * @param[out] gtable       (&) The 2D VBAP gain table energy normalised;
 *                          FLAT: S x L
 * @param[out] N_gtable     (&) number of points in the gain table, N_gtable=S
 * @param[out] nPairs       (&) number of loudspeaker pairs
 */
void generateVBAPgainTable2D_srcs(/* Input arguments */
                                  float* src_dirs_deg,
                                  int S,
                                  float* ls_dirs_deg,
                                  int L,
                                  /* Output arguments */
                                  float** gtable,
                                  int* N_gtable,
                                  int* nPairs);

/**
 * Generates a 2-D VBAP gain table based on specified loudspeaker directions
 *
 * This function generates the VBAP gains for a grid: -180:az_res_deg:180
 * azimuths, which should be accessed as:
 * \code{.c}
 *    aziIndex = (int)(matlab_fmodf(AZI + 180.0f, 360.0f)/az_res_deg + 0.5f);
 *    idx2d = aziIndex;
 *    for (ls = 0; ls < L; ls++){
 *        gains2D[ls] =  gtable[idx2d*L+ls];}
 * \endcode
 *
 * 'gains2D' are then the loudspeaker gains to pan the source to [AZI 0], using
 * the nearest grid point.
 *
 * @note The VBAP gains are ENERGY normalised; i.e. sum(gains^2) = 1
 *
 * @param[in]  ls_dirs_deg Loudspeaker directions in degrees (elev assumed to
 *                         be 0 for all); FLAT: L x 2
 * @param[in]  L           Number of loudspeakers
 * @param[in]  az_res_deg  Azimuthal resolution in degrees
 * @param[out] gtable      (&) The 2D VBAP gain table energy normalised;
 *                         FLAT: S x L
 * @param[out] N_gtable    (&) number of points in the gain table, N_gtable=S
 * @param[out] nPairs      (&) number of loudspeaker pairs
 */
void generateVBAPgainTable2D(/* Input arguments */
                             float* ls_dirs_deg,
                             int L,
                             int az_res_deg,
                             /* Output arguments */
                             float** gtable,
                             int* N_gtable,
                             int* nPairs);

/**
 * Calculates the frequency dependent pValues, which can be applied to ENERGY
 * normalised VBAP gains, to compensate for the room effect on the perceived
 * loudness fluctuations of sources when panning between loudspeakers
 *
 * This should be applied as:
 * \code{.c}
 *   if(pValues[band] != 2.0f){
 *      gains3D_sum_pvf = 0.0f;
 *      for (i = 0; i < nLoudspeakers; i++){
 *          gains3D_sum_pvf += powf(MAX(gains[i], 0.0f), pValues[band]);}
 *      gains3D_sum_pvf = powf(gains3D_sum_pvf, 1.0f/(pValues[band]+2.23e-13f));
 *      for (i = 0; i < nLoudspeakers; i++){
 *          gains_p[i] = gains[i] / (gains3D_sum_pvf+2.23e-13f);}
 *   }
 * \endcode
 *
 * Where "gains" are the original energy normalised VBAP gains and "gains_p"
 * have amplitude normalisation for the low frequencies, and energy
 * normalisation at the high frequencies [1].
 *
 * @param[in]  DTT     0..1 '0': for normal room, '1' for anechoic room, '0.5'
 *                     for listening room
 * @param[in]  freq    Frequency vector in Hz; nFreq x 1
 * @param[in]  nFreq   Number of frequencies in the frequency vector
 * @param[out] pValues pValues for each frequency; nFreq x 1
 *
 * @see [1] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V.
 *          (2014). Gain normalisation in amplitude panning as a function of
 *          frequency and room reverberance. 55th International Conference of
 *          the AES. Helsinki, Finland.
 */
void getPvalues(/* Input arguments */
                float DTT,
                float* freq,
                int nFreq,
                /* Output arguments */
                float* pValues);


/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */
/*
 * Note that the main functions are largely based on the VBAP Matlab library
 * (BSD-3-Clause License) found here:
 * https://github.com/polarch/Vector-Base-Amplitude-Panning
 */

/**
 * Computes the 3D convex-hull of a spherical grid of loudspeaker directions
 *
 * @note Compared with sphDelaunay(), this function also omits triangles where
 *       the normals and the centroid to the triangles have an angle larger than
 *       pi/2. Trianges which have an aperture larger than #APERTURE_LIMIT_DEG
 *       are also ommited.
 *
 * @param[in]  ls_dirs_deg        Loudspeaker directions in degrees; FLAT: L x 2
 * @param[in]  L                  Number of loudspeakers
 * @param[in]  omitLargeTriangles '0' normal triangulation, '1' remove large
 *                                triangles
 * @param[out] out_vertices       (&) loudspeaker directions in cartesian
 *                                coordinates; FLAT: L x 3
 * @param[out] numOutVertices     (&) number of loudspeakers
 * @param[out] out_faces          (&) loudspeaker triangle indices;
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
 * Inverts a 3x3 loudspeaker matrix
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
 * Computes a set of points which surround the source direction given a specific
 * degree of spread
 *
 * @param[in]  src_azi_rad  Source azimuth, in radians
 * @param[in]  src_elev_rad Source elevation, in radians
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
 * Calculates 3D VBAP gains given pre-computed loudspeaker triangles for each
 * source direction
 *
 * @param[in]  src_dirs     Source directions in degrees; FLAT: src_num x 2
 * @param[in]  src_num      Number of sources
 * @param[in]  ls_num       Number of loudspeakers
 * @param[in]  ls_groups    Loudspeaker triangle indices, see findLsTriplets();
 *                          FLAT: nFaces x 3
 * @param[in]  nFaces       Number of true loudspeaker triangles
 * @param[in]  spread       Spreading in degrees, 0: VBAP, >0: MDAP
 * @param[in]  layoutInvMtx Inverted 3x3 loudspeaker matrix flattened, see
 *                          invertLsMtx3D(); FLAT: nFaces x 9
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
 * @param[in]  ls_dirs_deg Loudspeaker directions in degrees; FLAT: L x 1
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
 * Inverts a 2x2 loudspeaker matrix
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
 * @param[in]  src_dirs     Source directions in degrees; FLAT: src_num x 1
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

#endif /* __SAF_VBAP_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup VBAP */
