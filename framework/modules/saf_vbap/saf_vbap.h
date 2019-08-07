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
 * Filename: saf_vbap.h
 * --------------------
 * VBAP functions largely derived from the MATLAB library by Archontis Politis,
 * found here: https://github.com/polarch/Vector-Base-Amplitude-Panning
 * 
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 02.10.2017
 */

#ifndef __SAF_VBAP_H_INCLUDED__
#define __SAF_VBAP_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: generateVBAPgainTable3D_srcs
 * --------------------------------------
 * Generates a 3-D VBAP gain table based on specified source and loudspeaker
 * directions;
 * Note: gtable is returned as NULL if the triangulation failed
 * Further Note: The VBAP gains are ENERGY normalised; i.e. sum(gains^2) = 1
 *
 * Input Arguments:
 *     src_dirs_deg       - Source directions in DEGREES; FLAT: S x 2
 *     S                  - number of Sources
 *     ls_dirs_deg        - Loudspeaker directions in DEGREES; FLAT: L x 2
 *     L                  - number of loudspeakers
 *     omitLargeTriangles - 0: normal triangulation, 1: remove large triangles
 *     enableDummies      - 0: disabled, 1: enabled. Dummies are placed at +/-90
 *                          elevation if required
 *     spread             - spreading factor in DEGREES, 0: VBAP, >0: MDAP
 * Output Arguments:
 *     gtable             - & The 3D VBAP gain table ENERGY NORMALISED;
 *                          FLAT: N_gtable x L
 *     N_gtable           - & number of points in the gain table
 *     nTriangles         - & number of loudspeaker triangles
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

/*
 * Function: generateVBAPgainTable3D
 * ---------------------------------
 * Generates a 3-D VBAP gain table based on specified loudspeaker directions;
 * Note: gtable is returned as NULL if the triangulation failed
 * This function generates the VBAP gains for a grid: -180:az_res_deg:180
 * azimuths and -90:el_res_deg:90 elevations, which should be accessed as:
 *      N_azi = (int)(360.0f / aziRes + 0.5f) + 1;
 *      aziIndex = (int)(matlab_fmodf(AZI + 180.0f, 360.0f)/az_res_deg + 0.5f);
 *      elevIndex = (int)((ELEV + 90.0f) / el_res_deg + 0.5f);
 *      idx3d = elevIndex * N_azi + aziIndex;
 *      for (ls = 0; ls < L; ls++)
 *          gains3D[ls] =  gtable[idx3d*L+ls];
 *
 * gains3D are the loudspeaker gains to pan the source to [AZI ELEV], using the
 * nearest grid point
 * Further Note: The VBAP gains are ENERGY normalised; i.e. sum(gains^2) = 1
 *
 * Input Arguments:
 *     ls_dirs_deg        - Loudspeaker directions in DEGREES; FLAT: L x 2
 *     L                  - number of loudspeakers
 *     az_res_deg         - azimuthal resolution in DEGREES
 *     el_res_deg         - elevation resolution in DEGREES
 *     omitLargeTriangles - 0: normal triangulation, 1: remove large triangles
 *     enableDummies      - 0: disabled, 1: enabled. Dummies are placed at +/-90
 *                          elevation if required
 *     spread             - spreading factor in DEGREES, 0: VBAP, >0: MDAP
 * Output Arguments:
 *     gtable             - & The 3D VBAP gain table ENERGY NORMALISED;
 *                          FLAT: N_gtable x L
 *     N_gtable           - & number of points in the gain table
 *     nTriangles         - & number of loudspeaker triangles
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

/*
 * Function: compressVBAPgainTable3D
 * ---------------------------------
 * Compresses a VBAP gain table to use less memory and CPU (by removing the
 * elements that are zero). Handy for large grid sizes for interpolation
 * purposes. Therefore, the gains are also re-normalised to have the AMPLITUDE-
 * preserving property.
 * If 'vbap_gtable' is generated by generateVBAPgainTable3D, then the compressed
 * tables should be accessed as:
 *      N_azi = (int)(360.0f / aziRes + 0.5f) + 1;
 *      aziIndex = (int)(matlab_fmodf(AZI + 180.0f, 360.0f)/az_res_deg + 0.5f);
 *      elevIndex = (int)((ELEVATION + 90.0f) / el_res_deg + 0.5f);
 *      idx3d = elevIndex * N_azi + aziIndex;
 *      for (i = 0; i < 3; i++){
 *          gains[i] =  vbap_gtableComp[idx3d*3+i];
 *          idx[i] =  vbap_gtableIdx[idx3d*3+i];
 *      }
 *
 * gains are then the gains for loudspeakers(idx) to pan the source to [AZI ELEV], using the nearest grid point
 * NOTE: The VBAP gains are AMPLITUDE normalised; i.e. sum(gains) = 1
 *
 * Input Arguments:
 *     vbap_gtable     - The 3D VBAP gain table; nTable x nDirs
 *     nTable          - number of points in the gain table
 *     nDirs           - number of loudspeakers
 * Output Arguments:
 *     vbap_gtableComp - & The compressed 3D VBAP gain table AMPLITUDE-
 *                       NORMALISED; FLAT: nTable x 3
 *     vbap_gtableIdx  - & The indices for the compressed 3D VBAP gain table;
 *                       FLAT: nTable x 3
 */
void compressVBAPgainTable3D(/* Input arguments */
                             float* vbap_gtable,
                             int nTable,
                             int nDirs,
                             /* Output arguments */
                             float** vbap_gtableComp,
                             int** vbap_gtableIdx);
    
/*
 * Function: VBAPgainTable2InterpTable
 * -----------------------------------
 * Renormalises a vbap gain table in place, so it may be utilised for
 * interpolation of data (e.g. powermaps or HRTFs).
 * Note: The VBAP gains are AMPLITUDE normalised; i.e. sum(gains) = 1
 * Further Note: vbap_gtable is renormalised "in-place"
 *
 * Input Arguments:
 *     vbap_gtable - The 3D VBAP gain table; nTable x nDirs
 *     nTable      - number of points in the gain table
 *     nDirs       - number of loudspeakers
 */
void VBAPgainTable2InterpTable(float* vbap_gtable,
                               int nTable,
                               int nDirs);
    
/*
 * Function: generateVBAPgainTable2D_srcs
 * --------------------------------------
 * Generates a 2-D VBAP gain table based on specified source and loudspeaker
 * directions;
 * Note: source and loudspeaker directions are required to be inter-leaved, i.e.
 *     [src_az1, 0; src_az2, 0; src_az3, 0;]
 * Further Note: The VBAP gains are ENERGY normalised; i.e. sum(gains^2) = 1
 *
 * Input Arguments:
 *     src_dirs_deg - Source directions in DEGREES (elev assumed to be 0
 *                    for all); FLAT: S x 2
 *     S            - number of Sources
 *     ls_dirs_deg  - Loudspeaker directions in DEGREES (elev assumed to
 *                    be 0 for all); FLAT: L x 2
 *     L            - number of loudspeakers
 * Output Arguments:
 *     gtable       - & The 2D VBAP gain table ENERGY NORMALISED;
 *                    FLAT: S x L
 *     N_gtable     - & number of points in the gain table, N_gtable=S
 *     nPairs       - & number of loudspeaker pairs
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
    
/*
 * Function: generateVBAPgainTable2D
 * ---------------------------------
 * Generates a 2-D VBAP gain table based on specified loudspeaker directions.
 * Note that this function generates the VBAP gains for a grid:
 * -180:az_res_deg:180 azimuths. Which should be accessed as:
 *      aziIndex = (int)(matlab_fmodf(AZI + 180.0f, 360.0f)/az_res_deg + 0.5f);
 *      idx2d = aziIndex;
 *      for (ls = 0; ls < L; ls++){
 *          gains2D[ls] =  gtable[idx2d*L+ls];}
 *
 * gains2D are then the loudspeaker gains to pan the source to [AZI 0],
 * using the nearest grid point
 * Note: The VBAP gains are ENERGY normalised; i.e. sum(gains^2) = 1
 *
 * Input Arguments:
 *     ls_dirs_deg - Loudspeaker directions in DEGREES (elev assumed to
 *                    be 0 for all); FLAT: L x 2
 *     L           - number of loudspeakers
 *     az_res_deg  - azimuthal resolution in DEGREES
 * Output Arguments:
 *     gtable      - & The 2D VBAP gain table ENERGY NORMALISED;
 *                   FLAT: S x L
 *     N_gtable    - & number of points in the gain table, N_gtable=S
 *     nPairs      - & number of loudspeaker pairs
 */
void generateVBAPgainTable2D(/* Input arguments */
                             float* ls_dirs_deg,
                             int L,
                             int az_res_deg,
                             /* Output arguments */
                             float** gtable,
                             int* N_gtable,
                             int* nPairs);
    
/*
 * Function: getPvalues
 * --------------------
 * Calculates the frequency dependent pValues, which can be applied to ENERGY,
 * normalised VBAP gains. This is performed as:
 * if(pValues[band] != 2.0f){
 *     gains3D_sum_pvf = 0.0f;
 *     for (i = 0; i < nLoudspeakers; i++){
 *          gains3D_sum_pvf += powf(MAX(gains[i], 0.0f), pValues[band]);}
 *     gains3D_sum_pvf = powf(gains3D_sum_pvf, 1.0f/(pValues[band]+2.23e-13f));
 *     for (i = 0; i < nLoudspeakers; i++){
 *         gains_p[i] = gains[i] / (gains3D_sum_pvf+2.23e-13f);}
 * }
 *
 * Where "gains" are the original energy normalised VBAP gains and "gains_p"
 * have amplitude normalisation for the low frequencies, and energy
 * normalisation at the high frequencies [1].
 *
 * Input Arguments:
 *     DTT     - 0..1 0: for normal room, 1: for anechoic room, 0.5: for
 *               listening room
 *     freq    - frequency vector in Hz; nFreq x 1
 *     nFreq   - number of frequencies in the frequency vector
 * Output Arguments:
 *     pValues - pValues for each frequency; nFreq x 1
 *
 * [1] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V. (2014).
 *     Gain normalisation in amplitude panning as a function of frequency and
 *     room reverberance. 55th International Conference of the AES. Helsinki,
 *     Finland.
 */
void getPvalues(/* Input arguments */
                float DTT,
                float* freq,
                int nFreq,
                /* Output arguments */
                float* pValues);
    
    
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __SAF_VBAP_H_INCLUDED__ */
