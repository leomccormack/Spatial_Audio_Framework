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
 * Filename: saf_sh.h (include header)
 * -----------------------------------
 * A collection of spherical harmonic related functions. Many of which have been
 * derived from Matlab libraries by Archontis Politis; found here:
 *     https://github.com/polarch/Spherical-Harmonic-Transform
 *     https://github.com/polarch/Array-Response-Simulator
 *     https://github.com/polarch/Spherical-Array-Processing
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 22.05.2016
 */

#ifndef __SAF_SH_H_INCLUDED__
#define __SAF_SH_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
#include "../saf_utilities/saf_complex.h"

/* ========================================================================== */
/*                                    Enums                                   */
/* ========================================================================== */
    
/*
 * Enum: ARRAY_CONSTRUCTION_TYPES
 * ------------------------------
 * Microphone/Hydrophone array contruction types
 *
 * Options:
 *     ARRAY_CONSTRUCTION_OPEN              - Open array, omni-directional
 *                                            sensors
 *     ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL  - Open array, directional sensors
 *     ARRAY_CONSTRUCTION_RIGID             - Rigid baffle, omni-directional
 *                                            sensors
 *     ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL - Rigid baffle, directional sensors
 */
typedef enum _ARRAY_CONSTRUCTION_TYPES {
    ARRAY_CONSTRUCTION_OPEN,
    ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL,
    ARRAY_CONSTRUCTION_RIGID,
    ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL
    
}ARRAY_CONSTRUCTION_TYPES;
    
/*
 * Enum: SECTOR_PATTERNS
 * ---------------------
 * Sector pattern designs for directionally-contraining sound-fields [1].
 *
 * Options:
 *     SECTOR_PATTERN_PWD      - plane-wave decomposition/Hyper-cardioid
 *     SECTOR_PATTERN_MAXRE    - Spatially tapered hyper-cardioid, such that it
 *                               has maximum energy concentrated in the look-
 *                               direction
 *     SECTOR_PATTERN_CARDIOID - cardioid pattern
 *
 * [1] Politis, A., & Pulkki, V. (2016). Acoustic intensity, energy-density and
 *     diffuseness estimation in a directionally-constrained region. arXiv
 *     preprint arXiv:1609.03409
 */
typedef enum _SECTOR_PATTERNS{
    SECTOR_PATTERN_PWD,
    SECTOR_PATTERN_MAXRE,
    SECTOR_PATTERN_CARDIOID
}SECTOR_PATTERNS;


/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

/*
 * Function: yawPitchRoll2Rzyx
 * ---------------------------
 * Contructs a 3x3 rotation matrix from the Euler angles, using the
 * yaw-pitch-roll (zyx) convention
 *
 * Input Arguments:
 *     yaw              - yaw angle in radians
 *     pitch            - pitch angle in radians
 *     roll             - roll angle in radians
 *     rollPitchYawFLAG - 1: use Rxyz, i.e. apply roll, pitch and then yaw,
 *                        0: Rzyx / y-p-r
 * Output Arguments:
 *     R                - zyx rotation matrix; 3 x 3
 */
void yawPitchRoll2Rzyx (/* Input Arguments */
                        float yaw,
                        float pitch,
                        float roll,
                        int rollPitchYawFLAG,
                        /* Output Arguments */
                        float R[3][3]);

/*
 * Function: unitSph2Cart
 * ----------------------
 * Converts spherical coordinates to cartesian coordinates of unit length
 *
 * Input Arguments:
 *     azi_rad  - azimuth in radians
 *     elev_rad - elevation in radians
 * Output Arguments:
 *     xyz      - unit cartesian coords, xyz; 3 x 1
 */
void unitSph2Cart(/* Input Arguments */
                  float azi_rad,
                  float elev_rad,
                  /* Output Arguments */
                  float xyz[3]);

/*
 * Function: unitCart2Sph
 * ----------------------
 * Converts cartesian coordinates of unit length to spherical coordinates
 *
 * Input Arguments:
 *     xyz         - unit cartesian coords, xyz
 * Output Arguments:
 *     AziElev_rad - azimuth and elevation in radians
 */
void unitCart2Sph(/* Input Arguments */
                  float xyz[3],
                  /* Output Arguments */
                  float AziElev_rad[2]);

/*
 * Function: unitCart2Sph_aziElev
 * ------------------------------
 * Converts cartesian coordinates of unit length to spherical coordinates
 *
 * Input Arguments:
 *     xyz      - unit cartesian coords, xyz
 * Output Arguments:
 *     azi_rad  - & azimuth in radians
 *     elev_rad - & elevation in radians
 */
void unitCart2Sph_aziElev(/* Input Arguments */
                          float xyz[3],
                          /* Output Arguments */
                          float* azi_rad,
                          float* elev_rad);
    
    
/* ========================================================================== */
/*                    SH and Beamforming related Functions                    */
/* ========================================================================== */
    
/*
 * Function: unnorm_legendreP
 * --------------------------
 * Calculates unnormalised legendre polynomials up to order N, for all values
 * in vector x [1].
 * Note: This INCLUDES the Condon-Shortley phase term. It is functionally
 * identical to MatLab's legendre function, 'unnorm' (with default settings).
 *
 * Input Arguments:
 *     n    - order of  legendre polynomial
 *     x    - vector of input values; lenX x 1
 *     lenX - number of input values
 * Output Arguments:
 *     y  - resulting unnormalised legendre values for each x value;
 *          FLAT: (n+1) x lenX
 *
 * [1] M, Abramowitz., I.A. Stegun. (1965). "Handbook of Mathematical
 *     Functions: Chapter 8", Dover Publications.
 */
void unnorm_legendreP(/* Input Arguments */
                      int n,
                      double* x,
                      int lenX,
                      /* Output Arguments */
                      double* y);

/*
 * Function: unnorm_legendreP
 * --------------------------
 * Calculates unnormalised legendre polynomials up to order N, for all values
 * in vector x. It uses a recursive approach, which makes it more suitable for
 * computing the legendre values in a real-time loop.
 * Note: This does NOT INCLUDE the Condon-Shortley phase term.
 *
 * Input Arguments:
 *     n          - order of  legendre polynomial
 *     x          - vector of input values; lenX x 1
 *     lenX       - number of input values
 *     Pnm_minus1 - previous Pnm, (not used for n=1); FLAT: (n+1) x lenX
 *     Pnm_minus2 - previous previous Pnm, (not used for n=0);
 *                  FLAT: (n+1) x lenX
 * Output Arguments:
 *     Pnm        - resulting unnormalised legendre values for each x value;
 *                  FLAT: (n+1) x lenX
 */
void unnorm_legendreP_recur(/* Input Arguments */
                            int n,
                            float* x,
                            int lenX,
                            float* Pnm_minus1,
                            float* Pnm_minus2,
                            /* Output Arguments */
                            float* Pnm);
    
/*
 * Function: getRSH
 * ----------------
 * This function returns REAL spherical harmonics [1] for multiple directions on
 * the sphere. WITHOUT the 1/sqrt(4*pi) term. i.e. max(omni) = 1
 * Note: Compared to 'getRSH_recur', this function uses 'unnorm_legendreP' and
 * double precision, so is more suitable for determining 'Y' in an
 * initialisation stage. This version is indeed slower, but more precise;
 * especially for high orders.
 * Further Note: this function is mainly INTENDED FOR AMBISONICS, due to the
 * omission of the 1/sqrt(4*pi) scaling, and the directions are given in
 * [azimuth elevation] (degrees).
 * In Ambisonics literature, the format convention of 'Y' is referred to as
 * ACN/N3D
 *
 * Input Arguments:
 *     order    - order of spherical harmonic expansion
 *     dirs_deg - directions on the sphere [azi, ELEVATION] convention, in
 *                DEGREES; FLAT: nDirs x 2
 *     nDirs    - number of directions
 * Output Arguments:
 *     Y        - & the SH weights [WITHOUT the 1/sqrt(4*pi)];
 *                FLAT: (order+1)^2 x nDirs
 *
 * [1] Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8).
 *     Berlin: Springer.
 */
void getRSH(/* Input Arguments */
            int order,
            float* dirs_deg,
            int nDirs,
            /* Output Arguments */
            float** Y);

/*
 * Function: getRSH_recur
 * ----------------------
 * This function returns REAL spherical harmonics [1] for multiple directions on
 * the sphere. WITHOUT the 1/sqrt(4*pi) term. i.e. max(omni) = 1
 * Note: Compared to 'getRSH', this function uses 'unnorm_legendreP_recur' and
 * single precision, so is more suitable for determining 'Y' in a real-time
 * loop. It sacrifices some precision, as numerical error propogates through
 * the recursion, but it is faster.
 * Further Note: this function is mainly INTENDED FOR AMBISONICS, due to the
 * omission of the 1/sqrt(4*pi) scaling, and the directions are given in
 * [azimuth elevation] (degrees).
 * In Ambisonics literature, the format convention of 'Y' is referred to as
 * ACN/N3D
 *
 * Input Arguments:
 *     order    - order of spherical harmonic expansion
 *     dirs_deg - directions on the sphere [azi, ELEVATION] convention, in
 *                DEGREES; FLAT: nDirs x 2
 *     nDirs    - number of directions
 * Output Arguments:
 *     Y        - & the SH weights [WITHOUT the 1/sqrt(4*pi)];
 *                FLAT: (order+1)^2 x nDirs
 *
 * [1] Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8).
 *     Berlin: Springer.
 */
void getRSH_recur(/* Input Arguments */
                  int order,
                  float* dirs_deg,
                  int nDirs,
                  /* Output Arguments */
                  float** Y);
    
/*
 * Function: getSHreal
 * -------------------
 * This function returns REAL spherical harmonics [1] for each direction on the
 * sphere. WITH the 1/sqrt(4*pi) term.  i.e. max(omni) = 1/sqrt(4*pi)
 * Note: Compared to 'getSHreal_recur', this function employs 'unnorm_legendreP'
 * and double precision, which is slower but more precise.
 *
 * Input Arguments:
 *     order    - order of spherical harmonic expansion
 *     dirs_rad - directions on the sphere [azi, INCLINATION] convention, in
 *                RADIANS; FLAT: nDirs x 2
 *     nDirs    - number of directions
 * Output Arguments:
 *     Y        - & the SH weights [WITH the 1/sqrt(4*pi)];
 *                FLAT: (order+1)^2 x nDirs
 *
 * [1] Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8).
 *     Berlin: Springer.
 */
void getSHreal(/* Input Arguments */
               int order,
               float* dirs_rad,
               int nDirs,
               /* Output Arguments */
               float* Y);
    
/*
 * Function: getSHreal_recur
 * -------------------------
 * This function returns REAL spherical harmonics [1] for each direction on the
 * sphere. WITH the 1/sqrt(4*pi) term.  i.e. max(omni) = 1/sqrt(4*pi)
 * Note: Compared to 'getSHreal', this function employs 'unnorm_legendreP_recur'
 * and single precision, which is faster but less precise.
 *
 * Input Arguments:
 *     order    - order of spherical harmonic expansion
 *     dirs_rad - directions on the sphere [azi, INCLINATION] convention,
 *                in RADIANS; FLAT: nDirs x 2
 *     nDirs    - number of directions
 * Output Arguments:
 *     Y        - & the SH weights [WITH the 1/sqrt(4*pi)];
 *                FLAT: (order+1)^2 x nDirs
 *
 * [1] Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8).
 *     Berlin: Springer.
 */
void getSHreal_recur(/* Input Arguments */
                     int order,
                     float* dirs_rad,
                     int nDirs,
                     /* Output Arguments */
                     float* Y);
    
/*
 * Function: getSHcomplex
 * ----------------------
 * This function returns COMPLEX spherical harmonics [1]for each direction on
 * the sphere. WITH the 1/sqrt(4*pi) term.  i.e. max(omni) = 1/sqrt(4*pi) + i0
 * Note: This function employs 'unnorm_legendreP' and double precision.
 *
 * Input Arguments:
 *     order    - order of spherical harmonic expansion
 *     dirs_rad - directions on the sphere [azi, INCLINATION] convention, in
 *                RADIANS; FLAT: nDirs x 2
 *     nDirs    - number of directions
 * Output Arguments:
 *     Y        - & the SH weights [WITH the 1/sqrt(4*pi)];
 *                FLAT: (order+1)^2 x nDirs
 *
 * [1] Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8).
 *     Berlin: Springer.
 */
void getSHcomplex(/* Input Arguments */
                  int order,
                  float* dirs_rad,
                  int nDirs,
                  /* Output Arguments */
                  float_complex* Y);
    
/*
 * Function: complex2realSHMtx
 * ---------------------------
 * Returns the unitary transformation matrix T_c2r. It expresses the real
 * spherical harmonics with respect to the complex ones, so that
 * r_N = T_c2r * y_N, where r_N and y_N is are the real and complex SH vectors,
 * respectively.
 
 * Input Arguments:
 *     order - order of spherical harmonic expansion
 * Output Arguments:
 *     T_c2r - transformation matrix for complex->real;
 *             FLAT: (order+1)^2 x (order+1)^2
 */
void complex2realSHMtx(/* Input Arguments */
                       int order,
                       /* Output Arguments */
                       float_complex* T_c2r);
    
/*
 * Function: real2complexSHMtx
 * ---------------------------
 * Returns the unitary transformation matrix T_r2c the expresses the complex
 * spherical harmonics with respect to the real ones, so that y_N = T_r2c * r_N,
 * where r_N and y_N are the real and complex SH vectors, respectively.
 
 * Input Arguments:
 *     order - order of spherical harmonic expansion
 * Output Arguments:
 *     T_c2r - transformation matrix for real->complex;
 *             FLAT: (order+1)^2 x (order+1)^2
 */
void real2complexSHMtx(/* Input Arguments */
                       int order,
                       /* Output Arguments */
                       float_complex* T_r2c);
    
/*
 * Function: complex2realCoeffs
 * ----------------------------
 * Converts SH coeffs from the complex to real basis
 *
 * Input Arguments:
 *     order - order of spherical harmonic expansion
 *     C_N   - complex coeffients; FLAT: (order+1)^2 x K
 *     K     - number of columns
 * Output Arguments:
 *     R_N  - real coefficients; FLAT: (order+1)^2 x K
 */
void complex2realCoeffs(/* Input Arguments */
                        int order,
                        float_complex* C_N,
                        int K,
                        /* Output Arguments */
                        float* R_N);
    
/*
 * Function: getSHrotMtxReal
 * -------------------------
 * Generates a real-valued spherical harmonic rotation matrix [1]
 * (assumes ACN channel ordering convention).
 * Note: the normalisation convention does not matter, as e.g. only dipoles
 * are used to rotated dipoles, quadrapoles to rotate quadrapoles etc.
 *
 * Input Arguments:
 *     R      - zyx rotation matrix; 3 x 3
 *     L      - order of spherical harmonic expansion
 * Output Arguments:
 *     RotMtx - SH domain rotation matrix; FLAT: (L+1)^2 x (L+1)^2
 *
 * [1] Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical
 *     Harmonics. Direct Determination by Recursion Page: Additions and
 *     Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100.
 */
void getSHrotMtxReal(float R[3][3],
                     float* RotMtx,
                     int L);

/*
 * Function: computeVelCoeffsMtx
 * -----------------------------
 * Computes the matrices that generate the coefficients of the beampattern of
 * order (sectorOrder+1) that is essentially the product of a pattern of
 * order=sectorOrder and a dipole. It is used in "beamWeightsVelocityPatterns".
 * For the derivation of the matrices see [1]
 *
 * Input Arguments:
 *     sectorOrder - order of patterns
 * Output Arguments:
 *     A_xyz       - Velocity coefficients;
 *                   FLAT: (sectorOrder+2)^2  x (sectorOrder+1)^2 x 3
 *
 * [1] Politis, A. and Pulkki, V., 2016. Acoustic intensity, energy-density and
 *     diffuseness estimation in a directionally-constrained region.  arXiv
 *     preprint arXiv:1609.03409
 */ 
void computeVelCoeffsMtx(/* Input Arguments */
                         int sectorOrder,
                         /* Output Arguments */
                         float_complex* A_xyz);
    
/*
 * Function: computeSectorCoeffsEP
 * -------------------------------
 * Computes the beamforming matrices of sector and velocity coefficients for
 * ENERGY-preserving (EP) sectors for real SH.
 * This partitioning of the sound-field into spatially-localised sectors has
 * been used for parametric reproduction in [1] and also for sound-field
 * visualision [2].
 *
 * Input Arguments:
 *     sectorOrder  - order of sector patterns
 *     A_xyz        - Velocity coefficients (see "computeVelCoeffsMtx");
 *                    FLAT: (sectorOrder+2)^2 x (sectorOrder+1)^2 x 3
 *     pattern      - see "SECTOR_PATTERNS" enum for the options
 *     sec_dirs_deg - sector directions [azi elev], in DEGREES;
 *                    FLAT: nSecDirs x 2
 *     nSecDirs     - number of sectors
 * Output Arguments:
 *     sectorCoeffs - the sector coefficients;
 *                    FLAT: (nSecDirs*4) x (orderSec+2)^2
 *
 * [1] Politis, A., Vilkamo, J., & Pulkki, V. (2015). Sector-based parametric
 *     sound field reproduction in the spherical harmonic domain. IEEE Journal
 *     of Selected Topics in Signal Processing, 9(5), 852-866.
 * [2] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of angular
 *     spectra based on a directional re-assignment approach for ambisonic
 *     sound-field visualisation". IEEE International Conference on Acoustics,
 *     Speech and Signal Processing (ICASSP).
 */
void computeSectorCoeffsEP(/* Input Arguments */
                           int orderSec,
                           float_complex*  A_xyz,
                           SECTOR_PATTERNS pattern,
                           float* sec_dirs_deg,
                           int nSecDirs,
                           /* Output Arguments */
                           float* sectorCoeffs);
    
/*
 * Function: computeSectorCoeffsAP
 * -------------------------------
 * Computes the beamforming matrices of sector and velocity coefficients for
 * AMPLITUDE-preserving (AP) sectors for real SH.
 * This partitioning of the sound-field into spatially-localised sectors has
 * been used for parametric reproduction in [1] and also for sound-field
 * visualision [2].
 *
 * Input Arguments:
 *     sectorOrder  - order of sector patterns
 *     A_xyz        - Velocity coefficients (see "computeVelCoeffsMtx");
 *                    FLAT: (sectorOrder+2)^2 x (sectorOrder+1)^2 x 3
 *     pattern      - see "SECTOR_PATTERNS" enum for the options
 *     sec_dirs_deg - sector directions [azi elev], in DEGREES;
 *                    FLAT: nSecDirs x 2
 *     nSecDirs     - number of sectors
 * Output Arguments:
 *     sectorCoeffs - the sector coefficients;
 *                    FLAT: (nSecDirs*4) x (orderSec+2)^2
 *
 * [1] Politis, A., Vilkamo, J., & Pulkki, V. (2015). Sector-based parametric
 *     sound field reproduction in the spherical harmonic domain. IEEE Journal
 *     of Selected Topics in Signal Processing, 9(5), 852-866.
 * [2] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of angular
 *     spectra based on a directional re-assignment approach for ambisonic
 *     sound-field visualisation". IEEE International Conference on Acoustics,
 *     Speech and Signal Processing (ICASSP).
 */
void computeSectorCoeffsAP(/* Input Arguments */
                           int orderSec,
                           float_complex* A_xyz,
                           SECTOR_PATTERNS pattern,
                           float* sec_dirs_deg,
                           int nSecDirs,
                           /* Output Arguments */
                           float* sectorCoeffs);
    
/*
 * Function: beamWeightsCardioid2Spherical
 * ---------------------------------------
 * Generate spherical coefficients for cardioids. For a specific order N of a
 * higher order cardioid of the form D(theta)=(1/2)^N * (1+cos(theta))^N,
 * generate the beamweights for the same pattern in the SHD. Because the pattern
 * is axisymmetric only the N+1 coefficients of m=0 are returned.
 *
 * Input Arguments:
 *     N   - order of spherical harmonic expansion
 * Output Arguments:
 *     b_n - beamformer weights; (N+1) x 1
 */
void beamWeightsCardioid2Spherical(/* Input Arguments */
                                   int N,
                                   /* Output Arguments */
                                   float* b_n);
    
    /**** NOT IMPLEMENTED YET ****/
/*
 * Function: beamWeightsDolphChebyshev2Spherical
 * ---------------------------------------------
 * Generate beamweights in the SHD for Dolph-Chebyshev beampatterns, with
 * mainlobe and sidelobe control [1]. Because the pattern is axisymmetric only
 * the N+1 coefficients of m=0 are returned.
 *
 * Input Arguments:
 *     N          - order of spherical harmonic expansion
 *     paramType  - 0: side-lobe level control, 1: mainlobe width control
 *     arrayParam - sidelobe level 1/R or mainlobe with 2*a0
 * Output Arguments:
 *     b_n        - beamformer weights; (N+1) x 1
 *
 * [1] Koretz, A. and Rafaely, B., 2009. Dolph-Chebyshev beampattern design for
 *     spherical arrays. IEEE Transactions on Signal Processing, 57(6),
 *     pp.2417-2420.
 */
void beamWeightsDolphChebyshev2Spherical(/* Input Arguments */
                                         int N,
                                         int paramType,
                                         float arrayParam,
                                         /* Output Arguments */
                                         float* b_n);
    
/*
 * Function: beamWeightsHypercardioid2Spherical
 * --------------------------------------------
 * The hypercardioid is the pattern that maximises the directivity-factor for a
 * certain SH order N. The hypercardioid is also the plane-wave decomposition
 * beamformer in the SHD, also called 'regular' because the beamweights are just
 * the SH values on the beam-direction. Since the pattern is axisymmetric only
 * the N+1 coefficients of m=0 are returned.
 *
 * Input Arguments:
 *     N  - order of spherical harmonic expansion
 * Output Arguments:
 *     b_n - beamformer weights; (N+1) x 1
 */
void beamWeightsHypercardioid2Spherical(/* Input Arguments */
                                        int N,
                                        /* Output Arguments */
                                        float* b_n);

/*
 * Function: beamWeightsMaxEV
 * --------------------------
 * Generate the beamweights for the a maximum energy-vector beampattern in the
 * SHD. This pattern originates from ambisonic-related research and it maximises
 * the ambisonic energy-vector, which is essentially the directional centroid of
 * the squared pattern. IT can also be seen as the pattern that maximizes the
 * acoustic intensity vector of a diffuse field weighted with this pattern.
 * In practice it is almost the same as a supercardioid that maximizes front-
 * back power ratio for a certain order, and it can be used as such. Because
 * the pattern is axisymmetric only the N+1 coefficients of m=0 are returned.
 * Details for their theory can be found e.g. in [1].
 *
 * Input Arguments:
 *     N  - order of spherical harmonic expansion
 * Output Arguments:
 *     b_n - beamformer weights; (N+1) x 1
 *
 * [1] Zotter, F., Pomberger, H. and Noisternig, M., 2012. Energy-preserving
 *     ambisonic decoding. Acta Acustica united with Acustica, 98(1), pp.37-47.
 */
void beamWeightsMaxEV(/* Input Arguments */
                      int N,
                      /* Output Arguments */
                      float* b_n);
    
/*
 * Function: beamWeightsVelocityPatternsReal
 * -----------------------------------------
 * If the sound-field is weighted with an axisymmetric spatial distribution
 * described by the N+1 SH coefficients b_n, then the beamweights capturing the
 * velocity signals for the weighted sound-field are of an order one higher
 * than the weighting pattern, and can be derived from it. This type of
 * beamforming has some applications for spatial sound reproduction and
 * acoustic analysis, see [1].
 *
 * Input Arguments:
 *     order     - order of spherical harmonic expansion
 *     b_n       - axisymmetric beamformer weights; (order+1) x 1
 *     azi_rad   - orientation, azimuth in RADIANS
 *     elev_rad  - orientation, ELEVATION in RADIANS
 *     A_xyz     - Velocity coefficients (see "computeVelCoeffsMtx");
 *                 FLAT: (sectorOrder+2)^2 x (sectorOrder+1)^2 x 3
 * Output Arguments:
 *     velCoeffs - beamforming coefficients for velocity patterns;
 *                 FLAT: (order+2)^2 x 3
 *
 * [1] Politis, A. and Pulkki, V., 2016. Acoustic intensity, energy-density and
 *     diffuseness estimation in a directionally-constrained region. arXiv
 *     preprint arXiv:1609.03409.
 */
void beamWeightsVelocityPatternsReal(/* Input Arguments */
                                     int order,
                                     float* b_n,
                                     float azi_rad,
                                     float elev_rad,
                                     float_complex* A_xyz,
                                     /* Output Arguments */
                                     float* velCoeffs);

/*
 * Function: beamWeightsVelocityPatternsComplex
 * --------------------------------------------
 * If the sound-field is weighted with an axisymmetric spatial distribution
 * described by the N+1 SH coefficients b_n, then the beamweights capturing the
 * velocity signals for the weighted sound-field are of an order one higher
 * than the weighting pattern, and can be derived from it. This type of
 * beamforming has some applications for spatial sound reproduction and
 * acoustic analysis, see [1].
 *
 * Input Arguments:
 *     order     - order of spherical harmonic expansion
 *     b_n       - axisymmetric beamformer weights; (order+1) x 1
 *     azi_rad   - orientation, azimuth in RADIANS
 *     elev_rad  - orientation, ELEVATION in RADIANS
 *     A_xyz     - Velocity coefficients (see "computeVelCoeffsMtx");
 *                 FLAT: (sectorOrder+2)^2 x (sectorOrder+1)^2 x 3
 * Output Arguments:
 *     velCoeffs - beamforming coefficients for velocity patterns;
 *                 FLAT: (order+2)^2 x 3
 *
 * [1] Politis, A. and Pulkki, V., 2016. Acoustic intensity, energy-density and
 *     diffuseness estimation in a directionally-constrained region. arXiv
 *     preprint arXiv:1609.03409.
 */
void beamWeightsVelocityPatternsComplex(/* Input Arguments */
                                        int order,
                                        float* b_n,
                                        float azi_rad,
                                        float elev_rad,
                                        float_complex* A_xyz,
                                        /* Output Arguments */
                                        float_complex* velCoeffs);
    
/*
 * Function: rotateAxisCoeffsReal
 * ------------------------------
 * Returns spherical coefficients for a rotated axisymmetric pattern.
 *
 * Input Arguments:
 *     order   - order of spherical harmonic expansion
 *     c_n     - coefficients describing a rotationally symmetric pattern order
 *               N, expressed as a sum of spherical harmonics of degree m=0;
 *               (N+1) x 1
 *     theta_0 - POLAR rotation for the pattern, in RADIANS
 *     phi_0   - azimuthal rotation for the pattern, in RADIANS
 * Output Arguments:
 *     c_nm    - coefficients of rotated pattern expressed as a sum of SHs;
 *               (N+1)^2 x 1
 */
void rotateAxisCoeffsReal(/* Input arguments */
                          int order,
                          float* c_n,
                          float theta_0,
                          float phi_0,
                          /* Output arguments */
                          float* c_nm);          
    
/*
 * Function: rotateAxisCoeffsComplex
 * ---------------------------------
 * Returns spherical coefficients for a rotated axisymmetric pattern.
 *
 * Input Arguments:
 *     order   - order of spherical harmonic expansion
 *     c_n     - coefficients describing a rotationally symmetric pattern order
 *               N, expressed as a sum of spherical harmonics of degree m=0;
 *               (N+1) x 1
 *     theta_0 - POLAR rotation for the pattern, in RADIANS
 *     phi_0   - azimuthal rotation for the pattern, in RADIANS
 * Output Arguments:
 *     c_nm    - coefficients of rotated pattern expressed as a sum of SHs;
 *               (N+1)^2 x 1
 */
void rotateAxisCoeffsComplex(/* Input arguments */
                             int order,
                             float* c_n,
                             float theta_0,
                             float phi_0,
                             /* Output arguments */
                             float_complex* c_nm);

/*
 * Function: checkCondNumberSHTReal
 * --------------------------------
 * Computes the condition numbers for a least-squares SHT
 *
 * Input Arguments:
 *     order    - order of spherical harmonic expansion
 *     dirs_rad - directions on the sphere [azi, INCLINATION] convention, in
 *                RADIANS; FLAT: nDirs x 2
 *     nDirs    - number of directions
 *     w        - integration weights; nDirs x 1
 * Output Arguments:
 *     cond_N   - condition numbers; (order+1) x 1
 */
void checkCondNumberSHTReal(/* Input arguments */
                            int order,
                            float* dirs_rad,
                            int nDirs,
                            float* w,
                            /* Output arguments */
                            float* cond_N);


/* ========================================================================== */
/*                     Localisation Functions in the  SHD                     */
/* ========================================================================== */

/*
 * Function: generatePWDmap
 * ------------------------
 * Generates a powermap based on the energy of plane-wave decomposition (PWD)/
 * hyper-cardioid beamformers.
 *
 * Input Arguments:
 *     order      - analysis order
 *     Cx         - correlation/covarience matrix;
 *                  FLAT: (order+1)^2 x (order+1)^2
 *     Y_grid     - steering vectors for each grid direcionts;
 *                  FLAT: (order+1)^2 x nGrid_dirs
 *     nGrid_dirs - number of grid directions
 * Output Arguments:
 *     pmap       - resulting PWD powermap; nGrid_dirs x 1
 */
void generatePWDmap(/* Input arguments */
                    int order,
                    float_complex* Cx,
                    float_complex* Y_grid,
                    int nGrid_dirs,
                    /* Output arguments */
                    float* pmap);
/*
 * Function: generateMVDRmap
 * -------------------------
 * Generates a powermap based on the energy of adaptive minimum variance
 * distortionless response (MVDR) beamformers.
 *
 * Input Arguments:
 *     order      - analysis order
 *     Cx         - correlation/covarience matrix;
 *                  FLAT: (order+1)^2 x (order+1)^2
 *     Y_grid     - steering vectors for each grid direcionts;
 *                  FLAT: (order+1)^2 x nGrid_dirs
 *     nGrid_dirs - number of grid directions
 *     regPar     - regularisation parameter, for diagonal loading of Cx
 * Output Arguments:
 *     pmap       - resulting MVDR powermap; nGrid_dirs x 1
 *     w_MVDR     - optional. weights will be copied to this, unless it's NULL;
 *                  FLAT: nSH x nGrid_dirs || NULL
 */
void generateMVDRmap(/* Input arguments */
                     int order,
                     float_complex* Cx,
                     float_complex* Y_grid,
                     int nGrid_dirs,
                     float regPar,
                     /* Output arguments */
                     float* pmap,
                     float_complex* w_MVDR);

/*
 * Function: generateCroPaCLCMVmap
 * -------------------------------
 * EXPERIMENTAL! Generates a powermap utilising the CroPaC LCMV post-filter
 * described in [1].
 * The spatial post-filter is estimated for all directions on the grid, and is
 * used to supress reverb/noise interference that may be present in an MVDR map.
 * Unlike in the paper, the second column for the contraints 'A', is
 * Y.*diag(Cx), rather than utilising a maximum energy beamformer. The post-
 * filters are then applied to the MVDR powermap map derived in the sherical
 * harmonic domain, rather than an MVDR beamformer generated directly in the
 * microphone array signal domain, like in the paper. Otherwise, the algorithm
 * is the same.
 *
 * Input Arguments:
 *     order      - analysis order
 *     Cx         - correlation/covarience matrix;
 *                  FLAT: (order+1)^2 x (order+1)^2
 *     Y_grid     - steering vectors for each grid direcionts;
 *                  FLAT: (order+1)^2 x nGrid_dirs
 *     nGrid_dirs - number of grid directions
 *     regPar     - regularisation parameter, for diagonal loading of Cx
 *     lambda     - parameter controlling how harsh CroPaC is applied, 0..1;
 *                  0: fully cropac, 1: fully mvdr
 * Output Arguments:
 *     pmap       - resulting CroPaC LCMV powermap; nGrid_dirs x 1
 *
 * [1] Delikaris-Manias, S., Vilkamo, J., & Pulkki, V. (2016). Signal-dependent
 *     spatial filtering based on weighted-orthogonal beamformers in the
 *     spherical harmonic domain. IEEE/ACM Transactions on Audio, Speech and
 *     Language Processing (TASLP), 24(9), 1507-1519.
 */
void generateCroPaCLCMVmap(/* Input arguments */
                           int order,
                           float_complex* Cx,
                           float_complex* Y_grid,
                           int nGrid_dirs,
                           float regPar,
                           float lambda,
                           /* Output arguments */
                           float* pmap);
    
/*
 * Function: generateMUSICmap
 * --------------------------
 * Generates an activity-map based on the sub-space multiple-signal
 * classification (MUSIC) method.
 *
 * Input Arguments:
 *     order        - analysis order
 *     Cx           - correlation/covarience matrix;
 *                    FLAT: (order+1)^2 x (order+1)^2
 *     Y_grid       - steering vectors for each grid direcionts;
 *                    FLAT: (order+1)^2 x nGrid_dirs
 *     nSources     - number of sources present in sound scene
 *     nGrid_dirs   - number of grid directions
 *     logScaleFlag - 1: log(pmap), 0: pmap.
 * Output Arguments:
 *     pmap         - resulting MUSIC pseudo-spectrum; nGrid_dirs x 1
 */
void generateMUSICmap(/* Input arguments */
                      int order,
                      float_complex* Cx,
                      float_complex* Y_grid,
                      int nSources,
                      int nGrid_dirs,
                      int logScaleFlag,
                      /* Output arguments */
                      float* pmap);

/*
 * Function: generateMinNormMap
 * ----------------------------
 * Generates an activity-map based on the sub-space minimum-norm (MinNorm)
 * method.
 *
 * Input Arguments:
 *     order        - analysis order
 *     Cx           - correlation/covarience matrix;
 *                    FLAT: (order+1)^2 x (order+1)^2
 *     Y_grid       - steering vectors for each grid direcionts;
 *                    FLAT: (order+1)^2 x nGrid_dirs
 *     nSources     - number of sources present in sound scene
 *     nGrid_dirs   - number of grid directions
 *     logScaleFlag - 1: log(pmap), 0: pmap.
 * Output Arguments:
 *     pmap         - resulting MinNorm pseudo-spectrum; nGrid_dirs x 1
 */
void generateMinNormMap(/* Input arguments */
                        int order,
                        float_complex* Cx,
                        float_complex* Y_grid,
                        int nSources,
                        int nGrid_dirs,
                        int logScaleFlag,
                        /* Output arguments */
                        float* pmap);


/* ========================================================================== */
/*                   Cylindrical/Spherical Bessel Functions                   */
/* ========================================================================== */

/*
 * Function: bessel_Jn
 * -------------------
 * Computes the (cylindrical) Bessel function of the first kind: Jn
 * returns the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * Input Arguments:
 *     N    - function order (highest is ~30 given numerical precision)
 *     z    - input values; nZ x 1
 *     nZ   - number of input values
 * Output Arguments:
 *     J_n  - Bessel values (set as NULL if not required); FLAT: nZ x (N+1)
 *     dJ_n - Bessel derivative values (set as NULL if not required);
 *            FLAT: nZ x (N+1)
 */
void bessel_Jn(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               double* J_n,
               double* dJ_n);

/*
 * Function: bessel_Yn
 * -------------------
 * Computes the (cylindrical) Bessel function of the second kind: Yn
 * returns the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * Input Arguments:
 *     N    - function order (highest is ~30 given numerical precision)
 *     z    - input values; nZ x 1
 *     nZ   - number of input values
 * Output Arguments:
 *     Y_n  - Bessel values (set as NULL if not required); FLAT: nZ x (N+1)
 *     dY_n - Bessel derivative values (set as NULL if not required);
 *            FLAT: nZ x (N+1)
 */
void bessel_Yn(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               double* Y_n,
               double* dY_n);

/*
 * Function: hankel_Hn1
 * --------------------
 * Computes the (cylindrical) Hankel function of the first kind: Hn1
 * returns the Hankel values and their derivatives up to order N for all values
 * in vector z
 *
 * Input Arguments:
 *     N     - function order (highest is ~30 given numerical precision)
 *     z     - input values; nZ x 1
 *     nZ    - number of input values
 * Output Arguments:
 *     Hn1_n  - Hankel values (set as NULL if not required); FLAT: nZ x (N+1)
 *     dHn1_n - Hankel derivative values (set as NULL if not required);
 *              FLAT: nZ x (N+1)
 */
void hankel_Hn1(/* Input arguments */
                int N,
                double* z,
                int nZ,
                /* Output arguments */
                double_complex* Hn1_n,
                double_complex* dHn1_n);

/*
 * Function: hankel_Hn2
 * --------------------
 * Computes the (cylindrical) Hankel function of the second kind: Hn2
 * returns the Hankel values and their derivatives up to order N for all values
 * in vector z
 *
 * Input Arguments:
 *     N     - function order (highest is ~30 given numerical precision)
 *     z     - input values; nZ x 1
 *     nZ    - number of input values
 * Output Arguments:
 *     Hn2_n  - Hankel values (set as NULL if not required); FLAT: nZ x (N+1)
 *     dHn2_n - Hankel derivative values (set as NULL if not required);
 *              FLAT: nZ x (N+1)
 */
void hankel_Hn2(/* Input arguments */
                int N,
                double* z,
                int nZ,
                /* Output arguments */
                double_complex* Hn2_n,
                double_complex* dHn2_n);

/*
 * Function: bessel_jn
 * -------------------
 * Computes the spherical Bessel function of the first kind: jn
 * returns the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * Input Arguments:
 *     N    - function order (highest is ~30 given numerical precision)
 *     z    - input values; nZ x 1
 *     nZ   - number of input values
 * Output Arguments:
 *     maxN - & maximum function order that could be computed <=N
 *     j_n  - Bessel values (set as NULL if not required); FLAT: nZ x (N+1)
 *     dj_n - Bessel derivative values (set as NULL if not required);
 *            FLAT: nZ x (N+1)
 */
void bessel_jn(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               int* maxN,
               double* j_n,
               double* dj_n);
    
/*
 * Function: bessel_in
 * -------------------
 * Computes the modified spherical Bessel function of the first kind: in
 * returns the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * Input Arguments:
 *     N    - function order (highest is ~30 given numerical precision)
 *     z    - input values; nZ x 1
 *     nZ   - number of input values
 * Output Arguments:
 *     maxN - & maximum function order that could be computed <=N
 *     i_n  - Bessel values (set as NULL if not required); FLAT: nZ x (N+1)
 *     di_n - Bessel derivative values (set as NULL if not required);
 *            FLAT: nZ x (N+1)
 */
void bessel_in(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               int* maxN,
               double* i_n,
               double* di_n);

/*
 * Function: bessel_yn
 * -------------------
 * Computes the spherical Bessel function of the second kind (Neumann): yn
 * returns the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * Input Arguments:
 *     N    - function order (highest is ~30 given numerical precision)
 *     z    - input values; nZ x 1
 *     nZ   - number of input values
 * Output Arguments:
 *     maxN - & maximum function order that could be computed <=N
 *     y_n  - Bessel values (set as NULL if not required); FLAT: nZ x (N+1)
 *     dy_n - Bessel derivative values (set as NULL if not required);
 *            FLAT: nZ x (N+1)
 */
void bessel_yn(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               int* maxN,
               double* y_n,
               double* dy_n);
    
/*
 * Function: bessel_kn
 * -------------------
 * Computes the modified spherical Bessel function of the second kind: kn
 * returns the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * Input Arguments:
 *     N    - function order (highest is ~30 given numerical precision)
 *     z    - input values; nZ x 1
 *     nZ   - number of input values
 * Output Arguments:
 *     maxN - & maximum function order that could be computed <=N
 *     k_n  - Bessel values (set as NULL if not required); FLAT: nZ x (N+1)
 *     dk_n - Bessel derivative values (set as NULL if not required);
 *            FLAT: nZ x (N+1)
 */
void bessel_kn(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               int* maxN,
               double* k_n,
               double* dk_n);

/*
 * Function: hankel_hn1
 * --------------------
 * Computes the spherical Hankel function of the first kind: hn1
 * returns the Hankel values and their derivatives up to order N for all values
 * in vector z
 *
 * Input Arguments:
 *     N     - function order (highest is ~30 given numerical precision)
 *     z     - input values; nZ x 1
 *     nZ    - number of input values
 * Output Arguments:
 *     maxN  - & maximum function order that could be computed <=N
 *     h_n1  - Hankel values (set as NULL if not required); FLAT: nZ x (N+1)
 *     dh_n1 - Hankel derivative values (set as NULL if not required);
 *             FLAT: nZ x (N+1)
 */
void hankel_hn1(/* Input arguments */
                int N,
                double* z,
                int nZ,
                /* Output arguments */
                int* maxN,
                double_complex* h_n1,
                double_complex* dh_n1);

/*
 * Function: hankel_hn2
 * --------------------
 * Computes the spherical Hankel function of the second kind: hn2
 * returns the Hankel values and their derivatives up to order N for all values
 * in vector z
 *
 * Input Arguments:
 *     N     - function order (highest is ~30 given numerical precision)
 *     z     - input values; nZ x 1
 *     nZ    - number of input values
 * Output Arguments:
 *     maxN  - & maximum function order that could be computed <=N
 *     h_n2  - Hankel values (set as NULL if not required); FLAT: nZ x (N+1)
 *     dh_n2 - Hankel derivative values (set as NULL if not required);
 *             FLAT: nZ x (N+1)
 */
void hankel_hn2(/* Input arguments */
                int N,
                double* z,
                int nZ,
                /* Output arguments */
                int* maxN,
                double_complex* h_n2,
                double_complex* dh_n2);
    
    
/* ========================================================================== */
/*              Microphone/Hydrophone array processing functions              */
/* ========================================================================== */

/*
 * Function: cylModalCoeffs
 * ------------------------
 * Calculates the modal coefficients for open/rigid cylindrical arrays
 *
 * Input Arguments:
 *     order     - max order (highest is ~30 given numerical precision)
 *     kr        - wavenumber*radius; nBands x 1
 *     nBands    - number of frequency bands/bins
 *     arrayType - see 'ARRAY_CONSTRUCTION_TYPES' enum
 * Output Arguments:
 *     b_N       - modal coefficients per kr and 0:order;
 *                 FLAT: nBands x (order+1)
 */
void cylModalCoeffs(/* Input arguments */
                    int order,
                    double* kr,
                    int nBands,
                    ARRAY_CONSTRUCTION_TYPES arrayType,
                    /* Output arguments */
                    double_complex* b_N);

/*
 * Function: sphArrayAliasLim
 * --------------------------
 * Returns a simple estimate of the spatial aliasing limit (the kR = maxN rule)
 *
 * Input Arguments:
 *     r    - array radius, meters
 *     c    - speed of sound, m/s
 *     maxN - order
 * Returns:
 *     spatial aliasing limit estimate
 */
float sphArrayAliasLim(/* Input arguments */
                       float r,
                       float c,
                       int maxN);

/*
 * Function: sphArrayNoiseThreshold
 * --------------------------------
 * Computes the frequencies that the noise in the output channels of a spherical
 * microphone array (SMA), after performing the spherical harmonic transform
 * (SHT) and equalisation of the output signals, reaches a certain user-defined
 * threshold maxG_db. The frequencies are computed only at the lower range of
 * each order, where its response decays rapidly, ignoring for example the nulls
 * of an open array at the higher frequencies. The estimation of the limits are
 * based on a linear approximation of the log-log response found e.g. in [1]
 *
 * Input Arguments:
 *     maxN      - maximum order of the array
 *     Nsensors  - number of sensors
 *     r         - mic radius, meters
 *     c         - speed of sound, m/s
 *     arrayType - see 'ARRAY_CONSTRUCTION_TYPES' enum
 *     dirCoeff  - only for directional (open) arrays, 1: omni, 0.5: card,
 *                 0:dipole
 *     maxG_db   - max allowed amplification for the noise level,
 *                 maxG_db = 20*log10(maxG)
 * Output Arguments:
 *     f_lim     - spatial aliasing limit estimate
 *
 * [1] Sector-based Parametric Sound Field Reproduction in the Spherical
 *     Harmonic Domain A Politis, J Vilkamo, V Pulkki. 2015. IEEE Journal of
 *     Selected Topics in Signal Processing 9 (5), 852 - 866
 */
void sphArrayNoiseThreshold(/* Input arguments */
                            int maxN,
                            int Nsensors,
                            float r,
                            float c,
                            ARRAY_CONSTRUCTION_TYPES arrayType,
                            double dirCoeff,
                            float maxG_db,
                            /* Output arguments */
                            float* f_lim);

/*
 * Function: sphModalCoeffs
 * ------------------------
 * Calculates the modal coefficients for open/rigid spherical arrays
 *
 * Input Arguments:
 *     order     - max order (highest is ~30 given numerical precision)
 *     kr        - wavenumber*radius; nBands x 1
 *     nBands    - number of frequency bands/bins
 *     arrayType - see 'ARRAY_CONSTRUCTION_TYPES' enum
 *     dirCoeff  - only for directional (open) arrays, 1: omni, 0.5: card,
 *                 0:dipole
 * Output Arguments:
 *     b_N       - modal coefficients per kr and 0:order;
 *                 FLAT: nBands x (order+1)
 */
void sphModalCoeffs(/* Input arguments */
                    int order,
                    double* kr,
                    int nBands,
                    ARRAY_CONSTRUCTION_TYPES arrayType,
                    double dirCoeff,
                    /* Output arguments */
                    double_complex* b_N);

/*
 * Function: sphScattererModalCoeffs
 * ---------------------------------
 * Calculates the modal coefficients for a rigid spherical scatterer with omni-
 * directional sensors (Assumes all sensors are placed the same distance from
 * the scatterer, w.r.t. the origin)
 *
 * Input Arguments:
 *     order  - max order (highest is ~30 given numerical precision)
 *     kr     - wavenumber*array_radius; nBands x 1
 *     kR     - wavenumber*scatterer_radius; nBands x 1
 *     nBands - number of frequency bands/bins
 * Output Arguments:
 *     b_N    - modal coefficients per kr and 0:order; FLAT: nBands x (order+1)
 */
void sphScattererModalCoeffs(/* Input arguments */
                             int order,
                             double* kr,
                             double* kR,
                             int nBands,
                             /* Output arguments */
                             double_complex* b_N);

/*
 * Function: sphScattererDirModalCoeffs
 * ------------------------------------
 * Calculates the modal coefficients for a rigid spherical scatterer with
 * directional sensors (Assumes all sensors are placed the same distance from
 * the scatterer, w.r.t. the origin)
 *
 * Input Arguments:
 *     order    - max order (highest is ~30 given numerical precision)
 *     kr       - wavenumber*array_radius; nBands x 1
 *     kR       - wavenumber*scatterer_radius; nBands x 1
 *     nBands   - number of frequency bands/bins
 *     dirCoeff - directivity coefficient, 1: omni, 0.5: card, 0:dipole
 * Output Arguments:
 *     b_N      - modal coefficients per kr and 0:order; FLAT: nBands x (order+1)
 */
void sphScattererDirModalCoeffs(/* Input arguments */
                                int order,
                                double* kr,
                                double* kR,
                                int nBands,
                                double dirCoeff,
                                /* Output arguments */
                                double_complex* b_N);
    
/*
 * Function: sphDiffCohMtxTheory
 * -----------------------------
 * Calculates the theoretical diffuse coherence matrix for a spherical array
 *
 * Input Arguments:
 *     order           - max order (highest is ~30 given numerical precision)
 *     sensor_dirs_rad - spherical coords of the sensors in RADIANS, [azi ELEV];
 *                       FLAT: N_sensors x 2
 *     N_sensors       - number of sensors
 *     arrayType       - see 'ARRAY_CONSTRUCTION_TYPES' enum
 *     dirCoeff        - only for directional (open) arrays, 1: omni, 0.5: card,
 *                       0:dipole
 *     kr              - wavenumber*sensor_radius; nBands x 1
 *     kR              - wavenumber*scatterer_radius, set to NULL if not
 *                       applicable; nBands x 1
 *     nBands          - number of frequency bands/bins
 * Output Arguments:
 *     M_diffcoh       - theoretical diffuse coherence matrix per frequency;
 *                       FLAT: N_sensors x N_sensors x nBands
 */
void sphDiffCohMtxTheory(/* Input arguments */
                         int order,
                         float* sensor_dirs_rad,
                         int N_sensors,
                         ARRAY_CONSTRUCTION_TYPES arrayType,
                         double dirCoeff,
                         double* kr,
                         double* kR,
                         int nBands,
                         /* Output arguments */
                         double* M_diffcoh);

/*
 * Function: simulateCylArray
 * --------------------------
 * Simulates a cylindrical microphone array, returning the transfer functions
 * for each (plane wave) source direction on the surface of the cylinder
 *
 * Input Arguments:
 *     order           - max order (highest is ~30 given numerical precision)
 *     kr              - wavenumber*radius; nBands x 1
 *     nBands          - number of frequency bands/bins
 *     sensor_dirs_rad - spherical coords of the sensors in RADIANS, [azi ELEV];
 *                       FLAT: N_sensors x 2
 *     N_sensors       - number of sensors
 *     src_dirs_deg    - spherical coords of the plane waves in DEGREES,
 *                       [azi ELEV]; FLAT: N_srcs x 2
 *     N_srcs          - number sources (DoAs of plane waves)
 *     arrayType       - see 'ARRAY_CONSTRUCTION_TYPES' enum
 * Output Arguments:
 *     H_array         - simulated array response for each plane wave;
 *                       FLAT: nBands x N_sensors x N_srcs
 */
void simulateCylArray(/* Input arguments */
                      int order,
                      double* kr,
                      int nBands,
                      float* sensor_dirs_rad,
                      int N_sensors,
                      float* src_dirs_deg,
                      int N_srcs,
                      ARRAY_CONSTRUCTION_TYPES arrayType,
                      /* Output arguments */
                      float_complex* H_array);

/*
 * Function: simulateSphArray
 * --------------------------
 * Simulates a spherical microphone array, returning the transfer functions
 * for each (plane wave) source direction on the surface of the sphere.
 *
 * Input Arguments:
 *     order           - max order (highest is ~30 given numerical precision)
 *     kr              - wavenumber*array_radius; nBands x 1
 *     kR              - wavenumber*scatterer_radius, set to NULL if not needed
 *     nBands          - number of frequency bands/bins
 *     sensor_dirs_rad - spherical coords of the sensors in RADIANS, [azi ELEV];
 *                       FLAT: N_sensors x 2
 *     N_sensors       - number of sensors
 *     src_dirs_deg    - spherical coords of the plane waves in DEGREES,
 *                       [azi ELEV]; FLAT: N_srcs x 2
 *     N_srcs          - number sources (DoAs of plane waves)
 *     arrayType       - see 'ARRAY_CONSTRUCTION_TYPES' enum
 *     dirCoeff        - only for directional (open) arrays, 1: omni, 0.5: card,
 *                       0:dipole
 * Output Arguments:
 *     H_array         - simulated array response for each plane wave;
 *                       FLAT: nBands x N_sensors x N_srcs
 */
void simulateSphArray(/* Input arguments */
                      int order,
                      double* kr,
                      double* kR,
                      int nBands,
                      float* sensor_dirs_rad,
                      int N_sensors,
                      float* src_dirs_deg,
                      int N_srcs,
                      ARRAY_CONSTRUCTION_TYPES arrayType,
                      double dirCoeff,
                      /* Output arguments */
                      float_complex* H_array);

/*
 * Function: evaluateSHTfilters
 * ----------------------------
 * generates some objective measures, which evaluate the performance of the
 * spatial encoding filters. This analysis is performed by comparing the spatial
 * resolution of the spherical harmonic components generated by the encoding
 * filters, with the ideal SH components. For more information, the reader is
 * directed to [1,2].
 *
 * Input Arguments:
 *     order      - transform/encoding order
 *     M_array2SH - encoding matrix per frequency;
 *                  FLAT: nBands x (order+1)^2 x nSensors
 *     nSensors   - number of sensors
 *     nBands     - number of frequency bands/bins
 *     H_array    - measured/modelled array responses for many directions;
 *                  FLAT: nBands x nSensors x nDirs
 *     nDirs      - number of directions the array was measured/modelled
 *     Y_grid     - spherical harmonics weights for each grid direction;
 *                  FLAT: nDirs x (order+1)^2
 * Output Arguments:
 *     cSH        - absolute values of the spatial correlation per band and
 *                  order; FLAT: nBands x (order+1)
 *     lSH        - level difference per band and order;
 *                  FLAT: nBands x (order+1)
 *
 * [1] Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording with
 *     higher order ambisonics-objective measurements and validation of
 *     spherical microphone. In Audio Engineering Society Convention 120.
 * [2] Politis, A., Gamper, H. (2017). "Comparing Modelled And Measurement-Based
 *     Spherical Harmonic Encoding Filters For Spherical Microphone Arrays. In
 *     IEEE Workshop on Applications of Signal Processing to Audio and Acoustics
 *     (WASPAA).
 */
void evaluateSHTfilters(/* Input arguments */
                        int order,
                        float_complex* M_array2SH,
                        int nSensors,
                        int nBands,
                        float_complex* H_array,
                        int nDirs,
                        float_complex* Y_grid,
                        /* Output arguments */
                        float* cSH,
                        float* lSH);
    
    
#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_SH_H_INCLUDED__ */
