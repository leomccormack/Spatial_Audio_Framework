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
 * @file saf_sh.h
 * @brief Public part of the "saf_sh" module
 *
 * A collection of spherical harmonic related functions. Many of which have been
 * derived from Matlab libraries by Archontis Politis [1-3].
 *
 * @see [1] https://github.com/polarch/Spherical-Harmonic-Transform
 * @see [2] https://github.com/polarch/Array-Response-Simulator
 * @see [3] https://github.com/polarch/Spherical-Array-Processing
 *
 * @author Leo McCormack
 * @date 22.05.2016
 */

#ifndef __SAF_SH_H_INCLUDED__
#define __SAF_SH_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "../saf_utilities/saf_complex.h"

#define ORDER2NSH(order) ((order+1)*(order+1))

/* ========================================================================== */
/*                                    Enums                                   */
/* ========================================================================== */

/**
 * Microphone/Hydrophone array contruction types
 */
typedef enum _ARRAY_CONSTRUCTION_TYPES {
    ARRAY_CONSTRUCTION_OPEN,             /**< Open array, omni-directional
                                          *   sensors */
    ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL, /**< Open array, directional sensors */
    ARRAY_CONSTRUCTION_RIGID,            /**< Rigid baffle, omni-directional
                                          *   sensors */
    ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL /**< Rigid baffle, directional sensors
                                          */
}ARRAY_CONSTRUCTION_TYPES;
    
/**
 * Sector pattern designs for directionally-contraining sound-fields [1].
 *
 * @see [1] Politis, A., & Pulkki, V. (2016). Acoustic intensity, energy-density
 *          and diffuseness estimation in a directionally-constrained region.
 *          arXiv preprint arXiv:1609.03409
 */
typedef enum _SECTOR_PATTERNS{
    SECTOR_PATTERN_PWD,     /**< Plane-wave decomposition/Hyper-cardioid */
    SECTOR_PATTERN_MAXRE,   /**< Spatially tapered hyper-cardioid, such that it
                             *   has maximum energy concentrated in the look-
                             *   direction */
    SECTOR_PATTERN_CARDIOID /**< Cardioid pattern */
    
}SECTOR_PATTERNS;


/* ========================================================================== */
/*                               Misc. Functions                              */
/* ========================================================================== */

/**
 * Constructs a 3x3 rotation matrix from the Euler angles, using the
 * yaw-pitch-roll (zyx) convention.
 *
 * @param[in]  yaw              Yaw angle in radians
 * @param[in]  pitch            Pitch angle in radians
 * @param[in]  roll             Roll angle in radians
 * @param[in]  rollPitchYawFLAG '1' to use Rxyz, i.e. apply roll, pitch and then
 *                              yaw, '0' Rzyx / y-p-r
 * @param[out] R                zyx rotation matrix; 3 x 3
 */
void yawPitchRoll2Rzyx (/* Input Arguments */
                        float yaw,
                        float pitch,
                        float roll,
                        int rollPitchYawFLAG,
                        /* Output Arguments */
                        float R[3][3]);

/**
 * Converts spherical coordinates to cartesian coordinates of unit length.
 *
 * @param[in]  azi_rad  Azimuth in radians
 * @param[in]  elev_rad Elevation in radians
 * @param[out] xyz      Unit cartesian coords, xyz; 3 x 1
 */
void unitSph2Cart(/* Input Arguments */
                  float azi_rad,
                  float elev_rad,
                  /* Output Arguments */
                  float xyz[3]);

/**
 * Converts cartesian coordinates of unit length to spherical coordinates.
 *
 * @param[in]  xyz         Unit cartesian coords, xyz
 * @param[out] AziElev_rad Azimuth and elevation in radians
 */
void unitCart2Sph(/* Input Arguments */
                  float xyz[3],
                  /* Output Arguments */
                  float AziElev_rad[2]);

/**
 * Converts cartesian coordinates of unit length to spherical coordinates.
 *
 * @param[in]  xyz      Unit cartesian coords, xyz
 * @param[out] azi_rad  (&) azimuth in radians
 * @param[out] elev_rad (&) elevation in radians
 */
void unitCart2Sph_aziElev(/* Input Arguments */
                          float xyz[3],
                          /* Output Arguments */
                          float* azi_rad,
                          float* elev_rad);

/**
 * Calculates unnormalised legendre polynomials up to order N, for all values in
 * vector x [1].
 *
 * @note This INCLUDES the Condon-Shortley phase term. It is functionally
 *       identical to Matlab's legendre function (with default settings
 *       ['unnorm']).
 *
 * @param[in]  n    Order of  legendre polynomial
 * @param[in]  x    Vector of input values; lenX x 1
 * @param[in]  lenX Number of input values
 * @param[out] y    Resulting unnormalised legendre values for each x value;
 *                  FLAT: (n+1) x lenX
 *
 * @see [1] M, Abramowitz., I.A. Stegun. (1965). "Handbook of Mathematical
 *          Functions: Chapter 8", Dover Publications.
 */
void unnorm_legendreP(/* Input Arguments */
                      int n,
                      double* x,
                      int lenX,
                      /* Output Arguments */
                      double* y);

/**
 * Calculates unnormalised legendre polynomials up to order N, for all values in
 * vector x.
 *
 * It uses a recursive approach, which makes it more suitable for computing the
 * legendre values in a real-time loop.
 *
 * @note This does NOT INCLUDE the Condon-Shortley phase term.
 *
 * @param[in]  n          Order of  legendre polynomial
 * @param[in]  x          Vector of input values; lenX x 1
 * @param[in]  lenX       Number of input values
 * @param[in]  Pnm_minus1 Previous Pnm, (not used for n=1); FLAT: (n+1) x lenX
 * @param[in]  Pnm_minus2 Previous previous Pnm, (not used for n=0);
 *                        FLAT: (n+1) x lenX
 * @param[out] Pnm        Resulting unnormalised legendre values for each x
 *                        value; FLAT: (n+1) x lenX
 */
void unnorm_legendreP_recur(/* Input Arguments */
                            int n,
                            float* x,
                            int lenX,
                            float* Pnm_minus1,
                            float* Pnm_minus2,
                            /* Output Arguments */
                            float* Pnm);


/* ========================================================================== */
/*                    SH and Beamforming related Functions                    */
/* ========================================================================== */

/**
 * Computes REAL spherical harmonics [1] for each direction on the sphere.
 *
 * The real spherical harmonics are computed WITH the 1/sqrt(4*pi) term.
 * i.e. max(omni) = 1/sqrt(4*pi). Compared to getSHreal_recur(), this function
 * employs unnorm_legendreP() and double precision, which is slower but more
 * precise.
 *
 * @param[in]  order    Order of spherical harmonic expansion
 * @param[in]  dirs_rad Directions on the sphere [azi, INCLINATION] convention,
 *                      in RADIANS; FLAT: nDirs x 2
 * @param[in]  nDirs    Number of directions
 * @param[out] Y        The SH weights [WITH the 1/sqrt(4*pi)];
 *                      FLAT: (order+1)^2 x nDirs
 *
 * @see [1] Rafaely, B. (2015). Fundamentals of spherical array processing
 *          (Vol. 8). Berlin: Springer.
 */
void getSHreal(/* Input Arguments */
               int order,
               float* dirs_rad,
               int nDirs,
               /* Output Arguments */
               float* Y);

/**
 * Computes REAL spherical harmonics [1] for each direction on the sphere.
 *
 * The real spherical harmonics are computed WITH the 1/sqrt(4*pi) term.
 * i.e. max(omni) = 1/sqrt(4*pi). Compared to getSHreal(), this function employs
 * unnorm_legendreP_recur() and single precision, which is faster but less
 * precise.
 *
 * @param[in]  order    Order of spherical harmonic expansion
 * @param[in]  dirs_rad Directions on the sphere [azi, INCLINATION] convention,
 *                      in RADIANS; FLAT: nDirs x 2
 * @param[in]  nDirs    Number of directions
 * @param[out] Y        The SH weights [WITH the 1/sqrt(4*pi)];
 *                      FLAT: (order+1)^2 x nDirs
 *
 * @see [1] Rafaely, B. (2015). Fundamentals of spherical array processing
 *          (Vol. 8). Berlin: Springer.
 */
void getSHreal_recur(/* Input Arguments */
                     int order,
                     float* dirs_rad,
                     int nDirs,
                     /* Output Arguments */
                     float* Y);

/**
 * Computes COMPLEX spherical harmonics [1] for each direction on the sphere.
 *
 * The real spherical harmonics are computed WITH the 1/sqrt(4*pi) term.
 * i.e. max(omni) = 1/sqrt(4*pi) + i0. This function employs unnorm_legendreP()
 * and double precision.
 *
 * @param[in]  order    Order of spherical harmonic expansion
 * @param[in]  dirs_rad Directions on the sphere [azi, INCLINATION] convention,
 *                      in RADIANS; FLAT: nDirs x 2
 * @param[in]  nDirs    Number of directions
 * @param[out] Y        The SH weights [WITH the 1/sqrt(4*pi)];
 *                      FLAT: (order+1)^2 x nDirs
 *
 * @see [1] Rafaely, B. (2015). Fundamentals of spherical array processing
 *          (Vol. 8). Berlin: Springer.
 */
void getSHcomplex(/* Input Arguments */
                  int order,
                  float* dirs_rad,
                  int nDirs,
                  /* Output Arguments */
                  float_complex* Y);

/**
 * Computes a complex to real spherical harmonic transform matrix
 *
 * Computes the unitary transformation matrix T_c2r. It expresses the real
 * spherical harmonics with respect to the complex ones, so that
 * r_N = T_c2r * y_N, where r_N and y_N is are the real and complex SH vectors,
 * respectively.
 *
 * @param[in]  order Order of spherical harmonic expansion
 * @param[out] T_c2r Transformation matrix for complex->real;
 *                   FLAT: (order+1)^2 x (order+1)^2
 */
void complex2realSHMtx(/* Input Arguments */
                       int order,
                       /* Output Arguments */
                       float_complex* T_c2r);

/**
 * Computes a real to complex spherical harmonic transform matrix
 *
 * Computes the unitary transformation matrix T_r2c the expresses the complex
 * spherical harmonics with respect to the real ones, so that y_N = T_r2c * r_N,
 * where r_N and y_N are the real and complex SH vectors, respectively.
 *
 * @param[in]  order Order of spherical harmonic expansion
 * @param[out] T_r2c Transformation matrix for real->complex;
 *                   FLAT: (order+1)^2 x (order+1)^2
 */
void real2complexSHMtx(/* Input Arguments */
                       int order,
                       /* Output Arguments */
                       float_complex* T_r2c);

/**
 * Converts SH coefficients from the complex to real basis
 *
 * @param[in]  order Order of spherical harmonic expansion
 * @param[in]  C_N   Complex coeffients; FLAT: (order+1)^2 x K
 * @param[in]  K     Number of columns
 * @param[out] R_N   Real coefficients; FLAT: (order+1)^2 x K
 */
void complex2realCoeffs(/* Input Arguments */
                        int order,
                        float_complex* C_N,
                        int K,
                        /* Output Arguments */
                        float* R_N);

/**
 * Generates a real-valued spherical harmonic rotation matrix [1] (assumes ACN
 * channel ordering convention).
 *
 * Note that the normalisation convention does not matter, as e.g. only dipoles
 * are used to rotated dipoles, quadrapoles to rotate quadrapoles etc.
 *
 * @param[in]  R      zyx rotation matrix; 3 x 3
 * @param[in]  L      Order of spherical harmonic expansion
 * @param[out] RotMtx SH domain rotation matrix; FLAT: (L+1)^2 x (L+1)^2
 *
 * @see [1] Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real
 *          Spherical Harmonics. Direct Determination by Recursion Page:
 *          Additions and Corrections. Journal of Physical Chemistry A, 102(45),
 *          9099?9100.
 */
void getSHrotMtxReal(float R[3][3],
                     float* RotMtx,
                     int L);

/**
 * Computes the matrices that generate the coefficients of the beampattern of
 * order (sectorOrder+1) that is essentially the product of a pattern of
 * order=sectorOrder and a dipole.
 *
 * It is used in "beamWeightsVelocityPatterns". For the derivation of the
 * matrices see [1].
 *
 * @param[in]  sectorOrder Order of patterns
 * @param[out] A_xyz       Velocity coefficients;
 *                         FLAT: (sectorOrder+2)^2  x (sectorOrder+1)^2 x 3
 *
 * @see [1] Politis, A. and Pulkki, V., 2016. Acoustic intensity, energy-density
 *          and diffuseness estimation in a directionally-constrained region.
 *          arXiv preprint arXiv:1609.03409
 */ 
void computeVelCoeffsMtx(/* Input Arguments */
                         int sectorOrder,
                         /* Output Arguments */
                         float_complex* A_xyz);

/**
 * Computes the beamforming matrices of sector and velocity coefficients for
 * ENERGY-preserving (EP) sectors for real SH.
 *
 * This partitioning of the sound-field into spatially-localised sectors has
 * been used for parametric sound-field reproduction in [1] and visualisation in
 * [2,3].
 *
 * @param[in]  orderSec     Order of sector patterns
 * @param[in]  A_xyz        Velocity coefficients (see "computeVelCoeffsMtx");
 *                          FLAT: (sectorOrder+2)^2 x (sectorOrder+1)^2 x 3
 * @param[in]  pattern      See "SECTOR_PATTERNS" enum for the options
 * @param[in]  sec_dirs_deg Sector directions [azi elev], in DEGREES;
 *                          FLAT: nSecDirs x 2
 * @param[in]  nSecDirs     Number of sectors
 * @param[out] sectorCoeffs The sector coefficients;
 *                          FLAT: (nSecDirs*4) x (orderSec+2)^2
 * @returns                 Normalisation coefficient
 *
 * @see [1] Politis, A., Vilkamo, J., & Pulkki, V. (2015). Sector-based
 *          parametric sound field reproduction in the spherical harmonic
 *          domain. IEEE Journal of Selected Topics in Signal Processing, 9(5),
 *          852-866.
 * @see [2] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of
 *          angular spectra based on a directional re-assignment approach for
 *          ambisonic sound-field visualisation". IEEE International Conference '
 *          on Acoustics, Speech and Signal Processing (ICASSP).
 * @see [3] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 */
float computeSectorCoeffsEP(/* Input Arguments */
                            int orderSec,
                            float_complex*  A_xyz,
                            SECTOR_PATTERNS pattern,
                            float* sec_dirs_deg,
                            int nSecDirs,
                            /* Output Arguments */
                            float* sectorCoeffs);

/**
 * Computes the beamforming matrices of sector and velocity coefficients for
 * AMPLITUDE-preserving (AP) sectors for real SH.
 *
 * This partitioning of the sound-field into spatially-localised sectors has
 * been used for parametric sound-field reproduction in [1] and visualision in
 * [2,3].
 *
 * @param[in]  orderSec     Order of sector patterns
 * @param[in]  A_xyz        Velocity coefficients (see "computeVelCoeffsMtx");
 *                          FLAT: (sectorOrder+2)^2 x (sectorOrder+1)^2 x 3
 * @param[in]  pattern      See "SECTOR_PATTERNS" enum for the options
 * @param[in]  sec_dirs_deg Sector directions [azi elev], in DEGREES;
 *                          FLAT: nSecDirs x 2
 * @param[in]  nSecDirs     Number of sectors
 * @param[out] sectorCoeffs The sector coefficients;
 *                          FLAT: (nSecDirs*4) x (orderSec+2)^2
 * @returns                 Normalisation coefficient
 *
 * @see [1] Politis, A., Vilkamo, J., & Pulkki, V. (2015). Sector-based
 *          parametric sound field reproduction in the spherical harmonic
 *          domain. IEEE Journal of Selected Topics in Signal Processing, 9(5),
 *          852-866.
 * @see [2] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of
 *          angular spectra based on a directional re-assignment approach for
 *          ambisonic sound-field visualisation". IEEE International Conference
 *          on Acoustics, Speech and Signal Processing (ICASSP).
 * @see [3] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 */
float computeSectorCoeffsAP(/* Input Arguments */
                            int orderSec,
                            float_complex* A_xyz,
                            SECTOR_PATTERNS pattern,
                            float* sec_dirs_deg,
                            int nSecDirs,
                            /* Output Arguments */
                            float* sectorCoeffs);

/**
 * Generates spherical coefficients for cardioids.
 *
 * For a specific order N of a higher order cardioid of the form
 * D(theta)=(1/2)^N * (1+cos(theta))^N, generate the beamweights for the same
 * pattern in the SHD. Because the pattern is axisymmetric only the N+1
 * coefficients of m=0 are returned.
 *
 * @param[in]  N   Order of spherical harmonic expansion
 * @param[out] b_n Beamformer weights; (N+1) x 1
 */
void beamWeightsCardioid2Spherical(/* Input Arguments */
                                   int N,
                                   /* Output Arguments */
                                   float* b_n);

/**
 * Generates beamweights in the SHD for Dolph-Chebyshev beampatterns, with
 * mainlobe and sidelobe control [1].
 *
 * Because the pattern is axisymmetric only the N+1 coefficients of m=0 are
 * returned.
 *
 * @note NOT IMPLEMENTED YET
 *
 * @param[in]  N          Order of spherical harmonic expansion
 * @param[in]  paramType  '0' side-lobe level control, '1' mainlobe width
 *                        control
 * @param[in]  arrayParam Sidelobe level 1/R or mainlobe with 2*a0
 * @param[out] b_n        Beamformer weights; (N+1) x 1
 *
 * @see [1] Koretz, A. and Rafaely, B., 2009. Dolph-Chebyshev beampattern design
 *          for spherical arrays. IEEE Transactions on Signal Processing, 57(6),
 *          pp.2417-2420.
 */
void beamWeightsDolphChebyshev2Spherical(/* Input Arguments */
                                         int N,
                                         int paramType,
                                         float arrayParam,
                                         /* Output Arguments */
                                         float* b_n);

/**
 * Generates beamweights in the SHD for hypercardioid beampatterns.
 *
 * The hypercardioid is the pattern that maximises the directivity-factor for a
 * certain SH order N. The hypercardioid is also the plane-wave decomposition
 * beamformer in the SHD, also called 'regular' because the beamweights are just
 * the SH values on the beam-direction. Since the pattern is axisymmetric only
 * the N+1 coefficients of m=0 are returned.
 *
 * @param[in]  N   Order of spherical harmonic expansion
 * @param[out] b_n Beamformer weights; (N+1) x 1
 */
void beamWeightsHypercardioid2Spherical(/* Input Arguments */
                                        int N,
                                        /* Output Arguments */
                                        float* b_n);

/**
 * Generates beamweights in the SHD for maximum energy-vector beampatterns
 *
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
 * @param[in]  N   Order of spherical harmonic expansion
 * @param[out] b_n Beamformer weights; (N+1) x 1
 *
 * @see [1] Zotter, F., Pomberger, H. and Noisternig, M., 2012. Energy-
 *          preserving ambisonic decoding. Acta Acustica united with Acustica,
 *          98(1), pp.37-47.
 */
void beamWeightsMaxEV(/* Input Arguments */
                      int N,
                      /* Output Arguments */
                      float* b_n);

/**
 * Generates beamforming coefficients for velocity patterns (REAL)
 *
 * If the sound-field is weighted with an axisymmetric spatial distribution
 * described by the N+1 SH coefficients b_n, then the beamweights capturing the
 * velocity signals for the weighted sound-field are of an order one higher
 * than the weighting pattern, and can be derived from it. This type of
 * beamforming has some applications for spatial sound reproduction and
 * acoustic analysis, see [1].
 *
 * @param[in]  order     Order of spherical harmonic expansion
 * @param[in]  b_n       Axisymmetric beamformer weights; (order+1) x 1
 * @param[in]  azi_rad   Orientation, azimuth in RADIANS
 * @param[in]  elev_rad  Orientation, ELEVATION in RADIANS
 * @param[in]  A_xyz     Velocity coefficients; see computeVelCoeffsMtx();
 *                       FLAT: (sectorOrder+2)^2 x (sectorOrder+1)^2 x 3
 * @param[out] velCoeffs Beamforming coefficients for velocity patterns;
 *                       FLAT: (order+2)^2 x 3
 *
 * @see [1] Politis, A. and Pulkki, V., 2016. Acoustic intensity, energy-density
 *          and diffuseness estimation in a directionally-constrained region.
 *          arXiv preprint arXiv:1609.03409.
 */
void beamWeightsVelocityPatternsReal(/* Input Arguments */
                                     int order,
                                     float* b_n,
                                     float azi_rad,
                                     float elev_rad,
                                     float_complex* A_xyz,
                                     /* Output Arguments */
                                     float* velCoeffs);

/**
 * Generates beamforming coefficients for velocity patterns (COMPLEX)
 *
 * If the sound-field is weighted with an axisymmetric spatial distribution
 * described by the N+1 SH coefficients b_n, then the beamweights capturing the
 * velocity signals for the weighted sound-field are of an order one higher
 * than the weighting pattern, and can be derived from it. This type of
 * beamforming has some applications for spatial sound reproduction and
 * acoustic analysis, see [1].
 *
 * @param[in]  order     Order of spherical harmonic expansion
 * @param[in]  b_n       Axisymmetric beamformer weights; (order+1) x 1
 * @param[in]  azi_rad   Orientation, azimuth in RADIANS
 * @param[in]  elev_rad  Orientation, ELEVATION in RADIANS
 * @param[in]  A_xyz     Velocity coefficients; see computeVelCoeffsMtx();
 *                       FLAT: (sectorOrder+2)^2 x (sectorOrder+1)^2 x 3
 * @param[out] velCoeffs Beamforming coefficients for velocity patterns;
 *                       FLAT: (order+2)^2 x 3
 *
 * @see [1] Politis, A. and Pulkki, V., 2016. Acoustic intensity, energy-density
 *          and diffuseness estimation in a directionally-constrained region.
 *          arXiv preprint arXiv:1609.03409.
 */
void beamWeightsVelocityPatternsComplex(/* Input Arguments */
                                        int order,
                                        float* b_n,
                                        float azi_rad,
                                        float elev_rad,
                                        float_complex* A_xyz,
                                        /* Output Arguments */
                                        float_complex* velCoeffs);

/**
 * Generates spherical coefficients for a rotated axisymmetric pattern (REAL)
 *
 * @param[in]  order   Order of spherical harmonic expansion
 * @param[in]  c_n     Coefficients describing a rotationally symmetric pattern
 *                     order N, expressed as a sum of spherical harmonics of
 *                     degree m=0; (N+1) x 1
 * @param[in]  theta_0 POLAR rotation for the pattern, in RADIANS
 * @param[in]  phi_0   Azimuthal rotation for the pattern, in RADIANS
 * @param[out] c_nm    Coefficients of rotated pattern expressed as a sum of
 *                     SHs; (N+1)^2 x 1
 */
void rotateAxisCoeffsReal(/* Input arguments */
                          int order,
                          float* c_n,
                          float theta_0,
                          float phi_0,
                          /* Output arguments */
                          float* c_nm);          

/**
 * Generates spherical coefficients for a rotated axisymmetric pattern (COMPLEX)
 *
 * @param[in]  order   Order of spherical harmonic expansion
 * @param[in]  c_n     Coefficients describing a rotationally symmetric pattern
 *                     order N, expressed as a sum of spherical harmonics of
 *                     degree m=0; (N+1) x 1
 * @param[in]  theta_0 POLAR rotation for the pattern, in RADIANS
 * @param[in]  phi_0   Azimuthal rotation for the pattern, in RADIANS
 * @param[out] c_nm    Coefficients of rotated pattern expressed as a sum of
 *                     SHs; (N+1)^2 x 1
 */
void rotateAxisCoeffsComplex(/* Input arguments */
                             int order,
                             float* c_n,
                             float theta_0,
                             float phi_0,
                             /* Output arguments */
                             float_complex* c_nm);

/**
 * Computes the condition numbers for a least-squares SHT
 *
 * @param[in]  order    Order of spherical harmonic expansion
 * @param[in]  dirs_rad Directions on the sphere [azi, INCLINATION] convention,
 *                      in RADIANS; FLAT: nDirs x 2
 * @param[in]  nDirs    Number of directions
 * @param[in]  w        Integration weights; nDirs x 1
 * @param[out] cond_N   Condition numbers; (order+1) x 1
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

/**
 * Generates a powermap based on the energy of plane-wave decomposition (PWD)/
 * hyper-cardioid beamformers.
 *
 * @param[in]  order      Analysis order
 * @param[in]  Cx         Correlation/covarience matrix;
 *                        FLAT: (order+1)^2 x (order+1)^2
 * @param[in]  Y_grid     Steering vectors for each grid direcionts;
 *                        FLAT: (order+1)^2 x nGrid_dirs
 * @param[in]  nGrid_dirs Number of grid directions
 * @param[out] pmap       Resulting PWD powermap; nGrid_dirs x 1
 */
void generatePWDmap(/* Input arguments */
                    int order,
                    float_complex* Cx,
                    float_complex* Y_grid,
                    int nGrid_dirs,
                    /* Output arguments */
                    float* pmap);

/**
 * Generates a powermap based on the energy of adaptive minimum variance
 * distortion-less response (MVDR) beamformers.
 *
 * @param[in]  order      Analysis order
 * @param[in]  Cx         Correlation/covarience matrix;
 *                        FLAT: (order+1)^2 x (order+1)^2
 * @param[in]  Y_grid     Steering vectors for each grid direcionts;
 *                        FLAT: (order+1)^2 x nGrid_dirs
 * @param[in]  nGrid_dirs Number of grid directions
 * @param[in]  regPar     Regularisation parameter, for diagonal loading of Cx
 * @param[out] pmap       Resulting MVDR powermap; nGrid_dirs x 1
 * @param[out] w_MVDR     (Optional) weights will be copied to this, unless
 *                        it's NULL; FLAT: nSH x nGrid_dirs || NULL
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

/**
 * EXPERIMENTAL! Generates a powermap utilising the CroPaC LCMV post-filter
 * described in [1].
 *
 * The spatial post-filter is estimated for all directions on the grid, and is
 * used to supress reverb/noise interference that may be present in an MVDR map.
 * Unlike in the paper, the second column for the contraints 'A', is
 * Y.*diag(Cx), rather than utilising a maximum energy beamformer. The post-
 * filters are then applied to the MVDR powermap map derived in the sherical
 * harmonic domain, rather than an MVDR beamformer generated directly in the
 * microphone array signal domain, like in the paper. Otherwise, the algorithm
 * is the same.
 *
 * @param[in]  order      Analysis order
 * @param[in]  Cx         Correlation/covarience matrix;
 *                        FLAT: (order+1)^2 x (order+1)^2
 * @param[in]  Y_grid     Steering vectors for each grid direcionts;
 *                        FLAT: (order+1)^2 x nGrid_dirs
 * @param[in]  nGrid_dirs Number of grid directions
 * @param[in]  regPar     Regularisation parameter, for diagonal loading of Cx
 * @param[in]  lambda     Parameter controlling how harsh CroPaC is applied,
 *                        0..1; 0: fully CroPaC, 1: fully MVDR
 * @param[out] pmap       Resulting CroPaC LCMV powermap; nGrid_dirs x 1
 *
 * @see [1] Delikaris-Manias, S., Vilkamo, J., & Pulkki, V. (2016). Signal-
 *          dependent spatial filtering based on weighted-orthogonal beamformers
 *          in the spherical harmonic domain. IEEE/ACM Transactions on Audio,
 *          Speech and Language Processing (TASLP), 24(9), 1507-1519.
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

/**
 * Generates an activity-map based on the sub-space multiple-signal
 * classification (MUSIC) method.
 *
 * @param[in]  order        Analysis order
 * @param[in]  Cx           Correlation/covarience matrix;
 *                          FLAT: (order+1)^2 x (order+1)^2
 * @param[in]  Y_grid       Steering vectors for each grid direcionts;
 *                          FLAT: (order+1)^2 x nGrid_dirs
 * @param[in]  nSources     Number of sources present in sound scene
 * @param[in]  nGrid_dirs   Number of grid directions
 * @param[in]  logScaleFlag '1' log(pmap), '0' pmap.
 * @param[out] pmap         Resulting MUSIC pseudo-spectrum; nGrid_dirs x 1
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

/**
 * Generates an activity-map based on the sub-space minimum-norm (MinNorm)
 * method.
 *
 * @param[in]  order        Analysis order
 * @param[in]  Cx           Correlation/covarience matrix;
 *                          FLAT: (order+1)^2 x (order+1)^2
 * @param[in]  Y_grid       Steering vectors for each grid direcionts;
 *                          FLAT: (order+1)^2 x nGrid_dirs
 * @param[in]  nSources     Number of sources present in sound scene
 * @param[in]  nGrid_dirs   Number of grid directions
 * @param[in]  logScaleFlag '1' log(pmap), '0' pmap.
 * @param[out] pmap         Resulting MinNorm pseudo-spectrum; nGrid_dirs x 1
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

/**
 * Computes the (cylindrical) Bessel function of the first kind: Jn
 *
 * Computes the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * @param[in]  N    Function order (highest is ~30 given numerical precision)
 * @param[in]  z    Input values; nZ x 1
 * @param[in]  nZ   Number of input values
 * @param[out] J_n  Bessel values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 * @param[out] dJ_n Bessel derivative values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 */
void bessel_Jn(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               double* J_n,
               double* dJ_n);

/**
 * Computes the (cylindrical) Bessel function of the second kind: Yn
 *
 * Computes the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * @param[in]  N    Function order (highest is ~30 given numerical precision)
 * @param[in]  z    Input values; nZ x 1
 * @param[in]  nZ   Number of input values
 * @param[out] Y_n  Bessel values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 * @param[out] dY_n Bessel derivative values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 */
void bessel_Yn(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               double* Y_n,
               double* dY_n);

/**
 * Computes the (cylindrical) Hankel function of the first kind: Hn1
 *
 * Computes the Hankel values and their derivatives up to order N for all values
 * in vector z
 *
 * @param[in]  N      Function order (highest is ~30 given numerical precision)
 * @param[in]  z      Input values; nZ x 1
 * @param[in]  nZ     Number of input values
 * @param[out] Hn1_n  Hankel values (set as NULL if not required);
 *                    FLAT: nZ x (N+1)
 * @param[out] dHn1_n Hankel derivative values (set as NULL if not required);
 *                    FLAT: nZ x (N+1)
 */
void hankel_Hn1(/* Input arguments */
                int N,
                double* z,
                int nZ,
                /* Output arguments */
                double_complex* Hn1_n,
                double_complex* dHn1_n);

/**
 * Computes the (cylindrical) Hankel function of the second kind: Hn2
 *
 * Computes the Hankel values and their derivatives up to order N for all values
 * in vector z
 *
 * @param[in]  N      Function order (highest is ~30 given numerical precision)
 * @param[in]  z      Input values; nZ x 1
 * @param[in]  nZ     Number of input values
 * @param[out] Hn2_n  Hankel values (set as NULL if not required);
 *                    FLAT: nZ x (N+1)
 * @param[out] dHn2_n Hankel derivative values (set as NULL if not required);
 *                    FLAT: nZ x (N+1)
 */
void hankel_Hn2(/* Input arguments */
                int N,
                double* z,
                int nZ,
                /* Output arguments */
                double_complex* Hn2_n,
                double_complex* dHn2_n);

/**
 * Computes the spherical Bessel function of the first kind: jn
 *
 * Computes the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * @param[in]  N    Function order (highest is ~30 given numerical precision)
 * @param[in]  z    Input values; nZ x 1
 * @param[in]  nZ   Number of input values
 * @param[out] maxN (&) maximum function order that could be computed <=N
 * @param[out] j_n  Bessel values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 * @param[out] dj_n Bessel derivative values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 */
void bessel_jn(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               int* maxN,
               double* j_n,
               double* dj_n);

/**
 * Computes the modified spherical Bessel function of the first kind: in
 *
 * Computes the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * @param[in]  N    Function order (highest is ~30 given numerical precision)
 * @param[in]  z    Input values; nZ x 1
 * @param[in]  nZ   Number of input values
 * @param[out] maxN (&) maximum function order that could be computed <=N
 * @param[out] i_n  Bessel values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 * @param[out] di_n Bessel derivative values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 */
void bessel_in(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               int* maxN,
               double* i_n,
               double* di_n);

/**
 * Computes the spherical Bessel function of the second kind (Neumann): yn
 *
 * Computes the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * @param[in]  N    Function order (highest is ~30 given numerical precision)
 * @param[in]  z    Input values; nZ x 1
 * @param[in]  nZ   Number of input values
 * @param[out] maxN (&) maximum function order that could be computed <=N
 * @param[out] y_n  Bessel values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 * @param[out] dy_n Bessel derivative values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 */
void bessel_yn(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               int* maxN,
               double* y_n,
               double* dy_n);

/**
 * Computes the modified spherical Bessel function of the second kind: kn
 *
 * Computes the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * @param[in]  N    Function order (highest is ~30 given numerical precision)
 * @param[in]  z    Input values; nZ x 1
 * @param[in]  nZ   Number of input values
 * @param[out] maxN (&) maximum function order that could be computed <=N
 * @param[out] k_n  Bessel values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 * @param[out] dk_n Bessel derivative values (set as NULL if not required);
 *                  FLAT: nZ x (N+1)
 */
void bessel_kn(/* Input arguments */
               int N,
               double* z,
               int nZ,
               /* Output arguments */
               int* maxN,
               double* k_n,
               double* dk_n);

/**
 * Computes the spherical Hankel function of the first kind: hn1
 *
 * Computes the Hankel values and their derivatives up to order N for all values
 * in vector z
 *
 * @param[in]  N     Function order (highest is ~30 given numerical precision)
 * @param[in]  z     Input values; nZ x 1
 * @param[in]  nZ    Number of input values
 * @param[out] maxN  (&) maximum function order that could be computed <=N
 * @param[out] h_n1  Hankel values (set as NULL if not required);
 *                   FLAT: nZ x (N+1)
 * @param[out] dh_n1 Hankel derivative values (set as NULL if not required);
 *                   FLAT: nZ x (N+1)
 */
void hankel_hn1(/* Input arguments */
                int N,
                double* z,
                int nZ,
                /* Output arguments */
                int* maxN,
                double_complex* h_n1,
                double_complex* dh_n1);

/**
 * Computes the spherical Hankel function of the second kind: hn2
 *
 * Computes the Hankel values and their derivatives up to order N for all values
 * in vector z
 *
 * @param[in]  N     Function order (highest is ~30 given numerical precision)
 * @param[in]  z     Input values; nZ x 1
 * @param[in]  nZ    Number of input values
 * @param[out] maxN  (&) maximum function order that could be computed <=N
 * @param[out] h_n2  Hankel values (set as NULL if not required);
 *                   FLAT: nZ x (N+1)
 * @param[out] dh_n2 Hankel derivative values (set as NULL if not required);
 *                   FLAT: nZ x (N+1)
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

/**
 * Calculates the modal coefficients for open/rigid cylindrical arrays
 *
 * @param[in]  order     Max order (highest is ~30 given numerical precision)
 * @param[in]  kr        wavenumber*radius; nBands x 1
 * @param[in]  nBands    Number of frequency bands/bins
 * @param[in]  arrayType See 'ARRAY_CONSTRUCTION_TYPES' enum
 * @param[out] b_N       Modal coefficients per kr and 0:order;
 *                       FLAT: nBands x (order+1)
 */
void cylModalCoeffs(/* Input arguments */
                    int order,
                    double* kr,
                    int nBands,
                    ARRAY_CONSTRUCTION_TYPES arrayType,
                    /* Output arguments */
                    double_complex* b_N);

/**
 * Returns a simple estimate of the spatial aliasing limit (the kR = maxN rule)
 *
 * @param[in] r    Array radius, meters
 * @param[in] c    Speed of sound, m/s
 * @param[in] maxN Order
 * @returns        Spatial aliasing limit estimate, in Hz
 */
float sphArrayAliasLim(/* Input arguments */
                       float r,
                       float c,
                       int maxN);

/**
 * Computes the frequncies (per order), at which the noise of a SHT of a SMA
 * exceeds a specified maximum level.
 *
 * Computes the frequencies that the noise in the output channels of a spherical
 * microphone array (SMA), after performing the spherical harmonic transform
 * (SHT) and equalisation of the output signals, reaches a certain user-defined
 * threshold maxG_db. The frequencies are computed only at the lower range of
 * each order, where its response decays rapidly, ignoring for example the nulls
 * of an open array at the higher frequencies. The estimation of the limits are
 * based on a linear approximation of the log-log response found e.g. in [1]
 *
 * @param[in]  maxN      Maximum order of the array
 * @param[in]  Nsensors  Number of sensors
 * @param[in]  r         Mic radius, meters
 * @param[in]  c         Speed of sound, m/s
 * @param[in]  arrayType See 'ARRAY_CONSTRUCTION_TYPES' enum
 * @param[in]  dirCoeff  Only for directional (open) arrays, 1: omni, 0.5: card,
 *                       0:dipole
 * @param[in]  maxG_db   Max allowed amplification for the noise level,
 *                       maxG_db = 20*log10(maxG)
 * @param[out] f_lim     Noise limit estimate; (maxN+1) x 1
 *
 * @see [1] Politis, A., Vilkamo, J., & Pulkki, V. (2015). Sector-based
 *          parametric sound field reproduction in the spherical harmonic
 *          domain. IEEE Journal of Selected Topics in Signal Processing, 9(5),
 *          852-866.
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

/**
 * Calculates the modal coefficients for open/rigid spherical arrays
 *
 * @param[in]  order     Max order (highest is ~30 given numerical precision)
 * @param[in]  kr        wavenumber*radius; nBands x 1
 * @param[in]  nBands    Number of frequency bands/bins
 * @param[in]  arrayType See 'ARRAY_CONSTRUCTION_TYPES' enum
 * @param[in]  dirCoeff  Only for directional (open) arrays, 1: omni, 0.5: card,
 *                       0:dipole
 * @param[out] b_N       Modal coefficients per kr and 0:order;
 *                       FLAT: nBands x (order+1)
 */
void sphModalCoeffs(/* Input arguments */
                    int order,
                    double* kr,
                    int nBands,
                    ARRAY_CONSTRUCTION_TYPES arrayType,
                    double dirCoeff,
                    /* Output arguments */
                    double_complex* b_N);

/**
 * Calculates the modal coefficients for a rigid spherical scatterer with
 * omni-directional sensors
 *
 * Assumes all sensors are placed the same distance from the scatterer, w.r.t.
 * the origin
 *
 * @param[in]  order  Max order (highest is ~30 given numerical precision)
 * @param[in]  kr     wavenumber*array_radius; nBands x 1
 * @param[in]  kR     wavenumber*scatterer_radius; nBands x 1
 * @param[in]  nBands Number of frequency bands/bins
 * @param[out] b_N    Modal coefficients per kr and 0:order;
 *                    FLAT: nBands x (order+1)
 */
void sphScattererModalCoeffs(/* Input arguments */
                             int order,
                             double* kr,
                             double* kR,
                             int nBands,
                             /* Output arguments */
                             double_complex* b_N);

/**
 * Calculates the modal coefficients for a rigid spherical scatterer with
 * directional sensors
 *
 * Assumes all sensors are placed the same distance from the scatterer, w.r.t.
 * the origin
 *
 * @param[in]  order    Max order (highest is ~30 given numerical precision)
 * @param[in]  kr       wavenumber*array_radius; nBands x 1
 * @param[in]  kR       wavenumber*scatterer_radius; nBands x 1
 * @param[in]  nBands   Number of frequency bands/bins
 * @param[in]  dirCoeff Directivity coefficient, 1: omni, 0.5: card, 0:dipole
 * @param[out] b_N      Modal coefficients per kr and 0:order;
 *                      FLAT: nBands x (order+1)
 */
void sphScattererDirModalCoeffs(/* Input arguments */
                                int order,
                                double* kr,
                                double* kR,
                                int nBands,
                                double dirCoeff,
                                /* Output arguments */
                                double_complex* b_N);
    
/**
 * Calculates the theoretical diffuse coherence matrix for a spherical array
 *
 * @param[in]  order           Max order (highest is ~30 given numerical
 *                             precision)
 * @param[in]  sensor_dirs_rad Spherical coords of the sensors in RADIANS,
 *                             [azi ELEV]; FLAT: N_sensors x 2
 * @param[in]  N_sensors       Number of sensors
 * @param[in]  arrayType       See 'ARRAY_CONSTRUCTION_TYPES' enum
 * @param[in]  dirCoeff        Only for directional (open) arrays, 1: omni,
 *                             0.5: card, 0:dipole
 * @param[in]  kr              wavenumber*sensor_radius; nBands x 1
 * @param[in]  kR              wavenumber*scatterer_radius, set to NULL if not
 *                             applicable; nBands x 1
 * @param[in]  nBands          Number of frequency bands/bins
 * @param[out] M_diffcoh       Theoretical diffuse coherence matrix per
 *                             frequency; FLAT: N_sensors x N_sensors x nBands
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

/**
 * Simulates a cylindrical microphone array, returning the transfer functions
 * for each (plane wave) source direction on the surface of the cylinder
 *
 * @param[in]  order           Max order (highest is ~30 given numerical
 *                             precision)
 * @param[in]  kr              wavenumber*radius; nBands x 1
 * @param[in]  nBands          Number of frequency bands/bins
 * @param[in]  sensor_dirs_rad Spherical coords of the sensors in RADIANS,
 *                             [azi ELEV]; FLAT: N_sensors x 2
 * @param[in]  N_sensors       Number of sensors
 * @param[in]  src_dirs_deg    Spherical coords of the plane waves in DEGREES,
 *                             [azi ELEV]; FLAT: N_srcs x 2
 * @param[in]  N_srcs          Number sources (DoAs of plane waves)
 * @param[in]  arrayType       See 'ARRAY_CONSTRUCTION_TYPES' enum
 * @param[out] H_array         Simulated array response for each plane wave;
 *                             FLAT: nBands x N_sensors x N_srcs
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

/**
 * Simulates a spherical microphone array, returning the transfer functions for
 * each (plane wave) source direction on the surface of the sphere.
 *
 * @param[in]  order           Max order (highest is ~30 given numerical
 *                             precision)
 * @param[in]  kr              wavenumber*array_radius; nBands x 1
 * @param[in]  kR              wavenumber*scatterer_radius, set to NULL if not
 *                             needed
 * @param[in]  nBands          Number of frequency bands/bins
 * @param[in]  sensor_dirs_rad Spherical coords of the sensors in RADIANS,
 *                             [azi ELEV]; FLAT: N_sensors x 2
 * @param[in]  N_sensors       Number of sensors
 * @param[in]  src_dirs_deg    Spherical coords of the plane waves in DEGREES,
 *                             [azi ELEV]; FLAT: N_srcs x 2
 * @param[in]  N_srcs          Number sources (DoAs of plane waves)
 * @param[in]  arrayType       See 'ARRAY_CONSTRUCTION_TYPES' enum
 * @param[in]  dirCoeff        Only for directional (open) arrays, 1: omni,
 *                             0.5: card, 0:dipole
 * @param[out]  H_array        Simulated array response for each plane wave;
 *                             FLAT: nBands x N_sensors x N_srcs
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

/**
 * Generates some objective measures, which evaluate the performance of the
 * spatial encoding filters.
 *
 * This analysis is performed by comparing the spatial resolution of the
 * spherical harmonic components generated by the encoding filters, with the
 * ideal SH components. For more information, the reader is directed to [1,2].
 *
 * @param[in]  order      Transform/encoding order
 * @param[in]  M_array2SH Encoding matrix per frequency;
 *                        FLAT: nBands x (order+1)^2 x nSensors
 * @param[in]  nSensors   Number of sensors
 * @param[in]  nBands     Number of frequency bands/bins
 * @param[in]  H_array    Measured/modelled array responses for many directions;
 *                        FLAT: nBands x nSensors x nDirs
 * @param[in]  nDirs      Number of directions the array was measured/modelled
 * @param[in]  Y_grid     Spherical harmonics weights for each grid direction;
 *                        FLAT: nDirs x (order+1)^2
 * @param[out] cSH        Absolute values of the spatial correlation per band
 *                        and order; FLAT: nBands x (order+1)
 * @param[out] lSH        Level difference per band and order;
 *                        FLAT: nBands x (order+1)
 *
 * @see [1] Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording
 *          with higher order ambisonics--objective measurements and validation
 *          of spherical microphone. In Audio Engineering Society Convention
 *          120.
 * @see [2] Politis, A., Gamper, H. (2017). "Comparing Modelled And Measurement-
 *          Based Spherical Harmonic Encoding Filters For Spherical Microphone
 *          Arrays. In IEEE Workshop on Applications of Signal Processing to
 *          Audio and Acoustics (WASPAA).
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
