/*
 Copyright 2016-2018 Leo McCormack
 
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
 *     saf_sh.h (include header)
 * Description:
 *     A collection of spherical harmonic related functions. Many of which have been
 *     derived from Matlab libraries by Archontis Politis; found here:
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
#endif
    
#include "../saf_utilities/saf_complex.h"
    
/****************/
/* Enum options */
/****************/

typedef enum _ARRAY_CONSTRUCTION_TYPES {
    ARRAY_CONSTRUCTION_OPEN,
    ARRAY_CONSTRUCTION_OPEN_DIRECTIONAL,
    ARRAY_CONSTRUCTION_RIGID,
    ARRAY_CONSTRUCTION_RIGID_DIRECTIONAL
}ARRAY_CONSTRUCTION_TYPES;
    

/*******************/
/* Misc. Functions */
/*******************/
    
/* Contructs a 3x3 rotation matrix from the Euler angles, using the yaw-pitch-roll (zyx) convention */
void yawPitchRoll2Rzyx (/* Input arguments */
                        float yaw,                /* yaw angle in radians */
                        float pitch,              /* pitch angle in radians */
                        float roll,               /* roll angle in radians */
                        int rollPitchYawFLAG,     /* 1: use Rxyz, i.e. apply roll, pitch and then yaw, 0: Rzyx / y-p-r */
                        /* Output arguments */
                        float R[3][3]);           /* zyx rotation matrix */

/* converts spherical coordinates to cartesian coordinates of unit length */
void unitSph2Cart(/* Input arguments */
                  float azi_rad,                  /* azimuth in radians */
                  float elev_rad,                 /* elevation in radians */
                  /* Output arguments */
                  float xyz[3]);                  /* unit cartesian coords, xyz */

/* converts cartesian coordinates of unit length to spherical coordinates */
void unitCart2Sph(/* Input arguments */
                  float xyz[3],                   /* unit cartesian coords, xyz */
                  /* Output arguments */
                  float AziElev_rad[2]);          /* azimuth and elevation in radians */

/* converts cartesian coordinates of unit length to spherical coordinates */
void unitCart2Sph_aziElev(/* Input arguments */
                          float xyz[3],           /* unit cartesian coords, xyz */
                          /* Output arguments */
                          float* azi_rad,         /* & azimuth in radians */
                          float* elev_rad);       /* & elevation in radians */
    
    
/****************************************/
/* SH and Beamforming related Functions */
/****************************************/
    
/* calculates unnormalised legendre polynomials up to order N, for all values in vector x
 * This INCLUDES the Condon-Shortley phase term. It is functionally identical to MatLab's legendre function, 'unnorm' (default).
 * M, Abramowitz., I.A. Stegun. (1965). "Handbook of Mathematical Functions: Chapter 8", Dover Publications.  */
void unnorm_legendreP(/* Input arguments */
                      int n,                      /* order of  legendre polynomial */
                      double* x,                  /* vector of input values; lenX x 1 */
                      int lenX,                   /* number of input values */
                      /* Output arguments */
                      double* y);                 /* resulting unnormalised legendre values for each x value; FLAT: (n+1) x lenX */

/* calculates unnormalised legendre polynomials values up to order N, for all values in vector x
 * does NOT INCLUDE the Condon-Shortley phase term
 * It uses a recursive approach, which makes it more suitable for computing the legendre values in a real-time loop */
void unnorm_legendreP_recur(/* Input arguments */
                            int n,                /* order of  legendre polynomial */
                            float* x,             /* vector of input values; lenX x 1 */
                            int lenX,             /* number of input values */
                            float* Pnm_minus1,    /* previous Pnm, (not used for n=1); FLAT: (n+1) x lenX */
                            float* Pnm_minus2,    /* previous previous Pnm, (not used for n=0); FLAT: (n+1) x lenX */
                            /* Output arguments */
                            float* Pnm);          /* resulting unnormalised legendre values for each x value; FLAT: (n+1) x lenX */
    
/* INTENDED FOR AMBISONICS, due to the omission of the 1/sqrt(4*pi) scaling, and the directions are given in
 * [azimuth elevation] (degrees). In Ambisonics literature, the format convention of 'Y' is referred to as ACN/N3D
 *
 * returns REAL spherical harmonics for multiple directions on the sphere. WITHOUT the 1/sqrt(4*pi) term. i.e. max(omni)=1
 * Compared to 'getRSH_recur', this function uses 'unnorm_legendreP' and double precision, so is more suitable for determining 'Y' in an
 * initialisation stage. This version is indeed slower, but more precise, especially for high orders.
 * For more information, the reader is  directed to:
 * Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8). Berlin: Springer. */
void getRSH(/* Input arguments */
            int order,                            /* order of spherical harmonic expansion */
            float* dirs_deg,                      /* directions on the sphere [azi, ELEVATION] convention, degrees; FLAT: nDirs x 2 */
            int nDirs,                            /* number of directions */
            /* Output arguments */
            float** Y);                           /* & the SH weights: FLAT: (order+1)^2 x nDirs */
    
/* INTENDED FOR AMBISONICS, due to the omission of the 1/sqrt(4*pi) term, and the directions are given in
 * [azimuth elevation] (degrees). In Ambisonics literature, the format convention of 'Y' is referred to as ACN/N3D
 *
 * returns REAL spherical harmonics for multiple directions on the sphere. WITHOUT the 1/sqrt(4*pi) scaling. i.e. max(omni)=1
 * Compared to 'getRSH', this function uses 'unnorm_legendreP_recur' and single precision, so is more suitable for determining 'Y' in a real-time
 * loop. It sacrifices some precision, as numerical error propogates through the recursion, but it is faster.
 * For more information, the reader is  directed to:
 * Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8). Berlin: Springer. */
void getRSH_recur(/* Input arguments */
                  int order,                      /* order of spherical harmonic expansion */
                  float* dirs_deg,                /* directions on the sphere [azi, ELEVATION] convention, degrees; FLAT: nDirs x 2 */
                  int nDirs,                      /* number of directions */
                  /* Output arguments */
                  float** Y);                     /* & the SH weights; FLAT: (order+1)^2 x nDirs */
    
/* returns real spherical harmonics for each direction on the sphere. WITH the 1/sqrt(4*pi) term.  i.e. max(omni)= 1/sqrt(4*pi)
 * compared to 'getSHreal_recur', this function employs 'unnorm_legendreP' and double precision, which is slower but more precise.
 * For more information, the reader is  directed to:
 * Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8). Berlin: Springer. */
void getSHreal(/* Input arguments */
               int order,                         /* order of spherical harmonic expansion */
               float* dirs_rad,                   /* directions on the sphere [azi, INCLINATION] convention, radians; FLAT: nDirs x 2 */
               int nDirs,                         /* number of directions */
               /* Output arguments */
               float* Y);                         /* the SH weights; FLAT:  (order+1)^2 x nDirs */
    
/* returns real spherical harmonics for each direction on the sphere. WITH the 1/sqrt(4*pi) term.  i.e. max(omni)= 1/sqrt(4*pi)
 * compared to 'getSHreal', this function employs 'unnorm_legendreP_recur' and single precision, which is faster but less precise.
 * For more information, the reader is  directed to:
 * Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8). Berlin: Springer. */
void getSHreal_recur(/* Input arguments */
                     int order,                   /* order of spherical harmonic expansion */
                     float* dirs_rad,             /* directions on the sphere [azi, INCLINATION] convention, radians; FLAT: nDirs x 2 */
                     int nDirs,                   /* number of directions */
                     /* Output arguments */
                     float* Y);                   /* the SH weights; FLAT: (order+1)^2 x nDirs */
    
/* returns complex spherical harmonics for each direction on the sphere. WITH the 1/sqrt(4*pi) term.  i.e. max(cabs(omni))= 1/sqrt(4*pi)
 * this function employs 'unnorm_legendreP' and double precision.
 * For more information, the reader is  directed to:
 * Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8). Berlin: Springer. */
void getSHcomplex(/* Input arguments */
                  int order,                      /* order of spherical harmonic expansion */
                  float* dirs_rad,                /* directions on the sphere [azi, INCLINATION] convention, radians; FLAT: nDirs x 2 */
                  int nDirs,                      /* number of directions */
                  /* Output arguments */
                  float_complex* Y);              /* the SH weights: (order+1)^2 x nDirs */
    
/* Returns the unitary transformation matrix T_c2r. It expresses the real spherical harmonics with respect to the complex ones,
 * so that r_N = T_c2r * y_N, where r_N and y_N is are the real and complex SH vectors, respectively */
void complex2realSHMtx(int order,                 /* order */
                       float_complex* T_c2r);     /* transformation matrix for complex->real; FLAT: (order+1)^2 x (order+1)^2  */
    
/* Returns the unitary transformation matrix T_r2c the expresses the complex spherical harmonics with respect to the real ones,
 * so that y_N = T_r2c * r_N, where r_N and y_N are the real and complex SH vectors respectively. */
void real2complexSHMtx(int order,                 /* order */
                       float_complex* T_r2c);     /* transformation matrix for real->complex; FLAT: (order+1)^2 x (order+1)^2  */
    
/* Convert SH coeffs from the complex to real basis */
void complex2realCoeffs(int order,                 /* order */
                        float_complex* C_N,        /* complex coeffients; FLAT: (order+1)^2 x K */
                        int K,                     /* number of columns */
                        float* R_N);               /* real coefficients; FLAT: (order+1)^2 x K */
    
/* generates a real-valued spherical harmonic rotation matrix (assumes ACN/N3D convention)
 * For more information, the reader is referred to:
 *     Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 *     by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
void getSHrotMtxReal(float R[3][3],               /* zyx rotation matrix */
                     float* RotMtx,               /* the rotation matrix; FLAT: (L+1)^2 x (L+1)^2 */
                     int L);                      /* order */
    
/* this routine computes the matrices that generate the coefficients of the beampattern of order (sectorOrder+1) that is
 * essentially the product of a pattern of order=sectorOrder and a dipole. It is used in "beamWeightsVelocityPatterns".
 * For the derivation of the matrices see:
 *     Politis, A. and Pulkki, V., 2016. Acoustic intensity, energy-density and diffuseness estimation in a
 *     directionally-constrained region.  arXiv preprint arXiv:1609.03409 */
void computeVelCoeffsMtx(/* Input arguments */
                         int sectorOrder,         /* order of patterns; */
                         /* Output arguments */
                         float_complex* A_xyz);   /* Velocity coefficients; FLAT: (sectorOrder+2)^2  x (sectorOrder+1)^2 x 3 */
    
/* Generate spherical coefficients for cardioids. For a specific order N of a higher order cardioid of the form
 * D(theta)=(1/2)^N * (1+cos(theta))^N, generate the beamweights for the same pattern in the SHD. Because the
 * pattern is axisymmetric only the N+1 coefficients of m=0 are returned. */
void beamWeightsCardioid2Spherical(/* Input arguments */
                                   int N,         /* order */
                                   /* Output arguments */
                                   float* b_n);   /* beamformer weights; (N+1) x 1 */
    
    // NOT IMPLEMENTED YET
/* Generate beamweights in the SHD for Dolph-Chebyshev beampatterns, with mainlobe and sidelobe control.
 * Because the pattern is axisymmetric only the N+1 coefficients of m=0 are returned.
 *     Koretz, A. and Rafaely, B., 2009. Dolph-Chebyshev beampattern design for spherical arrays.
 *     IEEE Transactions on Signal Processing, 57(6), pp.2417-2420. */
void beamWeightsDolphChebyshev2Spherical(/* Input arguments */
                                         int N,         /* order */
                                         int paramType, /* 0: side-lobe level control, 1: mainlobe width control */
                                         float arrayParam, /* sidelobe level 1/R or mainlobe with 2*a0 */
                                         /* Output arguments */
                                         float* b_n);   /* beamformer weights; (N+1) x 1 */
    
/* The hypercardioid is the pattern that maximises the directivity-factor for a certain SH order N. The
 * hypercardioid is also the plane-wave decomposition beamformer in the SHD, also called 'regular' because
 * the beamweights are just the SH values on the beam-direction. Since the pattern is axisymmetric only
 * the N+1 coefficients of m=0 are  returned. */
void beamWeightsHypercardioid2Spherical(/* Input arguments */
                                        int N,    /* order */
                                        /* Output arguments */
                                        float* b_n); /* beamformer weights; (N+1) x 1 */

/* Generate the beamweights for the a maximum energy-vector beampattern in the SHD. This pattern originates
 * from ambisonic-related research and it maximises the ambisonic energy-vector, which is essentially the
 * directional centroid of the squared pattern. IT can also be seen as the pattern that maximizes the acoustic
 * intensity vector of a diffuse field weighted with this pattern. In practice it is almost the same as
 * a supercardioid that maximizes front-back power ratio for a certain order, and it can be used as such.
 * Because the pattern is axisymmetric only the N+1 coefficients of m=0 are returned. Details for their theory
 * can be found e.g. in
 *     Zotter, F., Pomberger, H. and Noisternig, M., 2012.
 *     Energy-preserving ambisonic decoding. Acta Acustica united with Acustica, 98(1), pp.37-47. */
void beamWeightsMaxEV(/* Input arguments */
                      int N,                      /* order */
                      /* Output arguments */
                      float* b_n);                /* beamformer weights; (N+1) x 1 */
    
/* If the sound-field is weighted with an axisymmetric spatial distribution described by the N+1 SH coefficients
 * b_n, then the beamweights capturing the velocity signals for the weighted sound-field are of an order one higher
 * than the weighting pattern, and can be derived from it. This type of beamforming has some applications for spatial
 * sound reproduction and acoustic analysis, see
 *     Politis, A. and Pulkki, V., 2016. Acoustic intensity, energy-density and diffuseness estimation in a
 *     directionally-constrained region. arXiv preprint arXiv:1609.03409. */
void beamWeightsVelocityPatternsReal(/* Input arguments */
                                     int order,                    /* order */
                                     float* b_n,                   /* axisymmetric beamformer weights; (order+1) x 1 */
                                     float azi_rad,                /* orientation, azimuth in radius */
                                     float elev_rad,               /* orientation, elevation in radius */
                                     float_complex* A_xyz,         /* FLAT: (order+2)^2 x (order+1)^2 x 3 */
                                     /* Output arguments */
                                     float* velCoeffs);            /* FLAT: (order+2)^2 x 3 */
    
/* If the sound-field is weighted with an axisymmetric spatial distribution described by the N+1 SH coefficients
 * b_n, then the beamweights capturing the velocity signals for the weighted sound-field are of an order one higher
 * than the weighting pattern, and can be derived from it. This type of beamforming has some applications for spatial
 * sound reproduction and acoustic analysis, see
 *     Politis, A. and Pulkki, V., 2016. Acoustic intensity, energy-density and diffuseness estimation in a
 *     directionally-constrained region. arXiv preprint arXiv:1609.03409. */
void beamWeightsVelocityPatternsComplex(/* Input arguments */
                                        int order,                 /* order */
                                        float* b_n,                /* axisymmetric beamformer weights; (order+1) x 1 */
                                        float azi_rad,             /* orientation, azimuth in radius */
                                        float elev_rad,            /* orientation, elevation in radius */
                                        float_complex* A_xyz,      /* FLAT: (order+2)^2 x (order+1)^2 x 3 */
                                        /* Output arguments */
                                        float_complex* velCoeffs); /* FLAT: (order+2)^2 x 3 */
    
/* returns spherical coefficients for a rotated axisymmetric pattern */
void rotateAxisCoeffsReal(/* Input arguments */
                          int order,              /* order */
                          float* c_n,             /* coefficients describing a rotationally symmetric pattern order N,
                                                   * expressed as a sum of spherical harmonics of degree m=0; (N+1) x 1 */
                          float theta_0,          /* polar rotation for the pattern, radians */
                          float phi_0,            /* azimuthal rotation for the pattern, radians */
                          /* Output arguments */
                          float* c_nm);           /* coefficients of rotated pattern expressed as a sum of SHs; (N+1)^2 x 1 */
    
/* returns spherical coefficients for a rotated axisymmetric pattern */
void rotateAxisCoeffsComplex(/* Input arguments */
                             int order,           /* order */
                             float* c_n,          /* coefficients describing a rotationally symmetric pattern order N,
                                                   * expressed as a sum of spherical harmonics of degree m=0; (N+1) x 1 */
                             float theta_0,       /* polar rotation for the pattern, radians */
                             float phi_0,         /* azimuthal rotation for the pattern, radians */
                             /* Output arguments */
                             float_complex* c_nm); /* coefficients of rotated pattern expressed as a sum of SHs; (N+1)^2 x 1 */
    
/* computes the condition numbers for a least-squares SHT */
void checkCondNumberSHTReal(/* Input arguments */
                            int order,            /* order */
                            float* dirs_rad,      /* directions on the sphere [azi, INCLINATION] convention, radians; FLAT: nDirs x 2 */
                            int nDirs,            /* number of directions */
                            float* w,             /* weights; nDirs x 1 */
                            /* Output arguments */
                            float* cond_N);       /* condition numbers; (order+1) x 1 */
    
    
/**************************************/
/* Localisation Functions in the  SHD */
/**************************************/
    
/* generates a powermap utilising the PWD method */
void generatePWDmap(/* Input arguments */
                    int order,                    /* analysis order */
                    float_complex* Cx,            /* covarience matrix; FLAT: (order+1)^2 x (order+1)^2 */
                    float_complex* Y_grid,        /* steering vectors for grid direcionts; FLAT: (order+1)^2 x nGrid_dirs  */
                    int nGrid_dirs,               /* number of grid directions */
                    /* Output arguments */
                    float* pmap);                 /* resulting PWD powermap; nGrid_dirs x 1 */

/* generates a powermap utilising the MVDR method*/
void generateMVDRmap(/* Input arguments */
                     int order,                   /* analysis order */
                     float_complex* Cx,           /* covarience matrix; FLAT: (order+1)^2 x (order+1)^2 */
                     float_complex* Y_grid,       /* steering vectors for grid direcionts; FLAT: (order+1)^2 x nGrid_dirs  */
                     int nGrid_dirs,              /* number of grid directions */
                     float regPar,                /* regularisation parameter, for diagonal loading of Cx */
                     /* Output arguments */
                     float* pmap,                 /* resulting MVDR powermap; nGrid_dirs x 1 */
                     float_complex* w_MVDR);      /* optional. weights will be copied to this, unless it's NULL; FLAT: nSH x nGrid_dirs || NULL */

/* EXPERIMENTAL! Generates a powermap utilising the CroPaC LCMV post-filter described in:
 * Delikaris-Manias, S., Vilkamo, J., & Pulkki, V. (2016). Signal-dependent spatial filtering based on
 * weighted-orthogonal beamformers in the spherical harmonic domain. IEEE/ACM Transactions on Audio,
 * Speech and Language Processing (TASLP), 24(9), 1507-1519.
 *
 * The spatial post-filter is estimated for all directions on the grid, and is used to supress reverb/noise
 * interference that may be present in an MVDR map. Unlike in the paper, the second column for the contraints
 * 'A', is Y.*diag(Cx), rather than utilising a maximum energy beamformer. The post-filters are then applied
 * to the MVDR powermap map derived in the sherical harmonic domain, rather than an MVDR beamformer generated
 * directly in the microphone array signal domain, like in the paper. Otherwise, the algorithm is the same. */
void generateCroPaCLCMVmap(/* Input arguments */
                           int order,             /* analysis order */
                           float_complex* Cx,     /* covarience matrix; FLAT: (order+1)^2 x (order+1)^2 */
                           float_complex* Y_grid, /* steering vectors for grid direcionts; FLAT: (order+1)^2 x nGrid_dirs  */
                           int nGrid_dirs,        /* number of grid directions */
                           float regPar,          /* regularisation parameter, for diagonal loading of Cx */
                           float lambda,          /* parameter controlling how harsh CroPaC is applied, 0..1; 0: fully cropac, 1: fully mvdr */
                           /* Output arguments */
                           float* pmap);          /* resulting CroPaC LCMV powermap; nGrid_dirs x 1 */
    
/* generates a powermap utilising the subspace-based MUSIC method */
void generateMUSICmap(/* Input arguments */
                      int order,                  /* analysis order */
                      float_complex* Cx,          /* covarience matrix; FLAT: (order+1)^2 x (order+1)^2 */
                      float_complex* Y_grid,      /* steering vectors for grid direcionts; FLAT: (order+1)^2 x nGrid_dirs  */
                      int nSources,               /* number of sources present in sound scene */
                      int nGrid_dirs,             /* number of grid directions */
                      int logScaleFlag,           /* 1: log(pmap), 0: pmap. */
                      /* Output arguments */
                      float* pmap);               /* resulting MUSIC pseudo-spectrum; nGrid_dirs x 1 */

/* generates a powermap utilising the subspace-based MinNorm method */
void generateMinNormMap(/* Input arguments */
                        int order,                /* analysis order */
                        float_complex* Cx,        /* covarience matrix; FLAT: (order+1)^2 x (order+1)^2 */
                        float_complex* Y_grid,    /* steering vectors for grid direcionts; FLAT: (order+1)^2 x nGrid_dirs  */
                        int nSources,             /* number of sources present in sound scene */
                        int nGrid_dirs,           /* number of grid directions */
                        int logScaleFlag,         /* 1: log(pmap), 0: pmap. */
                        /* Output arguments */
                        float* pmap);             /* resulting MinNorm pseudo-spectrum; nGrid_dirs x 1 */

    
/******************************************/
/* Cylindrical/Spherical Bessel Functions */
/******************************************/
    
/* (cylindrical) Bessel function of the first kind: Jn
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_Jn(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical precision) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               double* J_n,                       /* Bessel values (set as NULL if not required); FLAT: nZ x (N+1) */
               double* dJ_n);                     /* Bessel derivative values (set as NULL if not required); FLAT: nZ x (N+1) */
    
/* (cylindrical) Bessel function of the second kind: Yn
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_Yn(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical precision) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               double* Y_n,                       /* Bessel values (set as NULL if not required); FLAT: nZ x (N+1) */
               double* dY_n);                     /* Bessel derivative values (set as NULL if not required); FLAT: nZ x (N+1) */
    
/* (cylindrical) Hankel function of the first kind: Hn1
 * returns the Hankel values and their derivatives up to order N for all values in vector z  */
void hankel_Hn1(/* Input arguments */
                int N,                            /* function order (highest is ~30 given numerical precision) */
                double* z,                        /* input values; nZ x 1 */
                int nZ,                           /* number of input values */
                /* Output arguments */
                double_complex* Hn1_n,            /* Hankel values (set as NULL if not required); FLAT: nZ x (N+1) */
                double_complex* dHn1_n);          /* Hankel derivative values (set as NULL if not required); FLAT: nZ x (N+1) */
    
/* (cylindrical) Hankel function of the second kind: Hn2
 * returns the Hankel values and their derivatives up to order N for all values in vector z  */
void hankel_Hn2(/* Input arguments */
                int N,                            /* function order (highest is ~30 given numerical precision) */
                double* z,                        /* input values; nZ x 1 */
                int nZ,                           /* number of input values */
                /* Output arguments */
                double_complex* Hn2_n,            /* Hankel values (set as NULL if not required); FLAT: nZ x (N+1) */
                double_complex* dHn2_n);          /* Hankel derivative values (set as NULL if not required); FLAT: nZ x (N+1) */
    
/* spherical Bessel function of the first kind: jn
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_jn(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical precision) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               int* maxN,                         /* & maximum function order that could be computed <=N */
               double* j_n,                       /* Bessel values (set as NULL if not required); FLAT: nZ x (N+1) */
               double* dj_n);                     /* Bessel derivative values (set as NULL if not required); FLAT: nZ x (N+1) */
    
/* modified spherical Bessel function of the first kind: in
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_in(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical precision) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               int* maxN,                         /* & maximum function order that could be computed <=N */
               double* i_n,                       /* Bessel values (set as NULL if not required); FLAT: nZ x (N+1) */
               double* di_n);                     /* Bessel derivative values (set as NULL if not required); FLAT: nZ x (N+1) */

/* spherical Bessel function of the second kind (Neumann): yn
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_yn(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical precision) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               int* maxN,                         /* & maximum function order that could be computed <=N */
               double* y_n,                       /* Bessel values (set as NULL if not required); FLAT: nZ x (N+1) */
               double* dy_n);                     /* Bessel derivative values (set as NULL if not required); FLAT: nZ x (N+1) */
    
/* modified spherical Bessel function of the second kind: kn
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_kn(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical precision) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               int* maxN,                         /* & maximum function order that could be computed <=N */
               double* k_n,                       /* Bessel values (set as NULL if not required); FLAT: nZ x (N+1) */
               double* dk_n);                     /* Bessel derivative values (set as NULL if not required); FLAT: nZ x (N+1) */

/* spherical Hankel function of the first kind: hn1
 * returns the Hankel values and their derivatives up to order N for all values in vector z */
void hankel_hn1(/* Input arguments */
                int N,                            /* function order (highest is ~30 given numerical precision) */
                double* z,                        /* input values; nZ x 1 */
                int nZ,                           /* number of input values */
                /* Output arguments */
                int* maxN,                        /* & maximum function order that could be computed <=N */
                double_complex* h_n1,             /* Hankel values (set as NULL if not required); FLAT: nZ x (N+1) */
                double_complex* dh_n1);           /* Hankel derivative values (set as NULL if not required); FLAT: nZ x (N+1) */

/* spherical Hankel function of the second kind: hn2
 * returns the Hankel values and their derivatives up to order N for all values in vector z */
void hankel_hn2(/* Input arguments */
                int N,                            /* function order (highest is ~30 given numerical precision) */
                double* z,                        /* input values; nZ x 1 */
                int nZ,                           /* number of input values */
                /* Output arguments */
                int* maxN,                        /* & maximum function order that could be computed <=N */
                double_complex* h_n2,             /* Hankel values (set as NULL if not required); FLAT: nZ x (N+1) */
                double_complex* dh_n2);           /* Hankel derivative values (set as NULL if not required); FLAT: nZ x (N+1) */
    
    
/*****************************************/
/* Microphone array processing functions */
/*****************************************/

/* calculates the modal coefficients for open/rigid cylindrical arrays */
void cylModalCoeffs(/* Input arguments */
                    int order,                    /* max order (highest is ~30 given numerical precision) */
                    double* kr,                   /* wavenumber*radius; nBands x 1 */
                    int nBands,                   /* number of frequency bands/bins */
                    ARRAY_CONSTRUCTION_TYPES arrayType, /* see 'ARRAY_CONSTRUCTION_TYPES' enum */
                    /* Output arguments */
                    double_complex* b_N);         /* modal coefficients per kr and 0:order; FLAT: nBands x (order+1) */
    
/* returns the simplest estimate of the spatial aliasing limit (the kR = maxN rule) */
float sphArrayAliasLim(/* Input arguments */
                       float r,                   /* mic radius, meters */
                       float c,                   /* speed of sound, m/s */
                       int maxN);                 /* order */
    
/* returns the frequencies that the noise in the output channels of a SMA, after performing the SHT and equalization of
 the output signals, reaches a certain user-defined threshold maxG_db. The frequencies are computed only at the lower range of each order,
 where its response decays rapidly, ignoring for example the nulls of an open array at the higher frequencies. The estimation of the limits are
 based on a linear approximation of the log-log response found e.g. in
     Sector-based Parametric Sound Field Reproduction in the Spherical Harmonic Domain A Politis, J Vilkamo, V Pulkki
     IEEE Journal of Selected Topics in Signal Processing 9 (5), 852 - 866 */
void sphArrayNoiseThreshold(/* Input arguments */
                            int maxN,             /* maximum order of the array */
                            int Nsensors,         /* number of sensors */
                            float r,              /* mic radius, meters */
                            float c,              /* speed of sound, m/s */
                            ARRAY_CONSTRUCTION_TYPES arrayType, /* see 'ARRAY_CONSTRUCTION_TYPES' enum */
                            double dirCoeff,      /* only for directional (open) arrays, 1: omni, 0.5: card, 0:dipole */
                            float maxG_db,        /* max allowed amplification for the noise level, maxG_db = 20*log10(maxG) */
                            /* Output arguments */
                            float* f_lim);        /* frequency points that the threhsold is reached; (max_N+1) x 1 */

/* calculates the modal coefficients for open/rigid spherical arrays */
void sphModalCoeffs(/* Input arguments */
                    int order,                    /* max order (highest is ~30 given numerical precision) */
                    double* kr,                   /* wavenumber*radius; nBands x 1 */
                    int nBands,                   /* number of frequency bands/bins */
                    ARRAY_CONSTRUCTION_TYPES arrayType, /* see 'ARRAY_CONSTRUCTION_TYPES' enum */
                    double dirCoeff,              /* only for directional (open) arrays, 1: omni, 0.5: card, 0:dipole */
                    /* Output arguments */
                    double_complex* b_N);         /* modal coefficients per kr and 0:order; FLAT: nBands x (order+1) */

/* calculates the modal coefficients for a rigid spherical scatterer with omnidirectional sensors
   (Assumes all sensors are placed the same distance from the scatterer, w.r.t. the origin) */
void sphScattererModalCoeffs(/* Input arguments */
                             int order,           /* max order (highest is ~30 given numerical precision) */
                             double* kr,          /* wavenumber*sensor_radius; nBands x 1 */
                             double* kR,          /* wavenumber*scatterer_radius; nBands x 1 */
                             int nBands,          /* number of frequency bands/bins */
                             /* Output arguments */
                             double_complex* b_N);/* modal coefficients per kr and 0:order; FLAT: nBands x (order+1) */
    
/* calculates the modal coefficients for a rigid spherical scatterer with directional sensors
   (Assumes all sensors are placed the same distance from the scatterer, w.r.t. the origin) */
void sphScattererDirModalCoeffs(/* Input arguments */
                                int order,        /* max order (highest is ~30 given numerical precision) */
                                double* kr,       /* wavenumber*sensor_radius; nBands x 1 */
                                double* kR,       /* wavenumber*scatterer_radius; nBands x 1 */
                                int nBands,       /* number of frequency bands/bins */
                                double dirCoeff,  /* directivity coefficient, 1: omni, 0.5: card, 0:dipole */
                                /* Output arguments */
                                double_complex* b_N);/* modal coefficients per kr and 0:order; FLAT: nBands x (order+1) */
    
/* calculates the theoretical diffuse coherence matrix for a spherical array */
void sphDiffCohMtxTheory(/* Input arguments */
                         int order,               /* max order (highest is ~30 given numerical precision) */
                         float* sensor_dirs_rad,  /* spherical coords of the sensors in RADIANS, [azi elev]; FLAT: N_sensors x 2 */
                         int N_sensors,           /* number of sensors */
                         ARRAY_CONSTRUCTION_TYPES arrayType, /* see 'ARRAY_CONSTRUCTION_TYPES' enum */
                         double dirCoeff,         /* only for directional (open) arrays, 1: omni, 0.5: card, 0:dipole */
                         double* kr,              /* wavenumber*sensor_radius; nBands x 1 */
                         double* kR,              /* wavenumber*scatterer_radius, set to NULL if not applicable; nBands x 1 */
                         int nBands,              /* number of frequency bands/bins */
                         /* Output arguments */
                         double* M_diffcoh);      /* theoretical diffuse coherence matrix per frequency; FLAT: N_sensors x N_sensors x nBands */

/* simulates a cylindrical microphone array, returning the transfer functions for each (plane wave) source direction
 * on the surface of the cylinder. */
void simulateCylArray(/* Input arguments */
                      int order,                  /* max order (highest is ~30 given numerical precision) */
                      double* kr,                 /* wavenumber*radius; nBands x 1 */
                      int nBands,                 /* number of frequency bands/bins */
                      float* sensor_dirs_rad,     /* spherical coords of the sensors in RADIANS, [azi elev]; FLAT: N_sensors x 2 */
                      int N_sensors,              /* number of sensors */
                      float* src_dirs_deg,        /* spherical coords of the plane waves in DEGREES, [azi elev]; FLAT: N_srcs x 2 */
                      int N_srcs,                 /* number sources (DoAs of plane waves) */
                      ARRAY_CONSTRUCTION_TYPES arrayType, /* see 'ARRAY_CONSTRUCTION_TYPES' enum */ 
                      /* Output arguments */
                      float_complex* H_array);    /* simulated array response for each plane wave; FLAT: nBands x N_sensors x N_srcs */

/* simulates a spherical microphone array, returning the transfer functions for each (plane wave) source direction
 * on the surface of the sphere. */
void simulateSphArray(/* Input arguments */
                      int order,                  /* max order (highest is ~30 given numerical precision) */
                      double* kr,                 /* wavenumber*radius; nBands x 1 */
                      double* kR,                 /* wavenumber*scatterer_radius, set to NULL if not applicable */
                      int nBands,                 /* number of frequency bands/bins */
                      float* sensor_dirs_rad,     /* spherical coords of the sensors in RADIANS, [azi elev]; FLAT: N_sensors x 2 */
                      int N_sensors,              /* number of sensors */
                      float* src_dirs_deg,        /* spherical coords of the plane waves in DEGREES, [azi elev]; FLAT: N_srcs x 2 */
                      int N_srcs,                 /* number sources (DoAs of plane waves) */
                      ARRAY_CONSTRUCTION_TYPES arrayType, /* see 'ARRAY_CONSTRUCTION_TYPES' enum */
                      double dirCoeff,            /* only for directional (open) arrays, 1: omni, 0.5: card, 0:dipole */
                      /* Output arguments */
                      float_complex* H_array);    /* simulated array response for each plane wave; FLAT: nBands x N_sensors x N_srcs */

/* generates some objective measures, which evaluate the performance of the spatial encoding filters. This analysis is performed by comparing the
 * the spatial resolution of the spherical harmonic components generated by the encoding filters, with the ideal SH components.
 * For more information, the reader is directed to:
 * Moreau, S., Daniel, J., Bertet, S., 2006, 3D sound field recording with higher order ambisonics-objective measurements and
 * validation of spherical microphone. In Audio Engineering Society Convention 120.
 * and:
 * Politis, A., Gamper, H. (2017). "Comparing Modelled And Measurement-Based Spherical Harmonic Encoding Filters For Spherical
 * Microphone Arrays. In IEEE Workshop on Applications of Signal Processing to Audio and Acoustics (WASPAA).
 */
void evaluateSHTfilters(/* Input arguments */
                        int order,                /* transform order */
                        float_complex* M_array2SH,/* encoding matrices; FLAT: nBands x (order+1)^2 x nSensors */
                        int nSensors,             /* number of sensors */
                        int nBands,               /* number of frequency bands/bins */
                        float_complex* H_array,   /* measured/modelled array responses for many directions; FLAT: nBands x nSensors x nDirs */
                        int nDirs,                /* number of directions the array was measured/modelled */
                        float_complex* Y_grid,    /* spherical harmonics weights for each grid direction; FLAT: nDirs x (order+1)^2 */
                        /* Output arguments */
                        float* cSH,               /* absolute values of the spatial correlation per band and order; FLAT: nBands x (order+1) */
                        float* lSH);              /* level difference per band and order; FLAT: nBands x (order+1) */

#ifdef __cplusplus
}
#endif

#endif /* __SAF_SH_H_INCLUDED__ */




