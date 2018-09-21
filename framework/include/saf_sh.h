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
 *     A collection of spherical harmonic related functions. Some of which have been
 *     derived from the Matlab library by Archontis Politis; found here:
 *     https://github.com/polarch/Spherical-Harmonic-Transform
 *     and MATLAB code by Symeon Delikaris-Manias
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
    
#include "saf_utilities.h"
    
#ifndef M_PI
  #define M_PI ( 3.14159265359f )
#endif
    
/*************************/
/* Processing modes tags */
/*************************/

typedef enum _BEAMFORMING_WEIGHT_TYPES {
    BFW_BASIC,                 /* beamforming weights = spherical harmonic weights for one direction on the sphere */
    BFW_MAX_RE,                /* maximum-energy beamformer */
    BFW_DOLPH_CHEBY_MAIN,      /* Dolph-Chebyshev beamfomer */
    BFW_DOLPH_CHEBY_DESIRED    /* Dolph-Chebyshev beamfomer */

} BEAMFORMING_WEIGHT_TYPES;
    
typedef enum _ARRAY_CONSTRUCTION_TYPES {
    ARRAY_CONSTRUCTION_OPEN,
    ARRAY_CONSTRUCTION_RIGID,
    ARRAY_CONSTRUCTION_DIRECTIONAL
}ARRAY_CONSTRUCTION_TYPES;
    
/******************/
/* Main Functions */
/******************/
    
/* NOTE: legendreP be removed in a future version, use "unnorm_legendreP" */
/* Computes unnormalised legendre polynomial of order 0 to L, at position x
 * see: http://mathworld.wolfram.com/LegendrePolynomial.html */
void legendreP(/* Input arguments */
               int L,                             /* maximum order of legendre polynomial */
               float x,                           /* position */
               /* Output arguments */
               float* ppm);                       /* the polynomials for orders 0 to L */
    
/* calculates unnormalised legendre values up to order N, for all values in vector x */
/* M, Abramowitz., I.A. Stegun. (1965). "Handbook of Mathematical Functions: Chapter 8", Dover Publications.  */
void unnorm_legendreP(/* Input arguments */
                      int n,                      /* order of  legendre polynomial */
                      double* x,                  /* vector of input values; lenX x 1 */
                      int lenX,                   /* number of input values */
                      /* Output arguments */
                      double* y);                 /* resulting unnormalised legendre values for each x value; FLAT: (n+1) x lenX */
    
/* returns real spherical harmonics for multiple directions on the sphere. WITHOUT the 1/sqrt(4*pi) scaling
 * For more information, the reader is  directed to:
 * Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8). Berlin: Springer. */
void getRSH(/* Input arguments */
            int N,                                /* order of spherical harmonic expansion */
            float* dirs_deg,                      /* directions on the sphere [azi, elev] convention; FLAT: nDirs x 2 */
            int nDirs,                            /* number of directions */
            /* Output arguments */
            float** Y);                           /* & the SH weights: FLAT: (N+1)^2 x nDirs */
    
/* returns real spherical harmonics for a direction on the sphere. WITH the 1/sqrt(4*pi) scaling
 * For more information, the reader is  directed to:
 * Rafaely, B. (2015). Fundamentals of spherical array processing (Vol. 8). Berlin: Springer. */
void getSHreal(/* Input arguments */
               int L,                             /* order of spherical harmonic expansion */
               float azi_rad,                     /* azimuth in radians */
               float incl_rad,                    /* pi/2-elevation (inclination) in radians */
               /* Output arguments */
               float* Y);                         /* the SH weights: (L+1)^2 x 1 */

/* Contructs a 3x3 rotation matrix from the Euler angles, using the yaw-pitch-roll (zyx) convention */
void yawPitchRoll2Rzyx (/* Input arguments */
                        float yaw,                /* yaw angle in radians */
                        float pitch,              /* pitch angle in radians */
                        float roll,               /* roll angle in radians */
                        int rollPitchYawFLAG,     /* 1: use Rxyz, i.e. apply roll, pitch and then yaw, 0: Rzyx / y-p-r */
                        /* Output arguments */
                        float R[3][3]);           /* zyx rotation matrix */
    
/* generates a real-valued spherical harmonic rotation matrix (assumes ACN/N3D convention)
 * For more information, the reader is referred to:
 * Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 * by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
void getSHrotMtxReal (float R[3][3],              /* zyx rotation matrix */
                      float* RotMtx,              /* the rotation matrix; FLAT: (L+1)^2 x (L+1)^2 */
                      int L);                     /* order */
    
/* generates beamforming weights for a direction on the sphere */
void calcBFweights(/* Input arguments */
                   BEAMFORMING_WEIGHT_TYPES BFW_type, /* see BEAMFORMING_WEIGHT_TYPES enum */
                   int order,                     /* order of spherical harmonic expansion */
                   float azi_rad,                 /* azimuth in radians */
                   float elev_rad,                /* elevation in radians */
                   /* Output arguments */
                   float* weights);               /* the resulting beamforming weights; (L+1)^2 x 1 */
    
/* converts spherical coordinates (with r=1) to cartesian coordinates of unit length */
void unitSph2Cart(/* Input arguments */
                  float azi_rad,                  /* azimuth in radians */
                  float elev_rad,                 /* elevation in radians */
                  /* Output arguments */
                  float xyz[3]);                  /* unit cartesian coords, xyz */

/* converts cartesian coordinates (unit length) to spherical coordinates (r=1) */
void unitCart2Sph(/* Input arguments */
                  float xyz[3],                   /* unit cartesian coords, xyz */
                  /* Output arguments */
                  float AziElev_rad[2]);          /* azimuth and elevation in radians */
   
/* converts cartesian coordinates (unit length) to spherical coordinates (r=1) */
void unitCart2Sph_aziElev(/* Input arguments */
                          float xyz[3],           /* unit cartesian coords, xyz */
                          /* Output arguments */
                          float* azi_rad,         /* & azimuth in radians */
                          float* elev_rad);       /* & elevation in radians */
    
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
    
/* generates a powermap utilising the subspace-based MUSIC method*/
void generateMUSICmap(/* Input arguments */
                      int order,                  /* analysis order */
                      float_complex* Cx,          /* covarience matrix; FLAT: (order+1)^2 x (order+1)^2 */
                      float_complex* Y_grid,      /* steering vectors for grid direcionts; FLAT: (order+1)^2 x nGrid_dirs  */
                      int nSources,               /* number of sources present in sound scene */
                      int nGrid_dirs,             /* number of grid directions */
                      int logScaleFlag,           /* 1: log(pmap), 0: pmap. */
                      /* Output arguments */
                      float* pmap);               /* resulting MUSIC pseudo-spectrum; nGrid_dirs x 1 */

/* generates a powermap utilising the subspace-based MinNorm method*/
void generateMinNormMap(/* Input arguments */
                        int order,                /* analysis order */
                        float_complex* Cx,        /* covarience matrix; FLAT: (order+1)^2 x (order+1)^2 */
                        float_complex* Y_grid,    /* steering vectors for grid direcionts; FLAT: (order+1)^2 x nGrid_dirs  */
                        int nSources,             /* number of sources present in sound scene */
                        int nGrid_dirs,           /* number of grid directions */
                        int logScaleFlag,         /* 1: log(pmap), 0: pmap. */
                        /* Output arguments */
                        float* pmap);             /* resulting MinNorm pseudo-spectrum; nGrid_dirs x 1 */

/* (cylindrical) Bessel function of the first kind: Jn
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_Jn(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical error) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               double* J_n,                       /* Bessel values; nZ x (N+1) */
               double* dJ_n);                     /* Bessel derivative values; nZ x (N+1) */
    
/* (cylindrical) Bessel function of the second kind: Yn
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_Yn(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical error) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               double* Y_n,                       /* Bessel values; nZ x (N+1) */
               double* dY_n);                     /* Bessel derivative values; nZ x (N+1) */
    
/* (cylindrical) Hankel function of the first kind: Hn1
 * returns the Hankel values and their derivatives up to order N for all values in vector z  */
void hankel_Hn1(/* Input arguments */
                int N,                            /* function order (highest is ~30 given numerical error) */
                double* z,                        /* input values; nZ x 1 */
                int nZ,                           /* number of input values */
                /* Output arguments */
                double_complex* Hn1_n,            /* Hankel values; nZ x (N+1) */
                double_complex* dHn1_n);          /* Hankel derivative values; nZ x (N+1) */
    
/* (cylindrical) Hankel function of the second kind: Hn2
 * returns the Hankel values and their derivatives up to order N for all values in vector z  */
void hankel_Hn2(/* Input arguments */
                int N,                            /* function order (highest is ~30 given numerical error) */
                double* z,                        /* input values; nZ x 1 */
                int nZ,                           /* number of input values */
                /* Output arguments */
                double_complex* Hn2_n,            /* Hankel values; nZ x (N+1) */
                double_complex* dHn2_n);          /* Hankel derivative values; nZ x (N+1) */
    
/* spherical Bessel function of the first kind: jn
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_jn(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical error) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               int* maxN,                         /* & maximum function order that could be computed <=N */
               double* j_n,                       /* Bessel values; nZ x (N+1) */
               double* dj_n);                     /* Bessel derivative values; nZ x (N+1) */
    
/* modified spherical Bessel function of the first kind: in
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_in(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical error) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               int* maxN,                         /* & maximum function order that could be computed <=N */
               double* i_n,                       /* Bessel values; nZ x (N+1) */
               double* di_n);                     /* Bessel derivative values; nZ x (N+1) */

/* spherical Bessel function of the second kind (Neumann): yn
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_yn(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical error) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               int* maxN,                         /* & maximum function order that could be computed <=N */
               double* y_n,                       /* Bessel values; nZ x (N+1) */
               double* dy_n);                     /* Bessel derivative values; nZ x (N+1) */
    
/* modified spherical Bessel function of the second kind: kn
 * returns the Bessel values and their derivatives up to order N for all values in vector z  */
void bessel_kn(/* Input arguments */
               int N,                             /* function order (highest is ~30 given numerical error) */
               double* z,                         /* input values; nZ x 1 */
               int nZ,                            /* number of input values */
               /* Output arguments */
               int* maxN,                         /* & maximum function order that could be computed <=N */
               double* k_n,                       /* Bessel values; nZ x (N+1) */
               double* dk_n);                     /* Bessel derivative values; nZ x (N+1) */

/* spherical Hankel function of the first kind: hn1
 * returns the Hankel values and their derivatives up to order N for all values in vector z */
void hankel_hn1(/* Input arguments */
                int N,                            /* function order (highest is ~30 given numerical error) */
                double* z,                        /* input values; nZ x 1 */
                int nZ,                           /* number of input values */
                /* Output arguments */
                int* maxN,                        /* & maximum function order that could be computed <=N */
                double_complex* h_n1,             /* Hankel values; nZ x (N+1) */
                double_complex* dh_n1);           /* Hankel derivative values; nZ x (N+1) */

/* spherical Hankel function of the second kind: hn2
 * returns the Hankel values and their derivatives up to order N for all values in vector z */
void hankel_hn2(/* Input arguments */
                int N,                            /* function order (highest is ~30 given numerical error) */
                double* z,                        /* input values; nZ x 1 */
                int nZ,                           /* number of input values */
                /* Output arguments */
                int* maxN,                        /* & maximum function order that could be computed <=N */
                double_complex* h_n2,             /* Hankel values; nZ x (N+1) */
                double_complex* dh_n2);           /* Hankel derivative values; nZ x (N+1) */
    
/* calculates the modal coefficients for open/rigid cylindrical arrays */
void cylModalCoeffs(/* Input arguments */
                    int order,                    /* max order (highest is ~30 given numerical error) */
                    double* kr,                   /* wavenumber*radius; nBands x 1 */
                    int nBands,                   /* number of frequency bands/bins */
                    ARRAY_CONSTRUCTION_TYPES arrayType, /* see 'ARRAY_CONSTRUCTION_TYPES' enum */
                    /* Output arguments */
                    double_complex* b_N);         /* modal coefficients per kr and 0:order; FLAT: nBands x (order+1) */
    
/* calculates the modal coefficients for open/rigid spherical arrays */
void sphModalCoeffs(/* Input arguments */
                    int order,                    /* max order (highest is ~30 given numerical error) */
                    double* kr,                   /* wavenumber*radius; nBands x 1 */
                    int nBands,                   /* number of frequency bands/bins */
                    ARRAY_CONSTRUCTION_TYPES arrayType, /* see 'ARRAY_CONSTRUCTION_TYPES' enum */
                    double dirCoeff,              /* only for directional (open) arrays, 0: omni, 0.5: card, 1:dipole */
                    /* Output arguments */
                    double_complex* b_N);         /* modal coefficients per kr and 0:order; FLAT: nBands x (order+1) */

/* simulates a cylindrical microphone array, returning the transfer functions for each (plane wave) source direction
 * on the surface of the cylinder. */
void simulateCylArray(/* Input arguments */
                      int order,                  /* max order (highest is ~30 given numerical error) */
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
                      int order,                  /* max order (highest is ~30 given numerical error) */
                      double* kr,                 /* wavenumber*radius; nBands x 1 */
                      int nBands,                 /* number of frequency bands/bins */
                      float* sensor_dirs_rad,     /* spherical coords of the sensors in RADIANS, [azi elev]; FLAT: N_sensors x 2 */
                      int N_sensors,              /* number of sensors */
                      float* src_dirs_deg,        /* spherical coords of the plane waves in DEGREES, [azi elev]; FLAT: N_srcs x 2 */
                      int N_srcs,                 /* number sources (DoAs of plane waves) */
                      ARRAY_CONSTRUCTION_TYPES arrayType, /* see 'ARRAY_CONSTRUCTION_TYPES' enum */
                      double dirCoeff,            /* only for directional (open) arrays, 0: omni, 0.5: card, 1:dipole */
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




