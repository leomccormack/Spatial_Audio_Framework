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
 *     saf_sh_internal.h
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

#ifndef __SHT_INTERNAL_H_INCLUDED__
#define __SHT_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <assert.h>
#include "saf_sh.h" 

#ifdef __cplusplus
extern "C" {
#endif
 
    
#ifndef M_PI
  #define M_PI ( 3.14159265359f )
#endif
    
    
/****************************/
/* Misc. Internal Functions */
/****************************/
    
/* classic factorial algorithm */
unsigned long factorial(unsigned long f);
    
/* computes the Wigner 3j symbol through the Racah formula found in
 * http://mathworld.wolfram.com/Wigner3j-Symbol.html, Eq.7. */
float wigner_3j(/* Input arguments */
                int j1,                           /* Wigner 3 j-symbol, j1 */
                int j2,                           /* Wigner 3 j-symbol, j2 */
                int j3,                           /* Wigner 3 j-symbol, j3 */
                int m1,                           /* Wigner 3 j-symbol, m1 */
                int m2,                           /* Wigner 3 j-symbol, m2 */
                int m3);                          /* Wigner 3 j-symbol, m3 */

/* constructs the (N1+1)^2x(N2+1)^2x(N+1)^2 matrix of Gaunt coefficients which represent the integral of three
 * spherical harmonics such as:
 * G^q_{q',q''} = \int_\Omega Y_{q'}Y_{q''}Y^*_{q} \mathrm{d}\Omega.
 *
 * With Gaunt coefficients, the spherical harmonic coefficients of the product of two spherical functions can
 * be given directly as a linear relationship between the harmonic coefficients of the two functions. */
void gaunt_mtx(/* Input arguments */
               int N1,                            /* order of first harmonic coeffient */
               int N2,                            /* order of second harmonic coefficient */
               int N,                             /* target order */
               /* Output arguments */
               float* A);                         /* gaunt matrix; flat: (N1+1)^2 x (N2+1)^2 x (N+1)^2 */
    
    
/****************************************************/
/* Internal Functions for spherical Hankels/Bessels */
/****************************************************/
    
/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
void SPHI(int N, double X, int *NM, double *SI, double *DI);
    
/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
void SPHK(int N, double X, int *NM, double *SK, double *DK);

/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
void SPHJ(int N, double X, int *NM, double *SJ, double *DJ);

/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
int MSTA1(double X, int MP);
    
/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
int MSTA2(double X, int N, int MP);
    
/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
double ENVJ(int N, double X);

/* Original Fortran code: "Fortran Routines for Computation of Special Functions":
 * jin.ece.uiuc.edu/routines/routines.html.
 * C implementation by J-P Moreau, Paris (www.jpmoreau.fr) */
void SPHY(int N, double X, int *NM, double *SY, double *DY);
    
    
/*******************************************************/
/* Internal Functions for spherical harmonic rotations */
/*******************************************************/
    
/* Used in the calculation of spherical harmonic rotation matrices
 * Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 * by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
float getP(int i, int l, int a, int b, float** R_1, float** R_lm1);

/* Used in the calculation of spherical harmonic rotation matrices
 * Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 * by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
float getU(int l, int m, int n, float** R_1, float** R_lm1);

/* Used in the calculation of spherical harmonic rotation matrices
 * Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 * by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
float getV(int l, int m, int n, float** R_1, float** R_lm1);

/* Used in the calculation of spherical harmonic rotation matrices
 * Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 * by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
float getW(int l, int m, int n, float** R_1, float** R_lm1);

/* Used in the calculation of spherical harmonic rotation matrices
 * Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real Spherical Harmonics. Direct Determination
 * by Recursion Page: Additions and Corrections. Journal of Physical Chemistry A, 102(45), 9099?9100. */
float getW(int l, int m, int n, float** R_1, float** R_lm1);
    
    
#ifdef __cplusplus
}
#endif

#endif /* __SHT_INTERNAL_H_INCLUDED__ */




















