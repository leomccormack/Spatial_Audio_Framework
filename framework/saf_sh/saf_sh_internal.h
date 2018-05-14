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
 *     A collection of spherical harmonic related functions. Some of which have been
 *     derived from the Matlab library by Archontis Politis; found here:
 *     https://github.com/polarch/Spherical-Harmonic-Transform
 *     and MATLAB code by Symeon Delikaris-Manias
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
 
/* Calculates Chebyshev Polynomial Coefficients */
void ChebyshevPolyCoeff (int n,              /* order of spherical harmonic expansion */
                         float* t_coeff);    /* resulting Chebyshev Polynomial Coefficients */

/* Calculates Legendre Polynomial Coefficients */
void LegendrePolyCoeff(int n,                /* order of spherical harmonic expansion */
                       float* p_coeff);      /* resulting Legendre Polynomial Coefficients */

/* Calculates Dolph-chebyshev weights */
void dolph_chebyshev(int M,                  /* order of spherical harmonic expansion */
                     float* d,               /* resulting weights */
                     int type );             /* 0: 1: */
    
/* Calculates max_rE weights */
void maxre3d(int M,                          /* order of spherical harmonic expansion */
             float* gm );                    /* resulting weights */

#ifdef __cplusplus
}
#endif

#endif /* __SHT_INTERNAL_H_INCLUDED__ */




















