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
 * @file saf_sh_internal.h
 * @ingroup SH
 * @brief Internal header for the Spherical Harmonic Transform and Spherical
 *        Array Processing module (#SAF_SH_MODULE)
 *
 * A collection of spherical harmonic related functions. Many of which have been
 * derived from the MATLAB libraries found in [1-3].
 *
 * @see [1] https://github.com/polarch/Spherical-Harmonic-Transform
 *          Copyright (c) 2015, Archontis Politis, BSD-3-Clause License
 * @see [2] https://github.com/polarch/Array-Response-Simulator
 *          Copyright (c) 2015, Archontis Politis, BSD-3-Clause License
 * @see [3] https://github.com/polarch/Spherical-Array-Processing
 *          Copyright (c) 2016, Archontis Politis, BSD-3-Clause License
 *
 * @author Leo McCormack
 * @date 22.05.2016
 * @license ISC
 */

#ifndef __SAF_SH_INTERNAL_H_INCLUDED__
#define __SAF_SH_INTERNAL_H_INCLUDED__

#include "saf_sh.h" 
#include "../saf_utilities/saf_utilities.h" /* for linear algebra speed-ups */
#include "saf_externals.h" 

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/** Internal data structure for sphPWD */
typedef struct _sphPWD_data {
    int order, nSH, nDirs;
    float_complex* grid_svecs;
    float_complex* A_Cx; 
    float* grid_dirs_xyz; 
    float* pSpec;
    float* P_minus_peak;
    float* P_tmp;
    float* VM_mask;

}sphPWD_data;

/** Internal data structure for sphMUSIC */
typedef struct _sphMUSIC_data {
    int order, nSH, nDirs;
    float_complex* grid_svecs;
    float_complex* VnA;
    float* grid_dirs_xyz;
    float* abs_VnA;
    float* pSpec;
    float* pSpecInv;
    float* P_minus_peak;
    float* P_tmp;
    float* VM_mask;

}sphMUSIC_data;

/** Internal data structure for sphESPRIT */
typedef struct _sphESPRIT_data {
    int N, NN, maxK;

    /* matrices and indices to be pre-computed */
    double* rWVnimu[6];             /* [6] x FLAT: N^2 x N^2 */
    double_complex* WVnimu[6];      /* [6] x FLAT: N^2 x N^2 */
    int* idx_from_Ynm2Ynimu[12];    /* [12] x nIdx[i] x 1 */
    int nIdx[12];                   /* number of indices in each "idx_from_Ynm2Ynimu" vector */

    /* matrices for run-time */
    void* hZpinv, *hZeigmp, *hZglslv;
    double_complex* Us_1m1, *Us_m1m1, *Us_11, *Us_m11, *Us_m10, *Us_10, *Us_00;
    double_complex* WVnimu0_Us1m1, *WVnimu1_Usm1m1, *WVnimu2_Us11;
    double_complex* WVnimu3_Usm11, *WVnimu4_Usm10,  *WVnimu5_Us10;
    double_complex* LambdaXYp, *LambdaXYm, *LambdaZ;
    double_complex* pinvUs, *PsiXYp, *PsiXYm, *PsiZ;
    double_complex* tmp_KK, *V, *PhiXYp, *PhiXYm, *PhiZ;

}sphESPRIT_data;


/* ========================================================================== */
/*                          Misc. Internal Functions                          */
/* ========================================================================== */

/**
 * Computes the Wigner 3j symbol through the Racah formula found in
 * http://mathworld.wolfram.com/Wigner3j-Symbol.html, Eq.7
 *
 * @param[in] j1 Wigner 3 j-symbol, j1
 * @param[in] j2 Wigner 3 j-symbol, j2
 * @param[in] j3 Wigner 3 j-symbol, j3
 * @param[in] m1 Wigner 3 j-symbol, m1
 * @param[in] m2 Wigner 3 j-symbol, m2
 * @param[in] m3 Wigner 3 j-symbol, m3
 * @returns wigner_3j symbol
 */
float wigner_3j(/* Input arguments */
                int j1,
                int j2,
                int j3,
                int m1,
                int m2,
                int m3);

/**
 * Constructs a matrix of Guant coefficients
 *
 * Constructs the (N1+1)^2 x (N2+1)^2 x (N+1)^2 matrix of Gaunt coefficients
 * which represent the integral of three spherical harmonics, such as
 * G^q_{q',q''} = int_Omega Y_{q'}Y_{q''}Y^*_{q} mathrm{d} Omega.
 *
 * With Gaunt coefficients, the spherical harmonic coefficients of the product
 * of two spherical functions can be given directly as a linear relationship
 * between the harmonic coefficients of the two functions.
 *
 * @param[in]  N1 Order of first harmonic coeffient
 * @param[in]  N2 Order of second harmonic coefficient
 * @param[in]  N  Target order
 * @param[out] A  Gaunt matrix; FLAT: (N1+1)^2 x (N2+1)^2 x (N+1)^2
 */
void gaunt_mtx(/* Input arguments */
               int N1,
               int N2,
               int N,
               /* Output arguments */
               float* A);


/* ========================================================================== */
/*             Internal functions for spherical harmonic rotations            */
/* ========================================================================== */

/** Helper function for getSHrotMtxReal() */
float getP(int M, int i, int l, int a, int b, float R_1[3][3], float* R_lm1);

/** Helper function for getSHrotMtxReal() */
float getU(int M, int l, int m, int n, float R_1[3][3], float* R_lm1);

/** Helper function for getSHrotMtxReal() */
float getV(int M, int l, int m, int n, float R_1[3][3], float* R_lm1);

/** Helper function for getSHrotMtxReal() */
float getW(int M, int l, int m, int n, float R_1[3][3], float* R_lm1);

/** Helper function for getSHrotMtxReal() */
float getW(int M, int l, int m, int n, float R_1[3][3], float* R_lm1);


/* ========================================================================== */
/*                        Internal functions for sphESPRIT                    */
/* ========================================================================== */

/** Helper function for sphESPRIT_create() */
void getWnimu(int order, int mm, int ni, int mu, double* Wnimu);

/** Helper function for sphESPRIT_create() */
void getVnimu(int order, int ni, int mu, double* Vnimu);

/** Helper function for sphESPRIT_create() */
void muni2q(int order, int ni, int mu, int* idx_nimu, int* idx_nm);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_SH_INTERNAL_H_INCLUDED__ */
