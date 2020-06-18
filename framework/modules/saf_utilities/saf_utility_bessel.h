/*
 * Copyright 2020 Leo McCormack
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
 * @file saf_utility_bessel.h
 * @brief Utility: A collection of routines for computing spherical and
 *        cylindrical Bessel and Hankel functions, including their derivatives
 * @author Leo McCormack
 * @date 26.05.2020
 */

#ifndef SAF_BESSEL_H_INCLUDED
#define SAF_BESSEL_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "saf_utility_complex.h"

/* ========================================================================== */
/*                        Cylindrical Bessel Functions                        */
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


/* ========================================================================== */
/*                         Spherical Bessel Functions                         */
/* ========================================================================== */

/**
 * Computes the spherical Bessel function of the first kind: jn
 *
 * Computes the Bessel values and their derivatives up to order N for all values
 * in vector z
 *
 * @note If the function fails to compute the function up to the specified
 *       order 'N', then the function will compute up to maximum order possible,
 *       and let the user know via the 'maxN' parameter. (i.e., always check if
 *       N=maxN, and handle things accordingly if maxN is lower).
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
 * in vector z.
 *
 * @note If the function fails to compute the function up to the specified
 *       order 'N', then the function will compute up to maximum order possible,
 *       and let the user know via the 'maxN' parameter. (i.e., always check if
 *       N=maxN, and handle things accordingly if maxN is lower).
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
 * @note If the function fails to compute the function up to the specified
 *       order 'N', then the function will compute up to maximum order possible,
 *       and let the user know via the 'maxN' parameter. (i.e., always check if
 *       N=maxN, and handle things accordingly if maxN is lower).
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
 * @note If the function fails to compute the function up to the specified
 *       order 'N', then the function will compute up to maximum order possible,
 *       and let the user know via the 'maxN' parameter. (i.e., always check if
 *       N=maxN, and handle things accordingly if maxN is lower).
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
 * in vector z.
 *
 * @note If the function fails to compute the function up to the specified
 *       order 'N', then the function will compute up to maximum order possible,
 *       and let the user know via the 'maxN' parameter. (i.e., always check if
 *       N=maxN, and handle things accordingly if maxN is lower).
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
 * @note If the function fails to compute the function up to the specified
 *       order 'N', then the function will compute up to maximum order possible,
 *       and let the user know via the 'maxN' parameter. (i.e., always check if
 *       N=maxN, and handle things accordingly if maxN is lower).
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


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_BESSEL_H_INCLUDED */
