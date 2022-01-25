/*
 * Copyright 2018-2019 Leo McCormack
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
 *@addtogroup CDF4SAP
 *@{
 * @file saf_cdf4sap.h
 * @brief Main header for the Covariance Domain Framework module
 *        (#SAF_CDF4SAP_MODULE)
 *
 * Covariance Domain Framework for Spatial Audio Processing (CDF4SAP). This is a
 * direct C port of the MATLAB function given in [1], which was originally
 * written by Juha Vilkamo. It is explained in further detail in [1,2].
 *
 * @see [1] Vilkamo, J., Ba"ckstro"m, T., & Kuntz, A. (2013). Optimized
 *          covariance domain framework for time--frequency processing of
 *          spatial audio. Journal of the Audio Engineering Society, 61(6),
 *          403-411.
 * @see [2] Vilkamo, J., & Ba"ckstro"m, T. (2018). Time--Frequency Processing:
 *          Methods and Tools. In Parametric Time--Frequency Domain Spatial
 *          Audio. John Wiley & Sons.
 *
 * @author Leo McCormack
 * @date 25.11.2016
 * @license ISC
 */

#ifndef __SAF_CDF4SAP_H_INCLUDED__
#define __SAF_CDF4SAP_H_INCLUDED__

#include "../saf_utilities/saf_utility_complex.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the Covariance Domain Framework
 *
 * @note Use this function for real-valued input/output matrices. For
 *       complex-valued input/output matrices use cdf4sap_cmplx_create().
 *
 * @param[in] phCdf  The address (&) of the CDF4SAP handle
 * @param[in] nXcols Number of columns/rows in square input matrix 'Cx'
 * @param[in] nYcols Number of columns/rows in square input matrix 'Cy'
 */
void cdf4sap_create(/* Input Arguments */
                    void ** const phCdf,
                    int nXcols,
                    int nYcols);

/**
 * Creates an instance of the Covariance Domain Framework
 *
 * @note Use this function for complex-valued input/output matrices. For
 *       real-valued input/output matrices use cdf4sap_create().
 *
 * @param[in] phCdf  The address (&) of the CDF4SAP handle
 * @param[in] nXcols Number of columns/rows in square input matrix 'Cx'
 * @param[in] nYcols Number of columns/rows in square input matrix 'Cy'
 */
void cdf4sap_cmplx_create(/* Input Arguments */
                          void ** const phCdf,
                          int nXcols,
                          int nYcols);

/**
 * Destroys an instance of the Covariance Domain Framework
 *
 * @note Use this function for real-valued input/output matrices. For
 *       complex-valued input/output matrices use cdf4sap_cmplx_destroy().
 *
 * @param[in] phCdf The address (&) of the CDF4SAP handle
 */
void cdf4sap_destroy(/* Input Arguments */
                     void ** const phCdf);

/**
 * Destroys an instance of the Covariance Domain Framework
 *
 * @note Use this function for complex-valued input/output matrices. For
 *       real-valued input/output matrices use cdf4sap_destroy().
 *
 * @param[in] phCdf The address (&) of the CDF4SAP handle
 */
void cdf4sap_cmplx_destroy(/* Input Arguments */
                           void ** const phCdf);

/**
 * Computes the optimal mixing matrices
 *
 * This function solves the problem of determining the optimal mixing matrices
 * 'M' and 'Mr', such that the covariance matrix of the output:
 * y_out = M*x + Mr*decorrelated(x), is aligned with the target matrix 'Cy',
 * given the covariance matrix of input x, Cx=x*x^H, and a prototype mixing
 * matrix Q.
 * For the derivation and a more detailed description, the reader is referred
 * to [1,2].
 *
 * @note Use this function for real-valued input/output matrices. For
 *       complex-valued input/output use formulate_M_and_Cr_cmplx().
 *
 * @note For an example of how to use this function, one may refer to the
 *       implementation of the parametric binaural Ambisonic decoder (described
 *       in [3]) found here: https://github.com/leomccormack/CroPaC-Binaural,
 *       or the relevant unit tests.
 *
 * @test test__formulate_M_and_Cr()
 *
 * @param[in]  hCdf          Covariance Domain Framework handle
 * @param[in]  Cx            Covariance matrix of input 'x';
 *                           FLAT: nXcols x nXcols
 * @param[in]  Cy            Target covariance matrix; FLAT: nYcols x nYcols
 * @param[in]  Q             Prototype matrix; FLAT: nYcols x nXcols
 * @param[in]  useEnergyFLAG Set to '0' to apply energy compensation to 'M'
 *                           instead of outputing 'Cr'. Set to '1' to output 'Cr'
 * @param[in]  reg           Regularisation term (suggested: 0.2f)
 * @param[out] M             Mixing matrix; FLAT: nYcols x nXcols
 * @param[out] Cr            Mixing matrix residual, set to NULL if not needed;
 *                           FLAT: nYcols x nYcols
 *
 * @see [1] Vilkamo, J., Ba"ckstro"m, T., & Kuntz, A. (2013). Optimized
 *          covariance domain framework for time--frequency processing of
 *          spatial audio. Journal of the Audio Engineering Society, 61(6),
 *          403-411.
 * @see [2] Vilkamo, J., & Ba"ckstro"m, T. (2018). Time--Frequency Processing:
 *          Methods and Tools. In Parametric Time--Frequency Domain Spatial
 *          Audio. John Wiley & Sons.
 * @see [3] McCormack, L., Delikaris-Manias, S. (2019). "Parametric first-order
 *          ambisonic decoding for headphones utilising the Cross-Pattern
 *          Coherence algorithm". inProc 1st EAA Spatial Audio Signal Processing
 *          Symposium, Paris, France.
 */
void formulate_M_and_Cr(/* Input Arguments */
                        void * const hCdf,
                        float* Cx,
                        float* Cy,
                        float* Q,
                        int useEnergyFLAG,
                        float reg,
                        /* Output Arguments */
                        float* M,
                        float* Cr);

/**
 * Computes the optimal mixing matrices
 *
 * This function solves the problem of determining the optimal mixing matrices
 * 'M' and 'Mr', such that the covariance matrix of the output:
 * y_out = M*x + Mr*decorrelated(x), is aligned with the target matrix 'Cy',
 * given the covariance matrix of input x, Cx=x*x^H, and a prototype mixing
 * matrix Q.
 * For the derivation and a more detailed description, the reader is referred
 * to [1,2].
 *
 * @note Use this function for complex-valued input/output matrices. For
 *       real-valued input/output use formulate_M_and_Cr().
 *
 * @note For an example of how to use this function, one may refer to the
 *       implementation of the parametric binaural Ambisonic decoder (described
 *       in [3]) found here: https://github.com/leomccormack/CroPaC-Binaural,
 *       or the relevant unit tests.
 *
 * @test test__formulate_M_and_Cr_cmplx()
 *
 * @param[in]  hCdf          Covariance Domain Framework handle
 * @param[in]  Cx            Covariance matrix of input 'x';
 *                           FLAT: nXcols x nXcols
 * @param[in]  Cy            Target covariance matrix; FLAT: nYcols x nYcols
 * @param[in]  Q             Prototype matrix; FLAT: nYcols x nXcols
 * @param[in]  useEnergyFLAG Set to '0' to apply energy compensation to 'M'
 *                           instead of outputing 'Cr'. Set to '1' to output
 *                           'Cr'
 * @param[in]  reg           Regularisation term (suggested: 0.2f)
 * @param[out] M             Mixing matrix; FLAT: nYcols x nXcols
 * @param[out] Cr            Mixing matrix residual, set to NULL if not needed;
 *                           FLAT: nYcols x nYcols
 *
 * @see [1] Vilkamo, J., Ba"ckstro"m, T., & Kuntz, A. (2013). Optimized
 *          covariance domain framework for time--frequency processing of
 *          spatial audio. Journal of the Audio Engineering Society, 61(6),
 *          403-411.
 * @see [2] Vilkamo, J., & Ba"ckstro"m, T. (2018). Time--Frequency Processing:
 *          Methods and Tools. In Parametric Time--Frequency Domain Spatial
 *          Audio. John Wiley & Sons.
 * @see [3] McCormack, L., Delikaris-Manias, S. (2019). "Parametric first-order
 *          ambisonic decoding for headphones utilising the Cross-Pattern
 *          Coherence algorithm". inProc 1st EAA Spatial Audio Signal Processing
 *          Symposium, Paris, France.
 */
void formulate_M_and_Cr_cmplx(/* Input Arguments */
                              void * const hCdf,
                              float_complex* Cx,
                              float_complex* Cy,
                              float_complex* Q,
                              int useEnergyFLAG,
                              float reg,
                              /* Output Arguments */
                              float_complex* M,
                              float_complex* Cr);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif  /* __SAF_CDF4SAP_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup CDF4SAP */
