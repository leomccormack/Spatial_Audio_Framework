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

/*
 * Filename: saf_cdf4sap.h (include header)
 * ----------------------------------------
 * Covariance Domain Framework for Spatial Audio Processing (CDF4SAP) [1,2]. A C
 * implementation of the MatLab code, originally written by Dr. Juha Vilkamo.
 *
 * [1] Vilkamo, J., Bäckström, T., & Kuntz, A. (2013). Optimized covariance
 *     domain framework for time–frequency processing of spatial audio. Journal
 *     of the Audio Engineering Society, 61(6), 403-411.
 * [2] Vilkamo, J., & Backstrom, T. (2018). Time--Frequency Processing: Methods
 *     and Tools. In Parametric Time-Frequency Domain Spatial Audio. John Wiley
 *     & Sons.
 *
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 25.11.2016
 */

#ifndef __SAF_CDF4SAP_H_INCLUDED__
#define __SAF_CDF4SAP_H_INCLUDED__

#include "../saf_utilities/saf_complex.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: cdf4sap_create
 * ------------------------
 * Creates an instance of the Covariance Domain Framework.
 * Note: Use this function if using REAL-VALUED input/output matrices.
 *
 * Input Arguments:
 *     phCdf  - address of Covariance Domain Framework handle
 *     nXcols - number of columns/rows in square input matrix 'Cx'
 *     nYcols - number of columns/rows in square input matrix 'Cy'
 */
void cdf4sap_create(/* Input Arguments */
                    void ** const phCdf,
                    int nXcols,
                    int nYcols);
    
/*
 * Function: cdf4sap_cmplx_create
 * ------------------------------
 * Creates an instance of the Covariance Domain Framework.
 * Note: Use this function if using COMPLEX-VALUED input/output matrices.
 *
 * Input Arguments:
 *     phCdf  - address of Covariance Domain Framework handle
 *     nXcols - number of columns/rows in square input matrix 'Cx'
 *     nYcols - number of columns/rows in square input matrix 'Cy'
 */
void cdf4sap_cmplx_create(/* Input Arguments */
                          void ** const phCdf,
                          int nXcols,
                          int nYcols);
    
/*
 * Function: cdf4sap_destroy
 * ------------------------
 * Destroys an instance of the Covariance Domain Framework.
 * Note: Use this function if using REAL-VALUED input/output matrices.
 *
 * Input Arguments:
 *     phCdf - address of Covariance Domain Framework handle
 */
void cdf4sap_destroy(/* Input Arguments */
                     void ** const phCdf);
    
/*
 * Function: cdf4sap_cmplx_destroy
 * -------------------------------
 * Destroys an instance of the Covariance Domain Framework.
 * Note: Use this function if using COMPLEX-VALUED input/output matrices.
 *
 * Input Arguments:
 *     phCdf - address of Covariance Domain Framework handle
 */
void cdf4sap_cmplx_destroy(/* Input Arguments */
                           void ** const phCdf);
    
/*
 * Function: formulate_M_and_Cr
 * ----------------------------
 * This function solves the problem of determining the optimal mixing matrices
 * 'M' and 'Cr', such that the covariance matrix of the output:
 * y_out = M*x + Mr*decorrelated(x), is aligned with the target matrix 'Cy',
 * given the covariance matrix of input x, Cx=x*x^H, and a prototype decoding
 * matrix Q.
 * For the derivation and a more detailed description, the reader is referred
 * to [1,2].
 * Note: Use this function if using REAL-VALUED input/output matrices.
 *
 * Input Arguments:
 *     hCdf          - Covariance Domain Framework handle
 *     Cx            - Covariance matrix of input 'x'; FLAT: nXcols x nXcols
 *     Cy            - Target covariance matrix; FLAT: nXcols x nXcols
 *     Q             - Prototype matrix; FLAT: nYcols x nXcols
 *     useEnergyFLAG - 0: Apply energy compensation to 'M' instead of outputing
 *                     'Cr', 1: output 'Cr' too
 *     reg           - regularisation term (suggested: 0.2f)
 * Output Arguments:
 *     M             - Mixing matrix; FLAT: nYcols x nXcols
 *     Cr            - Mixing matrix for residual; FLAT: nYcols x nYcols
 *
 * [1] Vilkamo, J., Bäckström, T., & Kuntz, A. (2013). Optimized covariance
 *     domain framework for time–frequency processing of spatial audio. Journal
 *     of the Audio Engineering Society, 61(6), 403-411.
 * [2] Vilkamo, J., & Backstrom, T. (2018). Time--Frequency Processing: Methods
 *     and Tools. In Parametric Time-Frequency Domain Spatial Audio. John Wiley
 *     & Sons.
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

/*
 * Function: formulate_M_and_Cr_cmplx
 * ----------------------------------
 * This function solves the problem of determining the optimal mixing matrices
 * 'M' and 'Cr', such that the covariance matrix of the output:
 * y_out = M*x + B*decorrelated(x), is aligned with the target matrix 'Cy',
 * given the covariance matrix of input x, Cx=x*x^H, and a prototype decoding
 * matrix Q.
 * For the derivation and a more detailed description, the reader is referred
 * to [1,2].
 * Note: Use this function if using COMPLEX-VALUED input/output matrices.
 *
 * Input Arguments:
 *     hCdf          - Covariance Domain Framework handle
 *     Cx            - Covariance matrix of input 'x'; FLAT: nXcols x nXcols
 *     Cy            - Target covariance matrix; FLAT: nXcols x nXcols
 *     Q             - Prototype matrix; FLAT: nYcols x nXcols
 *     useEnergyFLAG - 0: Apply energy compensation to 'M' instead of outputing
 *                     'Cr', 1: output 'Cr' too
 *     reg           - regularisation term (suggested: 0.2f)
 * Output Arguments:
 *     M             - Mixing matrix; FLAT: nYcols x nXcols
 *     Cr            - Mixing matrix for residual; FLAT: nYcols x nYcols
 *
 * [1] Vilkamo, J., Bäckström, T., & Kuntz, A. (2013). Optimized covariance
 *     domain framework for time–frequency processing of spatial audio. Journal
 *     of the Audio Engineering Society, 61(6), 403-411.
 * [2] Vilkamo, J., & Backstrom, T. (2018). Time--Frequency Processing: Methods
 *     and Tools. In Parametric Time-Frequency Domain Spatial Audio. John Wiley
 *     & Sons.
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
                              float* Cr);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif  /* __SAF_CDF4SAP_H_INCLUDED__ */


