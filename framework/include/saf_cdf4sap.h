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
 *     saf_cdf4sap.h (include header)
 * Description:
 *     Covariance Domain Framework for Spatial Audio Processing (CDF4SAP). A C implementation
 *     of the MatLab code, originally written by Dr. Juha Vilkamo.
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 25.11.2016
 */

#ifndef __SAF_CDF4SAP_H_INCLUDED__
#define __SAF_CDF4SAP_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h> 
#include "saf_utilities.h"   /* for blas/lapack and complex number support */

#ifdef __cplusplus
extern "C" {
#endif
    
/* creates an instance of the Covariance Domain Framework, REAL-VALUED MATRICES */
void cdf4sap_create(void ** const phCdf,              /* address of Covariance Domain Framework handle */
                    int nXcols,                       /* number of columns/rows in Cx */
                    int nYcols);                      /* number of columns/rows in Cy */
    
/* creates an instance of the Covariance Domain Framework, COMPLEX-VALUED MATRICES */
void cdf4sap_cmplx_create(void ** const phCdf,        /* address of Covariance Domain Framework handle */
                          int nXcols,                 /* number of columns/rows in Cx */
                          int nYcols);                /* number of columns/rows in Cy */
    
/* destroys an instance of the Covariance Domain Framework, for REAL-VALUED MATRICES */
void cdf4sap_destroy(void ** const phCdf);            /* address of Covariance Domain Framework handle */
    
/* destroys an instance of the Covariance Domain Framework, for COMPLEX-VALUED MATRICES */
void cdf4sap_cmplx_destroy(void ** const phCdf);      /* address of Covariance Domain Framework handle */
    
/* This function solves the problem of determining the optimal mixing matrices 'M' and 'Cr',
 * such that the covariance matrix of the output: y_out = M*x + B*decorrelated(x), is aligned
 * with the target matrix 'Cy', given the covariance matrix of input x, Cx=x*x^H, and a prototype
 * decoding matrix Q.
 *
 * For the derivation and a more detailed description, the reader is referred to:
 *     Vilkamo, J., Bäckström, T., & Kuntz, A. (2013). Optimized covariance domain framework
 *     for time–frequency processing of spatial audio. Journal of the Audio Engineering Society,
 *     61(6), 403-411.
 * and:
 *     Vilkamo, J., & Backstrom, T. (2018). Time--Frequency Processing: Methods and Tools. In
 *     Parametric Time-Frequency Domain Spatial Audio. John Wiley & Sons. */
void formulate_M_and_Cr(void * const hCdf,            /* Covariance Domain Framework handle */
                        float* Cx,                    /* covariance matrix of input, x; FLAT: nXcols x nXcols */
                        float* Cy,                    /* target covariance matrix; FLAT: nYcols x nYcols */
                        float* Q,                     /* prototype matrix; FLAT: nYcols x nXcols */
                        int useEnergyFLAG,            /* 0: apply energy compensation to 'M' instead of outputing 'Cr', 1: output 'Cr' too */
                        float reg,                    /* regularisation term */
                        float* M,                     /* Mixing matrix; FLAT: nYcols x nXcols */
                        float* Cr);                   /* Mixing matrix for residual; FLAT: nYcols x nYcols */

/* This is the same as 'formulate_M_and_Cr', except it employs COMPLEX-VALUED MATRICES */
void formulate_M_and_Cr_cmplx(void * const hCdf,      /* Covariance Domain Framework handle */
                              float_complex* Cx,      /* covariance matrix of input, x; FLAT: nXcols x nXcols */
                              float_complex* Cy,      /* target covariance matrix; FLAT: nYcols x nYcols */
                              float_complex* Q,       /* prototype matrix; FLAT: nYcols x nXcols */
                              int useEnergyFLAG,      /* 0: apply energy compensation to 'M' instead of outputing 'Cr', 1: output 'Cr' too */
                              float reg,              /* regularisation term */
                              float_complex* M,       /* Mixing matrix; FLAT: nYcols x nXcols */
                              float_complex* Cr);     /* Mixing matrix for residual; FLAT: nYcols x nYcols */
    

#ifdef __cplusplus
}
#endif


#endif  /* __SAF_CDF4SAP_H_INCLUDED__ */
