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
 *     of the matlab code by J. Vilkamo.
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

typedef struct _cdf4sap {
    int xsize, ysize, qsize1, qsize2;
    float* S_Cy, *S_Cy_Mtx, *U_Cy, *VT_Cy, *A_Cy, *Ky;
    float* S_Cx, *S_Cx_Mtx, *U_Cx, *VT_Cx, *A_Cx, *Kx, *Sx_reg_diag, *Sx_reg_invMtx, *Kx_reg_inverse;
    float* Q, *QCx, *QCxQ, *KxQ, *KxQG, *KxQGKy, *G_hat;
    float* S, *U, *VT, *Vlamb, *P, *KyP, *M_CM, *MCx, *Cy_tilde;
    float* G, *Mtmp;
    float* lambda;

}cdf4sap;

    
typedef struct _cdf4sap_cmplx {
    int xsize, ysize, qsize1, qsize2;
    float* S_Cy_real, *S_Cx_real, *S_real;
    float_complex* S_Cy, *S_Cy_Mtx, *U_Cy, *VH_Cy, *A_Cy, *Ky;
    float_complex* S_Cx, *S_Cx_Mtx, *U_Cx, *VH_Cx, *A_Cx, *Kx, *Sx_reg_diag, *Sx_reg_invMtx, *Kx_reg_inverse;
    float_complex* Q, *QCx, *QCxQ, *KxQ, *KxQG, *KxQGKy, *G_hat;
    float_complex* S, *U, *VH, *Vlamb, *P, *KyP, *M_CM, *MCx, *Cy_tilde;
    float_complex* G, *Mtmp;
    float_complex* lambda;

}cdf4sap_cmplx;


/* creates an instance of the Covariance Domain Framework */
void cdf4sap_alloc(void ** const phCdf,                             /* address of Covariance Domain Framework handle */
                   int CxMN,                                        /* dimension m=n length of Cx */
                   int CyMN,                                        /* dimension m=n length of Cy */
                   int QmM,                                         /* dimension m of Q */
                   int QmN);                                        /* dimension n of Q */
    
    
/* creates an instance of the Covariance Domain Framework */
void cdf4sap_alloc_cmplx(void ** const phCdf,                       /* address of Covariance Domain Framework handle */
                         int CxMN,                                  /* dimension m=n length of Cx */
                         int CyMN,                                  /* dimension m=n length of Cy */
                         int QmM,                                   /* dimension m of Q */
                         int QmN);                                  /* dimension n of Q */
    
    
/* destroys an instance of the Covariance Domain Framework */
void cdf4sap_free(void ** const phCdf);                             /* address of Covariance Domain Framework handle */
    
    
/* destroys an instance of the Covariance Domain Framework */
void cdf4sap_free_cmplx(void ** const phCdf);                       /* address of Covariance Domain Framework handle */
    

/* Formulate mixing and residual matrices
 * note: input arguments must be in row-major order, despite the processing using column-major;
 * input matricies do not need to be contiguously allocated for this
 *
 * For the derivation and a detailed description, the reader is referred to:
 * Vilkamo, J., Bäckström, T., & Kuntz, A. (2013). Optimized covariance domain
 * framework for time–frequency processing of spatial audio. Journal of the Audio
 * Engineering Society, 61(6), 403-411. */
void cdf4sap_formulate_M_and_Cr(/* Input arguments */
                                void * const hCdf,                  /* Covariance Domain Framework handle */
                                float** Cx,                         /* [xsize][xsize], Input Covarience matrix */
                                float** Cy,                         /* [ysize][ysize], Target Covarience matrix */
                                float** Qm,                         /* [qsize1][qsize2], decoding matrix */
                                int flag,                           /* 0 or 1; 1=Use energy compensation instead of residual */
                                float regularisationConstant,       /* suggested value, 0.2f */
                                /* Output arguments */
                                float** M,                          /* [xsize][ysize], Output Mixing matrix */
                                float** Cr);                        /* [ysize][ysize], Output Residual matrix */


/* Formulate mixing and residual matrices
 * note: input arguments must be in row-major order, despite the processing using column-major
 *
 * For the derivation and a detailed description, the reader is referred to:
 * Vilkamo, J., Bäckström, T., & Kuntz, A. (2013). Optimized covariance domain
 * framework for time–frequency processing of spatial audio. Journal of the Audio
 * Engineering Society, 61(6), 403-411. */
void cdf4sap_formulate_M_and_Cr_cmplx(/* Input arguments */
                                      void * const hCdf,            /* Covariance Domain Framework handle */
                                      float_complex** Cx,           /* [xsize][xsize], Input Covarience matrix */
                                      float_complex** Cy,           /* [ysize][ysize], Target Covarience matrix */
                                      float_complex** Qm,           /* [qsize1][qsize2], decoding matrix */
                                      int flag,                     /* 0 or 1; 1=Use energy compensation instead of residual */
                                      float regularisationConstant, /* suggested value, 0.2f */
                                      /* Output arguments */
                                      float_complex** M,            /* [xsize][ysize], Output Mixing matrix */
                                      float_complex** Cr);          /* [ysize][ysize], Output Residual matrix */

#ifdef __cplusplus
}
#endif


#endif  /* __SAF_CDF4SAP_H_INCLUDED__ */
