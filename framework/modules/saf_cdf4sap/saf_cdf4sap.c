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
 * @file saf_cdf4sap.c
 * @ingroup CDF4SAP
 * @brief Public source for the Covariance Domain Framework module
 *        (#SAF_CDF4SAP_MODULE)
 *
 * Covariance Domain Framework for Spatial Audio Processing (CDF4SAP). This is a
 * direct C port of the Matlab function given in [1], which was originally
 * written by Juha Vilkamo. The algorithm is explained in further detail in [2].
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

#include "saf_cdf4sap.h"
#include "saf_externals.h" 
#include "../saf_utilities/saf_utilities.h"

/**
 * Main data structure for the Covariance Domain Framework for Spatial Audio
 * Processing (CDF4SAP), for real-valued matrices.
 */
typedef struct _cdf4sap_data {
    /* Dimensions of Cx and Cy */
    int nXcols, nYcols;
    
    /* intermediate vectors & matrices */
    void* hSVD;
    float* lambda, *U_Cy, *S_Cy, *Ky, *U_Cx, *S_Cx, *s_Cx, *Kx, *Kx_reg_inverse, *U, *V, *P;
    float* G_hat, *Cx_QH;
    float* GhatH_Ky, *QH_GhatH_Ky, *KxH_QH_GhatH_Ky, *lambda_UH;
    float* P_Kxreginverse;
    float* Cx_MH, *Cy_tilde;
    float* G_M;
    
}cdf4sap_data;

/**
* Main data structure for the Covariance Domain Framework for Spatial Audio
* Processing (CDF4SAP), for complex-valued matrices.
*/
typedef struct _cdf4sap_cmplx_data {
    /* Dimensions of Cx and Cy */
    int nXcols, nYcols;
    
    /* intermediate vectors & matrices */
    void* hSVD;
    float_complex* Cr_cmplx;
    float_complex* lambda, *U_Cy, *S_Cy, *S_Cx, *Ky, *U_Cx, *Kx, *Kx_reg_inverse, *U, *V, *P;
    float* s_Cx, *G_hat_diag;
    float_complex* G_hat, *Cx_QH;
    float_complex* GhatH_Ky, *QH_GhatH_Ky, *KxH_QH_GhatH_Ky, *lambda_UH;
    float_complex *P_Kxreginverse;
    float_complex *Cx_MH, *Cy_tilde;
    float_complex* G_M;
    
}cdf4sap_cmplx_data;

void cdf4sap_create
(
    void ** const phCdf,
    int nXcols,
    int nYcols
)
{
    *phCdf = malloc1d(sizeof(cdf4sap_data));
    cdf4sap_data *h = (cdf4sap_data*)(*phCdf);
    
    h->nXcols = nXcols;
    h->nYcols = nYcols;
    h->lambda = malloc1d(nYcols * nXcols * sizeof(float));

    /* For the SVD */
    utility_ssvd_create(&h->hSVD, SAF_MAX(nXcols, nYcols), SAF_MAX(nXcols, nYcols));
    
    /* For the decomposition of Cy */
    h->U_Cy = malloc1d(nYcols*nYcols*sizeof(float));
    h->S_Cy = malloc1d(nYcols*nYcols*sizeof(float));
    h->Ky = malloc1d(nYcols*nYcols*sizeof(float));
    
    /* For the decomposition of Cx */
    h->U_Cx = malloc1d(nXcols*nXcols*sizeof(float));
    h->S_Cx = malloc1d(nXcols*nXcols*sizeof(float));
    h->s_Cx = malloc1d(nXcols*sizeof(float));
    h->Kx = malloc1d(nXcols*nXcols*sizeof(float));
    
    /* For the formulation of regularised Kx^-1 */
    h->Kx_reg_inverse = malloc1d(nXcols*nXcols*sizeof(float));
    
    /* For the formulation of normalisation matrix G_hat */
    h->G_hat = malloc1d(nYcols*nYcols*sizeof(float));
    h->Cx_QH = malloc1d(nXcols*nYcols*sizeof(float));
    
    /* For the formulation of optimal P */
    h->GhatH_Ky = malloc1d(nYcols*nYcols*sizeof(float));
    h->QH_GhatH_Ky = malloc1d(nXcols*nYcols*sizeof(float));
    h->KxH_QH_GhatH_Ky = malloc1d(nXcols*nYcols*sizeof(float));
    h->U = malloc1d(nXcols*nXcols*sizeof(float));
    h->V = malloc1d(nYcols*nYcols*sizeof(float));
    h->lambda_UH = malloc1d(nYcols*nXcols*sizeof(float));
    h->P = malloc1d(nYcols*nXcols*sizeof(float));
    
    /* For the formulation of M */
    h->P_Kxreginverse = malloc1d(nYcols*nXcols*sizeof(float));
    
    /* For the formulation of the residual covariance matrix */
    h->Cx_MH = malloc1d(nXcols*nYcols*sizeof(float));
    h->Cy_tilde = malloc1d(nYcols*nYcols*sizeof(float));
    
    /* For using energy compensation instead of residuals */
    h->G_M = malloc1d(nYcols*nXcols*sizeof(float));
}

void cdf4sap_cmplx_create
(
    void ** const phCdf,
    int nXcols,
    int nYcols
)
{
    *phCdf = malloc1d(sizeof(cdf4sap_cmplx_data));
    cdf4sap_cmplx_data *h = (cdf4sap_cmplx_data*)(*phCdf);
    
    h->nXcols = nXcols;
    h->nYcols = nYcols;
    h->lambda = malloc1d(nYcols * nXcols * sizeof(float_complex));
    h->Cr_cmplx = malloc1d(nYcols * nYcols * sizeof(float_complex));

    /* For the SVD */
    utility_csvd_create(&h->hSVD, SAF_MAX(nXcols, nYcols), SAF_MAX(nXcols, nYcols));

    /* For the decomposition of Cy */
    h->U_Cy = malloc1d(nYcols*nYcols*sizeof(float_complex));
    h->S_Cy = malloc1d(nYcols*nYcols*sizeof(float_complex));
    h->Ky = malloc1d(nYcols*nYcols*sizeof(float_complex));
    
    /* For the decomposition of Cx */
    h->U_Cx = malloc1d(nXcols*nXcols*sizeof(float_complex));
    h->S_Cx = malloc1d(nXcols*nXcols*sizeof(float_complex));
    h->s_Cx = malloc1d(nXcols*sizeof(float));
    h->Kx = malloc1d(nXcols*nXcols*sizeof(float_complex));
    
    /* For the formulation of regularised Kx^-1 */
    h->Kx_reg_inverse = malloc1d(nXcols*nXcols*sizeof(float_complex));
    
    /* For the formulation of normalisation matrix G_hat */
    h->G_hat_diag = malloc1d(nYcols*sizeof(float));
    h->G_hat = malloc1d(nYcols*nYcols*sizeof(float_complex));
    h->Cx_QH = malloc1d(nXcols*nYcols*sizeof(float_complex));
    
    /* For the formulation of optimal P */
    h->GhatH_Ky = malloc1d(nYcols*nYcols*sizeof(float_complex));
    h->QH_GhatH_Ky = malloc1d(nXcols*nYcols*sizeof(float_complex));
    h->KxH_QH_GhatH_Ky = malloc1d(nXcols*nYcols*sizeof(float_complex));
    h->U = malloc1d(nXcols*nXcols*sizeof(float_complex));
    h->V = malloc1d(nYcols*nYcols*sizeof(float_complex));
    h->lambda_UH = malloc1d(nYcols*nXcols*sizeof(float_complex));
    h->P = malloc1d(nYcols*nXcols*sizeof(float_complex));
    
    /* For the formulation of M */
    h->P_Kxreginverse = malloc1d(nYcols*nXcols*sizeof(float_complex));
    
    /* For the formulation of the residual covariance matrix */
    h->Cx_MH = malloc1d(nXcols*nYcols*sizeof(float_complex));
    h->Cy_tilde = malloc1d(nYcols*nYcols*sizeof(float_complex));
    
    /* For using energy compensation instead of residuals */
    h->G_M = malloc1d(nYcols*nXcols*sizeof(float_complex));
}

void cdf4sap_destroy
(
    void ** const phCdf
)
{
    cdf4sap_data *h = (cdf4sap_data*)(*phCdf);
    
    if(h!=NULL){
        utility_ssvd_destroy(&h->hSVD);
        free(h->lambda);
        free(h->U_Cy);
        free(h->S_Cy);
        free(h->Ky);
        free(h->U_Cx);
        free(h->S_Cx);
        free(h->s_Cx);
        free(h->Kx);
        free(h->Kx_reg_inverse);
        free(h->G_hat);
        free(h->Cx_QH);
        free(h->GhatH_Ky);
        free(h->QH_GhatH_Ky);
        free(h->KxH_QH_GhatH_Ky);
        free(h->U);
        free(h->V);
        free(h->lambda_UH);
        free(h->P);
        free(h->P_Kxreginverse);
        free(h->Cx_MH);
        free(h->Cy_tilde);
        free(h->G_M);
        free(h);
        h = NULL;
    }
}

void cdf4sap_cmplx_destroy
(
    void ** const phCdf
)
{
    cdf4sap_cmplx_data *h = (cdf4sap_cmplx_data*)(*phCdf);
    
    if(h!=NULL){
        utility_csvd_destroy(&h->hSVD);
        free(h->lambda);
        free(h->Cr_cmplx);
        free(h->U_Cy);
        free(h->S_Cy);
        free(h->Ky);
        free(h->U_Cx);
        free(h->S_Cx);
        free(h->s_Cx);
        free(h->Kx);
        free(h->Kx_reg_inverse);
        free(h->G_hat_diag);
        free(h->G_hat);
        free(h->Cx_QH);
        free(h->GhatH_Ky);
        free(h->QH_GhatH_Ky);
        free(h->KxH_QH_GhatH_Ky);
        free(h->U);
        free(h->V);
        free(h->lambda_UH);
        free(h->P);
        free(h->P_Kxreginverse);
        free(h->Cx_MH);
        free(h->Cy_tilde);
        free(h->G_M);
        free(h);
        h = NULL;
    }
}

void formulate_M_and_Cr
(
    void * const hCdf,
    float* Cx,
    float* Cy,
    float* Q,
    int useEnergyFLAG,
    float reg,
    float* M,
    float* Cr
)
{
    cdf4sap_data *h = (cdf4sap_data*)(hCdf);
    int i, j, nXcols, nYcols;
    
    nXcols = h->nXcols;
    nYcols = h->nYcols;
    
    memset(h->lambda, 0, nYcols * nXcols * sizeof(float));
    for(i = 0; i<SAF_MIN(nXcols,nYcols); i++)
        h->lambda[i*nXcols + i] = 1.0f;

    /* Decomposition of Cy */
    utility_ssvd(h->hSVD, Cy, nYcols, nYcols, h->U_Cy, h->S_Cy, NULL, NULL);
    for(i=0; i< nYcols; i++)
        h->S_Cy[i*nYcols+i] = sqrtf(SAF_MAX(h->S_Cy[i*nYcols+i], 2.23e-20f));
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nYcols, nYcols, 1.0f,
                h->U_Cy, nYcols,
                h->S_Cy, nYcols, 0.0f,
                h->Ky, nYcols);

    /* Decomposition of Cx */
    utility_ssvd(h->hSVD, Cx, nXcols, nXcols, h->U_Cx, h->S_Cx, NULL, h->s_Cx);
    for(i=0; i< nXcols; i++){
        h->S_Cx[i*nXcols+i] = sqrtf(SAF_MAX(h->S_Cx[i*nXcols+i], 2.23e-20f));
        h->s_Cx[i] = sqrtf(SAF_MAX(h->s_Cx[i], 2.23e-20f));
    }
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nXcols, nXcols, nXcols, 1.0f,
                h->U_Cx, nXcols,
                h->S_Cx, nXcols, 0.0f,
                h->Kx, nXcols);

    /* Regularisation of S_Cx */
    int ind;
    float limit, maxVal;
    utility_simaxv(h->s_Cx, nXcols, &ind);
    limit = h->s_Cx[ind] * reg + 2.23e-13f;
    for(i=0; i < nXcols; i++)
        h->S_Cx[i*nXcols+i] = 1.0f / SAF_MAX(h->S_Cx[i*nXcols+i], limit);

    /* Formulate regularised Kx^-1 */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nXcols, nXcols, nXcols, 1.0f,
                h->S_Cx, nXcols,
                h->U_Cx, nXcols, 0.0f,
                h->Kx_reg_inverse, nXcols);

    /* Formulate normalisation matrix G_hat */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nXcols, nYcols, nXcols, 1.0f,
                Cx, nXcols,
                Q, nXcols, 0.0f,
                h->Cx_QH, nYcols);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nYcols, nXcols, 1.0f,
                Q, nXcols,
                h->Cx_QH, nYcols, 0.0f,
                h->G_hat, nYcols);
    maxVal = -2.23e13f;
    for(i=0; i< nYcols; i++)
        maxVal = h->G_hat[i*nYcols+i]  > maxVal ? h->G_hat[i*nYcols+i] : maxVal;
    limit = maxVal * 0.001f + 2.23e-13f;
    for(i=0; i < nYcols; i++)
        for(j=0; j < nYcols; j++)
            h->G_hat[i*nYcols+j] = i==j ? sqrtf(SAF_MAX(Cy[i*nYcols+j], 2.23e-13f) / SAF_MAX(h->G_hat[i*nYcols+j], limit)) : 0.0f;

    /* Formulate optimal P */
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nYcols, nYcols, nYcols, 1.0f,
                h->G_hat, nYcols,
                h->Ky, nYcols, 0.0f,
                h->GhatH_Ky, nYcols);
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nXcols, nYcols, nYcols, 1.0f,
                Q, nXcols,
                h->GhatH_Ky, nYcols, 0.0f,
                h->QH_GhatH_Ky, nYcols);
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, nXcols, nYcols, nXcols, 1.0f,
                h->Kx, nXcols,
                h->QH_GhatH_Ky, nYcols, 0.0f,
                h->KxH_QH_GhatH_Ky, nYcols);
    utility_ssvd(h->hSVD, h->KxH_QH_GhatH_Ky, nXcols, nYcols, h->U, NULL, h->V, NULL);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nYcols, nXcols, nXcols, 1.0f,
                h->lambda, nXcols,
                h->U, nXcols, 0.0f,
                h->lambda_UH, nXcols);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nXcols, nYcols, 1.0f,
                h->V, nYcols,
                h->lambda_UH, nXcols, 0.0f,
                h->P, nXcols);

    /* Formulate M */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nXcols, nXcols, 1.0f,
                h->P, nXcols,
                h->Kx_reg_inverse, nXcols, 0.0f,
                h->P_Kxreginverse, nXcols);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nXcols, nYcols, 1.0f,
                h->Ky, nYcols,
                h->P_Kxreginverse, nXcols, 0.0f,
                M, nXcols);

    /* Formulate residual covariance matrix */
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nXcols, nYcols, nXcols, 1.0f,
                Cx, nXcols,
                M, nXcols, 0.0f,
                h->Cx_MH, nYcols);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nYcols, nXcols, 1.0f,
                M, nXcols,
                h->Cx_MH, nYcols, 0.0f,
                h->Cy_tilde, nYcols);
    if(Cr != NULL)
        for(i=0; i < nYcols*nYcols; i++)
            Cr[i] = Cy[i] - h->Cy_tilde[i];

    /* Use energy compensation instead of residuals */
    if(useEnergyFLAG){
        for(i=0; i < nYcols; i++)
            for(j=0; j < nYcols; j++)
                h->G_hat[i*nYcols+j] = i==j ? sqrtf( SAF_MAX(Cy[i*nYcols+j], 2.23e-20f) / (h->Cy_tilde[i*nYcols+j]+2.23e-7f)) : 0.0f;
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nXcols, nYcols, 1.0f,
                    h->G_hat, nYcols,
                    M, nXcols, 0.0f,
                    h->G_M, nXcols);
        memcpy(M, h->G_M, nYcols*nXcols*sizeof(float));
        if(Cr != NULL)
            memset(Cr, 0, nYcols*nYcols*sizeof(float));
    }
}

void formulate_M_and_Cr_cmplx
(
    void * const hCdf,
    float_complex* Cx,
    float_complex* Cy,
    float_complex* Q,
    int useEnergyFLAG,
    float reg,
    float_complex* M,
    float_complex* Cr
)
{
    cdf4sap_cmplx_data *h = (cdf4sap_cmplx_data*)(hCdf);
    int i, j, nXcols, nYcols;
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta = cmplxf(0.0f, 0.0f);

    nXcols = h->nXcols;
    nYcols = h->nYcols;
    
    memset(h->lambda, 0, nYcols * nXcols * sizeof(float_complex));
    for(i = 0; i<SAF_MIN(nXcols,nYcols); i++)
        h->lambda[i*nXcols + i] = cmplxf(1.0f, 0.0f);
    
    /* Decomposition of Cy */
    utility_csvd(h->hSVD, Cy, nYcols, nYcols, h->U_Cy, h->S_Cy, NULL, NULL);
    for(i=0; i< nYcols; i++)
        h->S_Cy[i*nYcols+i] = cmplxf(sqrtf(SAF_MAX(crealf(h->S_Cy[i*nYcols+i]), 2.23e-20f)), 0.0f);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nYcols, nYcols, &calpha,
                h->U_Cy, nYcols,
                h->S_Cy, nYcols, &cbeta,
                h->Ky, nYcols);
    
    /* Decomposition of Cx */
    utility_csvd(h->hSVD, Cx, nXcols, nXcols, h->U_Cx, h->S_Cx, NULL, h->s_Cx);
    for(i=0; i< nXcols; i++){
        h->s_Cx[i] = sqrtf(SAF_MAX(h->s_Cx[i], 2.23e-13f));
        h->S_Cx[i*nXcols+i] = cmplxf(h->s_Cx[i], 0.0f);
    }
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nXcols, nXcols, nXcols, &calpha,
                h->U_Cx, nXcols,
                h->S_Cx, nXcols, &cbeta,
                h->Kx, nXcols);
    
    /* Regularisation of S_Cx */
    int ind;
    float limit, maxVal;
    //utility_simaxv(h->s_Cx, nXcols, &ind);
    ind = 0; /* utility_csvd returns the singular values in decending order */
    limit = h->s_Cx[ind] * reg + 2.23e-13f;
    for(i=0; i < nXcols; i++)
        h->S_Cx[i*nXcols+i] = cmplxf(1.0f / SAF_MAX(h->s_Cx[i], limit), 0.0f);
    
    /* Formulate regularised Kx^-1 */
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nXcols, nXcols, nXcols, &calpha,
                h->S_Cx, nXcols,
                h->U_Cx, nXcols, &cbeta,
                h->Kx_reg_inverse, nXcols);
    
    /* Formulate normalisation matrix G_hat */
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nXcols, nYcols, nXcols, &calpha,
                Cx, nXcols,
                Q, nXcols, &cbeta,
                h->Cx_QH, nYcols);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nYcols, nXcols, &calpha,
                Q, nXcols,
                h->Cx_QH, nYcols, &cbeta,
                h->G_hat, nYcols);
    /* h->G_hat: imaginary parts along the diagonal are ~0, so it's OK to take
     * the real instead below: */
    maxVal = -2.23e13f;
    for(i=0; i< nYcols; i++)
        maxVal = cabsf(h->G_hat[i*nYcols+i]) > maxVal ? cabsf(h->G_hat[i*nYcols+i]) : maxVal; // crealf->cabsf
    limit = maxVal * 0.001f + 2.23e-13f;
#if 0 /* DOES NOT PASS UNIT TEST: */
    cblas_scopy(nYcols, (float*)h->G_hat, 2*(nYcols+1), h->G_hat_diag, 1); /* take diagonal (real) */
    memset(h->G_hat, 0, nYcols*nYcols*sizeof(float_complex));
    for(i=0; i<nYcols; i++)
        h->G_hat_diag[i] = sqrtf(((float*)Cy)[2*(i*nYcols+i)]/SAF_MAX(h->G_hat_diag[i], limit));
    cblas_scopy(nYcols, (float*)h->G_hat_diag, 1, (float*)h->G_hat, /*re+im*/2*(nYcols+1)); /* load the diagonal */
#else
    for(i=0; i < nYcols; i++)
        for(j=0; j < nYcols; j++)
            h->G_hat[i*nYcols+j] = i==j ? cmplxf(crealf(csqrtf( ccdivf(Cy[i*nYcols+j], cmplxf(SAF_MAX(cabsf(h->G_hat[i*nYcols+j]), limit), 0.0f)))), 0.0f) : cmplxf(0.0f, 0.0f);  // changed crealf->cabsf
#endif

    /* Formulate optimal P */
    cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, nYcols, nYcols, nYcols, &calpha,
                h->G_hat, nYcols,
                h->Ky, nYcols, &cbeta,
                h->GhatH_Ky, nYcols);
    cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, nXcols, nYcols, nYcols, &calpha,
                Q, nXcols,
                h->GhatH_Ky, nYcols, &cbeta,
                h->QH_GhatH_Ky, nYcols);
    cblas_cgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, nXcols, nYcols, nXcols, &calpha,
                h->Kx, nXcols,
                h->QH_GhatH_Ky, nYcols, &cbeta,
                h->KxH_QH_GhatH_Ky, nYcols);
    utility_csvd(h->hSVD, h->KxH_QH_GhatH_Ky, nXcols, nYcols, h->U, NULL, h->V, NULL);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nYcols, nXcols, nXcols, &calpha,
                h->lambda, nXcols,
                h->U, nXcols, &cbeta,
                h->lambda_UH, nXcols);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nXcols, nYcols, &calpha,
                h->V, nYcols,
                h->lambda_UH, nXcols, &cbeta,
                h->P, nXcols);
    
    /* Formulate M */
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nXcols, nXcols, &calpha,
                h->P, nXcols,
                h->Kx_reg_inverse, nXcols, &cbeta,
                h->P_Kxreginverse, nXcols);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nXcols, nYcols, &calpha,
                h->Ky, nYcols,
                h->P_Kxreginverse, nXcols, &cbeta,
                M, nXcols);
    
    /* Formulate residual covariance matrix */
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nXcols, nYcols, nXcols, &calpha,
                Cx, nXcols,
                M, nXcols, &cbeta,
                h->Cx_MH, nYcols);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nYcols, nXcols, &calpha,
                M, nXcols,
                h->Cx_MH, nYcols, &cbeta,
                h->Cy_tilde, nYcols);
    if(Cr != NULL) {
#if 1
        cblas_sscal(nYcols*nYcols, 0.0f, ((float*)Cr)+1, 2);      /* set imag part to zero */
        cblas_scopy(nYcols*nYcols, (float*)Cy, 2, (float*)Cr, 2); /* copy real part */
        cblas_saxpy(nYcols*nYcols, -1.0f, (float*)h->Cy_tilde, 2, (float*)Cr, 2); /* subtract tilde to get residual (real only) */
#else
        for(i=0; i < nYcols*nYcols; i++){
            h->Cr_cmplx[i] = ccsubf(Cy[i], h->Cy_tilde[i]);
            Cr[i] = cmplxf(crealf(h->Cr_cmplx[i]), 0.0f);
            //Cr[i] = h->Cr_cmplx[i];
        }
#endif
    }

    /* Use energy compensation instead of residuals */
    if(useEnergyFLAG){
        for(i=0; i < nYcols; i++)
            for(j=0; j < nYcols; j++)
                h->G_hat[i*nYcols+j] = i==j ? csqrtf(ccdivf(Cy[i*nYcols+j], craddf(h->Cy_tilde[i*nYcols+j], 2.23e-13f))): cmplxf(0.0f, 0.0f);
        
        /* for(i=0; i < nYcols; i++)
              for(j=0; j < nYcols; j++)
                  h->G_hat[i*nYcols+j] = i==j ? cmplxf(crealf(csqrtf(ccdivf(Cy[i*nYcols+j], cmplxf(MAX(crealf(h->Cy_tilde[i*nYcols+j])+2.23e-7f, limit), 0.0f)))), 0.0f) : cmplxf(0.0f, 0.0f);
        */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nXcols, nYcols, &calpha,
                    h->G_hat, nYcols,
                    M, nXcols, &cbeta,
                    h->G_M, nXcols);
        memcpy(M, h->G_M, nYcols*nXcols*sizeof(float_complex));
        if(Cr != NULL)
            memset(Cr, 0, nYcols*nYcols*sizeof(float_complex));
    } 
}
