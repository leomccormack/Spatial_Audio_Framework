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
 *     saf_cdf4sap.c
 * Description:
 *     Covariance Domain Framework for Spatial Audio Processing (CDF4SAP). A C implementation
 *     of the MatLab code, originally written by Dr. Juha Vilkamo.
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 25.11.2016
 */

#include "saf_cdf4sap.h"

typedef struct _cdf4sap_data {
    /* Dimensions of Cx and Cy */
    int nXcols, nYcols;
    
    /* intermediate vectors & matrices */
    float* lambda, *U_Cy, *S_Cy, *Ky, *U_Cx, *S_Cx, *Kx, *Kx_reg_inverse, *U, *V, *P;
    float* G_hat, *Cx_QH;
    float* GhatH_Ky, *QH_GhatH_Ky, *KxH_QH_GhatH_Ky, *lambda_UH;
    float* P_Kxreginverse;
    float* Cx_MH, *Cy_tilde;
    float* G_M;
    
}cdf4sap_data;

typedef struct _cdf4sap_cmplx_data {
    /* Dimensions of Cx and Cy */
    int nXcols, nYcols;
    
    /* intermediate vectors & matrices */
    float_complex* lambda, *U_Cy, *S_Cy, *S_Cx, *Ky, *U_Cx, *Kx, *Kx_reg_inverse, *U, *V, *P;
    float* s_Cx;
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
    *phCdf = malloc(sizeof(cdf4sap_data));
    cdf4sap_data *h = (cdf4sap_data*)(*phCdf);
    
    h->nXcols = nXcols;
    h->nYcols = nYcols;
    h->lambda = malloc(nYcols * nXcols * sizeof(float));
    
    /* For the decomposition of Cy */
    h->U_Cy = malloc(nYcols*nYcols*sizeof(float));
    h->S_Cy = malloc(nYcols*nYcols*sizeof(float));
    h->Ky = malloc(nYcols*nYcols*sizeof(float));
    
    /* For the decomposition of Cx */
    h->U_Cx = malloc(nXcols*nXcols*sizeof(float));
    h->S_Cx = malloc(nXcols*nXcols*sizeof(float));
    h->Kx = malloc(nXcols*nXcols*sizeof(float));
    
    /* For the formulation of regularised Kx^-1 */
    h->Kx_reg_inverse = malloc(nXcols*nXcols*sizeof(float));
    
    /* For the formulation of normalisation matrix G_hat */
    h->G_hat = malloc(nYcols*nYcols*sizeof(float));
    h->Cx_QH = malloc(nXcols*nYcols*sizeof(float));
    
    /* For the formulation of optimal P */
    h->GhatH_Ky = malloc(nYcols*nYcols*sizeof(float));
    h->QH_GhatH_Ky = malloc(nXcols*nYcols*sizeof(float));
    h->KxH_QH_GhatH_Ky = malloc(nXcols*nYcols*sizeof(float));
    h->U = malloc(nXcols*nXcols*sizeof(float));
    h->V = malloc(nYcols*nYcols*sizeof(float));
    h->lambda_UH = malloc(nYcols*nXcols*sizeof(float));
    h->P = malloc(nYcols*nXcols*sizeof(float));
    
    /* For the formulation of M */
    h->P_Kxreginverse = malloc(nYcols*nXcols*sizeof(float));
    
    /* For the formulation of the residual covariance matrix */
    h->Cx_MH = malloc(nXcols*nYcols*sizeof(float));
    h->Cy_tilde = malloc(nYcols*nYcols*sizeof(float));
    
    /* For using energy compensation instead of residuals */
    h->G_M = malloc(nYcols*nXcols*sizeof(float));
}

void cdf4sap_cmplx_create
(
    void ** const phCdf,
    int nXcols,
    int nYcols
)
{
    *phCdf = malloc(sizeof(cdf4sap_cmplx_data));
    cdf4sap_cmplx_data *h = (cdf4sap_cmplx_data*)(*phCdf);
    
    h->nXcols = nXcols;
    h->nYcols = nYcols;
    h->lambda = malloc(nYcols * nXcols * sizeof(float_complex));
    
    /* For the decomposition of Cy */
    h->U_Cy = malloc(nYcols*nYcols*sizeof(float_complex));
    h->S_Cy = malloc(nYcols*nYcols*sizeof(float_complex));
    h->Ky = malloc(nYcols*nYcols*sizeof(float_complex));
    
    /* For the decomposition of Cx */
    h->U_Cx = malloc(nXcols*nXcols*sizeof(float_complex));
    h->S_Cx = malloc(nXcols*nXcols*sizeof(float_complex));
    h->s_Cx = malloc(nXcols*sizeof(float));
    h->Kx = malloc(nXcols*nXcols*sizeof(float_complex));
    
    /* For the formulation of regularised Kx^-1 */
    h->Kx_reg_inverse = malloc(nXcols*nXcols*sizeof(float_complex));
    
    /* For the formulation of normalisation matrix G_hat */
    h->G_hat = malloc(nYcols*nYcols*sizeof(float_complex));
    h->Cx_QH = malloc(nXcols*nYcols*sizeof(float_complex));
    
    /* For the formulation of optimal P */
    h->GhatH_Ky = malloc(nYcols*nYcols*sizeof(float_complex));
    h->QH_GhatH_Ky = malloc(nXcols*nYcols*sizeof(float_complex));
    h->KxH_QH_GhatH_Ky = malloc(nXcols*nYcols*sizeof(float_complex));
    h->U = malloc(nXcols*nXcols*sizeof(float_complex));
    h->V = malloc(nYcols*nYcols*sizeof(float_complex));
    h->lambda_UH = malloc(nYcols*nXcols*sizeof(float_complex));
    h->P = malloc(nYcols*nXcols*sizeof(float_complex));
    
    /* For the formulation of M */
    h->P_Kxreginverse = malloc(nYcols*nXcols*sizeof(float_complex));
    
    /* For the formulation of the residual covariance matrix */
    h->Cx_MH = malloc(nXcols*nYcols*sizeof(float_complex));
    h->Cy_tilde = malloc(nYcols*nYcols*sizeof(float_complex));
    
    /* For using energy compensation instead of residuals */
    h->G_M = malloc(nYcols*nXcols*sizeof(float_complex));
}

void cdf4sap_destroy
(
    void ** const phCdf
)
{
    cdf4sap_data *h = (cdf4sap_data*)(*phCdf);
    
    free(h->lambda);
    free(h->U_Cy);
    free(h->S_Cy);
    free(h->Ky);
    free(h->U_Cx);
    free(h->S_Cx);
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

void cdf4sap_cmplx_destroy
(
    void ** const phCdf
)
{
    cdf4sap_cmplx_data *h = (cdf4sap_cmplx_data*)(*phCdf);
    
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
    for(i = 0; i<nXcols; i++)
        h->lambda[i*nXcols + i] = 1.0f;
    
    /* Decomposition of Cy */
    utility_ssvd(Cy, nYcols, nYcols, h->U_Cy, h->S_Cy, NULL, NULL);
    for(i=0; i< nYcols; i++)
        h->S_Cy[i*nYcols+i] = sqrtf(h->S_Cy[i*nYcols+i]);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nYcols, nYcols, 1.0f,
                h->U_Cy, nYcols,
                h->S_Cy, nYcols, 0.0f,
                h->Ky, nYcols);
    
    /* Decomposition of Cx */
    utility_ssvd(Cx, nXcols, nXcols, h->U_Cx, h->S_Cx, NULL, NULL);
    for(i=0; i< nXcols; i++)
        h->S_Cx[i*nXcols+i] = sqrtf(h->S_Cx[i*nXcols+i]);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nXcols, nXcols, nXcols, 1.0f,
                h->U_Cx, nXcols,
                h->S_Cx, nXcols, 0.0f,
                h->Kx, nXcols);
    
    /* Regularisation of S_Cx */
    float maxVal, limit;
    maxVal = -2.23e7f;
    for(i=0; i< nXcols; i++)
        maxVal = h->S_Cx[i*nXcols+i]  > maxVal ? h->S_Cx[i*nXcols+i] : maxVal;
    limit = maxVal * reg + 2.23e-7f;
    for(i=0; i < nXcols; i++)
        h->S_Cx[i*nXcols+i] = 1.0f / MAX(h->S_Cx[i*nXcols+i], limit);
    
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
    maxVal = -2.23e7f;
    for(i=0; i< nYcols; i++)
        maxVal = h->G_hat[i*nYcols+i]  > maxVal ? h->G_hat[i*nYcols+i] : maxVal;
    limit = maxVal * 0.001f + 2.23e-7f;
    for(i=0; i < nYcols; i++)
        for(j=0; j < nYcols; j++)
            h->G_hat[i*nYcols+j] = i==j ? sqrtf(Cy[i*nYcols+j] / MAX(h->G_hat[i*nYcols+j], limit)) : 0.0f;
    
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
    utility_ssvd(h->KxH_QH_GhatH_Ky, nXcols, nYcols, h->U, NULL, h->V, NULL);
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
    for(i=0; i < nYcols*nYcols; i++)
        Cr[i] = Cy[i] - h->Cy_tilde[i];
    
    /* Use energy compensation instead of residuals */
    if(useEnergyFLAG){
        for(i=0; i < nYcols; i++)
            for(j=0; j < nYcols; j++)
                h->G_hat[i*nYcols+j] = i==j ? sqrtf(Cy[i*nYcols+j] / (h->Cy_tilde[i*nYcols+j]+2.23e-7f)) : 0.0f;
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nXcols, nYcols, 1.0f,
                    h->G_hat, nYcols,
                    M, nXcols, 0.0f,
                    h->G_M, nXcols);
        memcpy(M, h->G_M, nYcols*nXcols*sizeof(float));
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
    for(i = 0; i<nXcols; i++)
        h->lambda[i*nXcols + i] = cmplxf(1.0f, 0.0f);
    
    /* Decomposition of Cy */
    
    utility_csvd(Cy, nYcols, nYcols, h->U_Cy, h->S_Cy, NULL, NULL);
    for(i=0; i< nYcols; i++)
        h->S_Cy[i*nYcols+i] = sqrtf(h->S_Cy[i*nYcols+i]);
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nYcols, nYcols, &calpha,
                h->U_Cy, nYcols,
                h->S_Cy, nYcols, &cbeta,
                h->Ky, nYcols);
    
    /* Decomposition of Cx */
    utility_csvd(Cx, nXcols, nXcols, h->U_Cx, h->S_Cx, NULL, h->s_Cx);
    for(i=0; i< nXcols; i++){
        h->s_Cx[i] = sqrtf(h->s_Cx[i]);
        h->S_Cx[i*nXcols+i] = cmplxf(h->s_Cx[i], 0.0f);
    }
    cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nXcols, nXcols, nXcols, &calpha,
                h->U_Cx, nXcols,
                h->S_Cx, nXcols, &cbeta,
                h->Kx, nXcols);
    
    /* Regularisation of S_Cx */
    float maxVal, limit;
    maxVal = -2.23e7f;
    for(i=0; i< nXcols; i++)
        maxVal = h->s_Cx[i] > maxVal ? h->s_Cx[i] : maxVal;
    limit = maxVal * reg + 2.23e-7f;
    for(i=0; i < nXcols; i++)
        h->S_Cx[i*nXcols+i] = cmplxf(1.0f / MAX(h->s_Cx[i], limit), 0.0f);
    
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
                h->G_hat, nYcols); /* imaginary parts along the diagonal are ~0, so it's OK to take the real below: */
    maxVal = -2.23e7f;
    for(i=0; i< nYcols; i++)
        maxVal = crealf(h->G_hat[i*nYcols+i]) > maxVal ? crealf(h->G_hat[i*nYcols+i]) : maxVal;
    limit = maxVal * 0.001f + 2.23e-7f;
    for(i=0; i < nYcols; i++)
        for(j=0; j < nYcols; j++)
            h->G_hat[i*nYcols+j] = i==j ? cmplxf(sqrtf(crealf(Cy[i*nYcols+j]) / MAX(crealf(h->G_hat[i*nYcols+j]), limit)), 0.0f) : cmplxf(0.0f, 0.0f);
    
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
    utility_csvd(h->KxH_QH_GhatH_Ky, nXcols, nYcols, h->U, NULL, h->V, NULL);
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
    for(i=0; i < nYcols*nYcols; i++){
        Cr[i] = ccsubf(Cy[i], h->Cy_tilde[i]);
        Cr[i] = cmplxf(crealf(Cr[i]), 0.0f);
    }
  
    /* Use energy compensation instead of residuals */
    if(useEnergyFLAG){
        for(i=0; i < nYcols; i++)
            for(j=0; j < nYcols; j++)
                h->G_hat[i*nYcols+j] = i==j ? cmplxf(sqrtf(crealf(Cy[i*nYcols+j]) / (crealf(h->Cy_tilde[i*nYcols+j]+2.23e-7f))), 0.0f) : cmplxf(0.0f, 0.0f);
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nYcols, nXcols, nYcols, &calpha,
                    h->G_hat, nYcols,
                    M, nXcols, &cbeta,
                    h->G_M, nXcols);
        memcpy(M, h->G_M, nYcols*nXcols*sizeof(float_complex));
        memset(Cr, 0, nYcols*nYcols*sizeof(float_complex));
    }
}
 











