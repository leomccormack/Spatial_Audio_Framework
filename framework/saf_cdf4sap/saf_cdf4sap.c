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
 *     of the matlab code by J. Vilkamo.
 * Dependencies:
 *     saf_utilities
 * Author, date created:
 *     Leo McCormack, 25.11.2016
 */

#include "saf_cdf4sap.h"

void cdf4sap_alloc
(
    void ** const phCdf,
    int CxMN,
    int CyMN,
    int QmM,
    int QmN
)
{
    cdf4sap* pData = (cdf4sap*)malloc(sizeof(cdf4sap));
    if (pData == NULL) { return;/*error*/ }
    *phCdf = (void*)pData; 
    int xsize, ysize, qsize1, qsize2;

    xsize = CxMN;
    ysize = CyMN;
    qsize1 = QmM;
    qsize2 = QmN;
    pData->xsize = CxMN;
    pData->ysize = CyMN;
    pData->qsize1 = QmM;
    pData->qsize2 = QmN;

    pData->S_Cy = (float*)malloc(ysize * sizeof(float));
    pData->S_Cy_Mtx = (float*)malloc(ysize*ysize * sizeof(float)); //zeros
    pData->U_Cy = (float*)malloc(ysize*ysize * sizeof(float));
    pData->VT_Cy = (float*)malloc(ysize*ysize * sizeof(float));
    pData->A_Cy = (float*)malloc(ysize*ysize * sizeof(float));
    pData->Ky = (float*)malloc(ysize*ysize * sizeof(float));

    pData->S_Cx = (float*)malloc(xsize * sizeof(float));
    pData->S_Cx_Mtx = (float*)malloc(xsize*xsize * sizeof(float)); //zeros
    pData->U_Cx = (float*)malloc(xsize*xsize * sizeof(float));
    pData->VT_Cx = (float*)malloc(xsize*xsize * sizeof(float));
    pData->A_Cx = (float*)malloc(xsize*xsize * sizeof(float));
    pData->Kx = (float*)malloc(xsize*xsize * sizeof(float));
    pData->Sx_reg_diag = (float*)malloc(xsize * sizeof(float));
    pData->Sx_reg_invMtx = (float*)malloc(xsize*xsize * sizeof(float)); //zeros
    pData->Kx_reg_inverse = (float*)malloc(xsize*xsize * sizeof(float));

    pData->Q = (float*)malloc(qsize1*qsize2 * sizeof(float));
    pData->QCx = (float*)malloc(qsize1*qsize2 * sizeof(float));
    pData->QCxQ = (float*)malloc(qsize1*qsize1 * sizeof(float));
    pData->KxQ = (float*)malloc(qsize2*qsize1 * sizeof(float));
    pData->KxQG = (float*)malloc(qsize2*qsize1 * sizeof(float));
    pData->KxQGKy = (float*)malloc(qsize2*qsize1 * sizeof(float));
    pData->G_hat = (float*)malloc(qsize1*qsize1 * sizeof(float)); //zeros

    pData->S = (float*)malloc(MIN(qsize1, qsize2) * sizeof(float));
    pData->U = (float*)malloc(qsize2*qsize2 * sizeof(float));
    pData->VT = (float*)malloc(qsize1*qsize1 * sizeof(float));
    pData->Vlamb = (float*)malloc(qsize1*qsize2 * sizeof(float));
    pData->P = (float*)malloc(qsize1*qsize2 * sizeof(float));
    pData->KyP = (float*)malloc(qsize1*qsize2 * sizeof(float));
    pData->M_CM = (float*)malloc(qsize1*qsize2 * sizeof(float));
    pData->MCx = (float*)malloc(qsize1*qsize2 * sizeof(float));
    pData->Cy_tilde = (float*)malloc(qsize1*qsize1 * sizeof(float));

    pData->G = (float*)malloc(ysize*ysize * sizeof(float)); //zeros
    pData->Mtmp = (float*)malloc(qsize1*qsize2 * sizeof(float));

    pData->lambda = (float*)malloc(qsize1*qsize2 * sizeof(float)); //zeros
}

void cdf4sap_alloc_cmplx
(
    void ** const phCdf,
    int CxMN,
    int CyMN,
    int QmM,
    int QmN
)
{
    cdf4sap_cmplx* pData = (cdf4sap_cmplx*)malloc(sizeof(cdf4sap_cmplx));
    if (pData == NULL) { return;/*error*/ }
    *phCdf = (void*)pData;
    int xsize, ysize, qsize1, qsize2;
    
    xsize = CxMN;
    ysize = CyMN;
    qsize1 = QmM;
    qsize2 = QmN;
    pData->xsize = CxMN;
    pData->ysize = CyMN;
    pData->qsize1 = QmM;
    pData->qsize2 = QmN;
    
    pData->S_Cy_real = (float*)malloc(ysize * sizeof(float));
    pData->S_Cy = (float_complex*)malloc(ysize * sizeof(float_complex));
    pData->S_Cy_Mtx = (float_complex*)malloc(ysize*ysize * sizeof(float_complex)); //zeros
    pData->U_Cy = (float_complex*)malloc(ysize*ysize * sizeof(float_complex));
    pData->VH_Cy = (float_complex*)malloc(ysize*ysize * sizeof(float_complex));
    pData->A_Cy = (float_complex*)malloc(ysize*ysize * sizeof(float_complex));
    pData->Ky = (float_complex*)malloc(ysize*ysize * sizeof(float_complex));
    
    pData->S_Cx_real = (float*)malloc(xsize * sizeof(float));
    pData->S_Cx = (float_complex*)malloc(xsize * sizeof(float_complex));
    pData->S_Cx_Mtx = (float_complex*)malloc(xsize*xsize * sizeof(float_complex)); //zeros
    pData->U_Cx = (float_complex*)malloc(xsize*xsize * sizeof(float_complex));
    pData->VH_Cx = (float_complex*)malloc(xsize*xsize * sizeof(float_complex));
    pData->A_Cx = (float_complex*)malloc(xsize*xsize * sizeof(float_complex));
    pData->Kx = (float_complex*)malloc(xsize*xsize * sizeof(float_complex));
    pData->Sx_reg_diag = (float_complex*)malloc(xsize * sizeof(float_complex));
    pData->Sx_reg_invMtx = (float_complex*)malloc(xsize*xsize * sizeof(float_complex)); //zeros
    pData->Kx_reg_inverse = (float_complex*)malloc(xsize*xsize * sizeof(float_complex));
    
    pData->Q = (float_complex*)malloc(qsize1*qsize2 * sizeof(float_complex));
    pData->QCx = (float_complex*)malloc(qsize1*qsize2 * sizeof(float_complex));
    pData->QCxQ = (float_complex*)malloc(qsize1*qsize1 * sizeof(float_complex));
    pData->KxQ = (float_complex*)malloc(qsize2*qsize1 * sizeof(float_complex));
    pData->KxQG = (float_complex*)malloc(qsize2*qsize1 * sizeof(float_complex));
    pData->KxQGKy = (float_complex*)malloc(qsize2*qsize1 * sizeof(float_complex));
    pData->G_hat = (float_complex*)malloc(qsize1*qsize1 * sizeof(float_complex)); //zeros
    
    pData->S_real = (float*)malloc(MIN(qsize1, qsize2) * sizeof(float));
    pData->S = (float_complex*)malloc(MIN(qsize1, qsize2) * sizeof(float_complex));
    pData->U = (float_complex*)malloc(qsize2*qsize2 * sizeof(float_complex));
    pData->VH = (float_complex*)malloc(qsize1*qsize1 * sizeof(float_complex));
    pData->Vlamb = (float_complex*)malloc(qsize1*qsize2 * sizeof(float_complex));
    pData->P = (float_complex*)malloc(qsize1*qsize2 * sizeof(float_complex));
    pData->KyP = (float_complex*)malloc(qsize1*qsize2 * sizeof(float_complex));
    pData->M_CM = (float_complex*)malloc(qsize1*qsize2 * sizeof(float_complex));
    pData->MCx = (float_complex*)malloc(qsize1*qsize2 * sizeof(float_complex));
    pData->Cy_tilde = (float_complex*)malloc(qsize1*qsize1 * sizeof(float_complex));
    
    pData->G = (float_complex*)malloc(ysize*ysize * sizeof(float_complex)); //zeros
    pData->Mtmp = (float_complex*)malloc(qsize1*qsize2 * sizeof(float_complex));
    
    pData->lambda = (float_complex*)malloc(qsize1*qsize2 * sizeof(float_complex)); //zeros
}

void cdf4sap_free(void ** const phCdf)
{
    cdf4sap *pData = (cdf4sap*)(*phCdf);

    if (pData != NULL) {
        free((void*)pData->S_Cy);
        free((void*)pData->S_Cy_Mtx);
        free((void*)pData->U_Cy);
        free((void*)pData->VT_Cy);
        free((void*)pData->A_Cy);
        free((void*)pData->Ky);

        free((void*)pData->S_Cx);
        free((void*)pData->S_Cx_Mtx);
        free((void*)pData->U_Cx);
        free((void*)pData->VT_Cx);
        free((void*)pData->A_Cx);
        free((void*)pData->Kx);
        free((void*)pData->Sx_reg_diag);
        free((void*)pData->Sx_reg_invMtx);
        free((void*)pData->Kx_reg_inverse);

        free((void*)pData->Q);
        free((void*)pData->QCx);
        free((void*)pData->QCxQ);
        free((void*)pData->KxQ);
        free((void*)pData->KxQG);
        free((void*)pData->KxQGKy);
        free((void*)pData->G_hat);

        free((void*)pData->S);
        free((void*)pData->U);
        free((void*)pData->VT);
        free((void*)pData->Vlamb);
        free((void*)pData->P);
        free((void*)pData->KyP);
        free((void*)pData->M_CM);
        free((void*)pData->MCx);
        free((void*)pData->Cy_tilde);

        free((void*)pData->G);
        free((void*)pData->Mtmp);

        free((void*)pData->lambda);

        free(pData);
        pData = NULL;
    }
}

void cdf4sap_free_cmplx(void ** const phCdf)
{
    cdf4sap_cmplx *pData = (cdf4sap_cmplx*)(*phCdf);
    
    if (pData != NULL) {
        free((void*)pData->S_Cy_real);
        free((void*)pData->S_Cy);
        free((void*)pData->S_Cy_Mtx);
        free((void*)pData->U_Cy);
        free((void*)pData->VH_Cy);
        free((void*)pData->A_Cy);
        free((void*)pData->Ky);
        
        free((void*)pData->S_Cx_real);
        free((void*)pData->S_Cx);
        free((void*)pData->S_Cx_Mtx);
        free((void*)pData->U_Cx);
        free((void*)pData->VH_Cx);
        free((void*)pData->A_Cx);
        free((void*)pData->Kx);
        free((void*)pData->Sx_reg_diag);
        free((void*)pData->Sx_reg_invMtx);
        free((void*)pData->Kx_reg_inverse);
        
        free((void*)pData->Q);
        free((void*)pData->QCx);
        free((void*)pData->QCxQ);
        free((void*)pData->KxQ);
        free((void*)pData->KxQG);
        free((void*)pData->KxQGKy);
        free((void*)pData->G_hat);
        
        free((void*)pData->S_real);
        free((void*)pData->S);
        free((void*)pData->U);
        free((void*)pData->VH);
        free((void*)pData->Vlamb);
        free((void*)pData->P);
        free((void*)pData->KyP);
        free((void*)pData->M_CM);
        free((void*)pData->MCx);
        free((void*)pData->Cy_tilde);
        
        free((void*)pData->G);
        free((void*)pData->Mtmp);
        
        free((void*)pData->lambda);
        
        free(pData);
        pData = NULL;
    }
}

/* Vilkamo, J., Bäckström, T., & Kuntz, A. (2013). Optimized covariance domain
 * framework for time–frequency processing of spatial audio. Journal of the Audio
 * Engineering Society, 61(6), 403-411. */
void cdf4sap_formulate_M_and_Cr
(
    void * const hCdf,
    float** Cx, 
    float** Cy, 
    float** Qm, 
    int flag,
    float regularisationConstant,
    float** M,
    float** Cr
)
{
    cdf4sap *pData = (cdf4sap*)(hCdf);
    int i, j;
    int m, n, lda, ldu, ldvt, info, lwork;
    int xsize, ysize, qsize1, qsize2;
    float lim;
    float_complex tmp;

    /* prep */
    xsize = pData->xsize;
    ysize = pData->ysize;
    qsize1 = pData->qsize1;
    qsize2 = pData->qsize2;

    memset(pData->S_Cy_Mtx, 0, ysize*ysize * sizeof(float));
    memset(pData->S_Cx_Mtx, 0, xsize*xsize * sizeof(float));
    memset(pData->Sx_reg_invMtx, 0, xsize*xsize * sizeof(float));
    memset(pData->G_hat, 0, qsize1*qsize1 * sizeof(float));
    memset(pData->G, 0, ysize*ysize * sizeof(float));
    memset(pData->lambda, 0, qsize1*qsize2 * sizeof(float));

    float wkopt;
    float* work;
    for (i = 0; i<qsize1; i++)
        for (j = 0; j<qsize2; j++)
            pData->Q[j*qsize1 + i] = Qm[i][j];
    if (qsize1 < qsize2){
        for (i = 0; i < MIN(qsize2, qsize1); i++)
            pData->lambda[i*MIN(qsize2, qsize1) + i] = 1.0f;
    }
    else{
        for (i = 0; i < MIN(qsize2, qsize1); i++)
            pData->lambda[i*MAX(qsize2, qsize1) + i] = 1.0f;
    }
    
    /* Decomposition of Cy */
    m = ysize; n = ysize; lda = ysize; ldu = ysize; ldvt = ysize;
    for (i = 0; i<ysize; i++)
        for (j = 0; j<ysize; j++)
            pData->A_Cy[j*ysize + i] = Cy[i][j]; /* store in column-major order */
    lwork = -1;
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    sgesvd_("A", "N", (__CLPK_integer*)&m, (__CLPK_integer*)&n, pData->A_Cy, (__CLPK_integer*)&lda, pData->S_Cy, pData->U_Cy,
            (__CLPK_integer*)&ldu, pData->VT_Cy, (__CLPK_integer*)&ldvt, &wkopt, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
#else
    sgesvd_("A", "N", &m, &n, pData->A_Cy, &lda, pData->S_Cy, pData->U_Cy,
            &ldu, pData->VT_Cy, &ldvt, &wkopt, &lwork, &info);
#endif
    lwork = (int)wkopt;
    work = (float*)malloc(lwork * sizeof(float));
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    sgesvd_("A", "N", (__CLPK_integer*)&m, (__CLPK_integer*)&n, pData->A_Cy, (__CLPK_integer*)&lda, pData->S_Cy, pData->U_Cy,
            (__CLPK_integer*)&ldu, pData->VT_Cy, (__CLPK_integer*)&ldvt, work, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
#else
    sgesvd_("A", "N", &m, &n, pData->A_Cy, &lda, pData->S_Cy, pData->U_Cy,
            &ldu, pData->VT_Cy, &ldvt, work, &lwork, &info);
#endif
    free((void*)work);
    for (i = 0; i<ysize; i++)
        pData->S_Cy_Mtx[i*ysize + i] = sqrtf(MAX(2e-13f, pData->S_Cy[i]));
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ysize, ysize, ysize, 1,
                pData->U_Cy, ysize,
                pData->S_Cy_Mtx, ysize, 0,
                pData->Ky, ysize);

    /* Decomposition of Cx */
    m = xsize; n = xsize; lda = xsize; ldu = xsize; ldvt = xsize;
    for (i = 0; i<xsize; i++)
        for (j = 0; j<xsize; j++)
            pData->A_Cx[j*xsize + i] = Cx[i][j];
    lwork = -1;
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    sgesvd_("A", "N", (__CLPK_integer*)&m, (__CLPK_integer*)&n, pData->A_Cx, (__CLPK_integer*)&lda, pData->S_Cx, pData->U_Cx,
            (__CLPK_integer*)&ldu, pData->VT_Cx, (__CLPK_integer*)&ldvt, &wkopt, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
#else
    sgesvd_("A", "N", &m, &n, pData->A_Cx, &lda, pData->S_Cx, pData->U_Cx,
            &ldu, pData->VT_Cx, &ldvt, &wkopt, &lwork, &info);
#endif
    lwork = (int)wkopt;
    work = (float*)malloc(lwork * sizeof(float));
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    sgesvd_("A", "N", (__CLPK_integer*)&m, (__CLPK_integer*)&n, pData->A_Cx, (__CLPK_integer*)&lda, pData->S_Cx, pData->U_Cx,
            (__CLPK_integer*)&ldu, pData->VT_Cx, (__CLPK_integer*)&ldvt, work, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
#else
    sgesvd_("A", "N", &m, &n, pData->A_Cx, &lda, pData->S_Cx, pData->U_Cx,
            &ldu, pData->VT_Cx, &ldvt, work, &lwork, &info);
#endif
    free((void*)work);
    for (i = 0; i<xsize; i++)
        pData->S_Cx_Mtx[i*xsize + i] = sqrtf(MAX(2e-13f, pData->S_Cx[i]));
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, xsize, xsize, xsize, 1,
                pData->U_Cx, xsize,
                pData->S_Cx_Mtx, xsize, 0,
                pData->Kx, xsize); 

    /* Regularisation of Sx */
    lim = sqrtf(MAX(2e-13f, pData->S_Cx[0]))*regularisationConstant + 2e-13f; // [0] is max
    for (i = 0; i<xsize; i++)
        pData->Sx_reg_diag[i] = MAX(sqrtf(pData->S_Cx[i]), lim);

    /* Formulate regularised Kx^-1 */
    for (i = 0; i<xsize; i++)
        pData->Sx_reg_invMtx[i*xsize + i] = 1.0f / pData->Sx_reg_diag[i];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, xsize, xsize, xsize, 1,
                pData->Sx_reg_invMtx, xsize,
                pData->U_Cx, xsize, 0,
                pData->Kx_reg_inverse, xsize);

    /* Formulate normalisation matrix G_hat */
    for (i = 0; i<ysize; i++)
        for (j = 0; j<ysize; j++)
            pData->A_Cy[j*ysize + i] = Cy[i][j];
    for (i = 0; i<xsize; i++)
        for (j = 0; j<xsize; j++)
            pData->A_Cx[j*xsize + i] = Cx[i][j];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize1, qsize2, qsize2, 1,
                pData->Q, qsize1,
                pData->A_Cx, qsize2, 0,
                pData->QCx, qsize1);
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, qsize1, qsize1, qsize2, 1,
                pData->QCx, qsize1,
                pData->Q, qsize1, 0,
                pData->QCxQ, qsize1);
    lim = 2e-13f;
    for (i = 0; i<qsize1; i++)
        lim = MAX(lim, pData->QCxQ[i*qsize1 + i]);
    lim = lim*0.001f + 2e-13f;
    for (i = 0; i<qsize1; i++) {
        tmp = ccdivf(cmplxf(pData->A_Cy[i*qsize1 + i], 0.0f), cmplxf(MAX(pData->QCxQ[i*qsize1 + i], lim), 0.0f));
        pData->G_hat[i*qsize1 + i] = crealf(csqrtf(tmp));
    }

    /* Formulate optimal P */
    cblas_sgemm(CblasColMajor, CblasTrans, CblasTrans, qsize2, qsize1, qsize2, 1,
                pData->Kx, qsize2,
                pData->Q, qsize1, 0,
                pData->KxQ, qsize2);
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, qsize2, qsize1, qsize1, 1,
                pData->KxQ, qsize2,
                pData->G_hat, qsize1, 0,
                pData->KxQG, qsize2);
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize2, qsize1, qsize1, 1,
                pData->KxQG, qsize2,
                pData->Ky, qsize1, 0,
                pData->KxQGKy, qsize2);
    m = qsize2; n = qsize1; lda = qsize2; ldu = qsize2; ldvt = qsize1;
    lwork = -1;
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    sgesvd_("A", "A", (__CLPK_integer*)&m, (__CLPK_integer*)&n, pData->KxQGKy, (__CLPK_integer*)&lda, pData->S, pData->U,
            (__CLPK_integer*)&ldu, pData->VT, (__CLPK_integer*)&ldvt, &wkopt, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
#else
    sgesvd_("A", "A", &m, &n, pData->KxQGKy, &lda, pData->S, pData->U,
            &ldu, pData->VT, &ldvt, &wkopt, &lwork, &info);
#endif
    lwork = (int)wkopt;
    work = (float*)malloc(lwork * sizeof(float));
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    sgesvd_("A", "A", (__CLPK_integer*)&m, (__CLPK_integer*)&n, pData->KxQGKy, (__CLPK_integer*)&lda, pData->S, pData->U,
            (__CLPK_integer*)&ldu, pData->VT, (__CLPK_integer*)&ldvt, work, (__CLPK_integer*)&lwork, (__CLPK_integer*)&info);
#else
    sgesvd_("A", "A", &m, &n, pData->KxQGKy, &lda, pData->S, pData->U,
            &ldu, pData->VT, &ldvt, work, &lwork, &info);
#endif
    free((void*)work);
    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, qsize1, qsize2, qsize1, 1,
                pData->VT, qsize1, /* remember lapack returned VT not V! */
                pData->lambda, qsize1, 0,
                pData->Vlamb, qsize1);
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, qsize1, qsize2, qsize2, 1,
                pData->Vlamb, qsize1,
                pData->U, qsize2, 0,
                pData->P, qsize1);

    /* Formulate M */
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize1, qsize2, qsize1, 1,
                pData->Ky, qsize1,
                pData->P, qsize1, 0,
                pData->KyP, qsize1);
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize1, qsize2, qsize2, 1,
                pData->KyP, qsize1,
                pData->Kx_reg_inverse, qsize2, 0,
                pData->M_CM, qsize1);
     
    /* Formulate residual covariance matrix */
    for (i = 0; i<xsize; i++)
        for (j = 0; j<xsize; j++)
            pData->A_Cx[j*xsize + i] = Cx[i][j];
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize1, qsize2, qsize2, 1,
                pData->M_CM, qsize1,
                pData->A_Cx, qsize2, 0,
                pData->MCx, qsize1);
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, qsize1, qsize1, qsize2, 1,
                pData->MCx, qsize1,
                pData->M_CM, qsize1, 0,
                pData->Cy_tilde, qsize1); 

    /* Use energy compensation instead of residual */
    if (flag) {
        memcpy(pData->Mtmp, pData->M_CM, qsize2*qsize1 * sizeof(float));
        for (i = 0; i<qsize1; i++) {
            tmp = ccdivf(cmplxf(Cy[i][i], 0.0f), cmplxf(pData->Cy_tilde[i*qsize1 + i] + 2e-13f, 0.0f));
            pData->G[i*qsize1 + i] = crealf(csqrtf(tmp));
            cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize1, qsize2, qsize1, 1,
                        pData->G, qsize1,
                        pData->Mtmp, qsize1, 0,
                        pData->M_CM, qsize1);
        }
    }  

    /* output */
    for (i = 0; i<qsize1; i++)
        for (j = 0; j<qsize2; j++)
            M[i][j] = pData->M_CM[j*qsize1 + i]; /* transpose - back to row major */
    for (i = 0; i<qsize1; i++)
        for (j = 0; j<qsize1; j++)
            Cr[i][j] = Cy[i][j] - pData->Cy_tilde[j*qsize1 + i]; /* transpose - back to row major */
}

/* Vilkamo, J., Bäckström, T., & Kuntz, A. (2013). Optimized covariance domain
 * framework for time–frequency processing of spatial audio. Journal of the Audio
 * Engineering Society, 61(6), 403-411. */
void cdf4sap_formulate_M_and_Cr_cmplx
(
    void * const hCdf,
    float_complex** Cx,
    float_complex** Cy,
    float_complex** Qm,
    int flag,
    float regularisationConstant,
    float_complex** M,
    float_complex** Cr
)
{
    cdf4sap_cmplx *pData = (cdf4sap_cmplx*)(hCdf);
    int i, j;
    int m, n, lda, ldu, ldvt, info, lwork;
    int xsize, ysize, qsize1, qsize2;
    float lim;
    float_complex tmp;
    const float_complex calpha = cmplxf(1.0f, 0.0f); const float_complex cbeta = cmplxf(0.0f, 0.0f); /* blas */
    
    /* prep */
    xsize = pData->xsize;
    ysize = pData->ysize;
    qsize1 = pData->qsize1;
    qsize2 = pData->qsize2;
    
    memset(pData->S_Cy_Mtx, 0, ysize*ysize * sizeof(float_complex));
    memset(pData->S_Cx_Mtx, 0, xsize*xsize * sizeof(float_complex));
    memset(pData->Sx_reg_invMtx, 0, xsize*xsize * sizeof(float_complex));
    memset(pData->G_hat, 0, qsize1*qsize1 * sizeof(float_complex));
    memset(pData->G, 0, ysize*ysize * sizeof(float_complex));
    memset(pData->lambda, 0, qsize1*qsize2 * sizeof(float_complex));
    
    float* rwork;
    float_complex wkopt;
    float_complex* work;
    for (i = 0; i<qsize1; i++)
        for (j = 0; j<qsize2; j++)
            pData->Q[j*qsize1 + i] = Qm[i][j];
    if (qsize1 < qsize2) {
        for (i = 0; i < MIN(qsize2, qsize1); i++)
            pData->lambda[i*MIN(qsize2, qsize1) + i] = cmplxf(1.0f, 0.0f);
    }
    else {
        for (i = 0; i < MIN(qsize2, qsize1); i++)
            pData->lambda[i*MAX(qsize2, qsize1) + i] = cmplxf(1.0f, 0.0f);
    }
    
    /* Decomposition of Cy */
    m = ysize; n = ysize; lda = ysize; ldu = ysize; ldvt = ysize;
    rwork = (float*)malloc(MAX( 1, 5*MIN(m,n) )*sizeof(float));
    for (i = 0; i<ysize; i++)
        for (j = 0; j<ysize; j++)
            pData->A_Cy[j*ysize + i] = Cy[i][j]; /* store in column-major order */
    lwork = -1;
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgesvd_("A", "N", (__CLPK_integer*)&m, (__CLPK_integer*)&n, (__CLPK_complex*)pData->A_Cy, (__CLPK_integer*)&lda, pData->S_Cy_real, (__CLPK_complex*)pData->U_Cy,
            (__CLPK_integer*)&ldu, (__CLPK_complex*)pData->VH_Cy, (__CLPK_integer*)&ldvt, (__CLPK_complex*)&wkopt, (__CLPK_integer*)&lwork, rwork, (__CLPK_integer*)&info);
#elif INTEL_MKL_VERSION
    cgesvd_("A", "N", &m, &n, (MKL_Complex8*)pData->A_Cy, &lda, pData->S_Cy_real, (MKL_Complex8*)pData->U_Cy,
            &ldu, (MKL_Complex8*)pData->VH_Cy, &ldvt, (MKL_Complex8*)&wkopt, &lwork, rwork, &info);
#endif
    lwork = (int)crealf(wkopt);
    work = (float_complex*)malloc(lwork * sizeof(float_complex));
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgesvd_("A", "N", (__CLPK_integer*)&m, (__CLPK_integer*)&n, (__CLPK_complex*)pData->A_Cy, (__CLPK_integer*)&lda, pData->S_Cy_real, (__CLPK_complex*)pData->U_Cy,
            (__CLPK_integer*)&ldu, (__CLPK_complex*)pData->VH_Cy, (__CLPK_integer*)&ldvt, (__CLPK_complex*)work, (__CLPK_integer*)&lwork, rwork, (__CLPK_integer*)&info);
#elif INTEL_MKL_VERSION
    cgesvd_("A", "N", &m, &n, (MKL_Complex8*)pData->A_Cy, &lda, pData->S_Cy_real, (MKL_Complex8*)pData->U_Cy,
            &ldu, (MKL_Complex8*)pData->VH_Cy, &ldvt, (MKL_Complex8*)work, &lwork, rwork, &info);
#endif
    free((void*)work);
    free((void*)rwork);
    for (i = 0; i<ysize; i++)
        pData->S_Cy_Mtx[i*ysize + i] = cmplxf(sqrtf(MAX(pData->S_Cy_real[i], 2e-13f)), 0.0f);
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ysize, ysize, ysize, &calpha,
                pData->U_Cy, ysize,
                pData->S_Cy_Mtx, ysize, &cbeta,
                pData->Ky, ysize);
    
    /* Decomposition of Cx */
    m = xsize; n = xsize; lda = xsize; ldu = xsize; ldvt = xsize;
    rwork = (float*)malloc(MAX( 1, 5*MIN(m,n) )*sizeof(float));
    for (i = 0; i<xsize; i++)
        for (j = 0; j<xsize; j++)
            pData->A_Cx[j*xsize + i] = Cx[i][j];
    lwork = -1;
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgesvd_("A", "N", (__CLPK_integer*)&m, (__CLPK_integer*)&n, (__CLPK_complex*)pData->A_Cx, (__CLPK_integer*)&lda, pData->S_Cx_real, (__CLPK_complex*)pData->U_Cx,
            (__CLPK_integer*)&ldu, (__CLPK_complex*)pData->VH_Cx, (__CLPK_integer*)&ldvt, (__CLPK_complex*)&wkopt, (__CLPK_integer*)&lwork, rwork, (__CLPK_integer*)&info);
#elif INTEL_MKL_VERSION
    cgesvd_("A", "N", &m, &n, (MKL_Complex8*)pData->A_Cx, &lda, pData->S_Cx_real, (MKL_Complex8*)pData->U_Cx,
            &ldu, (MKL_Complex8*)pData->VH_Cx, &ldvt, (MKL_Complex8*)&wkopt, &lwork, rwork, &info);
#endif
    lwork = (int)crealf(wkopt);
    work = (float_complex*)malloc(lwork * sizeof(float_complex));
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgesvd_("A", "N", (__CLPK_integer*)&m, (__CLPK_integer*)&n, (__CLPK_complex*)pData->A_Cx, (__CLPK_integer*)&lda, pData->S_Cx_real, (__CLPK_complex*)pData->U_Cx,
            (__CLPK_integer*)&ldu, (__CLPK_complex*)pData->VH_Cx, (__CLPK_integer*)&ldvt, (__CLPK_complex*)work, (__CLPK_integer*)&lwork, rwork, (__CLPK_integer*)&info);
#elif INTEL_MKL_VERSION
    cgesvd_("A", "N", &m, &n, (MKL_Complex8*)pData->A_Cx, &lda, pData->S_Cx_real, (MKL_Complex8*)pData->U_Cx,
            &ldu, (MKL_Complex8*)pData->VH_Cx, &ldvt, (MKL_Complex8*)work, &lwork, rwork, &info);
#endif
    free((void*)work);
    free((void*)rwork);
    for (i = 0; i<xsize; i++)
        pData->S_Cx_Mtx[i*xsize + i] =  cmplxf(sqrtf(MAX(pData->S_Cx_real[i], 2e-13f)), 0.0f);
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, xsize, xsize, xsize, &calpha,
                pData->U_Cx, xsize,
                pData->S_Cx_Mtx, xsize, &cbeta,
                pData->Kx, xsize);
    
    /* Regularisation of Sx */
    lim = sqrtf(MAX(2e-13f, pData->S_Cx_real[0]))*regularisationConstant + 2e-13f; // [0] is max
    for (i = 0; i<xsize; i++)
        pData->Sx_reg_diag[i] = cmplxf(MAX(sqrtf(pData->S_Cx_real[i]), lim), 0.0f);
    
    /* Formulate regularised Kx^-1 */
    for (i = 0; i<xsize; i++)
        pData->Sx_reg_invMtx[i*xsize + i] = cmplxf(1.0f / crealf(pData->Sx_reg_diag[i]), 0.0f);
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, xsize, xsize, xsize, &calpha,
                pData->Sx_reg_invMtx, xsize,
                pData->U_Cx, xsize, &cbeta,
                pData->Kx_reg_inverse, xsize);
    
    /* Formulate normalisation matrix G_hat */
    for (i = 0; i<ysize; i++)
        for (j = 0; j<ysize; j++)
            pData->A_Cy[j*ysize + i] = Cy[i][j];
    for (i = 0; i<xsize; i++)
        for (j = 0; j<xsize; j++)
            pData->A_Cx[j*xsize + i] = Cx[i][j];
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize1, qsize2, qsize2, &calpha,
                pData->Q, qsize1,
                pData->A_Cx, qsize2, &cbeta,
                pData->QCx, qsize1);
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, qsize1, qsize1, qsize2, &calpha,
                pData->QCx, qsize1,
                pData->Q, qsize1, &cbeta,
                pData->QCxQ, qsize1);
    lim = 2e-13f;
    for (i = 0; i<qsize1; i++)
        lim = MAX(lim, cabsf(pData->QCxQ[i*qsize1 + i]));
    lim = lim*0.001f + 2e-13f;
    for (i = 0; i<qsize1; i++) {
        tmp = ccdivf(pData->A_Cy[i*qsize1 + i], cmplxf(MAX(cabsf(pData->QCxQ[i*qsize1 + i]), lim), 0.0f));
        pData->G_hat[i*qsize1 + i] = cmplxf(crealf(csqrtf(tmp)),0.0f);
    }
    
    /* Formulate optimal P */
    cblas_cgemm(CblasColMajor, CblasConjTrans, CblasConjTrans, qsize2, qsize1, qsize2, &calpha,
                pData->Kx, qsize2,
                pData->Q, qsize1, &cbeta,
                pData->KxQ, qsize2);
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, qsize2, qsize1, qsize1, &calpha,
                pData->KxQ, qsize2,
                pData->G_hat, qsize1, &cbeta,
                pData->KxQG, qsize2);
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize2, qsize1, qsize1, &calpha,
                pData->KxQG, qsize2,
                pData->Ky, qsize1, &cbeta,
                pData->KxQGKy, qsize2);
    m = qsize2; n = qsize1; lda = qsize2; ldu = qsize2; ldvt = qsize1;
    rwork = (float*)malloc(MAX( 1, 5*MIN(m,n) )*sizeof(float));
    lwork = -1;
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgesvd_("A", "A", (__CLPK_integer*)&m, (__CLPK_integer*)&n, (__CLPK_complex*)pData->KxQGKy, (__CLPK_integer*)&lda, pData->S_real, (__CLPK_complex*)pData->U,
            (__CLPK_integer*)&ldu, (__CLPK_complex*)pData->VH, (__CLPK_integer*)&ldvt, (__CLPK_complex*)&wkopt, (__CLPK_integer*)&lwork, rwork, (__CLPK_integer*)&info);
#elif INTEL_MKL_VERSION
    cgesvd_("A", "A", &m, &n, (MKL_Complex8*)pData->KxQGKy, &lda, pData->S_real, (MKL_Complex8*)pData->U,
            &ldu, (MKL_Complex8*)pData->VH, &ldvt, (MKL_Complex8*)&wkopt, &lwork, rwork, &info);
#endif
    lwork = (int)crealf(wkopt);
    work = (float_complex*)malloc(lwork * sizeof(float_complex));
#if defined(__APPLE__) && !defined(SAF_USE_INTEL_MKL)
    cgesvd_("A", "A", (__CLPK_integer*)&m, (__CLPK_integer*)&n, (__CLPK_complex*)pData->KxQGKy, (__CLPK_integer*)&lda, pData->S_real, (__CLPK_complex*)pData->U,
            (__CLPK_integer*)&ldu, (__CLPK_complex*)pData->VH, (__CLPK_integer*)&ldvt, (__CLPK_complex*)work, (__CLPK_integer*)&lwork, rwork, (__CLPK_integer*)&info);
#elif INTEL_MKL_VERSION
    cgesvd_("A", "A", &m, &n, (MKL_Complex8*)pData->KxQGKy, &lda, pData->S_real, (MKL_Complex8*)pData->U,
            &ldu, (MKL_Complex8*)pData->VH, &ldvt, (MKL_Complex8*)work, &lwork, rwork, &info);
#endif
    free((void*)work);
    free((void*)rwork);
    cblas_cgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, qsize1, qsize2, qsize1, &calpha,
                pData->VH, qsize1, /* remember lapack returned V^H not V! */
                pData->lambda, qsize1, &cbeta,
                pData->Vlamb, qsize1);
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, qsize1, qsize2, qsize2, &calpha,
                pData->Vlamb, qsize1,
                pData->U, qsize2, &cbeta,
                pData->P, qsize1);
    
    /* Formulate M */
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize1, qsize2, qsize1, &calpha,
                pData->Ky, qsize1,
                pData->P, qsize1, &cbeta,
                pData->KyP, qsize1);
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize1, qsize2, qsize2, &calpha,
                pData->KyP, qsize1,
                pData->Kx_reg_inverse, qsize2, &cbeta,
                pData->M_CM, qsize1);
    
    /* Formulate residual covariance matrix */
    for (i = 0; i<xsize; i++)
        for (j = 0; j<xsize; j++)
            pData->A_Cx[j*xsize + i] = Cx[i][j];
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize1, qsize2, qsize2, &calpha,
                pData->M_CM, qsize1,
                pData->A_Cx, qsize2, &cbeta,
                pData->MCx, qsize1);
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, qsize1, qsize1, qsize2, &calpha,
                pData->MCx, qsize1,
                pData->M_CM, qsize1, &cbeta,
                pData->Cy_tilde, qsize1); 
    
    /* Use energy compensation instead of residual */
    if (flag) {
        memcpy(pData->Mtmp, pData->M_CM, qsize2*qsize1 * sizeof(float_complex));
        for (i = 0; i<qsize1; i++) {
            tmp = ccdivf(Cy[i][i], craddf(pData->Cy_tilde[i*qsize1 + i], 2e-13f));
            pData->G[i*qsize1 + i] = csqrtf(tmp);
            cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, qsize1, qsize2, qsize1, &calpha,
                        pData->G, qsize1,
                        pData->Mtmp, qsize1, &cbeta,
                        pData->M_CM, qsize1);
        }
    }  
    
    /* output */
    for (i = 0; i<qsize1; i++)
        for (j = 0; j<qsize2; j++)
            M[i][j] = pData->M_CM[j*qsize1 + i]; /* transpose - back to row major */
    for (i = 0; i<qsize1; i++) {
        for (j = 0; j<qsize1; j++) {
            //Cr[i][j] = ccsubf(Cy[i][j], pData->Cy_tilde[j*qsize1 + i]); /* transpose - back to row major */
            Cr[i][j] = cmplxf(crealf(Cy[i][j]) - crealf(pData->Cy_tilde[j*qsize1 + i]), 0.0f); /* transpose - back to row major */
        }
    }
}













