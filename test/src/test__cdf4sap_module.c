/*
 * Copyright 2020-2021 Leo McCormack
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
 * @file test__cdf4sap_module.c
 * @brief Unit tests for the SAF cdf4sap module
 * @author Leo McCormack
 * @date 27.04.2020
 * @license ISC
 */

#include "saf_test.h"

void test__formulate_M_and_Cr(void){
    int i, j, it, nCHin, nCHout, lenSig;
    float reg, tmp;
    float** Q, **x, **y, **z, **Cx, **Cy, **Cz, **M, **Cr, **Mr, **Q_Cx, **Cp;
    float** decor, **z_r, **eye_nCHout;
    void* hCdf, *hCdf_res;

    /* Config */
    const float acceptedTolerance = 0.1f; /* Due to regularisation, the result will never be exact */
    /* However, this is a very generous tolerance value. If the number of input
     * and output channels are similar, then this tolerance can be much lower
     * (0.00001). The error is only ever high when there is a large discrepency
     * between the number of input and output channels. */
    const int nIterations = 1000;

    /* Loop through iterations */
    for(it=0; it<nIterations; it++){
        rand_0_1(&tmp, 1);
        nCHin = (int)(tmp*12.0f + 4.1f); /* random number between 4 and 16 */
        rand_0_1(&tmp, 1);
        nCHout = (int)(tmp*12.0f + 4.1f); /* random number between 4 and 16 */
        rand_0_1(&tmp, 1);
        lenSig = (int)(tmp*384.0f + 128.1f); /* random number between 128 and 512 */

        /* Define prototype decoder and compute input signal covariance matrix */
        Q = (float**)calloc2d(nCHout, nCHin, sizeof(float));
        for(i=0; i<SAF_MIN(nCHin, nCHout); i++)
            Q[i][i] = 1.0f; /* Identity */
        x = (float**)malloc2d(nCHin, lenSig, sizeof(float));
        rand_m1_1(FLATTEN2D(x), nCHin*lenSig);
        Cx = (float**)malloc2d(nCHin, nCHin, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nCHin, nCHin, lenSig, 1.0f,
                    FLATTEN2D(x), lenSig,
                    FLATTEN2D(x), lenSig, 0.0f,
                    FLATTEN2D(Cx), nCHin);

        /* Compute target covariance matrix */
        y = (float**)malloc2d(nCHout, lenSig, sizeof(float));
        rand_m1_1(FLATTEN2D(y), nCHout*lenSig);
        Cy = (float**)malloc2d(nCHout, nCHout, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nCHout, nCHout, lenSig, 1.0f,
                    FLATTEN2D(y), lenSig,
                    FLATTEN2D(y), lenSig, 0.0f,
                    FLATTEN2D(Cy), nCHout);

        /* Compute optimal mixing matrix - with energy compensation enabled */
        M = (float**)malloc2d(nCHout, nCHin, sizeof(float));
        reg = 0.2f;
        cdf4sap_create(&hCdf, nCHin, nCHout);
        formulate_M_and_Cr(hCdf, FLATTEN2D(Cx), FLATTEN2D(Cy), FLATTEN2D(Q), 1, reg, FLATTEN2D(M), NULL);

        /* Apply mixing matrix to 'x' and assert that it's covariance matrix matches
         * the target covariance matrix */
        z = (float**)malloc2d(nCHout, lenSig, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHin, 1.0f,
                    FLATTEN2D(M), nCHin,
                    FLATTEN2D(x), lenSig, 0.0f,
                    FLATTEN2D(z), lenSig);
        Cz = (float**)malloc2d(nCHout, nCHout, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nCHout, nCHout, lenSig, 1.0f,
                    FLATTEN2D(z), lenSig,
                    FLATTEN2D(z), lenSig, 0.0f,
                    FLATTEN2D(Cz), nCHout);
        if(nCHin>=nCHout){
            for(i=0; i<nCHout; i++)
                for(j=0; j<nCHout; j++)
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, Cy[i][j], Cz[i][j]);
        }
        else{ /* if nCHin<nCHout, then only the diagonal elements will match */
            for(i=0; i<nCHout; i++)
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, Cy[i][i], Cz[i][i]);
        }

        /* Determine prototype covariance matrix */
        Q_Cx = (float**)malloc2d(nCHout, nCHin, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, nCHin, nCHin, 1.0f,
                    FLATTEN2D(Q), nCHin,
                    FLATTEN2D(Cx), nCHin, 0.0f,
                    FLATTEN2D(Q_Cx), nCHin);
        Cp = (float**)malloc2d(nCHout, nCHout, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nCHout, nCHout, nCHin, 1.0f,
                    FLATTEN2D(Q_Cx), nCHin,
                    FLATTEN2D(Q), nCHin, 0.0f,
                    FLATTEN2D(Cp), nCHout);
        for(i=0; i<nCHout; i++)
            for(j=0; j<nCHout; j++)
                if(i!=j)
                    Cp[i][j] = 0.0f; /* Zero non-diagonal elements */

        /* Create perfectly incoherent frame. Note, in practice this would instead
         * be a decorrelated version of the prototype signals, [i.e.
         * decorrelate(Q*x) ]*/
        decor = (float**)malloc2d(nCHout, lenSig, sizeof(float));
        rand_m1_1(FLATTEN2D(decor), nCHout*lenSig);

        /* Now compute optimal mixing matrix, but this time also including the
         * residual mixing matrix */
        M = (float**)malloc2d(nCHout, nCHin, sizeof(float));
        reg = 0.2f;
        Cr = (float**)malloc2d(nCHout, nCHout, sizeof(float));
        formulate_M_and_Cr(hCdf, FLATTEN2D(Cx), FLATTEN2D(Cy), FLATTEN2D(Q), 0, reg, FLATTEN2D(M), FLATTEN2D(Cr));
        cdf4sap_create(&hCdf_res, nCHout, nCHout);
        Mr = (float**)calloc2d(nCHout, nCHout, sizeof(float));
        eye_nCHout = (float**)calloc2d(nCHout, nCHout, sizeof(float));
        for(i=0; i<nCHout; i++)
            eye_nCHout[i][i] = 1.0f;
        formulate_M_and_Cr(hCdf_res, FLATTEN2D(Cp), FLATTEN2D(Cr), FLATTEN2D(eye_nCHout), 0, reg, FLATTEN2D(Mr), NULL);

        /* Apply mixing matrix to x, and residual mixing matrix to the decorrelated
         * prototype signals, and sum */
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHin, 1.0f,
                    FLATTEN2D(M), nCHin,
                    FLATTEN2D(x), lenSig, 0.0f,
                    FLATTEN2D(z), lenSig);
        z_r = (float**)malloc2d(nCHout, lenSig, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHout, 1.0f,
                    FLATTEN2D(Mr), nCHout,
                    FLATTEN2D(decor), lenSig, 0.0f,
                    FLATTEN2D(z_r), lenSig);
        cblas_saxpy(nCHout*lenSig, 1.0f, FLATTEN2D(z_r), 1, FLATTEN2D(z), 1);

        /* Assert that the covariance matrix of 'z' matches the target covariance
         * matrix */
        if(nCHin>=nCHout){
            for(i=0; i<nCHout; i++)
                for(j=0; j<nCHout; j++)
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, Cy[i][j], Cz[i][j]);
        }
        else{ /* if nCHin<nCHout, then only the diagonal elements will match */
            for(i=0; i<nCHout; i++)
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, Cy[i][i], Cz[i][i]);
        }

        /* Clean-up */
        cdf4sap_destroy(&hCdf);
        cdf4sap_destroy(&hCdf_res);
        free(Q);
        free(x);
        free(y);
        free(z);
        free(Cx);
        free(Cy);
        free(Cz);
        free(M);
        free(Cr);
        free(Mr);
        free(Q_Cx);
        free(Cp);
        free(decor);
        free(z_r);
        free(eye_nCHout);
    }
}

void test__formulate_M_and_Cr_cmplx(void){
    int i, j, it, nCHin, nCHout, lenSig;
    float reg, tmp;
    float_complex** Q, **x, **y, **z, **Cx, **Cy, **Cz, **M, **Cr, **Mr, **Q_Cx, **Cp;
    float_complex** decor, **z_r, **eye_nCHout;
    void* hCdf, *hCdf_res;
    const float_complex calpha = cmplxf(1.0f, 0.0f), cbeta = cmplxf(0.0f, 0.0f);

    /* Config */
    const float acceptedTolerance = 0.1f; /* Due to regularisation, the result will never be exact */
    /* However, this is a very generous tolerance value. If the number of input
     * and output channels are similar, then this tolerance can be much lower
     * (0.00001). The error is only ever high when there is a large discrepency
     * between the number of input and output channels. */
    const int nIterations = 300;

    /* Loop through iterations */
    for(it=0; it<nIterations; it++){
        rand_0_1(&tmp, 1);
        nCHin = (int)(tmp*12.0f + 4.1f); /* random number between 4 and 16 */
        rand_0_1(&tmp, 1);
        nCHout = (int)(tmp*12.0f + 4.1f); /* random number between 4 and 16 */
        rand_0_1(&tmp, 1);
        lenSig = (int)(tmp*384.0f + 128.1f); /* random number between 128 and 512 */

        /* Define prototype decoder and compute input signal covariance matrix */
        Q = (float_complex**)calloc2d(nCHout, nCHin, sizeof(float_complex));
        for(i=0; i<SAF_MIN(nCHin, nCHout); i++)
            Q[i][i] = cmplxf(1.0f, 0.0f); /* Identity */
        x = (float_complex**)malloc2d(nCHin, lenSig, sizeof(float_complex));
        rand_cmplx_m1_1(FLATTEN2D(x), nCHin*lenSig);
        Cx = (float_complex**)malloc2d(nCHin, nCHin, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nCHin, nCHin, lenSig, &calpha,
                    FLATTEN2D(x), lenSig,
                    FLATTEN2D(x), lenSig, &cbeta,
                    FLATTEN2D(Cx), nCHin);

        /* Compute target covariance matrix */
        y = (float_complex**)malloc2d(nCHout, lenSig, sizeof(float_complex));
        rand_cmplx_m1_1(FLATTEN2D(y), nCHout*lenSig);
        Cy = (float_complex**)malloc2d(nCHout, nCHout, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nCHout, nCHout, lenSig, &calpha,
                    FLATTEN2D(y), lenSig,
                    FLATTEN2D(y), lenSig, &cbeta,
                    FLATTEN2D(Cy), nCHout);

        /* Compute optimal mixing matrix - with energy compensation enabled */
        M = (float_complex**)malloc2d(nCHout, nCHin, sizeof(float_complex));
        reg = 0.2f;
        cdf4sap_cmplx_create(&hCdf, nCHin, nCHout);
        formulate_M_and_Cr_cmplx(hCdf, FLATTEN2D(Cx), FLATTEN2D(Cy), FLATTEN2D(Q), 1, reg, FLATTEN2D(M), NULL);

        /* Apply mixing matrix to 'x' and assert that it's covariance matrix matches
         * the target covariance matrix */
        z = (float_complex**)malloc2d(nCHout, lenSig, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHin, &calpha,
                    FLATTEN2D(M), nCHin,
                    FLATTEN2D(x), lenSig, &cbeta,
                    FLATTEN2D(z), lenSig);
        Cz = (float_complex**)malloc2d(nCHout, nCHout, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nCHout, nCHout, lenSig, &calpha,
                    FLATTEN2D(z), lenSig,
                    FLATTEN2D(z), lenSig, &cbeta,
                    FLATTEN2D(Cz), nCHout);
        if(nCHin>=nCHout){
            for(i=0; i<nCHout; i++){
                for(j=0; j<nCHout; j++){
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, crealf(Cy[i][j]), crealf(Cz[i][j]));
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cimagf(Cy[i][j]), cimagf(Cz[i][j]));
                }
            }
        }
        else{ /* if nCHin<nCHout, then only the diagonal elements will match */
            for(i=0; i<nCHout; i++){
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, crealf(Cy[i][i]), crealf(Cz[i][i]));
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cimagf(Cy[i][i]), cimagf(Cz[i][i]));
            }
        }

        /* Determine prototype covariance matrix */
        Q_Cx = (float_complex**)malloc2d(nCHout, nCHin, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, nCHin, nCHin, &calpha,
                    FLATTEN2D(Q), nCHin,
                    FLATTEN2D(Cx), nCHin, &cbeta,
                    FLATTEN2D(Q_Cx), nCHin);
        Cp = (float_complex**)malloc2d(nCHout, nCHout, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, nCHout, nCHout, nCHin, &calpha,
                    FLATTEN2D(Q_Cx), nCHin,
                    FLATTEN2D(Q), nCHin, &cbeta,
                    FLATTEN2D(Cp), nCHout);
        for(i=0; i<nCHout; i++)
            for(j=0; j<nCHout; j++)
                if(i!=j)
                    Cp[i][j] = cmplxf(0.0f, 0.0f); /* Zero non-diagonal elements */

        /* Create perfectly incoherent frame. Note, in practice this would instead
         * be a decorrelated version of the prototype signals, [i.e.
         * decorrelate(Q*x) ]*/
        decor = (float_complex**)malloc2d(nCHout, lenSig, sizeof(float_complex));
        rand_cmplx_m1_1(FLATTEN2D(decor), nCHout*lenSig);

        /* Now compute optimal mixing matrix, but this time also including the
         * residual mixing matrix */
        M = (float_complex**)malloc2d(nCHout, nCHin, sizeof(float_complex));
        reg = 0.2f;
        Cr = (float_complex**)malloc2d(nCHout, nCHout, sizeof(float_complex));
        formulate_M_and_Cr_cmplx(hCdf, FLATTEN2D(Cx), FLATTEN2D(Cy), FLATTEN2D(Q), 0, reg, FLATTEN2D(M), FLATTEN2D(Cr));
        cdf4sap_cmplx_create(&hCdf_res, nCHout, nCHout);
        Mr = (float_complex**)calloc2d(nCHout, nCHout, sizeof(float_complex));
        eye_nCHout = (float_complex**)calloc2d(nCHout, nCHout, sizeof(float_complex));
        for(i=0; i<nCHout; i++)
            eye_nCHout[i][i] = cmplxf(1.0f, 0.0f);
        formulate_M_and_Cr_cmplx(hCdf_res, FLATTEN2D(Cp), FLATTEN2D(Cr), FLATTEN2D(eye_nCHout), 0, reg, FLATTEN2D(Mr), NULL);

        /* Apply mixing matrix to x, and residual mixing matrix to the decorrelated
         * prototype signals, and sum */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHin, &calpha,
                    FLATTEN2D(M), nCHin,
                    FLATTEN2D(x), lenSig, &cbeta,
                    FLATTEN2D(z), lenSig);
        z_r = (float_complex**)malloc2d(nCHout, lenSig, sizeof(float_complex));
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nCHout, lenSig, nCHout, &calpha,
                    FLATTEN2D(Mr), nCHout,
                    FLATTEN2D(decor), lenSig, &cbeta,
                    FLATTEN2D(z_r), lenSig);
        cblas_saxpy(/*re+im*/2*nCHout*lenSig, 1.0f, (const float*)FLATTEN2D(z_r), 1, (float*)FLATTEN2D(z), 1);

        /* Assert that the covariance matrix of 'z' matches the target covariance
         * matrix */
        if(nCHin>=nCHout){
            for(i=0; i<nCHout; i++){
                for(j=0; j<nCHout; j++){
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, crealf(Cy[i][j]), crealf(Cz[i][j]));
                    TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cimagf(Cy[i][j]), cimagf(Cz[i][j]));
                }
            }
        }
        else{ /* if nCHin<nCHout, then only the diagonal elements will match */
            for(i=0; i<nCHout; i++){
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, crealf(Cy[i][i]), crealf(Cz[i][i]));
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, cimagf(Cy[i][i]), cimagf(Cz[i][i]));
            }
        }

        /* Clean-up */
        cdf4sap_cmplx_destroy(&hCdf);
        cdf4sap_cmplx_destroy(&hCdf_res);
        free(Q);
        free(x);
        free(y);
        free(z);
        free(Cx);
        free(Cy);
        free(Cz);
        free(M);
        free(Cr);
        free(Mr);
        free(Q_Cx);
        free(Cp);
        free(decor);
        free(z_r);
        free(eye_nCHout);
    }
}

