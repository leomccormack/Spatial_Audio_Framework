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
 * @file test__hoa_module.c
 * @brief Unit tests for the SAF hoa module
 * @author Leo McCormack
 * @date 27.04.2020
 * @license ISC
 */

#include "saf_test.h"

void test__getLoudspeakerDecoderMtx(void){
    int i, j, k, nLS, order, nSH;
    float* ls_dirs_deg, *amp, *en;
    float** ls_dirs_rad, **decMtx_SAD, **decMtx_MMD, **decMtx_EPAD, **decMtx_AllRAD, **Ysrc, **LSout;

    /* Config */
    const float acceptedTolerance = 0.00001f;
    int nTestOrders = 10;
    int testOrders[10] = {1,2,3,4,5,6,7,8,9,10};

    /* Loop over orders */
    for(i=0; i<nTestOrders; i++) {
        order = testOrders[i];
        nSH = ORDER2NSH(order);

        /* Pull an appropriate t-design for this order */
        ls_dirs_deg = (float*)__HANDLES_Tdesign_dirs_deg[2 * order-1];
        nLS = __Tdesign_nPoints_per_degree[2 * order-1];
        ls_dirs_rad = (float**)malloc2d(nLS, 2, sizeof(float));
        for(j=0; j<nLS; j++){
            ls_dirs_rad[j][0] = ls_dirs_deg[j*2] * SAF_PI/180.0f;
            ls_dirs_rad[j][1] = SAF_PI/2.0f - ls_dirs_deg[j*2+1] * SAF_PI/180.0f; /* elevation->inclination */
        }

        /* Compute decoders */
        decMtx_SAD = (float**)malloc2d(nLS, nSH, sizeof(float));
        decMtx_MMD = (float**)malloc2d(nLS, nSH, sizeof(float));
        decMtx_EPAD = (float**)malloc2d(nLS, nSH, sizeof(float));
        decMtx_AllRAD = (float**)malloc2d(nLS, nSH, sizeof(float));
        getLoudspeakerDecoderMtx(ls_dirs_deg, nLS, LOUDSPEAKER_DECODER_SAD, order, 0, FLATTEN2D(decMtx_SAD));
        getLoudspeakerDecoderMtx(ls_dirs_deg, nLS, LOUDSPEAKER_DECODER_MMD, order, 0, FLATTEN2D(decMtx_MMD));
        getLoudspeakerDecoderMtx(ls_dirs_deg, nLS, LOUDSPEAKER_DECODER_EPAD, order, 0, FLATTEN2D(decMtx_EPAD));
        getLoudspeakerDecoderMtx(ls_dirs_deg, nLS, LOUDSPEAKER_DECODER_ALLRAD, order, 0, FLATTEN2D(decMtx_AllRAD));

        /* SAD/MMD/EPAD should all be equivalent in this special/uniform case */
        for(j=0; j<nLS; j++)
            for(k=0; k<nSH; k++)
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, decMtx_SAD[j][k], decMtx_MMD[j][k]);
        for(j=0; j<nLS; j++)
            for(k=0; k<nSH; k++)
                TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, decMtx_SAD[j][k], decMtx_EPAD[j][k]);

        /* Compute output for PWs in direction of Loudspeakers: */
        Ysrc = (float**)malloc2d(nSH, nLS, sizeof(float));
        getSHreal(order, FLATTEN2D(ls_dirs_rad), nLS, FLATTEN2D(Ysrc));
        LSout = (float**)malloc2d(nLS, nLS, sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nLS, nLS, nSH, 1.0f,
                    FLATTEN2D(decMtx_EPAD), nSH,
                    FLATTEN2D(Ysrc), nLS, 0.0f,
                    FLATTEN2D(LSout), nLS);

        /* Compute amplitude and energy for each source */
        amp = (float*)calloc1d(nLS, sizeof(float));
        en = (float*)calloc1d(nLS, sizeof(float));
        for (int idxSrc=0; idxSrc<nLS; idxSrc++) {
            for (int idxLS=0; idxLS<nLS; idxLS++) {
                amp[idxSrc] += LSout[idxLS][idxSrc];
                en[idxSrc] += LSout[idxLS][idxSrc] * LSout[idxLS][idxSrc];
            }
        }
        /* Check output amplitude and Energy */
        for (int idxSrc=0; idxSrc<nLS; idxSrc++) {
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, amp[idxSrc], 1.0);
            TEST_ASSERT_FLOAT_WITHIN(acceptedTolerance, en[idxSrc], (float)nSH / (float)nLS);
        }

        /* Clean-up */
        free(decMtx_SAD);
        free(decMtx_MMD);
        free(decMtx_EPAD);
        free(decMtx_AllRAD);
        free(ls_dirs_rad);
        free(amp);
        free(en);
        free(Ysrc);
        free(LSout);
    }
}

void test__truncationEQ(void)
{
    double *kr;
    float *w_n, *gain, *gainDB;

    /* Config */
    const int order_truncated = 4;
    const int order_target = 42;
    const float softThreshold = 12.0f;
    const int enableMaxRE = 1;
    const double fs = 48000;
    const int nBands = 128;
    kr = malloc1d(nBands * sizeof(double));
    const double r = 0.085;
    const double c = 343.;
    w_n = calloc1d((order_truncated+1), sizeof(float));
    gain = malloc1d(nBands * sizeof(float));

    /* Prep */
    double* freqVector = malloc1d(nBands*sizeof(double));
    for (int k=0; k<nBands; k++)
    {
        freqVector[k] = (double)k * fs/(2.0*((double)nBands-1));
        kr[k] = 2.0*SAF_PId / c * freqVector[k] * r;
    }
    if (enableMaxRE) {
        // maxRE as order weighting
        float *maxRECoeffs = malloc1d((order_truncated+1) * sizeof(float));
        beamWeightsMaxEV(order_truncated, maxRECoeffs);
        for (int idx_n=0; idx_n<order_truncated+1; idx_n++) {
            w_n[idx_n] = maxRECoeffs[idx_n];
            w_n[idx_n] /= sqrtf((float)(2*idx_n+1) / (FOURPI));
        }
        float w_0 = w_n[0];
        for (int idx_n=0; idx_n<order_truncated+1; idx_n++)
            w_n[idx_n] /= w_0;
        free(maxRECoeffs);
    }
    else {
        // just truncation, no tapering
        for (int idx_n=0; idx_n<order_truncated+1; idx_n++)
            w_n[idx_n] = 1.0f;
    }

    truncationEQ(w_n, order_truncated, order_target, kr, nBands, softThreshold, gain);

    /* Asserting gain offset */
    TEST_ASSERT_TRUE(gain[0]-1.0 < 2.0e-6);

    /* Asserting that gain within 0 and 12 (+6db soft clip) */
    gainDB = malloc1d(nBands * sizeof(double));
    for (int idxBand=0; idxBand<nBands; idxBand++){
        gainDB[idxBand] = 20.0f*log10f(gain[idxBand]);
        TEST_ASSERT_TRUE(gainDB[idxBand] > 0-2.0e-6);
        TEST_ASSERT_TRUE(gainDB[idxBand] < softThreshold + 6.0 + 0-2.0e-6);
    }

    /* clean-up */
    free(kr);
    free(w_n);
    free(freqVector);
    free(gain);
    free(gainDB);
}
