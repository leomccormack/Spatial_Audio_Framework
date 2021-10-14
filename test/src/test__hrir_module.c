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
 * @file test__hrir_module.c
 * @brief Unit tests for the SAF hrir module
 * @author Leo McCormack
 * @date 27.04.2020
 * @license ISC
 */

#include "saf_test.h"

void test__resampleHRIRs(void){
    float* hrirs_out, *hrirs_tmp;
    int i, j, target_fs, hrirs_out_len, hrirs_tmp_len, max_ind;

    /* The Speex resampler has generally quite a good compromise between quality and speed.
     * This tolerance is quite high, but ultimately, it's how it sounds that matters.
     * If SAF_USE_INTEL_IPP is defined, then the Intel IPP resampler is employed instead */
    const float acceptedTolerance = 0.08f;

    /* Test 1 - passing a unit impulse through, and :qasserting the peak is where it should be */
    float* ir;
    ir = calloc1d(NUM_EARS * 256, sizeof(float));
    ir[10] = 1.0f;
    ir[256+10] = 1.0f;
    hrirs_out = NULL;
    for(i=0; i<100; i++){
        resampleHRIRs((float*)ir, 1, 256, 48000, 48000, 0, &hrirs_out, &hrirs_out_len); /* 1x samplerate */
        utility_simaxv(hrirs_out, hrirs_out_len, &max_ind);
        TEST_ASSERT_TRUE(max_ind==10);
        utility_simaxv(hrirs_out+hrirs_out_len, hrirs_out_len, &max_ind);
        TEST_ASSERT_TRUE(max_ind==10);
        free(hrirs_out); hrirs_out = NULL;
        resampleHRIRs((float*)ir, 1, 256, 48000, 96000, 0, &hrirs_out, &hrirs_out_len); /* 2x samplerate */
        utility_simaxv(hrirs_out, hrirs_out_len, &max_ind);
        TEST_ASSERT_TRUE(max_ind==20);
        utility_simaxv(hrirs_out+hrirs_out_len, hrirs_out_len, &max_ind);
        TEST_ASSERT_TRUE(max_ind==20);
        free(hrirs_out); hrirs_out = NULL;
        resampleHRIRs((float*)ir, 1, 256, 48000, 24000, 0, &hrirs_out, &hrirs_out_len); /* 0.5x samplerate */
        utility_simaxv(hrirs_out, hrirs_out_len, &max_ind);
        TEST_ASSERT_TRUE(max_ind==5);
        utility_simaxv(hrirs_out+hrirs_out_len, hrirs_out_len, &max_ind);
        TEST_ASSERT_TRUE(max_ind==5);
        free(hrirs_out); hrirs_out = NULL;
    }
    free(ir);

    /* Test 2 - converting 48e3 to 48e3 (i.e., no actual resampling, but still passing through the filter) */
    target_fs = 48000;
    hrirs_out = NULL;
    resampleHRIRs((float*)__default_hrirs, __default_N_hrir_dirs, __default_hrir_len, __default_hrir_fs,
                  target_fs, 0 /*do not zero pad*/, &hrirs_out, &hrirs_out_len);
    for(i=0; i<__default_N_hrir_dirs*NUM_EARS; i++) /* Loop over IRs */
        for(j=0; j<SAF_MIN(__default_hrir_len,hrirs_out_len); j++) /* Loop over Samples */
            TEST_ASSERT_TRUE(fabsf(((float*)__default_hrirs)[i*__default_hrir_len+j] - hrirs_out[i*hrirs_out_len+j]) <= acceptedTolerance);
    TEST_ASSERT_TRUE(__default_hrir_len==hrirs_out_len);
    free(hrirs_out);

    /* Test 2 - converting 48e3 to 96e3 and then back to 48e3 */
    target_fs = 96000;
    hrirs_tmp = NULL;
    resampleHRIRs((float*)__default_hrirs, __default_N_hrir_dirs, __default_hrir_len, __default_hrir_fs,
                  target_fs, 0 /*do not zero pad*/, &hrirs_tmp, &hrirs_tmp_len);
    target_fs = 48000;
    hrirs_out = NULL;
    resampleHRIRs(hrirs_tmp, __default_N_hrir_dirs, hrirs_tmp_len, 96000,
                  target_fs, 0 /*do not zero pad*/, &hrirs_out, &hrirs_out_len);
    for(i=0; i<__default_N_hrir_dirs*NUM_EARS; i++) /* Loop over IRs */
        for(j=0; j<SAF_MIN(__default_hrir_len,hrirs_out_len); j++) /* Loop over Samples */
            TEST_ASSERT_TRUE(fabsf(((float*)__default_hrirs)[i*__default_hrir_len+j] - hrirs_out[i*hrirs_out_len+j]) <= acceptedTolerance);
    TEST_ASSERT_TRUE(__default_hrir_len==hrirs_out_len);
    free(hrirs_tmp);
    free(hrirs_out);
}
