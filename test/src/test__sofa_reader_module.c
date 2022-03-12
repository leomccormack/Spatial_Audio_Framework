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
 * @file test__sofa_reader_module.c
 * @brief Unit tests for the SAF sofa reader module
 * @author Leo McCormack
 * @date 27.04.2020
 * @license ISC
 */

#include "saf_test.h"

#ifndef SAF_TEST_SOFA_FILE_PATH
# define SAF_TEST_SOFA_FILE_PATH "/Users/mccorml1/Documents/FABIAN_HRTF_DATABASE_V1/1 HRIRs/SOFA/FABIAN_HRIR_measured_HATO_20.sofa"
#endif

#ifdef SAF_ENABLE_SOFA_READER_MODULE

void test__saf_sofa_open(void){
    SAF_SOFA_ERROR_CODES error;
    saf_sofa_container sofa;
    for(int i=0; i<1 /* increase if timing... */; i++){
        /* Note that saf_sofa_open() reverts to mysofa_load(), if SAF_ENABLE_NETCDF is not defined */
        error = saf_sofa_open(&sofa, SAF_TEST_SOFA_FILE_PATH, SAF_SOFA_READER_OPTION_DEFAULT);
        saf_sofa_close(&sofa);
    }
}

void test__mysofa_load(void){
    int err;
    struct MYSOFA_HRTF *hrtf;
    for(int i=0; i<1 /* increase if timing... */; i++){
        hrtf = mysofa_load(SAF_TEST_SOFA_FILE_PATH, &err);
        mysofa_free(hrtf);
    }
}

void test__sofa_comparison(void){
#ifdef SAF_ENABLE_NETCDF
    SAF_SOFA_ERROR_CODES error;
    saf_sofa_container sofa;
    int err;
    struct MYSOFA_HRTF *hrtf;

    /* Load the same SOFA file */
    error = saf_sofa_open(&sofa, SAF_TEST_SOFA_FILE_PATH, SAF_SOFA_READER_OPTION_NETCDF);
    hrtf = mysofa_load(SAF_TEST_SOFA_FILE_PATH, &err);

    /* If both SOFA loaders were successful */
    if(error==SAF_SOFA_OK && err==MYSOFA_OK){
        /* Check that data is equivalent */
        for(unsigned int i=0; i<hrtf->M * hrtf->R * hrtf->N; i++)
            TEST_ASSERT_TRUE(fabsf(sofa.DataIR[i]- hrtf->DataIR.values[i])<0.000001f);
        for(unsigned int i=0; i<hrtf->M * hrtf->C; i++)
            TEST_ASSERT_TRUE(fabsf(sofa.SourcePosition[i]- hrtf->SourcePosition.values[i])<0.000001f);
        for(unsigned int i=0; i<hrtf->M * hrtf->C; i++)
            TEST_ASSERT_TRUE(fabsf(sofa.SourcePosition[i]- hrtf->SourcePosition.values[i])<0.000001f);
        for(unsigned int i=0; i<hrtf->_I * hrtf->R; i++)
            TEST_ASSERT_TRUE(fabsf(sofa.DataDelay[i]- hrtf->DataDelay.values[i])<0.000001f);
    }

    saf_sofa_close(&sofa);
    mysofa_free(hrtf);
#endif /* SAF_ENABLE_NETCDF */
}

#endif /* SAF_ENABLE_SOFA_READER_MODULE */
