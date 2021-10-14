/**
 * @file test__sofa_reader_module.c
 * @brief Unit tests for the SAF sofa reader module
 * @author Leo McCormack
 * @date 27.04.2020
 * @license Mixed (module dependent)
 */

#include "saf_test.h"

#ifndef SAF_TEST_SOFA_FILE_PATH
# define SAF_TEST_SOFA_FILE_PATH "/Users/mccorml1/Documents/FABIAN_HRTF_DATABASE_V1/1 HRIRs/SOFA/FABIAN_HRIR_measured_HATO_20.sofa"
#endif

#if defined(SAF_ENABLE_SOFA_READER_MODULE)
void test__saf_sofa_open(void){
    SAF_SOFA_ERROR_CODES error;
    saf_sofa_container sofa;
    for(int i=0; i<1 /* increase if timing... */; i++){
        /* Note that saf_sofa_open() reverts to mysofa_load(), if SAF_ENABLE_NETCDF is not defined */
        error = saf_sofa_open(&sofa, SAF_TEST_SOFA_FILE_PATH);
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
    SAF_SOFA_ERROR_CODES error;
    saf_sofa_container sofa;
    int err;
    struct MYSOFA_HRTF *hrtf;

    /* Load the same SOFA file */
    error = saf_sofa_open(&sofa, SAF_TEST_SOFA_FILE_PATH);
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
}
#endif
