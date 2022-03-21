/*
 * Copyright 2020-2021 Leo McCormack
 *
 * This software is dual-licensed. Please refer to the LICENCE.md file for more
 * information.
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
 * @file saf_test.c
 * @brief Unit test program for the Spatial_Audio_Framework
 *
 * New unit tests may be added with the following steps:
 *
 *  1) add a unit test function prototype to the include header: saf_test.h,
 *     for example:
 * \code{.c}
 *     void test__descriptiveNameOfNewUnitTest(void);
 * \endcode
 *
 *  2) add the source code for the test in the appropriate source file. For
 *     example, if the unit test relates to the reverb module, then add the
 *     following to test__reverb_module.c
 * \code{.c}
 *     void test__descriptiveNameOfNewUnitTest(void)
 *     {
 *         // Compact usage of the SAF function(s) under test, making use of the
 *         // unity unit test framework. Refer to existing unit tests for
 *         // examples on how to interface with unity, or visit:
 *         //     https://github.com/ThrowTheSwitch/Unity
 *     }
 * \endcode
 *
 *  3) add a call for the new unit test in the main test source file: saf_test.c
 *     for example:
 * \code{.c}
 *     RUN_TEST(test__descriptiveNameOfNewUnitTest);
 * \endcode
 *
 * @author Leo McCormack
 * @date 27.04.2020
 * @license Mixed (module dependent)
 */

#include "saf_test.h" 

static tick_t start;      /**< Start time for whole test program */
static tick_t start_test; /**< Start time for the current unit test */
/** Called before each unit test is executed */
void setUp(void) { start_test = timer_current(); }
/** Called after each unit test is executed */
void tearDown(void) { }
/** Displays the time taken to run the current unit test */
static void timerResult(void) {
    printf("    (Time elapsed: %lfs) \n", (double)timer_elapsed(start_test));
}

#undef RUN_TEST
/** A custom Unity RUN_TEST, which calls timerResult() upon exiting each test */
#define RUN_TEST(testfunc)  UNITY_NEW_TEST(#testfunc) \
    if (TEST_PROTECT()) {  setUp();  testfunc();  } \
    if (TEST_PROTECT() && (!TEST_IS_IGNORED))  {tearDown(); } \
    UnityConcludeTest(); timerResult();

/* Main test program */
int main_test(void) {
    printf("%s\n", SAF_VERSION_BANNER);
    printf("%s\n", SAF_EXTERNALS_CONFIGURATION_STRING);
    printf("Executing the Spatial_Audio_Framework unit testing program");
#ifdef NDEBUG
    printf(" (Release):\n");
#else
    printf(" (Debug):\n");
#endif

    /* initialise */
    timer_lib_initialize();
    start = timer_current();
    UNITY_BEGIN();

    /* SAF utilities modules unit tests */
    RUN_TEST(test__cylindricalBesselFunctions);
    RUN_TEST(test__sphericalBesselFunctions);
    RUN_TEST(test__cart2sph);
    RUN_TEST(test__delaunaynd);
    RUN_TEST(test__quaternion);
    RUN_TEST(test__saf_stft_50pc_overlap);
    RUN_TEST(test__saf_stft_LTI);
    RUN_TEST(test__saf_matrixConv);
    RUN_TEST(test__saf_rfft);
    RUN_TEST(test__saf_fft);
    RUN_TEST(test__qmf);
    RUN_TEST(test__smb_pitchShifter);
    RUN_TEST(test__sortf);
    RUN_TEST(test__sortz);
    RUN_TEST(test__cmplxPairUp);
    RUN_TEST(test__getVoronoiWeights);
    RUN_TEST(test__unique_i);
    RUN_TEST(test__latticeDecorrelator);
    RUN_TEST(test__butterCoeffs);
    RUN_TEST(test__evalIIRTransferFunction);
    RUN_TEST(test__faf_IIRFilterbank);
    RUN_TEST(test__gexpm);
    RUN_TEST(test__dvf_calcDVFShelfParams);
    RUN_TEST(test__dvf_interpDVFShelfParams);
    RUN_TEST(test__dvf_dvfShelfCoeffs);

    /* SAF cdf4sap module unit tests */
    RUN_TEST(test__formulate_M_and_Cr);
    RUN_TEST(test__formulate_M_and_Cr_cmplx);

    /* SAF hoa module unit tests */
    RUN_TEST(test__getLoudspeakerDecoderMtx);
    RUN_TEST(test__truncationEQ);

    /* SAF sh module unit tests */
    RUN_TEST(test__getSHreal);
    RUN_TEST(test__getSHreal_recur);
    RUN_TEST(test__getSHcomplex);
    RUN_TEST(test__getSHrotMtxReal);
    RUN_TEST(test__real2complexSHMtx);
    RUN_TEST(test__complex2realSHMtx);
    RUN_TEST(test__computeSectorCoeffsEP);
    RUN_TEST(test__checkCondNumberSHTReal);
    RUN_TEST(test__calculateGridWeights);
    RUN_TEST(test__sphMUSIC);
    RUN_TEST(test__sphPWD);
    RUN_TEST(test__sphESPRIT);
    RUN_TEST(test__sphModalCoeffs);

    /* SAF hrir module unit tests */
    RUN_TEST(test__resampleHRIRs);

    /* SAF reverb modules unit tests */
    RUN_TEST(test__ims_shoebox_RIR);
    RUN_TEST(test__ims_shoebox_TD);

    /* SAF vbap modules unit tests */

    /* SAF sofa reader module unit tests */
#if defined(SAF_ENABLE_SOFA_READER_MODULE)
    RUN_TEST(test__saf_sofa_open);
    RUN_TEST(test__mysofa_load);
    RUN_TEST(test__sofa_comparison);
#endif /* SAF_ENABLE_SOFA_READER_MODULE */

    /* SAF tracker module unit tests */
#ifdef SAF_ENABLE_TRACKER_MODULE
    RUN_TEST(test__tracker3d);
#endif /* SAF_ENABLE_TRACKER_MODULE */
    
/* SAF HADES module unit tests */
#ifdef SAF_ENABLE_HADES_MODULE
    RUN_TEST(test__hades);
#endif /* SAF_ENABLE_HADES_MODULE */

    /* SAF resources unit tests */
    RUN_TEST(test__afSTFT);
    RUN_TEST(test__realloc2d_r);
    RUN_TEST(test__malloc4d);
    RUN_TEST(test__malloc5d);
    RUN_TEST(test__malloc6d);

    /* SAF examples unit tests */
#ifdef SAF_ENABLE_EXAMPLES_TESTS
    RUN_TEST(test__saf_example_ambi_bin);
    RUN_TEST(test__saf_example_ambi_dec);
    RUN_TEST(test__saf_example_ambi_enc);
    RUN_TEST(test__saf_example_array2sh);
    RUN_TEST(test__saf_example_rotator);
    RUN_TEST(test__saf_example_spreader);
#endif /* SAF_ENABLE_EXAMPLES_TESTS */

    /* close */
    timer_lib_shutdown();
    printf("\nTotal time elapsed: %lfs", (double)timer_elapsed(start));
    return UNITY_END();
}
