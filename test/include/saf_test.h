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
 * @file saf_test.h
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

#ifndef __SAF_TEST_H_INCLUDED__
#define __SAF_TEST_H_INCLUDED__

#include "resources/unity.h" /* unit testing suite (MIT license) */
#include "resources/timer.h" /* for timing the individual tests */
#include "saf.h"             /* master framework include header */
#include "saf_externals.h"   /* to also include saf dependencies (cblas etc.) */
#ifdef SAF_ENABLE_EXAMPLES_TESTS
# include "ambi_bin.h"
# include "ambi_dec.h"
# include "ambi_drc.h"
# include "ambi_enc.h"
# include "ambi_roomsim.h"
# include "array2sh.h"
# include "beamformer.h"
# include "binauraliser.h"
# include "binauraliser_nf.h"
# include "decorrelator.h"
# include "dirass.h"
# include "matrixconv.h"
# include "multiconv.h"
# include "panner.h"
# include "pitch_shifter.h"
# include "powermap.h"
# include "rotator.h"
# include "sldoa.h"
# include "spreader.h"
# include "tvconv.h"
#endif /* SAF_ENABLE_EXAMPLES_TESTS */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/** Main unit testing program */
int main_test(void);

/* ========================================================================== */
/*                      SAF utilities module unit tests                       */
/* ========================================================================== */

/**
 * Testing bessel_Jn(), bessel_Yn() */
void test__cylindricalBesselFunctions(void);
/**
 * Testing bessel_jn(), bessel_in(), bessel_yn(), and bessel_kn() */
void test__sphericalBesselFunctions(void);
/**
 * Testing cart2sph() and sph2cart() are reversible */
void test__cart2sph(void);
/**
 * Testing that the delaunaynd() function can triangulate basic shapes */
void test__delaunaynd(void);
/**
 * Testing that quaternion2rotationMatrix() and rotationMatrix2quaternion()
 * are reversible */
void test__quaternion(void);
/**
 * Testing for perfect reconstruction of the saf_stft (when configured for 50%
 * window overlap) */
void test__saf_stft_50pc_overlap(void);
/** 
 * Testing for perfect reconstruction of the saf_stft (when configured for
 * linear time-invariant (LTI) filtering applications) */
void test__saf_stft_LTI(void);
/**
 * Testing the forward and backward real-(half)complex FFT (saf_rfft) */
void test__saf_rfft(void);
/**
 * Testing the forward and backward complex-complex FFT (saf_fft) */
void test__saf_fft(void);
/**
 * Testing the saf_matrixConv */
void test__saf_matrixConv(void);
/**
 * Testing the (near)-perfect reconstruction performance of the QMF filterbank
 */
void test__qmf(void);
/**
 * Testing that the smb_pitchShifter can shift the energy of input spectra by
 * one octave down */
void test__smb_pitchShifter(void);
/**
 * Testing the sortf() function (sorting real floating point numbers) */
void test__sortf(void);
/**
 * Testing the sortz() function (sorting complex double-floating point numbers)
 */
void test__sortz(void);
/**
 * Testing the cmplxPairUp() function (grouping up conjugate symmetric values)
 */
void test__cmplxPairUp(void);
/**
 * Testing that the weights from the getVoronoiWeights() function sum to 4pi
 * and that the weights are all identical for a uniform arrangement of points
 */
void test__getVoronoiWeights(void);
/**
 * Testing that the unique_i() function operates correctly */
void test__unique_i(void);
/**
 * Testing the performance of the latticeDecorrelator, verifying that the inter-
 * channel cross-correlation coefficients are near 0 */
void test__latticeDecorrelator(void);
/**
 * Testing that the coefficients computed with butterCoeffs() are numerically
 * similar to the "butter" function in Matlab */
void test__butterCoeffs(void);
/**
 * Testing that the magnitudes and phases returned by evalIIRTransferFunction()
 * are numerically similar to the "freqz" function in Matlab */
void test__evalIIRTransferFunction(void);
/**
 * Testing that the faf_IIRFilterbank can reconstruct the original signal power
 */
void test__faf_IIRFilterbank(void);
/**
 * Testing computing the matrix exponential - comparing the output to that of
 * the "expm" function in Matlab */
void test__gexpm(void);
/**
 * Calculate high shelf parameters, g0, gInf, fc, from the lookup table
 * coefficients (10 degree steps) */
void test__dvf_calcDVFShelfParams(void);
/**
 * Test the interpolation of high shelf parameters based on distance and
 * incidence angle parameters */
void test__dvf_interpDVFShelfParams(void);
/**
 * Test the generation of high shelf coeffs based on shelf gain and fc
 * parameters */
void test__dvf_dvfShelfCoeffs(void);


/* ========================================================================== */
/*                       SAF cdf4sap module unit tests                        */
/* ========================================================================== */

/**
 * Testing the formulate_M_and_Cr() function, and verifying that the output
 * mixing matrices yield signals that have the target covariance */
void test__formulate_M_and_Cr(void);
/**
 * Testing the formulate_M_and_Cr_cmplx() function, and verifying that the
 * output mixing matrices yield signals that have the target covariance */
void test__formulate_M_and_Cr_cmplx(void);


/* ========================================================================== */
/*                         SAF hoa module unit tests                          */
/* ========================================================================== */

/**
 * Testing to assure that (given a uniform loudspeaker layout), the SAD, MMD and
 * EPAD decoders are all equivalent */
void test__getLoudspeakerDecoderMtx(void);
/**
 * Testing the truncation EQ */
void test__truncationEQ(void);


/* ========================================================================== */
/*                          SAF sh module unit tests                          */
/* ========================================================================== */

/**
 * Testing the orthogonality of the getSHreal() function */
void test__getSHreal(void);
/**
 * Testing that the getSHreal_recur() function is somewhat numerically identical
 * to the full-fat getSHreal() function */
void test__getSHreal_recur(void);
/**
 * Testing the orthogonality of the getSHcomplex() function */
void test__getSHcomplex(void);
/**
 * Testing the spherical harmonic rotation matrix function getSHrotMtxReal() */
void test__getSHrotMtxReal(void);
/**
 * Testing the real to complex spherical harmonic conversion, using
 * getSHcomplex() as the reference */
void test__real2complexSHMtx(void);
/**
 * Testing the complex to real spherical harmonic conversion, using getSHreal()
 * as the reference */
void test__complex2realSHMtx(void);
/**
 * Testing the computeSectorCoeffsEP() and computeVelCoeffsMtx() functions and
 * comparing their output to that of a reference Matlab function */
void test__computeSectorCoeffsEP(void);
/**
 * Testing that for T-designs, the condition numbers are all equal to 1 */
void test__checkCondNumberSHTReal(void);
/**
 * Test grid weight approximation */
void test__calculateGridWeights(void);
/**
 * Testing the DoA estimation performance of sphMUSIC() */
void test__sphMUSIC(void);
/**
 * Testing the DoA estimation performance of sphPWD() */
void test__sphPWD(void);
/**
 * Testing the DoA estimation performance of sphESPRIT() */
void test__sphESPRIT(void);
/**
 * Testing the sphModalCoeffs() function */
void test__sphModalCoeffs(void);


/* ========================================================================== */
/*                         SAF hrir module unit tests                         */
/* ========================================================================== */

/**
 * Testing that resampleHRIRs() is resampling adequately */
void test__resampleHRIRs(void);


/* ========================================================================== */
/*                       SAF reverb module unit tests                         */
/* ========================================================================== */

/**
 * Testing the ims shoebox simulator, when applying the echograms in the time-
 * domain */
void test__ims_shoebox_TD(void);
/**
 * Testing the ims shoebox simulator, when generating room impulse respones
 * (RIRs) from the computed echograms */
void test__ims_shoebox_RIR(void);


/* ========================================================================== */
/*                         SAF vbap module unit tests                         */
/* ========================================================================== */


/* ========================================================================== */
/*                     SAF sofa reader module unit tests                      */
/* ========================================================================== */

#ifdef SAF_ENABLE_SOFA_READER_MODULE

/**
 * Testing the SAF SOFA reader that uses the netcdf library */
void test__saf_sofa_open(void);
/**
 * Testing the dependency free mysofa SOFA reader */
void test__mysofa_load(void);
/**
 * Testing that the two SOFA readers produce the same results */
void test__sofa_comparison(void);

#endif /* SAF_ENABLE_SOFA_READER_MODULE */


/* ========================================================================== */
/*                       SAF tracker module unit tests                        */
/* ========================================================================== */

#ifdef SAF_ENABLE_TRACKER_MODULE

/**
 * Testing that the particle-filtering based tracker is able to correctly track
 * two simultaneous targets */
void test__tracker3d(void);

#endif /* SAF_ENABLE_TRACKER_MODULE */


/* ========================================================================== */
/*                        SAF HADES module unit tests                         */
/* ========================================================================== */

#ifdef SAF_ENABLE_HADES_MODULE

/** Test for hades */
void test__hades(void);

#endif /* SAF_ENABLE_HADES_MODULE */


/* ========================================================================== */
/*                         SAF resources unit tests                           */
/* ========================================================================== */

/**
 * Testing the alias-free STFT filterbank (near)-perfect reconstruction
 * performance */
void test__afSTFT(void);
/**
 * Testing the realloc2d_r() function (reallocating 2-D array, while retaining
 * the previous data order; except truncated or extended) */
void test__realloc2d_r(void);
/**
 * Testing that malloc4d() works, and is truely contiguously allocated */
void test__malloc4d(void);
/**
 * Testing that malloc5d() works, and is truely contiguously allocated */
void test__malloc5d(void);
/**
 * Testing that malloc6d() works, and is truely contiguously allocated */
void test__malloc6d(void);


/* ========================================================================== */
/*                           SAF examples unit tests                          */
/* ========================================================================== */

#ifdef SAF_ENABLE_EXAMPLES_TESTS

/**
 * Testing the SAF ambi_bin.h example (this may also serve as a tutorial on how
 * to use it) */
void test__saf_example_ambi_bin(void);
/**
 * Testing the SAF ambi_dec.h example (this may also serve as a tutorial on how
 * to use it) */
void test__saf_example_ambi_dec(void);
/**
 * Testing the SAF ambi_enc.h example (this may also serve as a tutorial on how
 * to use it) */
void test__saf_example_ambi_enc(void);
/**
 * Testing the SAF array2sh.h example (this may also serve as a tutorial on how
 * to use it) */
void test__saf_example_array2sh(void);
/**
 * Testing the SAF rotator.h example (this may also serve as a tutorial on how
 * to use it) */
void test__saf_example_rotator(void);
/**
 * Testing the SAF spreader.h example (this may also serve as a tutorial on how
 * to use it) */
void test__saf_example_spreader(void);

#endif /* SAF_ENABLE_EXAMPLES_TESTS */


/* ========================================================================== */
/*                           SAF DVF module unit tests                          */
/* ========================================================================== */



#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_TEST_H_INCLUDED__ */
