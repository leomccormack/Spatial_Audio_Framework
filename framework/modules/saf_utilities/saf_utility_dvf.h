/*
* Copyright 2020-2021 Michael McCrea, Leo McCormack
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
 *@addtogroup Utilities
 *@{
 * @file saf_utility_dvf.h
 * @brief Distance variation function filter coefficient data [1].
 *
 * @see [1] S. Spagnol, E. Tavazzi, and F. Avanzini, “Distance rendering and
 *          perception of nearby virtual sound sources with a near-field filter
 *          model,” Applied Acoustics, vol. 115, pp. 61–73, Jan. 2017,
 *          doi: 10.1016/j.apacoust.2016.08.015.
 *
 * @author Michael McCrea
 * @date 20.02.2021
 * @license ISC
 */

#ifndef SAF_DVF_H_INCLUDED
#define SAF_DVF_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <math.h>

/**
 * Calculate the Distance Variation Function (DVF) filter coefficients,
 * as described in [1].
 * 
 * @see [1] S. Spagnol, E. Tavazzi, and F. Avanzini, “Distance rendering and
 *          perception of nearby virtual sound sources with a near-field filter
 *          model,” Applied Acoustics, vol. 115, pp. 61–73, Jan. 2017,
 *          doi: 10.1016/j.apacoust.2016.08.015.
 *
 * @param[in]  alpha    Lateral angle, similar the interaural-polar convention,
 *                      but specified as an offset from the interaural axis,
 *                      [0, 180] (deg). See `doaToIpsiInteraural()`
 *                      to convert frontal azimuth/elevation to the expected
 *                      format.
 * @param[in]  rho      Source distance, normalized to head radius, >= 1.
 * @param[in]  fs       Sample rate.
 * @param[out] b        Numerator coefficients for the DVF shelving filter.
 * @param[out] a        Denominator coefficients for the DVF shelving filter
 */
void calcDVFCoeffs(
                   /* Input Arguments */
                   float alpha,
                   float rho,
                   float fs,
                   /* Output Arguments */
                   float * b,
                   float * a);

/**
 * Calculate the shelving filter parameters for the Distance Variation Function
 * filter from the source (ipsilateral) azimuth and distance.
 *
 * @param[in]  theta Lateral angle, on the inter-aural axis [0..180] (deg)
 * @param[in]  rho   Source distance, normalized to head radius, >= 1
 * @param[out] iG0   Interpolated DC gain
 * @param[out] iGInf Interpolated high shelf gain
 * @param[out] iFc   Interpolated high shelf cutoff frequency
 */

void interpDVFShelfParams(
                         /* Input Arguments */
                         float theta,
                         float rho,
                         /* Output Arguments */
                         float* iG0,
                         float* iGInf,
                         float* iFc);

/**
 * Calculate the DVF filter coefficients from shelving filter parameters.
 *
 * @param[in]  g0   High shelf gain at DC [dB].
 * @param[in]  gInf High shelf gain at Nyquist frequency [dB].
 * @param[in]  fc   Shelf cutoff frequency [Hz].
 * @param[out] fs   Sample rate.
 * @param[out] b0   Numerator coefficient 1.
 * @param[out] b1   Numerator coefficient 2.
 * @param[out] a1   Denominator coefficient 2.
*/
void dvfShelfCoeffs(/* Input Arguments */
                   float g0,
                   float gInf,
                   float fc,
                   float fs,
                   /* Output Arguments */
                   float* b0,
                   float* b1,
                   float* a1);

/**
 * Calculate the high shelf gains and cutoff parameters, given a azimuth index
 * `i` and distance `rho`. This will be called twice per parameter change and
 * the shelf filter parameters will be linearly interpolated according to the
 * azimuth.
 *
 * @param[in]   i    Coefficient table row index
 * @param[in]   rho  Normalized source distance (1 = head radius)
 * @param[out]  g0   High shelf gain at DC [dB]
 * @param[out]  gInf High shelf gain at Nyquist frequency [dB]
 * @param[out]  fc   Shelf cutoff frequency [Hz]
 */
void calcDVFShelfParams(int i,
                        float rho,
                        float* g0,
                        float* gInf,
                        float* fc);

/**
 * Convert a frontal azimuth/elevation to a modified Interaural-Polar
 * coordinate. Whereas Interaural-Polar coordinates are with reference to the
 * median plane, alpha [0, 90], beta [0, 180] this modification is with
 * reference to the transverse plane (ipsilateral ear direction), alpha
 * [0, 180], beta [0, 90]. This is intended for the input to
 * `interpDVFShelfParams()` for calculating DVF filter parameters, which
 * are framed as an offset from the interaural axis, and based on a spherical
 * head model (i.e. elevation translates to a change in lateral angle).
 *
 * @param[in]  azimuth      Source DoA, 0˚ is forward-facing, angle increases
 *                          counter-clockwise (deg, [-360, 360]).
 * @param[in]  elevation    Source elevation, angles increase upward from the
 *                          horizon (deg, [-180, 180]).
 * @param[out] alphaLR      2-element array of lateral angle alpha for left and
 *                          right ear (deg, [0,180]]).
 * @param[out] betaLR       2-element array of vertal angle beta for left and
 *                          right ear (deg, [0,90]]).
*/
void doaToIpsiInteraural(float azimuth,
                         float elevation,
                         float* alphaLR,
                         float* betaLR);

#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_DVF_H_INCLUDED */
