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
 * Apply the Distance Variation function to the input signal, as described in [1].
 * 
 * @see [1] S. Spagnol, E. Tavazzi, and F. Avanzini, “Distance rendering and
 *          perception of nearby virtual sound sources with a near-field filter
 *          model,” Applied Acoustics, vol. 115, pp. 61–73, Jan. 2017,
 *          doi: 10.1016/j.apacoust.2016.08.015.
 *
 * @param[in]  theta      Ipsilateral azimuth, on the inter-aural axis [0..180]
 * *                      (degrees)
 * @param[in]  rho        Source distance, normalized to head radius, >= 1
 * @param[in]  in_signal  (&) Input signal pointer.
 * @param[in]  nSamples   Number of samples to process
 * @param[in]  fs         Sample rate
 * @param[out] wz         (&) Filter coefficients to be passed to the next block
 * @param[out] out_signal (&) Output signal pointer
 */
void applyDVF(/* Input Arguments */
              float theta,
              float rho,
              float* in_signal,
              int nSamples,
              float fs,
              /* Output Arguments */
              float* wz,
              float* out_signal);

/**
 * Apply the Distance Variation function to the input signal
 *
 * @test test__dvf_calcHighShelfParams(), test__dvf_interpHighShelfParams()
 *
 * @param[in]  theta Ipsilateral azimuth, on the inter-aural axis [0..180] (deg)
 * @param[in]  rho   Source distance, normalized to head radius, >= 1
 * @param[out] iG0   (&) interpolated DC gain
 * @param[out] iGInf (&) interpolated high shelf gain
 * @param[out] iFc   (&) interpolated high shelf cutoff frequency
 */
void interpHighShelfParams(/* Input Arguments */
                           float theta,
                           float rho,
                           /* Output Arguments */
                           float* iG0,
                           float* iGInf,
                           float* iFc);

/**
 * Apply the Distance Variation function to the input signal
 *
 * @test test__dvf_calcIIRCoeffs()
 *
 * @param[in]  g0   High shelf gain at DC [dB]
 * @param[in]  gInf High shelf gain at Nyquist frequency [dB]
 * @param[in]  fc   Shelf cutoff frequency [Hz]
 * @param[out] fs   Sample rate
 * @param[out] b0   (&) IIR numerator coefficient 1
 * @param[out] b1   (&) IIR numerator coefficient 2
 * @param[out] a1   (&) IIR denominator coefficient 2
*/
void calcIIRCoeffs(/* Input Arguments */
                   float g0,
                   float gInf,
                   float fc,
                   float fs,
                   /* Output Arguments */
                   float* b0,
                   float* b1,
                   float* a1);

/**
 * Calculate the High shelf gains and cutoff parameters, given a azimuth index
 * 'i' and distance 'rho'. This will be called twice per parameter change and
 * the shelf filter parameters will be linearly interpolated according to the
 * precise azimuth.
 *
 * @param[in]   i    Coefficient table row index ()
 * @param[in]   rho  Normalized source distance (1 = head radius)
 * @param[out]  g0   High shelf gain at DC [dB]
 * @param[out]  gInf High shelf gain at Nyquist frequency [dB]
 * @param[out]  fc   Shelf cutoff frequency [Hz]
 */
void calcHighShelfParams(int i,
                         float rho,
                         float* g0,
                         float* gInf,
                         float* fc);

/**
 * Convert user-supplied DoA theta (0˚ forward) to a theta value understood by
 * the filter designer (DoA on the interaural axis, 0˚ is ipsilateral). Note
 * because it's based on a spherical head, the sign is always positive. Note
 * input is expected to be [-180, 180]
 *
 * @param[in]  thetaFront Source DoA, 0˚ is forward-facing, positive angles move
 *                        counter-clockwise [deg, (-180, 180)]
 * @param[out] ipsiDoaLR  High shelf gain at DC [dB]
*/
void convertFrontalDoAToIpsilateral(float thetaFront,
                                    float* ipsiDoaLR);


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_DVF_H_INCLUDED */
