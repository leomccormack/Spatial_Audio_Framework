/*
* Copyright ...
*
* Permission to use ...
*
* THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
* REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
* AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
* INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
* LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
* OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
* PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef saf_utility_dvf_h
#define saf_utility_dvf_h

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <math.h>

#endif /* saf_utility_dvf_h */

/**
 * Apply the Distance Variation function to the input signal.
 *
 * @param[in]   theta Ipsilateral azimuth, on the inter-aural axis [0, 180] (deg)
 * @param[in]   rho Source distance, normalized to head radius, >= 1
 * @param[in]   in_signal (&) Input signal pointer.
 * @param[in]   nSamples Number of samples to process
 * @param[in]   fs  Sample rate
 * @param[out]  wz  (&) Filter coefficients to be passed to the next block
 * @param[out]  out_signal  (&) Output signal pointer
 */
void applyDVF(/* Input Arguments */
              float theta,   /* ipsilateral azimuth, on the inter-aural axis [0, 180] (deg) */
              float rho,     /* distance, normalized to head radius, >= 1 */
              float* in_signal,
              int nSamples,
              float fs,
              /* Output Arguments */
              float* wz,
              float* out_signal);

int levelUp(
             int a
            );

#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */
