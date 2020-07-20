/*
 * Copyright 2020 Leo McCormack
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
 * @file samplerate_interface.h
 * @brief A custom minimal interface to libsamplerate, which exposes only the
 *        functionality which SAF requires (i.e. off-line, multi-channel
 *        conversions)
 *
 * Internal changes to "src.sinc.c" include:
 *   - Raising maximum supported channels from 128 to 65536
 *   - Employing the optimised routines from the saf_utilities module to
 *     achieve much improved performance.
 *
 * Test example: 3000 x 8192 channels of noise @48kHz, which is then converted
 * to @44.1kHz and then back to @48kHz (with SRC_SINC_BEST_QUALITY)
 *   - before: 24.6 seconds (debug), 10.4 seconds (release)
 *   - after:  2.7 seconds (debug),  2.5 seconds (release)
 *
 * @author Leo McCormack
 * @date 20.07.2020
 */

#ifndef SAMPLERATE_INTERFACE_H
#define SAMPLERATE_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif	/* __cplusplus */

/** Mirror of the enum options sent internally to src_set_converter()
 * (DO NOT CHANGE TAG IDs!) */
typedef enum _RESAMPLE_QUALITY_OPTIONS
{
    RESAMPLE_BEST_QUALITY   = 0,
    RESAMPLE_MEDIUM_QUALITY = 1,
    RESAMPLE_FASTEST        = 2,
    RESAMPLE_ORDER_HOLD     = 3,
    RESAMPLE_LINEAR         = 4

}RESAMPLE_QUALITY_OPTIONS;

/**
 * Converts the samplerate of a multichannel input signal using sampleratelib
 *
 * @note Processing is bypassed if input_fs==output_fs. Instead, the outsig
 *       is truncated if length_insig>length_outsig, or zero padded if
 *       length_insig<length_outsig.
 *
 * @param[in]  insig         Input signal; FLAT: nChannels x length_insig
 * @param[in]  length_insig  Length of input signal, in samples
 * @param[in]  length_outsig Length of output signal, in samples
 * @param[in]  input_fs      Input samplerate
 * @param[in]  output_fs     Target/output samplerate
 * @param[in]  nChannels     Number of input/output channels
 * @param[in]  quality       See #_RESAMPLE_QUALITY_OPTIONS enum
 * @param[out] outsig        Resampled/Output signal;
 *                           FLAT: nChannels x length_outsig
 */
void sampleratelib_resample(/* Input Arguments */
                            float* insig,
                            int length_insig,
                            int length_outsig,
                            int input_fs,
                            int output_fs,
                            int nChannels,
                            RESAMPLE_QUALITY_OPTIONS quality,
                            /* Output Arguments */
                            float* outsig);


#ifdef __cplusplus
}		/* extern "C" */
#endif	/* __cplusplus */

#endif	/* SAMPLERATE_INTERFACE_H */

