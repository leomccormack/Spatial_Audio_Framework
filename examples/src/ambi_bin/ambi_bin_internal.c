/*
 * Copyright 2018 Leo McCormack
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
 * @file ambi_bin_internal.c
 * @brief A binaural Ambisonic decoder for reproducing Ambisonic sound scenes
 *        over headphones
 *
 * The decoder offers choice over many different binaural decoding options [1-4]
 * It also supports sound-field rotation for head-tracking and can accomodate
 * loading custom HRIR sets via the SOFA standard.
 *
 * @test test__saf_example_ambi_bin()
 *
 * @see [1] Z. Ben-Hur, F. Brinkmann, J. Sheaffer, S. Weinzierl, and B. Rafaely,
 *          "Spectral equalization in binaural signals represented by order-
 *          truncated spherical harmonics" The Journal of the Acoustical
 *          Society of America, vol. 141, no. 6, pp. 4087--4096, 2017.
 * @see [2] B. Bernschutz, A. V. Giner, C. Po"rschmann, and J. Arend, "Binaural
 *          reproduction of plane waves with reduced modal order" Acta Acustica
 *          united with Acustica, vol. 100, no. 5, pp. 972--983, 2014.
 * @see [3] Zaunschirm M, Scho"rkhuber C, Ho"ldrich R. Binaural rendering of
 *          Ambisonic signals by head-related impulse response time alignment
 *          and a diffuseness constraint. The Journal of the Acoustical Society
 *          of America. 2018 Jun 19;143(6):3616-27
 * @see [4] Scho"rkhuber C, Zaunschirm M, Ho"ldrich R. Binaural Rendering of
 *          Ambisonic Signals via Magnitude Least Squares. InProceedings of the
 *          DAGA 2018 (Vol. 44, pp. 339-342).
 *
 * @author Leo McCormack
 * @date 14.04.2018
 * @license ISC
 */

#include "ambi_bin.h"
#include "ambi_bin_internal.h"

void ambi_bin_setCodecStatus(void* const hAmbi, CODEC_STATUS newStatus)
{
    ambi_bin_data *pData = (ambi_bin_data*)(hAmbi);
    if(newStatus==CODEC_STATUS_NOT_INITIALISED){
        /* Pause until current initialisation is complete */
        while(pData->codecStatus == CODEC_STATUS_INITIALISING)
            SAF_SLEEP(10);
    }
    pData->codecStatus = newStatus;
}
