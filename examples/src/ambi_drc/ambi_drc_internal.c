/*
 * Copyright 2017-2018 Leo McCormack
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
 * @file ambi_drc_internal.c
 * @brief A frequency-dependent Ambisonic sound scene dynamic range compressor
 *        (DRC)
 *
 * The implementation can also keep track of the frequency-dependent gain
 * factors for the omnidirectional component over time (for optional plotting).
 * The design is based on the algorithm presented in [1].
 *
 * The DRC gain factors per band are determined based on the omnidirectional
 * component, which are then applied to all of the higher-order components;
 * thus, the spatial information of the Ambisonic sound scene is retained
 * (although, your perception of them may change due to the DRC).
 *
 * @see [1] McCormack, L., & Välimäki, V. (2017). "FFT-Based Dynamic Range
 *          Compression". in Proceedings of the 14th Sound and Music Computing
 *          Conference, July 5-8, Espoo, Finland.
 *
 * @author Leo McCormack
 * @date 07.01.2017
 * @license ISC
 */

#include "ambi_drc.h"
#include "ambi_drc_internal.h"

/* Adapted from:
 * D. Giannoulis, M. Massberg, and J. D. Reiss, “Digital dynamic range compressor design: Tutorial and analysis,”
 * Journal of the Audio Engineering Society, vol. 60, no. 6, pp. 399–408, June 2012. */
float ambi_drc_gainComputer
(
    float xG,
    float T,
    float R,
    float W
)
{
    float yG;
    if (2.0f*(xG - T) < -W)
        yG = xG;
    else if (2.0f*(fabsf(xG - T)) <= W)
        yG = xG + (1.0f / R - 1.0f) * powf(xG - T + W / 2.0f, 2.0f) / (2.0f*W);
    else if (2.0f*(xG - T) > W)
        yG = T + (xG - T) / R;
    else
        yG = 0.0f;

    return yG;
}

/* Adapted from:
 * D. Giannoulis, M. Massberg, and J. D. Reiss, “Digital dynamic range compressor design: Tutorial and analysis,”
 * Journal of the Audio Engineering Society, vol. 60, no. 6, pp. 399–408, June 2012. */
float ambi_drc_smoothPeakDetector
(
    float xL,
    float yL_z1,
    float alpha_a,
    float alpha_r
)
{
    float yL;
    if (xL > yL_z1)
        yL = alpha_a*yL_z1 + (1.0f - alpha_a) * xL;
    else
        yL = alpha_r*yL_z1 + (1.0f - alpha_r) * xL;

    return yL;
}

void ambi_drc_initTFT
(
    void* const hAmbi
)
{
    ambi_drc_data *pData = (ambi_drc_data*)(hAmbi);

    /* Initialise afSTFT */
    if (pData->hSTFT == NULL)
        afSTFT_create(&(pData->hSTFT), pData->new_nSH, pData->new_nSH, HOP_SIZE, 0, 1, AFSTFT_BANDS_CH_TIME);
    else if(pData->nSH!=pData->new_nSH){/* Or change the number of channels */
        afSTFT_channelChange(pData->hSTFT, pData->new_nSH, pData->new_nSH);
        afSTFT_clearBuffers(pData->hSTFT);
    }
    pData->nSH = pData->new_nSH; 
}

void ambi_drc_setInputOrder(SH_ORDERS inOrder, int* nSH)
{
    switch(inOrder){ 
        case SH_ORDER_FIRST:
            (*nSH) = 4;
            break;
        case SH_ORDER_SECOND:
            (*nSH) = 9;
            break;
        case SH_ORDER_THIRD:
            (*nSH) = 16;
            break;
        case SH_ORDER_FOURTH:
            (*nSH) = 25;
            break;
        case SH_ORDER_FIFTH:
            (*nSH) = 36;
            break;
        case SH_ORDER_SIXTH:
            (*nSH) = 49;
            break;
        case SH_ORDER_SEVENTH:
            (*nSH) = 64;
            break;
    }
}
