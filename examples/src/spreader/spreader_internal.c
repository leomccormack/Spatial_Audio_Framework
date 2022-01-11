/*
 * Copyright 2021 Leo McCormack
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
 * @file: spreader_internal.c
 * @brief An arbitrary array panner (HRIRs, microphone array IRs, etc.) with
 *        coherent and incoherent spreading modes, as described in [1].
 *
 * @see [1] McCormack, L. Politis, A., and Pulkki, V., 2021, October. Rendering
 *          of source spread for arbitrary playback setups based on spatial
 *          covariance matching. In 2021 IEEE Workshop on Applications of Signal
 *          Processing to Audio and Acoustics (WASPAA). IEEE
 *
 * @author Leo McCormack
 * @date 07.04.2021
 * @license ISC
 */

#include "spreader_internal.h"

void spreader_setCodecStatus(void* const hSpr, CODEC_STATUS newStatus)
{
    spreader_data *pData = (spreader_data*)(hSpr);
    if(newStatus==CODEC_STATUS_NOT_INITIALISED){
        /* Pause until current initialisation is complete */
        while(pData->codecStatus == CODEC_STATUS_INITIALISING)
            SAF_SLEEP(10);
    }
    pData->codecStatus = newStatus;
}
  







 
