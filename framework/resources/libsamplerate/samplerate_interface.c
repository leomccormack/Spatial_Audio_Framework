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

#include "samplerate_interface.h"
#include "samplerate.h"

# include "../../modules/saf_utilities/saf_utilities.h"

#include <assert.h>
#include <string.h>

void sampleratelib_resample
(
    float* insig,
    int length_insig,
    int length_outsig,
    int input_fs,
    int output_fs,
    int nChannels,
    RESAMPLE_QUALITY_OPTIONS quality,
    float* outsig
)
{
    SRC_DATA data;
    int i, err;

    if(input_fs==output_fs){
        memset(outsig, 0, nChannels*length_outsig*sizeof(float));
        for(i = 0; i<nChannels; i++)
            memcpy(&outsig[i*length_outsig], &insig[i*length_insig], MIN(length_insig, length_outsig)*sizeof(float));
        return;
    }

    data.data_in = insig;
    data.data_out = outsig;
    data.input_frames = length_insig; /* frames refers to number of samples */
    data.output_frames = length_outsig;
    data.src_ratio = (double)output_fs/(double)input_fs;
    err = src_simple(&data, (int)quality, nChannels);

    assert(err==0);
}
