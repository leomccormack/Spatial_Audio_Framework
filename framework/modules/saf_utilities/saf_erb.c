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

/*
 * Filename: saf_erb.c
 * -------------------
 * A function to ascertain frequencies that fall within critical bands.
 * [Equivalent-Rectangular Bandwidth (ERB)]
 *
 * Dependencies:
 *     none
 * Author, date created:
 *     Leo McCormack, 30.07.2018
 */

#include "saf_utilities.h"

/* erb_idx start from 1 (matlab style), not 0 */
void findERBpartitions
(
    float* centerFreq,
    int nBands,
    float maxFreqLim,
    int** erb_idx,
    float** erb_freqs,
    int* nERBBands     
)
{
    int band, counter, next_erb_idx;
    float band_centreFreq, erb, erb_centre, tmp;
    
    band_centreFreq = (powf(2.0f, 1.0f/3.0f)+1.0f)/2.0f;
    free1d((void*)&(*erb_idx));
    free1d((void*)&(*erb_freqs)); 
    (*erb_idx) = malloc1d(sizeof(int));
    (*erb_freqs) = malloc1d(sizeof(float));
    (*erb_idx)[0] = 1;
    (*erb_freqs)[0] = centerFreq[0];
    counter = 0;
    next_erb_idx = 0;
    while((*erb_freqs)[counter]<maxFreqLim){
        erb = 24.7f + 0.108f * (*erb_freqs)[counter] * band_centreFreq;
        (*erb_idx) = realloc1d((*erb_idx), (counter+2)*sizeof(int));
        (*erb_freqs) = realloc1d((*erb_freqs), (counter+2)*sizeof(float));
        (*erb_freqs)[counter+1] = (*erb_freqs)[counter] + erb;
        erb_centre = FLT_MAX;
        /*  find closest band frequency as upper partition limit */
        for(band=0; band<nBands; band++){
            tmp =fabsf((*erb_freqs)[counter+1] - centerFreq[band]);
            if(tmp <erb_centre){
                erb_centre = tmp;
                next_erb_idx = band;
            }
        }
        (*erb_idx)[counter+1] = next_erb_idx + 1;
        if((*erb_idx)[counter+1] == (*erb_idx)[counter])
            (*erb_idx)[counter+1] = (*erb_idx)[counter+1]+1;
        (*erb_freqs)[counter+1] = centerFreq[(*erb_idx)[counter+1]-1];
        counter++;
    }
    /* last limit set at last band */
    (*erb_idx) = realloc1d((*erb_idx), (counter + 2) * sizeof(int));
    (*erb_freqs) = realloc1d((*erb_freqs), (counter + 2) * sizeof(float));
    (*erb_idx)[counter+1] = nBands;
    (*erb_freqs)[counter+1] = centerFreq[nBands-1];
    (*nERBBands) = counter+2;
}
