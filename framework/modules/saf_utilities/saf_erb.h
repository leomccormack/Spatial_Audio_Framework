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
 * @file saf_erb.h
 * @brief A function to ascertain frequencies that fall within critical bands
 *        [Equivalent-Rectangular Bandwidth (ERB)]
 *
 * @author Leo McCormack
 * @date 30.07.2018  
 */

#ifndef SAF_ERB_H_INCLUDED
#define SAF_ERB_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

/**
 * This function takes a frequency vector and groups its frequencies into
 * critical bands [Equivalent-Rectangular Bandwidth (ERB)].
 *
 * e.g.
 *   - centerFreq[erb_idx[0]-1] -> centerFreq[erb_idx[1]-1] is ERB band 1
 *   - centerFreq[erb_idx[1]-1] -> centerFreq[erb_idx[2]-1] is ERB band 2
 *
 * @warning erb indices start from 1!
 *
 * @param[in]  centerFreq Frequency vector; nBands x 1
 * @param[in]  nBands     Number of bins/bands in frequency vector
 * @param[in]  maxFreqLim Past this frequency the bands are grouped into 1 band
 * @param[out] erb_idx    (&) ERB indices (start from 1); nERBBands x 1
 * @param[out] erb_freqs  (&) ERB frequencies; nERBBands x 1
 * @param[out] nERBBands  (&) Number of ERB bands; 1 x 1
 */
void findERBpartitions(/* Input Arguments */
                       float* centerFreq,
                       int nBands,
                       float maxFreqLim,
                       /* Output Arguments */
                       int** erb_idx,
                       float** erb_freqs,
                       int* nERBBands);  


#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_ERB_H_INCLUDED */
