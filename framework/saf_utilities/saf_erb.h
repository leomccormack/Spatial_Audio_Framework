/*
 Copyright 2016-2018 Leo McCormack
 
 Permission to use, copy, modify, and/or distribute this software for any purpose with or
 without fee is hereby granted, provided that the above copyright notice and this permission
 notice appear in all copies.
 
 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
 THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT
 SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR
 ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 OR PERFORMANCE OF THIS SOFTWARE.
*/
/*
 * Filename:
 *     saf_erb.h
 * Description:
 *     Contains a collection of useful memory allocation functions and cross-platform
 *     complex number wrappers. Optimised linear algebra routines utilising BLAS and LAPACK
 *     are also included.
 * Dependencies:
 *     Windows users only: custom Intel MKL '.lib/.dll' files are required.
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library by default.
 *     However, Mac users may elect to use a custom Intel MKL '.dylib' instead.
 *     Further instructions for both Windows/Mac users can be found here:
 *     https://github.com/leomccormack/Spatial_Audio_Framework
 * Author, date created:
 *     Leo McCormack, 30.07.2018
 */

#ifndef SAF_ERB_H_INCLUDED
#define SAF_ERB_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif
    
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
    
void findERBpartitions(float* centerFreq,
                       int nBands,
                       float maxFreqLim,    /* past this frequency the bands are grouped into 1 */
                       int** erb_idx,       /* & erb indices (start from 1); nERBBands x 1 */
                       float** erb_freqs,   /* & erb frequencies; nERBBands x 1 */
                       int* nERBBands);     /* & number of bands; 1 x 1 */
 
 
#ifdef __cplusplus
}/* extern "C" */
#endif

#endif /* SAF_ERB_H_INCLUDED */




