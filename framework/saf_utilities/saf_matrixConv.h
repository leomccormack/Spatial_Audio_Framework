/*
 Copyright 2019 Leo McCormack
 
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
 *     saf_matrixConv.h
 * Description:
 *     Matrix convolver.
 * Dependencies:
 *     Windows users only: custom Intel MKL '.lib/.dll' files are required.
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library by default.
 *     However, Mac users may elect to use a custom Intel MKL '.dylib' instead.
 *     Further instructions for both Windows/Mac users can be found here:
 *     https://github.com/leomccormack/Spatial_Audio_Framework
 * Author, date created:
 *     Leo McCormack, 06.04.2019
 */

#ifndef SAF_MATRIXCOLV_H_INCLUDED
#define SAF_MATRIXCOLV_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif
    
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
    
#if defined(INTEL_MKL_VERSION)
    
/* Not tested yet */
 
/* creates an instance of matrixConv.   */
void matrixConv_create(void ** const phMC,       /* address of matrixConv handle */
                       int hopSize,              /* hop size in samples. */
                       float* H,                 /* time-domain filters; flat: nCHout x nCHin x length_h*/
                       int length_h,             /* length of the filters */
                       int nCHin,                /* number of input channels */
                       int nCHout);              /* number of output channels */

/* destroys an instance ofsafFFT */
void matrixConv_destroy(void ** const phMC);     /* address of matrixConv handle */
    
/* performs the matrix convolution.  */
void matrixConv_apply(void * const hMC,          /* matrixConv handle */
                      float* inputSigs,          /* input signals; nCHin x hopSize */
                      float* outputSigs);        /* output signals; nCHout x hopSize  */
    
#else
    /* PLEASE NOTE: matrixConv only works with Intel's MKL. Apple's vDSP does not support non-power-of-2's FFTs. */
#endif

#ifdef __cplusplus
}/* extern "C" */
#endif

#endif /* SAF_MATRIXCOLV_H_INCLUDED */






