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
 *     saf_calloc.h
 * Description:
 *     Contains a collection of useful memory allocation functions and cross-platform
 *     complex number wrappers. Optimised linear algebra routines utilising BLAS and LAPACK
 *     are also included.
 * Dependencies:
 *     Windows users only: Intel's MKL must be installed, which can be freely aquired via:
 *     https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library.
 * Author, date created:
 *     Leo McCormack, 11.07.2016
 */


#ifndef SAF_CALLOC_H_INCLUDED
#define SAF_CALLOC_H_INCLUDED

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* e.g: float*** array = (float***)calloc3d(2, 4, 6, sizeof(float));
        free3d(array, 2, 4); */

#ifdef __cplusplus
extern "C" {
#endif

void***** calloc5d(int dim1, int dim2, int dim3, int dim4, int dim5, size_t _Size);
void****  calloc4d(int dim1, int dim2, int dim3, int dim4, size_t _Size);
void***   calloc3d(int dim1, int dim2, int dim3, size_t _Size);
void**    calloc2d(int dim1, int dim2, size_t _Size);
void*     calloc1d(int dim1, size_t _Size);

#ifdef __cplusplus
}/* extern "C" */
#endif

#endif /* SAF_CALLOC_H_INCLUDED */
