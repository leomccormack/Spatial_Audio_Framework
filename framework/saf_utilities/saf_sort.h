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
 *     saf_sort.h
 * Description:
 *     Contains a collection of useful memory allocation functions and cross-platform
 *     complex number wrappers. Optimised linear algebra routines utilising BLAS and LAPACK
 *     are also included.
 * Dependencies:
 *     Windows users only: Intel's MKL must be installed, which can be freely aquired via:
 *     https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries
 *     Mac users only: saf_utilities will utilise Apple's Accelerate library.
 * Author, date created:
 *     Leo McCormack, 30.07.2018
 */

#ifndef SAF_SORT_H_INCLUDED
#define SAF_SORT_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif
    
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* sort a vector of integer values into ascending/decending order. Optionally returning the new indices. */
void sorti(int* in_vec,          /* vector[len] to be sorted */
           int* out_vec,         /* if NULL, then in_vec is sorted "in-place" */
           int* new_idices,      /* set to NULL if you don't need them */
           int len,              /* number of elements in vectors, must be consistent with the input data */
           int descendFLAG);     /* !1:ascending, 1:descending */
    
/* sort a vector of floating point values into ascending/decending order. Optionally returning the new indices. */
void sortf(float* in_vec,        /* vector[len] to be sorted */
           float* out_vec,       /* if NULL, then in_vec is sorted "in-place" */
           int* new_idices,      /* set to NULL if you don't need them */
           int len,              /* number of elements in vectors, must be consistent with the input data */
           int descendFLAG);     /* !1:ascending, 1:descending */
    
/* sort a vector of double floating point values into ascending/decending order. Optionally returning the new indices. */
void sortd(double* in_vec,       /* vector[len] to be sorted */
           double* out_vec,      /* if NULL, then in_vec is sorted "in-place" */
           int* new_idices,      /* set to NULL if you don't need them */
           int len,              /* number of elements in vectors, must be consistent with the input data */
           int descendFLAG);     /* !1:ascending, 1:descending */

#ifdef __cplusplus
}/* extern "C" */
#endif

#endif /* SAF_SORT_H_INCLUDED */
