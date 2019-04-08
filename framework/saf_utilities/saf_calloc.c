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
 *     saf_calloc.c
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
 *     Leo McCormack, 11.07.2016
 */

#include "saf_calloc.h"

void*** calloc3d(int dim1, int dim2, int dim3, size_t _Size)
{
    void*** arr;
    int i, j;
    arr = (void***)malloc(dim1 * sizeof(void**));
    for (i = 0; i < dim1; i++){
        arr[i] = (void**)malloc(dim2 * sizeof(void*));
        for (j = 0; j < dim2; j++){
            arr[i][j] = (void*)calloc(dim3, _Size);
        }
    }
    return arr;
}

void** calloc2d(int dim1, int dim2, size_t _Size)
{
    void** arr;
    int i;
    arr = (void**)malloc(dim1 * sizeof(void*));
    for (i = 0; i < dim1; i++){
        arr[i] = (void*)calloc(dim2, _Size);
    }
    return arr;
}

void* calloc1d(int dim1, size_t _Size)
{
    void* arr;
    arr = (void*)calloc(dim1, _Size);
    return arr;
}
