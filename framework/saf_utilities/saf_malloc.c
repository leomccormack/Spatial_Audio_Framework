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
 *     saf_malloc.c
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

#include "saf_malloc.h"

void ***malloc3d(int dim1, int  dim2, int dim3, size_t _Size)
{
    void*** arr;
    int i, j;
    arr = (void***)malloc(dim1 * sizeof(void**));
    for (i = 0; i < dim1; i++){
        arr[i] = (void**)malloc(dim2 * sizeof(void*));
        for (j = 0; j < dim2; j++){
            arr[i][j] = (void*)malloc(dim3 * _Size);
        }
    }
    return arr;
}

//void **malloc2d(int dim1, int  dim2, size_t _Size)
//{
//    int i;
//    void ** arr2d;
//
//    arr2d = malloc(dim1*sizeof(void*));
//    arr2d[0] = malloc(dim1*dim2*_Size);
//    for (i=1; i<dim1; i++)
//        arr2d[i] = arr2d[i-1] + dim2;
//
//    return arr2d;
//}

void **malloc2d(int dim1, int  dim2, size_t _Size)
{
#if __STDC_VERSION__ >= 199901L
    void **arr;
    long i;
    if(dim1 < 1 || dim2 < 1) return(NULL);
    i = dim1 * sizeof(void *);
    i += dim1*dim2*_Size;
    arr = (void **)malloc((size_t)i);
    if(arr!= NULL) {
        arr[0] = (void *)(arr+dim1);
        for(i = 1; i < dim1; i++){
            arr[i] = (void *)(arr[i-1])+dim2*_Size;
        }
    }
    return(arr);
#else
    void** arr;
    int i;
    arr = (void**)malloc(dim1 * sizeof(void*));
    for (i = 0; i < dim1; i++){
        arr[i] = (void*)malloc(dim2 * _Size);
    }
    return arr;
#endif
}

void* malloc1d(int dim1, size_t _Size)
{
    void* arr;
    arr = (void*)malloc(dim1 * _Size);
    return arr;
}





