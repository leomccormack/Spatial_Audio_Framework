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
 *     saf_free.c
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

#include "saf_free.h"

void free3d(void ***arr, int dim1, int dim2)
{
    int i, j;
    if (arr){
        for (i = 0; i < dim1; i++) {
            for (j = 0; j < dim2; j++){
                free(arr[i][j]);
            }
            free(arr[i]);
        }
        free(arr);
        arr = NULL;
    }
}

//void free2d(void **arr)
//{
//    if(arr!=NULL){
//        //free(arr[0]);
//        free(arr);
//        arr = NULL;
//    }
//}

void free2d(void **arr, int dim1)
{
//#if __STDC_VERSION__ >= 199901L
//    if (arr){
//        free(arr);
//        arr = NULL;
//    }
//#else
    int i;
    if (arr){
        for (i = 0; i < dim1; i++)    {
            free(arr[i]);
        }
        free(arr);
        arr = NULL;
    }
//#endif
}

void free1d(void *arr)
{
    if (arr){
        free(arr);
        arr = NULL;
    }
}

