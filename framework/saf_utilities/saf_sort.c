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
 *     saf_sort.c
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

#include "saf_sort.h"
#include "saf_complex.h" /* TODO: add sort complex numbers, abs/real */

typedef struct saf_sort_int {
    int val;
    int idx;
}saf_sort_int;

typedef struct saf_sort_float {
    float val;
    int idx;
}saf_sort_float;

typedef struct saf_sort_double {
    double val;
    int idx;
}saf_sort_double;

static int cmp_asc_int(const void *a,const void *b) {
    struct saf_sort_int *a1 = (struct saf_sort_int*)a;
    struct saf_sort_int *a2 = (struct saf_sort_int*)b;
    if((*a1).val<(*a2).val)return -1;
    else if((*a1).val>(*a2).val)return 1;
    else return 0;
}

static int cmp_desc_int(const void *a,const void *b) {
    struct saf_sort_int *a1 = (struct saf_sort_int*)a;
    struct saf_sort_int *a2 = (struct saf_sort_int*)b;
    if((*a1).val>(*a2).val)return -1;
    else if((*a1).val<(*a2).val)return 1;
    else return 0;
}

static int cmp_asc_float(const void *a,const void *b) {
    struct saf_sort_float *a1 = (struct saf_sort_float*)a;
    struct saf_sort_float *a2 = (struct saf_sort_float*)b;
    if((*a1).val<(*a2).val)return -1;
    else if((*a1).val>(*a2).val)return 1;
    else return 0;
}

static int cmp_desc_float(const void *a,const void *b) {
    struct saf_sort_float *a1 = (struct saf_sort_float*)a;
    struct saf_sort_float *a2 = (struct saf_sort_float*)b;
    if((*a1).val>(*a2).val)return -1;
    else if((*a1).val<(*a2).val)return 1;
    else return 0;
}

static int cmp_asc_double(const void *a,const void *b) {
    struct saf_sort_double *a1 = (struct saf_sort_double*)a;
    struct saf_sort_double *a2 = (struct saf_sort_double*)b;
    if((*a1).val<(*a2).val)return -1;
    else if((*a1).val>(*a2).val)return 1;
    else return 0;
}

static int cmp_desc_double(const void *a,const void *b) {
    struct saf_sort_double *a1 = (struct saf_sort_double*)a;
    struct saf_sort_double *a2 = (struct saf_sort_double*)b;
    if((*a1).val>(*a2).val)return -1;
    else if((*a1).val<(*a2).val)return 1;
    else return 0;
}

void sorti
(
    int* in_vec,      /* vector[len] to be sorted */
    int* out_vec,     /* if NULL, then in_vec is sorted "in-place" */
    int* new_idices,  /* set to NULL if you don't need them */
    int len,          /* number of elements in vectors, must be consistent with the input data */
    int descendFLAG   /* !1:ascending, 1:descending */
)
{
    int i;
    struct saf_sort_int *data;
    
    data = malloc(len*sizeof(saf_sort_int));
    for(i=0;i<len;i++) {
        data[i].val=in_vec[i];
        data[i].idx=i;
    }
    if(descendFLAG)
        qsort(data,len,sizeof(data[0]),cmp_desc_int);
    else
        qsort(data,len,sizeof(data[0]),cmp_asc_int);
    for(i=0;i<len;i++){
        if (out_vec!=NULL)
            out_vec[i] = data[i].val;
        else
            in_vec[i] = data[i].val; /* overwrite input vector */
        if(new_idices!=NULL)
            new_idices[i] = data[i].idx;
    }
    free(data);
}

void sortf
(
    float* in_vec,    /* vector[len] to be sorted */
    float* out_vec,   /* if NULL, then in_vec is sorted "in-place" */
    int* new_idices,  /* set to NULL if you don't need them */
    int len,          /* number of elements in vectors, must be consistent with the input data */
    int descendFLAG   /* !1:ascending, 1:descending */
)
{
    int i;
    struct saf_sort_float *data;
    
    data = malloc(len*sizeof(saf_sort_float));
    for(i=0;i<len;i++) {
        data[i].val=in_vec[i];
        data[i].idx=i;
    }
    if(descendFLAG)
        qsort(data,len,sizeof(data[0]),cmp_desc_float);
    else
        qsort(data,len,sizeof(data[0]),cmp_asc_float);
    for(i=0;i<len;i++){
        if (out_vec!=NULL)
            out_vec[i] = data[i].val;
        else
            in_vec[i] = data[i].val; /* overwrite input vector */
        if(new_idices!=NULL)
            new_idices[i] = data[i].idx;
    }
    free(data);
}

void sortd
(
    double* in_vec,   /* vector[len] to be sorted */
    double* out_vec,  /* if NULL, then in_vec is sorted "in-place" */
    int* new_idices,  /* set to NULL if you don't need them */
    int len,          /* number of elements in vectors, must be consistent with the input data */
    int descendFLAG   /* !1:ascending, 1:descending */
)
{
    int i;
    struct saf_sort_double *data;
    
    data = malloc(len*sizeof(saf_sort_double));
    for(i=0;i<len;i++) {
        data[i].val=in_vec[i];
        data[i].idx=i;
    }
    if(descendFLAG)
        qsort(data,len,sizeof(data[0]),cmp_desc_double);
    else
        qsort(data,len,sizeof(data[0]),cmp_asc_double);
    for(i=0;i<len;i++){
        if (out_vec!=NULL)
            out_vec[i] = data[i].val;
        else
            in_vec[i] = data[i].val; /* overwrite input vector */
        if(new_idices!=NULL)
            new_idices[i] = data[i].idx;
    }
    free(data);
}

