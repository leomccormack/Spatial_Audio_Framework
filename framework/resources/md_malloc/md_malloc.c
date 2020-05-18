/*
 Copyright (c) 2019 Leo McCormack
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

/**
 * @file md_malloc.c
 * @brief Implementations of dynamic memory allocation functions for
 *        contiguous multidimensional "arrays"
 *
 * Taken from: https://github.com/leomccormack/md_malloc
 *
 * An example of allocating, indexing and freeing a 3-D "array":
 * \code{.c}
 *   float*** example3D = (float***)malloc3d(10, 20, 5, sizeof(float);
 *   // Due to the contiguous nature of the allocation, this is possible:
 *   memset(ADR3D(example3D), 0, 10*20*5*sizeof(float));
 *   // And my still be indexed normally as:
 *   example3D[3][19][2] = 22.0f;
 *   // To free, simply call:
 *   free(example3D);
 * \endcode
 *
 * @author Leo McCormack
 * @date 11.06.2019
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "md_malloc.h"

#ifndef MIN
# define MIN(a,b) (( (a) < (b) ) ? (a) : (b))
#endif
#ifndef MAX
# define MAX(a,b) (( (a) > (b) ) ? (a) : (b))
#endif

void* malloc1d(size_t dim1_data_size)
{
    void *ptr = malloc(dim1_data_size);
#ifndef NDEBUG
    if (ptr == NULL)
        fprintf(stderr, "Error: 'malloc1d' failed to allocate %zu bytes.\n", dim1_data_size);
#endif
    return ptr;
}

void* calloc1d(size_t dim1, size_t data_size)
{
    void *ptr = calloc(dim1, data_size);
#ifndef NDEBUG
    if (ptr == NULL)
        fprintf(stderr, "Error: 'calloc1d' failed to allocate %zu bytes.\n", dim1*data_size);
#endif
    return ptr;
}

void* realloc1d(void* ptr, size_t dim1_data_size)
{
    ptr = realloc(ptr, dim1_data_size);
#ifndef NDEBUG
    if (ptr == NULL)
        fprintf(stderr, "Error: 'realloc1d' failed to allocate %zu bytes.\n", dim1_data_size);
#endif
    return ptr;
}

void** malloc2d(size_t dim1, size_t dim2, size_t data_size)
{
    size_t i, stride;
    void** ptr;
    unsigned char* p2;
    stride = dim2*data_size;
    ptr = malloc(dim1*sizeof(void*) + dim1*stride);
#ifndef NDEBUG
    if(ptr==NULL)
        fprintf(stderr, "Error: 'malloc2d' failed to allocate %zu bytes.\n", dim1*sizeof(void*) + dim1*stride);
#endif
    p2 = (unsigned char*)(ptr + dim1);
    for(i=0; i<dim1; i++)
        ptr[i] = &p2[i*stride];
    return ptr;
}

void** calloc2d(size_t dim1, size_t dim2, size_t data_size)
{
    size_t i, stride;
    void** ptr;
    unsigned char* p2;
    stride = dim2*data_size;
    ptr = calloc(dim1, sizeof(void*) + stride);
#ifndef NDEBUG
    if(ptr==NULL)
        fprintf(stderr, "Error: 'calloc2d' failed to allocate %zu bytes.\n", dim1*sizeof(void*) + dim1*stride);
#endif
    p2 = (unsigned char*)(ptr + dim1);
    for(i=0; i<dim1; i++)
        ptr[i] = &p2[i*stride];
    return ptr;
}

void** realloc2d(void** ptr, size_t dim1, size_t dim2, size_t data_size)
{
    size_t i, stride;
    unsigned char* p2;
    stride = dim2*data_size;
    ptr = realloc(ptr, dim1*sizeof(void*) + dim1*stride);
#ifndef NDEBUG
    if(ptr==NULL)
        fprintf(stderr, "Error: 'realloc2d' failed to allocate %zu bytes.\n", dim1*sizeof(void*) + dim1*stride);
#endif
    p2 = (unsigned char*)(ptr + dim1);
    for(i=0;i<dim1;i++)
        ptr[i] = &p2[i*stride];
    return ptr;
}

void** realloc2d_r(void** ptr, size_t new_dim1, size_t new_dim2, size_t prev_dim1, size_t prev_dim2, size_t data_size)
{
    size_t i, stride;
    void** prev_data;

    /* Copy previous data */
    prev_data = malloc2d(prev_dim1, prev_dim2, data_size);
    memcpy(ADR2D(prev_data), ADR2D(ptr), prev_dim1*prev_dim2*data_size);

    /* Resize */
    unsigned char* p2;
    stride = new_dim2*data_size;
    ptr = realloc(ptr, new_dim1*sizeof(void*) + new_dim1*stride);
#ifndef NDEBUG
    if(ptr==NULL)
        fprintf(stderr, "Error: 'realloc2d' failed to allocate %zu bytes.\n", new_dim1*sizeof(void*) + new_dim1*stride);
#endif
    p2 = (unsigned char*)(ptr + new_dim1);
    for(i=0;i<new_dim1;i++)
        ptr[i] = &p2[i*stride];

    /* Repopulate */
    for(i=0; i<MIN(new_dim1,prev_dim1); i++)
        memcpy(ptr[i], prev_data[i], MIN(new_dim2,prev_dim2)*data_size);
    free(prev_data);
    return ptr;
}

void*** malloc3d(size_t dim1, size_t dim2, size_t dim3, size_t data_size)
{
    size_t i, j, stride1, stride2;
    void*** ptr;
    void** p2;
    unsigned char* p3;
    stride1 = dim2*dim3*data_size;
    stride2 = dim3*data_size;
    ptr = malloc(dim1*sizeof(void**) + dim1*dim2*sizeof(void*) + dim1*stride1);
#ifndef NDEBUG
    if(ptr==NULL)
        fprintf(stderr, "Error: 'malloc3d' failed to allocate %zu bytes.\n", dim1*sizeof(void**) + dim1*dim2*sizeof(void*) + dim1*stride1);
#endif
    p2 = (void**)(ptr + dim1);
    p3 = (unsigned char*)(p2 + dim1*dim2);
    for(i=0;i<dim1;i++)
        ptr[i] = &p2[i*dim2];
    for(i=0;i<dim1;i++)
        for(j=0;j<dim2;j++)
            p2[i*dim2+j] = &p3[i*stride1 + j*stride2];
    return ptr;
}

void*** calloc3d(size_t dim1, size_t dim2, size_t dim3, size_t data_size)
{
    size_t i, j, stride1, stride2;
    void*** ptr;
    void** p2;
    unsigned char* p3;
    stride1 = dim2*dim3*data_size;
    stride2 = dim3*data_size;
    ptr = calloc(dim1, sizeof(void**) + dim2*sizeof(void*) + stride1);
#ifndef NDEBUG
    if(ptr==NULL)
        fprintf(stderr, "Error: 'calloc3d' failed to allocate %zu bytes.\n", dim1*sizeof(void**) + dim1*dim2*sizeof(void*) + dim1*stride1);
#endif
    p2 = (void**)(ptr + dim1);
    p3 = (unsigned char*)(p2 + dim1*dim2);
    for(i=0;i<dim1;i++)
        ptr[i] = &p2[i*dim2];
    for(i=0;i<dim1;i++)
        for(j=0;j<dim2;j++)
            p2[i*dim2+j] = &p3[i*stride1 + j*stride2];
    return ptr;
}

void*** realloc3d(void*** ptr, size_t new_dim1, size_t new_dim2, size_t new_dim3, size_t data_size)
{
    size_t i, j, stride1, stride2;
    void** p2;
    unsigned char* p3;
    stride1 = new_dim2*new_dim3*data_size;
    stride2 = new_dim3*data_size;
    ptr = realloc(ptr, new_dim1*sizeof(void**) + new_dim1*new_dim2*sizeof(void*) + new_dim1*stride1);
#ifndef NDEBUG
    if(ptr==NULL)
        fprintf(stderr, "Error: 'realloc3d' failed to allocate %zu bytes.\n", new_dim1*sizeof(void**) + new_dim1*new_dim2*sizeof(void*) + new_dim1*stride1);
#endif
    p2 = (void**)(ptr + new_dim1);
    p3 = (unsigned char*)(p2 + new_dim1*new_dim2);
    for(i=0;i<new_dim1;i++)
        ptr[i] = &p2[i*new_dim2];
    for(i=0;i<new_dim1;i++)
        for(j=0;j<new_dim2;j++)
            p2[i*new_dim2+j] = &p3[i*stride1 + j*stride2];
    return ptr;
}

void*** realloc3d_r(void*** ptr, size_t new_dim1, size_t new_dim2, size_t new_dim3, size_t prev_dim1, size_t prev_dim2, size_t prev_dim3, size_t data_size)
{
    size_t i, j, stride1, stride2;
    void** p2;
    unsigned char* p3;
    void*** prev_data;

    /* Copy previous data */
    prev_data = malloc3d(prev_dim1, prev_dim2, prev_dim3, data_size);
    memcpy(ADR3D(prev_data), ADR3D(ptr), prev_dim1*prev_dim2*prev_dim3*data_size);

    /* Resize */
    stride1 = new_dim2*new_dim3*data_size;
    stride2 = new_dim3*data_size;
    ptr = realloc(ptr, new_dim1*sizeof(void**) + new_dim1*new_dim2*sizeof(void*) + new_dim1*stride1);
#ifndef NDEBUG
    if(ptr==NULL)
        fprintf(stderr, "Error: 'realloc3d' failed to allocate %zu bytes.\n", new_dim1*sizeof(void**) + new_dim1*new_dim2*sizeof(void*) + new_dim1*stride1);
#endif
    p2 = (void**)(ptr + new_dim1);
    p3 = (unsigned char*)(p2 + new_dim1*new_dim2);
    for(i=0;i<new_dim1;i++)
        ptr[i] = &p2[i*new_dim2];
    for(i=0;i<new_dim1;i++)
        for(j=0;j<new_dim2;j++)
            p2[i*new_dim2+j] = &p3[i*stride1 + j*stride2];

    /* Repopulate */
    for(i=0; i<MIN(new_dim1,prev_dim1); i++)
        for(j=0; j<MIN(new_dim2,prev_dim2); j++)
            memcpy(ptr[i][j], prev_data[i][j], MIN(new_dim3,prev_dim3)*data_size);
    free(prev_data);
    return ptr;
}

