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
/*
 * Filename:
 *     md_malloc.h
 * Description:
 *     Implementations of dynamic memory allocation functions for
 *     contiguous multidimensional "arrays".
 *     Taken from: https://github.com/leomccormack/md_malloc
 * Example:
 *     An example of allocating, indexing and freeing a 3-D "array":
 *       float*** example3D = (float***)malloc3d(10, 20, 5, sizeof(float);
 *       // Due to the contiguous nature of the allocation, this is possible:
 *       memset(&example3D[0][0][0], 0, 10*20*5*sizeof(float));
 *       // And my still be indexed normally as:
 *       example3D[3][19][2] = 22.0f;
 *       // To free, simply call:
 *       free(example3D);
 *       // Or:
 *       free3d((void***)example3D);
 * Dependencies:
 *     (none)
 * Author, date created:
 *     Leo McCormack, 11.06.2019
 */

#ifndef MD_MALLOC_INCLUDED
#define MD_MALLOC_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif
    
/* These macros return a pointer to the address of the first element in the
 * array. Use them for passing arrays to memset/memcpy, or blas/lapack functions
 * etc. e.g.
 *   float** array2D = (float**)malloc2d(10, 40, sizeof(float));
 *   memset(ADR2D(array2D), 0, 10*40*sizeof(float)); */
#define ADR1D(A) (&A[0])
#define ADR2D(A) (&A[0][0])
#define ADR3D(A) (&A[0][0][0])
    
/* 1-D */
void* malloc1d(size_t dim1_data_size);
void* calloc1d(size_t dim1, size_t data_size);
void* realloc1d(void* ptr, size_t dim1_data_size);
void free1d(void** ptr);
    
/* 2-D */
void** malloc2d(size_t dim1, size_t dim2, size_t data_size);
void** calloc2d(size_t dim1, size_t dim2, size_t data_size);
void** realloc2d(void** ptr, size_t dim1, size_t dim2, size_t data_size);
void free2d(void*** ptr);

/* 3-D */
void*** malloc3d(size_t dim1, size_t dim2, size_t dim3, size_t data_size);
void*** calloc3d(size_t dim1, size_t dim2, size_t dim3, size_t data_size);
void*** realloc3d(void*** ptr, size_t dim1, size_t dim2, size_t dim3, size_t data_size);
void free3d(void**** ptr);

    
#ifdef __cplusplus
} /*extern "C"*/
#endif

#endif /* MD_MALLOC_INCLUDED */
