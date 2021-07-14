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
 * @file md_malloc.h
 * @brief Contiguous memory allocation functions for multi-dimensional arrays
 *
 * Adapted from: https://github.com/leomccormack/md_malloc
 *
 * An example of allocating, indexing and freeing a 3-D "array":
 * \code{.c}
 *   float*** example3D = (float***)malloc3d(10, 20, 5, sizeof(float);
 *   // Due to the contiguous nature of the allocation, this is possible:
 *   memset(FLATTEN3D(example3D), 0, 10*20*5*sizeof(float));
 *   // And it may also be indexed normally as:
 *   example3D[3][19][2] = 22.0f;
 *   // To free, simply call:
 *   free(example3D);
 * \endcode
 *
 * @author Leo McCormack
 * @date 11.06.2019
 * @license MIT
 */

#ifndef MD_MALLOC_INCLUDED
#define MD_MALLOC_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/**
 * Use this macro when passing a 2-D dynamic multi-dimensional array to
 * memset, memcpy or any other function that expects a flat contiguous 1-D block
 * of data
 *
 * e.g.
 * \code{.c}
 *   float** array2D = (float**)malloc2d(10, 40, sizeof(float));
 *   memset(FLATTEN2D(array2D), 0, 10*40*sizeof(float));
 *   // ...
 *   free(array2D);
 * \endcode
 */
#define FLATTEN2D(A) (*A)  /* || (&A[0][0]) */

/**
 * Use this macro when passing a 3-D dynamic multi-dimensional array to
 * memset, memcpy or any other function that expects a flat contiguous 1-D block
 * of data
 */
#define FLATTEN3D(A) (**A) /* || (&A[0][0][0]) */

/**
 * Use this macro when passing a 4-D dynamic multi-dimensional array to
 * memset, memcpy or any other function that expects a flat contiguous 1-D block
 * of data
 */
#define FLATTEN4D(A) (***A) /* || (&A[0][0][0][0]) */

/**
 * Use this macro when passing a 5-D dynamic multi-dimensional array to
 * memset, memcpy or any other function that expects a flat contiguous 1-D block
 * of data
 */
#define FLATTEN5D(A) (****A) /* || (&A[0][0][0][0][0]) */

/**
 * Use this macro when passing a 6-D dynamic multi-dimensional array to
 * memset, memcpy or any other function that expects a flat contiguous 1-D block
 * of data
 */
#define FLATTEN6D(A) (*****A) /* || (&A[0][0][0][0][0][0]) */
    
/** 1-D malloc (same as malloc, but with error checking) */
void* malloc1d(size_t dim1_data_size);

/** 1-D calloc (same as calloc, but with error checking) */
void* calloc1d(size_t dim1, size_t data_size);

/** 1-D realloc (same as realloc, but with error checking) */
void* realloc1d(void* ptr, size_t dim1_data_size);

/** 2-D malloc (contiguously allocated, so use free() as usual to deallocate) */
void** malloc2d(size_t dim1, size_t dim2, size_t data_size);

/** 2-D calloc (contiguously allocated, so use free() as usual to deallocate) */
void** calloc2d(size_t dim1, size_t dim2, size_t data_size);

/** 2-D realloc which does NOT retain previous data order */
void** realloc2d(void** ptr, size_t dim1, size_t dim2, size_t data_size);

/**
 * 2-D realloc which does retain previous data order
 *
 * @test test__realloc2d_r()
 */
void** realloc2d_r(void** ptr, size_t new_dim1, size_t new_dim2,
                   size_t prev_dim1, size_t prev_dim2, size_t data_size);
 
/** 3-D malloc (contiguously allocated, so use free() as usual to deallocate) */
void*** malloc3d(size_t dim1, size_t dim2, size_t dim3, size_t data_size);

/** 3-D calloc (contiguously allocated, so use free() as usual to deallocate) */
void*** calloc3d(size_t dim1, size_t dim2, size_t dim3, size_t data_size);

/** 3-D realloc which does NOT retain previous data order */
void*** realloc3d(void*** ptr, size_t dim1, size_t dim2, size_t dim3,
                  size_t data_size);

/** 3-D realloc which does retain previous data order */
void*** realloc3d_r(void*** ptr, size_t new_dim1, size_t new_dim2,
                    size_t new_dim3, size_t prev_dim1, size_t prev_dim2,
                    size_t prev_dim3, size_t data_size);

/** 4-D malloc (contiguously allocated, so use free() as usual to deallocate) */
void**** malloc4d(size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                  size_t data_size);

/** 4-D calloc (contiguously allocated, so use free() as usual to deallocate) */
void**** calloc4d(size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                  size_t data_size);

/** 4-D realloc which does NOT retain previous data order */
void**** realloc4d(void**** ptr, size_t dim1, size_t dim2, size_t dim3,
                   size_t dim4, size_t data_size);

/** 5-D malloc (contiguously allocated, so use free() as usual to deallocate) */
void***** malloc5d(size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                   size_t dim5, size_t data_size);

/** 5-D calloc (contiguously allocated, so use free() as usual to deallocate) */
void***** calloc5d(size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                   size_t dim5, size_t data_size);

/** 5-D realloc which does NOT retain previous data order */
void***** realloc5d(void***** ptr, size_t dim1, size_t dim2, size_t dim3,
                    size_t dim4, size_t dim5, size_t data_size);

/** 6-D malloc (contiguously allocated, so use free() as usual to deallocate) */
void****** malloc6d(size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                    size_t dim5, size_t dim6, size_t data_size);

/** 6-D malloc (contiguously allocated, so use free() as usual to deallocate) */
void****** calloc6d(size_t dim1, size_t dim2, size_t dim3, size_t dim4,
                    size_t dim5, size_t dim6, size_t data_size);

/** 6-D realloc which does NOT retain previous data order */
void****** realloc6d(void****** ptr, size_t dim1, size_t dim2, size_t dim3,
                     size_t dim4, size_t dim5, size_t dim6, size_t data_size);


#ifdef __cplusplus
} /*extern "C"*/
#endif /* __cplusplus */

#endif /* MD_MALLOC_INCLUDED */
