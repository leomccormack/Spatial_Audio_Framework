(https://github.com/leomccormack/md_malloc)

# md_malloc

A header file comprising functions for dynamically allocating contiguous multi-dimensional arrays. The code is written in C and is MSVC compliant.

## Motivation

The main reasons for releasing these contiguous multi-dimensional memory allocation functions as a convenient header-only file are threefold. 

The first is to avoid (and hopefully discourage) the wide-spread adoption of non-contiguous memory allocation:
```c
/* For example, commonly, allocating an A x B, 2-D "array" is written as: */
float** array2D = malloc(A*sizeof(float*));
for(int i=0; i<A; i++)
    array2D[i] = malloc(B*sizeof(float));
    
/* This 2-D array may indeed be easily indexed (i,j) array2D[i][j]. However, since it 
 * employs multiple malloc calls, the second dimension vectors are not guaranteed to be
 * adjacent to one-another in memory i.e. the memory is not necessarily "contiguous". 
 * Therefore, calls to memset/memcpy or other functions which rely on contiguous memory
 * allocation may result in crashes or behave unpredictably. */
 
memset(&array2D[0][0], 0, A*B*sizeof(float)); /* could be fine, but possibly not */

/* freeing this 2-D "array" is also slightly cumbersome: */
for(int i=0; i<A; i++)
    free(array2D[i]);
free(array2D);
```

Second, is to still retain easy array indexing:
```c
/* One may, allocate a contiguous A x B, 2-D "array" with as single malloc call: */
float* array2D = malloc(A*B*sizeof(float));

/* However, these are also not ideal, as they must be indexed (i,j) as array2D[i*B+j], 
 * which gets messy for high-dimensional arrays; 
 * e.g. (i,j,k,l,p), array5D[i*B*C*D*E + j*C*D*E + k*D*E + l*E +p] */
```

And thirdly, to maintain easy passing of N-D arrays to other functions:
```c
/* Since C99, multidimensional arrays can be declared as */
float (*array2D)[B] = malloc(sizeof(float[A][B]));

/* These can be indexed "array2D[i][j]" and freed "free(array2D)" easily. However, 
 * these cannot be easily and flexibly passed to generalised functions due to the 
 * variable type. */
 
/* Furthermore, for cross-platform project development, supporting the ancient MSVC 
 * compiler is (unfortunately) a common limitation. 
 * Since Microsoft's C-compiler is still based on C89/C90, this C99-style may not be 
 * considered as a viable option regardless. */
```


## Getting Started

To include md_malloc in your project, simply add the following:

```c
#define MD_MALLOC_ENABLE
#include "md_malloc.h"
```
Arrays may be allocated, resized, and freed as:

```c
/* To allocate an A x B, 2-D "array" */
float** array2D = (float**)malloc2d(A, B, sizeof(float));

/* To resize: */
array2D = (float**)realloc2d(array2D, C, D, sizeof(float));

/* Note that "array2D" is: 1) contiguous, 2) easily indexable, and 3) compiles under MSVC */

memset(&array2D[0][0], 0, C*D*sizeof(float)); 
array2D[i][j] = 66.6f;

/* To free: */
free(array2D);

```

## Testing

This project also includes a test/test.c file, which checks for array contiguity by using memcpy and CBLAS calls on md_malloc allocated arrays and subsequently comparing their results to that of their static memory counterparts.
The test file also compares the time taken when using md_array, mangled arrays (array2d[i*dim2+j]), and C-99 style arrays. 

TLDR: md_malloc has roughly the same speed performance for 2-D arrays. Whereas, 3-D arrays begin to show a slow-down. It is up to you and the nature of your project, as to whether this is an acceptable compromise. Here is the test output when allocating, indexing, populating, and freeing a 290 x 300 2-D array 30 000 times, and a 290 x 300 x 295 3-D array 100 times:

``` 
********** Malloc Speed Test - 2D DATA **********
- Mangled array time taken 13 seconds 803 milliseconds
- md_malloc time taken 13 seconds 684 milliseconds
- C99-style malloc time taken 13 seconds 603 milliseconds
md_malloc was 0 seconds 119 milliseconds FASTER than Mangled array
md_malloc was 0 seconds 81 milliseconds SLOWER than C99-style array

********** Malloc Speed Test - 3D DATA **********
- Mangled array time taken 13 seconds 478 milliseconds
- md_malloc time taken 14 seconds 578 milliseconds
- C99-style malloc time taken 13 seconds 426 milliseconds
md_malloc was 1 seconds 100 milliseconds SLOWER than Mangled array
md_malloc was 1 seconds 152 milliseconds SLOWER than C99-style array

Program ended with exit code: 0
```
(Using a 2,9 GHz Intel Core i7 Skylake 6920H CPU, Xcode 10.2 - Apple Clang compiler -Os)


## References

The creation of this mini-project was inspired by discussions on stackoverflow regarding the topic, which can be found here:
* https://stackoverflow.com/questions/42094465/correctly-allocating-multi-dimensional-arrays
* https://stackoverflow.com/questions/4016842/dynamic-contiguous-3d-arrays-in-c
* https://stackoverflow.com/questions/35385154/dynamically-allocate-contiguous-block-of-memory-for-2d-array-of-unknown-data-typ


## Contact

Comments, bug reports, and suggested improvements are very much welcomed: leo.mccormack(at)aalto.fi

## License

The code is distributed under the MIT license.

