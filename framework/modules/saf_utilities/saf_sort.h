/*
 * Copyright 2016-2018 Leo McCormack
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 * PERFORMANCE OF THIS SOFTWARE.
 */

/*
 * Filename: saf_sort.h
 * --------------------
 * Contains some useful sorting functions
 *
 * Dependencies:
 *     none
 * Author, date created:
 *     Leo McCormack, 30.07.2018
 */

#ifndef SAF_SORT_H_INCLUDED
#define SAF_SORT_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
    
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: sorti
 * ---------------
 * Sort a vector of integer values into ascending/decending order. Optionally
 * returning the new indices.
 *
 * Input/Output Arguments:
 *     in_vec      - vector to be sorted; len x 1
 *     out_vec     - output vector. If NULL, then in_vec is sorted "in-place"
 *     new_idices  - indices to sort in_vec (set to NULL if you don't need them)
 *     len         - number of elements in vectors
 *     descendFLAG - !1:ascending, 1:descending
 */
void sorti(int* in_vec,
           int* out_vec,
           int* new_idices,
           int len,
           int descendFLAG);
    
/*
 * Function: sortf
 * ---------------
 * Sort a vector of floating-point values into ascending/decending order.
 * Optionally returning the new indices.
 *
 * Input/Output Arguments:
 *     in_vec      - vector to be sorted; len x 1
 *     out_vec     - output vector. If NULL, then in_vec is sorted "in-place"
 *     new_idices  - indices to sort in_vec (set to NULL if you don't need them)
 *     len         - number of elements in vectors
 *     descendFLAG - !1:ascending, 1:descending
 */
void sortf(float* in_vec,
           float* out_vec,
           int* new_idices,
           int len,
           int descendFLAG);

/*
 * Function: sortd
 * ---------------
 * Sort a vector of double floating-point values into ascending/decending order.
 * Optionally returning the new indices.
 *
 * Input/Output Arguments:
 *     in_vec      - vector to be sorted; len x 1
 *     out_vec     - output vector. If NULL, then in_vec is sorted "in-place"
 *     new_idices  - indices to sort in_vec (set to NULL if you don't need them)
 *     len         - number of elements in vectors
 *     descendFLAG - !1:ascending, 1:descending
 */
void sortd(double* in_vec,
           double* out_vec,
           int* new_idices,
           int len,
           int descendFLAG);

/*
 * Function: findClosestGridPoints
 * -------------------------------
 * Finds indicies into "grid_dirs" that are the closest to "target dirs". e.g.
 *     grid_dirs[idx_closest[0]] will be the closest direction in "grid_dirs"
 *     to target_dirs[0].
 *
 * Input Arguments:
 *     grid_dirs    - sph coordinates of grid directions; FLAT: nGrid x 2
 *     nGrid        - number of directions in grid
 *     target_dirs  - sph coordinates of target directions; FLAT: nTarget x 2
 *     nTarget      - number of target directions to find
 *     degFLAG      - 0: coordinates are in RADIANS, 1: coords are in DEGREES
 * Output Arguments:
 *     idx_closest  - resulting indices (set to NULL to ignore); nTarget x 1
 *     dirs_closest - grid_dirs(idx_closest); (set to NULL to ignore);
 *                    nTarget x 1
 *     angle_diff   - angle diff btwn target + grid dir, in degrees (set to NULL
 *                    to ignore); nTarget x 1
 */
void findClosestGridPoints(float* grid_dirs,
                           int nGrid,
                           float* target_dirs,
                           int nTarget,
                           int degFLAG,
                           int* idx_closest,
                           float* dirs_closest,
                           float* angle_diff);
    
    
#ifdef __cplusplus
}/* extern "C" */
#endif /* __cplusplus */

#endif /* SAF_SORT_H_INCLUDED */
