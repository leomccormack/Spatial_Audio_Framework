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

/**
 *@addtogroup Utilities
 *@{
 * @file saf_utility_sort.h
 * @brief A collection of useful sorting functions
 *
 * @author Leo McCormack
 * @date 30.07.2018
 * @license ISC
 */

#ifndef SAF_SORT_H_INCLUDED
#define SAF_SORT_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Sort a vector of integer values into ascending/decending order (optionally
 * returning the new indices as well)
 *
 * @param [in]  in_vec      Vector to be sorted; len x 1
 * @param [out] out_vec     Output vector (set to NULL if you don't want it)
 * @param [out] new_indices Indices used to sort 'in_vec' (set to NULL if you
 *                          don't want them)
 * @param [in]  len         Number of elements in vectors
 * @param [in]  descendFLAG '0' ascending, '1' descending
 */
void sorti(int* in_vec,
           int* out_vec,
           int* new_indices,
           int len,
           int descendFLAG);
    
/**
 * Sort a vector of floating-point values into ascending/decending order
 * (optionally returning the new indices as well)
 *
 * @test test__sortf()
 *
 * @param [in]  in_vec      Vector to be sorted; len x 1
 * @param [out] out_vec     Output vector (set to NULL if you don't want it)
 * @param [out] new_indices Indices used to sort 'in_vec' (set to NULL if you
 *                          don't want them)
 * @param [in]  len         Number of elements in vectors
 * @param [in]  descendFLAG '0' ascending, '1' descending
 */
void sortf(float* in_vec,
           float* out_vec,
           int* new_indices,
           int len,
           int descendFLAG);

/**
 * Sort a vector of double floating-point values into ascending/decending order
 * (optionally returning the new indices as well)
 *
 * @param [in]  in_vec      Vector to be sorted; len x 1
 * @param [out] out_vec     Output vector (set to NULL if you don't want it)
 * @param [out] new_indices Indices used to sort 'in_vec' (set to NULL if you
 *                          don't want them)
 * @param [in]  len         Number of elements in vectors
 * @param [in]  descendFLAG '0' ascending, '1' descending
 */
void sortd(double* in_vec,
           double* out_vec,
           int* new_indices,
           int len,
           int descendFLAG);

/**
 * Sort a vector of complex floating-point values into ascending/decending order
 *
 * @note The values are first sorted based on their real parts. Values with
 *       identical real parts are then sorted based on their imaginary parts.
 *
 * @param [in]  in_vec      Vector to be sorted; len x 1
 * @param [out] out_vec     Output vector
 * @param [in]  len         Number of elements in vectors
 * @param [in]  descendFLAG '0' ascending, '1' descending
 */
void sortc(float_complex* in_vec,
           float_complex* out_vec,
           int len,
           int descendFLAG);

/**
 * Sort a vector of complex double floating-point values into ascending/
 * decending order
 *
 * @note The values are first sorted based on their real parts. Values with
 *       identical real parts are then sorted based on their imaginary parts.
 *
 * @test test__sortz()
 *
 * @param [in]  in_vec      Vector to be sorted; len x 1
 * @param [out] out_vec     Output vector
 * @param [in]  len         Number of elements in vectors
 * @param [in]  descendFLAG '0' ascending, '1' descending
 */
void sortz(double_complex* in_vec,
           double_complex* out_vec,
           int len,
           int descendFLAG);

/**
 * Pairs up complex numbers and sorts them in ascending order based on their
 * real parts first, and then on their imaginary parts
 *
 * @note This function is the same as sortz() except that any values that are
 *       purely real, are pushed to the end of the output vector (and also in
 *       ascending order).
 *
 * @test test__cmplxPairUp()
 *
 * @param [in]  in_vec      Vector to be sorted; len x 1
 * @param [out] out_vec     Output vector
 * @param [in]  len         Number of elements in vectors
 */
void cmplxPairUp(double_complex* in_vec,
                 double_complex* out_vec,
                 int len);

/**
 * Finds indicies into "grid_dirs" that are the closest to "target dirs".
 *
 * e.g. grid_dirs[idx_closest[0]] will be the closest direction in "grid_dirs"
 * to target_dirs[0].
 *
 * @param [in]  grid_dirs    Spherical coordinates of grid directions;
 *                           FLAT: nGrid x 2
 * @param [in]  nGrid        Number of directions in grid
 * @param [in]  target_dirs  Spherical coordinates of target directions;
 *                           FLAT: nTarget x 2
 * @param [in]  nTarget      Number of target directions to find
 * @param [in]  degFLAG      '0' coordinates are in RADIANS, '1' coords are in
 *                           DEGREES
 * @param [out] idx_closest  Resulting indices (set to NULL to ignore);
 *                           nTarget x 1
 * @param [out] dirs_closest grid_dirs(idx_closest); (set to NULL to ignore);
 *                           nTarget x 1
 * @param [out] angle_diff   Angle diff between target a d grid dir, in degrees
 *                           (set to NULL to ignore); nTarget x 1
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

/**@} */ /* doxygen addtogroup Utilities */
