/*
 * Copyright 2017-2018 Leo McCormack
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
 * @file rotator_internal.c
 * @brief A basic spherical harmonic/ Ambisonic signals rotator, based on the
 *        recursive approach detailed in [1]
 *
 * @test test__saf_example_rotator()
 *
 * @see [1] Ivanic, J., Ruedenberg, K. (1998). Rotation Matrices for Real
 *          Spherical Harmonics. Direct Determination by Recursion Page:
 *          Additions and Corrections. Journal of Physical Chemistry A, 102(45),
 *          9099?9100.
 *
 * @author Leo McCormack
 * @date 02.11.2017
 * @license ISC
 */

#include "rotator.h"
#include "rotator_internal.h"


