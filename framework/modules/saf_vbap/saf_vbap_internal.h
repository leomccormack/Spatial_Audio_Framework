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
 * @file saf_vbap_internal.h
 * @ingroup VBAP
 * @brief Internal header for the VBAP/MDAP module (#SAF_VBAP_MODULE)
 *
 * VBAP functions largely derived from the MATLAB library found in [1].
 *
 * @see [1] https://github.com/polarch/Vector-Base-Amplitude-Panning
 *          Copyright (c) 2015, Archontis Politis, BSD-3-Clause License
 *
 * @author Leo McCormack
 * @date 02.10.2017
 * @license ISC
 */

#ifndef __SAF_VBAP_INTERNAL_H_INCLUDED__
#define __SAF_VBAP_INTERNAL_H_INCLUDED__

#include "saf_vbap.h"
#include "../saf_utilities/saf_utilities.h"
#include "saf_externals.h" 

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 * In degrees, if no ls_dirs have elevation +/- this value, dummies are placed
 * at +/- 90 elevation.  */
#define ADD_DUMMY_LIMIT ( 60.0f )
/**
 * if omitLargeTriangles==1, triangles with an aperture larger than this are
 * discarded */
#define APERTURE_LIMIT_DEG ( 180.0f )

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Cross product between two 3-element floating point vectors
 */
void ccross(float a[3], float b[3], float c[3]);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_VBAP_INTERNAL_H_INCLUDED__ */
