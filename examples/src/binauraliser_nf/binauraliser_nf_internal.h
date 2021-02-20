/*
 * Copyright 2017-2018 Michael McCrea, Leo McCormack
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
 * @file: binauraliser_nf_internal.h
 * @brief Convolves input audio (up to 64 channels) with interpolated HRTFs in
 *        the time-frequency domain, and applies optional near-field binaural
 *        filtering, as described in [1].
 *
 * The HRTFs are interpolated by applying amplitude-preserving VBAP gains to the
 * HRTF magnitude responses and inter-aural time differences (ITDs)
 * individually, before being re-combined. The example also allows the user to
 * specify an external SOFA file for the convolution, and rotations of the
 * source directions to accomodate head-tracking.
 *
 * @see [1] S. Spagnol, E. Tavazzi, and F. Avanzini, “Distance rendering and
 *          perception of nearby virtual sound sources with a near-field filter
 *          model,” Applied Acoustics, vol. 115, pp. 61–73, Jan. 2017,
 *          doi: 10.1016/j.apacoust.2016.08.015.
 *
 * @author Michael McCrea, Leo McCormack
 * @date 22.02.2022
 * @license ISC
 */

#ifndef __BINAURALISER_NF_INTERNAL_H_INCLUDED__
#define __BINAURALISER_NF_INTERNAL_H_INCLUDED__

#include "../binauraliser/binauraliser_internal.h"
#include <binauraliser_nf.h>    /* Include header for this example */
#include "saf.h"                /* Main include header for SAF */
#include "saf_externals.h"      /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for binauraliserNF. Contains all variables from binauraliser
 * (audio buffers, afSTFT, HRTFs, internal variables, flags, user parameters) plus
 * those specific to the near field variant.
 */
typedef struct _binauraliserNF {
    struct _binauraliser;               /**< "inherit" member vars of binauraliser struct */

    /* audio buffers */
    float**          binsrcsTD;         /**< near field DVF-filtered sources frame; (#MAX_NUM_INPUTS * #NUM_EARS) x #BINAURALISER_FRAME_SIZE */
    float_complex*** binauralTF;        /**< time-frequency domain output frame; #HYBRID_BANDS x (#MAX_NUM_INPUTS * #NUM_EARS) x #TIME_SLOTS
                                         *   Note: (#MAX_NUM_INPUTS * #NUM_EARS) dimensions are combined because afSTFT requires type float_complex***
                                         */
    /* misc. */
    float src_dists_m[MAX_NUM_INPUTS];  /**< source distance,  meters */
    float farfield_thresh_m;            /**< distance considered to be far field (no near field filtering),  meters */
    float farfield_headroom;            /**< scale factor applied to farfield_thresh_m when resetting to the far field, and for UI range, meters */
    float nearfield_limit_m;            /**< minimum distance allowed for near-field filtering, from head _center_,  meters, def. 0.15 */
    float head_radius;                  /**< head radius, used calculate normalized source distance meters, def. 0.09096 */
    float head_radius_recip;            /**< reciprocal of head radius */

} binauraliserNF_data;

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Initialise the filterbank used by binauraliserNF.
 *
 * @note Call this function before binauraliser_initHRTFsAndGainTables()
 */
void binauraliserNF_initTFT(void* const hBin);

/**
 * Resets the source distances to the default far field distance.
 */
void binauraliserNF_resetSourceDistances(void* const hBin);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __BINAURALISER_INTERNAL_NF_H_INCLUDED__ */
