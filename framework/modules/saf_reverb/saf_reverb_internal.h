/*
 * Copyright 2020 Leo McCormack
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
 * @file saf_reverb_internal.h
 * @ingroup Reverb
 * @brief Internal header for the reverb processing module (#SAF_REVERB_MODULE)
 *
 * A collection of reverb and room simulation algorithms.
 *
 * @author Leo McCormack
 * @date 06.05.2020
 * @license ISC
 */

#ifndef __REVERB_INTERNAL_H_INCLUDED__
#define __REVERB_INTERNAL_H_INCLUDED__

#include "saf_reverb.h"
#include "saf_externals.h"
#include "../saf_utilities/saf_utilities.h"
#include "../saf_sh/saf_sh.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                         IMS Shoebox Room Simulator                         */
/* ========================================================================== */

/** Number of wall for a shoebox room */
#define IMS_NUM_WALLS_SHOEBOX ( 6 )
/** FIR filter order (must be even) */
#define IMS_FIR_FILTERBANK_ORDER ( 400 )
/** IIR filter order (1st or 3rd) */
#define IMS_IIR_FILTERBANK_ORDER ( 3 )
/** Circular buffer length */
#define IMS_CIRC_BUFFER_LENGTH ( 4*8192U )
/** Circular buffer length, minus 1 */
#define IMS_CIRC_BUFFER_LENGTH_MASK ( IMS_CIRC_BUFFER_LENGTH - 1U )
/** Maximum number of samples that ims should expect to process at a time */
#define IMS_MAX_NSAMPLES_PER_FRAME ( 20000 )
/** Order of lagrange interpolation filters */
#define IMS_LAGRANGE_ORDER ( 2 )
/** Lagrange interpolator look-up table size */
#define IMS_LAGRANGE_LOOKUP_TABLE_SIZE ( 100 )
/** Index for the current echogram slot */
#define IMS_EG_CURRENT ( 0 )
/** Index for the previous echogram slot */
#define IMS_EG_PREV ( 1 )
/** Number of echogram slots */
#define IMS_EG_NUM_SLOTS ( 2 )
/** While a source or receiver ID is yet to be active, it is "IMS_UNASSIGNED" */
#define IMS_UNASSIGNED ( -1 )

/** Void pointer (just to improve code readability when working with arrays of
 *  handles) */
typedef void* voidPtr;

/** Union struct for Cartesian coordinates (access as .x,.y,.z, or .v[3]) */
typedef struct _ims_pos_xyz {
    union {
        struct { float x; /**< x Cartesian coordinate, in metres */
                 float y; /**< y Cartesian coordinate, in metres */
                 float z; /**< z Cartesian coordinate, in metres */
        };
        float v[3]; /**< x,y,z Cartesian coordinates, in metres */
    };
} ims_pos_xyz;

/** Supported receiver types */
typedef enum {
    RECEIVER_SH   /**< Spherical harmonic receiver */
    // RECEIVER_ARRAY /**< Microphone array/HRIR measurements-based receiver */

} RECEIVER_TYPES;

/** Source object */
typedef struct _ims_src_obj
{
    float* sig;      /**< Source signal pointer */
    ims_pos_xyz pos; /**< Source position */
    int ID;          /**< Unique source ID */

} ims_src_obj;

/** Receiver object */
typedef struct _ims_rec_obj
{
    float** sigs;        /**< Receiver signal pointers (one per channel) */ 
    RECEIVER_TYPES type; /**< Receiver type (see #RECEIVER_TYPES enum) */
    int nChannels;       /**< Number of channels for receiver */
    ims_pos_xyz pos;     /**< Receiver position */
    int ID;              /**< Unique receiver ID */

} ims_rec_obj;

/** Echogram structure */
typedef struct _echogram_data
{
    /* The echogram data: */
    int numImageSources;     /**< Number of image sources in current echogram */
    int nChannels;           /**< Number of channels */
    float** value;           /**< Echogram magnitudes per image source and
                              *   channel; nChannels x numImageSources */
    float* time;             /**< Propagation time (in seconds) for each image
                              *   source; numImageSources x 1 */
    int** order;             /**< Reflection order for each image and dimension;
                              *   numImageSources x 3 */
    ims_pos_xyz* coords;     /**< Reflection coordinates (Cartesian);
                              *   numImageSources x 3 */
    int* sortedIdx;          /**< Indices that sort the echogram based on
                              *   propagation time in ascending order;
                              *   numImageSources x 1 */

    /* Optional helper variables for run-time speed-ups */
    int include_rt_vars;     /**< 0: the below vars are disabled, 1: enabled */
    float* tmp1;             /**< 1st temporary vector; numImageSources x 1 */
    float* tmp2;             /**< 2nd temporary vector; numImageSources x 1 */
    int* rIdx;               /**< Current circular buffer read indices;
                              *   numImageSources x 1 */
    int* rIdx_frac[IMS_LAGRANGE_ORDER]; /**< Current circular buffer read
                              *   indices for fractional buffers;
                              *   IMS_LAGRANGE_ORDER x numImageSources */
    float** h_frac;          /**< Current fractional delay coeffs;
                              *   (IMS_LAGRANGE_ORDER+1) x numImageSources x */
    float** cb_vals;         /**< Current circular buffer values (per channel &
                              *   image source); nChannels x numImageSources */
    float** contrib;         /**< Total contribution (i.e. cb_vals .* value);
                              *   nChannels x numImageSources */
    float* ones_dummy;       /**< Just a vector of ones, for the cblas_sdot
                              *   sum hack, and fmodf; numImageSources x 1 */

} echogram_data;

/**
 * Helper structure, comprising variables used when computing echograms and
 * rendering RIRs. The idea is that there should be one instance of this per
 * source/reciever combination.
 */
typedef struct _ims_core_workspace
{
    /* Locals */
    float room[3];        /**< Room dimensions, in meters */
    float d_max;          /**< Maximum distance, in meters */
    int N_max;            /**< Maximum reflection order */
    ims_pos_xyz src;      /**< Source position */
    ims_pos_xyz rec;      /**< Receiver position */
    int nBands;           /**< Number of bands */

    /* Internal */
    int Nx, Ny, Nz;
    int lengthVec, numImageSources;
    int* validIDs;
    int* s_ord;
    float* II, *JJ, *KK;
    int* iII, *iJJ, *iKK;
    float* s_x, *s_y, *s_z, *s_d, *s_t, *s_att;

    /* Echograms */
    int refreshEchogramFLAG;    /**< 1: Refresh needed, 0: refresh not needed */
    void* hEchogram;            /**< Pressure echogram (single-channel) */
    void* hEchogram_rec;        /**< Echogram with the receiver directivities
                                 *   applied (multi-channel) */
    voidPtr* hEchogram_abs;     /**< Echograms with the receiver directivities
                                 *   and also wall absorption applied (multi-
                                 *   channel); one echogram per band */
    voidPtr* hPrevEchogram_abs; /**< Previous echograms (hEchogram_abs), one
                                 *   per band, which can be used for cross-
                                 *   fading. */

    /* Room impulse responses (only used/allocated when a render function is
     * called) */
    int refreshRIRFLAG;
    int rir_len_samples;
    float rir_len_seconds;
    float*** rir_bands; /* nBands x nChannels x rir_len_samples */
 
} ims_core_workspace;

/**
 * Main structure for IMS. It comprises variables describing the room, and the
 * source and receiver objects within it. It also includes "core workspace"
 * handles for each source/receiver combination.
 */
typedef struct _ims_scene_data
{
    /* Locals */
    float room_dims[3];       /**< Room dimensions, in meters */
    float c_ms;               /**< Speed of sound, in ms^1 */
    float fs;                 /**< Sampling rate */
    int nBands;               /**< Number of frequency bands */
    float** abs_wall;         /**< Wall aborption coeffs per wall; nBands x 6 */

    /* Source and receiver positions */
    ims_src_obj srcs[IMS_MAX_NUM_SOURCES];   /**< Source positions */
    ims_rec_obj recs[IMS_MAX_NUM_RECEIVERS]; /**< Receiver positions */
    long nSources;            /**< Current number of sources */
    long nReceivers;          /**< Current number of receivers */

    /* Internal */
    voidPtr** hCoreWrkSpc;    /**< One per source/receiver combination */
    float* band_centerfreqs;  /**< Octave band CENTRE frequencies; nBands x 1 */
    float* band_cutofffreqs;  /**< Octave band CUTOFF frequencies;
                               *   (nBands-1) x 1 */
    float** H_filt;           /**< nBands x (#IMS_FIR_FILTERBANK_ORDER+1) */
    ims_rir** rirs;           /**< One per source/receiver combination */

    /* Circular buffers (only used/allocated when applyEchogramTD() function is
     * called for the first time) */
    unsigned long wIdx[IMS_EG_NUM_SLOTS][IMS_MAX_NUM_RECEIVERS][IMS_MAX_NUM_SOURCES];  /**< current write indices for circular buffers */
    float*** circ_buffer[IMS_EG_NUM_SLOTS];  /**< [IMS_EG_NUM_SLOTS] x (nChannels x nBands x #IMS_CIRC_BUFFER_LENGTH) */

    /* IIR filterbank (only used/allocated when applyEchogramTD() function is
     * called for the first time) */
    voidPtr* hFaFbank;          /**< One per source */
    float*** src_sigs_bands;    /**< nSources x nBands x nSamples */

    /* Temporary receiver frame used for cross-fading (only used/allocated when
     * applyEchogramTD() function is called for the first time) */
    float*** rec_sig_tmp[IMS_EG_NUM_SLOTS]; /**< [IMS_EG_NUM_SLOTS] x (nReceivers x nChannels x nSamples) */
    float* interpolator_fIn;  /**< framesize x 1 */
    float* interpolator_fOut; /**< framesize x 1 */
    float* tmp_frame;         /**< framesize x 1 */
    int applyCrossFadeFLAG[IMS_MAX_NUM_RECEIVERS][IMS_MAX_NUM_SOURCES];
    int framesize;            /**< Curent framesize in samples */

    /* Lagrange interpolator look-up table */
    float lookup_fractions[IMS_LAGRANGE_LOOKUP_TABLE_SIZE];
    float lookup_H_frac[IMS_LAGRANGE_ORDER+1][IMS_LAGRANGE_LOOKUP_TABLE_SIZE];

} ims_scene_data;


/* =========================== Internal Functions =========================== */

/**
 * Creates an instance of the core workspace
 *
 * The idea is that there is one core workspace instance per source/receiver
 * combination.
 *
 * @param[in] phWork (&) address of the workspace handle
 * @param[in] nBands Number of bands
 */
void ims_shoebox_coreWorkspaceCreate(void** phWork,
                                     int nBands);

/**
 * Destroys an instance of the core workspace
 *
 * @param[in] phWork  (&) address of the workspace handle
 */
void ims_shoebox_coreWorkspaceDestroy(void** phWork);

/**
 * Creates an instance of an echogram container
 *
 * @param[in] phEcho          (&) address of the echogram container
 * @param[in] include_rt_vars 1: optional run-time helper functions enabled,
 *                            0: disabled
 */
void ims_shoebox_echogramCreate(void** phEcho,
                                int include_rt_vars);

/**
 * Resizes an echogram container
 *
 * @note The container is only resized if the number of image sources or
 *       channels have changed.
 *
 * @param[in] hEcho           echogram container
 * @param[in] numImageSources New number of image sources
 * @param[in] nChannels       New number of channels
 */
void ims_shoebox_echogramResize(void* hEcho,
                                int numImageSources,
                                int nChannels);

/**
 * Copies echogram data from container 'X' into container 'Y' (also resizing 'Y'
 * as needed)
 *
 * @warning Helper variables are resized (if needed), but do not copy values!
 *
 * @param[in] hEchoX  echogram container X
 * @param[in] hEchoY  echogram container Y
 */
void ims_shoebox_echogramCopy(void* hEchoX,
                              void* hEchoY);

/**
 * Destroys an instance of an echogram container
 *
 * @param[in] phEcho (&) address of the echogram container
 */
void ims_shoebox_echogramDestroy(void** phEcho);

/**
 * Calculates an echogram of a rectangular space using the image source method,
 * for a specific source/reciever combination up to a maximum propagation time
 *
 * Note the coordinates of source/receiver are specified from the left ground
 * corner of the room:
 *
 * \verbatim
 *
 *                ^x
 *             __|__    _
 *             |  |  |   |
 *             |  |  |   |
 *          y<----.  |   | l
 *             |     |   |
 *             |     |   |
 *             o_____|   -
 *
 *             |-----|
 *                w
 *
 * \endverbatim
 *
 * @param[in] hWork   workspace handle
 * @param[in] room    Room dimensions, in meters
 * @param[in] src     Source position, in meters
 * @param[in] rec     Receiver position, in meters
 * @param[in] maxTime Maximum propagation time to compute the echogram, seconds
 * @param[in] c_ms    Speed of source, in meters per second
 */
void ims_shoebox_coreInitT(void* hWork,
                           float room[3],
                           ims_pos_xyz src,
                           ims_pos_xyz rec,
                           float maxTime,
                           float c_ms);

/**
 * Calculates an echogram of a rectangular space using the image source method,
 * for a specific source/reciever combination up to a maximum reflection order
 *
 * @param[in] hWork workspace handle
 * @param[in] room  Room dimensions, in meters
 * @param[in] src   Source position, in meters
 * @param[in] rec   Receiver position, in meters
 * @param[in] maxN  Maximum propagation time to compute the echogram, seconds
 * @param[in] c_ms  Speed of source, in meters per second
 */
void ims_shoebox_coreInitN(void* hWork,
                           float room[3],
                           ims_pos_xyz src,
                           ims_pos_xyz rec,
                           int maxN,
                           float c_ms);

/**
 * Imposes spherical harmonic directivies onto the echogram computed with
 * ims_shoebox_coreInit() for a specific source/reciever combination
 *
 * @note Call ims_shoebox_coreInit() before applying the directivities
 *
 * @param[in] hWork    workspace handle
 * @param[in] sh_order Spherical harmonic order
 */
void ims_shoebox_coreRecModuleSH(void* hWork,
                                 int sh_order);

/**
 * Applies boundary absoption per frequency band, onto the echogram computed
 * with ims_shoebox_coreRecModuleSH() for a specific source/reciever combination
 *
 * Absorption coefficients are given for each of the walls on the respective
 * planes [x+ y+ z+; x- y- z-].
 *
 * @note Call ims_shoebox_coreRecModuleX before applying the absoption
 *
 * @param[in] hWork    workspace handle
 * @param[in] abs_wall Absorption coefficients; nBands x 6
 */
void ims_shoebox_coreAbsorptionModule(void* hWork,
                                      float** abs_wall);

/**
 * Renders a room impulse response for a specific source/reciever combination
 *
 * @note Call ims_shoebox_coreAbsorptionModule() before rendering rir
 *
 * @param[in]  hWork               workspace handle
 * @param[in]  fractionalDelayFLAG 0: disabled, 1: use Lagrange interpolation
 * @param[in]  fs                  SampleRate, Hz
 * @param[in]  H_filt              filterbank; nBands x (filterOrder+1)
 * @param[out] rir                 Room impulse response
 */
void ims_shoebox_renderRIR(void* hWork,
                           int fractionalDelayFLAG,
                           float fs,
                           float** H_filt,
                           ims_rir* rir);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __REVERB_INTERNAL_H_INCLUDED__ */
