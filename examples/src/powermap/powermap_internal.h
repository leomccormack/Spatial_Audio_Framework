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
 * @file powermap_internal.h
 * @brief A sound-field visualiser, which utilises spherical harmonic signals as
 *        input; note this code is a remnant from the work conducted in [1]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S. and Pulkki, V., 2017. Parametric
 *          acoustic camera for real-time sound capture, analysis and tracking.
 *          In Proceedings of the 20th International Conference on Digital Audio
 *          Effects (DAFx-17) (pp. 412-419)
 *
 * @author Leo McCormack
 * @date 26.04.2016
 * @license ISC
 */

#ifndef __POWERMAP_INTERNAL_H_INCLUDED__
#define __POWERMAP_INTERNAL_H_INCLUDED__

#include "powermap.h"      /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#if !defined(POWERMAP_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define POWERMAP_FRAME_SIZE ( FRAME_SIZE )          /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define POWERMAP_FRAME_SIZE ( 1024 )                /**< Framesize, in time-domain samples */
# endif
#endif
#define HOP_SIZE ( 128 )                              /**< STFT hop size */
#define HYBRID_BANDS ( HOP_SIZE + 5 )                 /**< Number of frequency bands */
#define TIME_SLOTS ( POWERMAP_FRAME_SIZE / HOP_SIZE ) /**< Number of STFT timeslots */
#define NUM_DISP_SLOTS ( 2 )                          /**< Number of display slots */
#define MAX_COV_AVG_COEFF ( 0.45f )                   /**< Maximum supported covariance averaging coefficient  */

/* Checks: */
#if (POWERMAP_FRAME_SIZE % HOP_SIZE != 0)
# error "POWERMAP_FRAME_SIZE must be an integer multiple of HOP_SIZE"
#endif
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/** Contains variables for scanning grids, and beamforming */
typedef struct _powermap_codecPars
{
    float* grid_dirs_deg;   /**< Spherical scanning grid directions, in degrees; FLAT: grid_nDirs x 2 */
    int grid_nDirs;         /**< Number of scanning directions */
    float* interp_dirs_deg; /**< 2D rectangular window interpolation directions, in degrees; FLAT: interp_nDirs x 2 */
    float* interp_table;    /**< Spherical->2D interpolation table; FLAT: interp_nDirs x grid_nDirs */
    int interp_nDirs;       /**< Number of interpolation directions */
    int interp_nTri;        /**< Number of triangles in the spherical triangulared grid */
    float* Y_grid[MAX_SH_ORDER];                 /**< real SH basis (real datatype); MAX_NUM_SH_SIGNALS x grid_nDirs */
    float_complex* Y_grid_cmplx[MAX_SH_ORDER];   /**< real SH basis (complex datatype); MAX_NUM_SH_SIGNALS x grid_nDirs */
    
}powermap_codecPars;
    
/**
 * Main structure for powermap. Contains variables for audio buffers, internal
 * variables, flags, user parameters
 */
typedef struct _powermap
{
    /* FIFO buffers */
    int FIFO_idx;                   /**< FIFO buffer index */
    float inFIFO[MAX_NUM_SH_SIGNALS][POWERMAP_FRAME_SIZE]; /**< Input FIFO buffer */

    /* TFT */
    float** SHframeTD;              /**< time-domain SH input frame; #MAX_NUM_SH_SIGNALS x #POWERMAP_FRAME_SIZE */
    float_complex*** SHframeTF;     /**< time-frequency domain SH input frame; #HYBRID_BANDS x #MAX_NUM_SH_SIGNALS x #TIME_SLOTS */
    void* hSTFT;                    /**< afSTFT handle */
    float freqVector[HYBRID_BANDS]; /**< Frequency vector (filterbank centre frequencies) */
    float fs;                       /**< Host sample rate, in Hz*/
    
    /* internal */
    float_complex Cx[HYBRID_BANDS][MAX_NUM_SH_SIGNALS*MAX_NUM_SH_SIGNALS];     /**< covariance matrices per band */
    int new_masterOrder;            /**< New maximum/master SH analysis order (current value will be replaced by this after next re-init) */
    int dispWidth;                  /**< Number of pixels on the horizontal in the 2D interpolated powermap image */
    
    /* ana configuration */
    CODEC_STATUS codecStatus;       /**< see #CODEC_STATUS */
    PROC_STATUS procStatus;         /**< see #PROC_STATUS */
    float progressBar0_1;           /**< Current (re)initialisation progress, between [0..1] */
    char* progressBarText;          /**< Current (re)initialisation step, string */
    powermap_codecPars* pars;       /**< codec parameters */
    
    /* display */
    float* pmap;                    /**< grid_nDirs x 1 */
    float* prev_pmap;               /**< grid_nDirs x 1 */
    float* pmap_grid[NUM_DISP_SLOTS]; /**< powermap interpolated to grid; interp_nDirs x 1 */
    int dispSlotIdx;                /**< Current display slot */
    float pmap_grid_minVal;         /**< Current minimum value in pmap (used to normalise [0..1]) */
    float pmap_grid_maxVal;         /**< Current maximum value in pmap (used to normalise [0..1]) */
    int recalcPmap;                 /**< set this to 1 to generate a new powermap */
    int pmapReady;                  /**< 0: powermap not started yet, 1: powermap is ready for plotting*/
    
    /* User parameters */
    int masterOrder;                /**< Current maximum/master SH analysis order */
    int analysisOrderPerBand[HYBRID_BANDS]; /**< SH analysis order per frequency band */
    float pmapEQ[HYBRID_BANDS];     /**< Equalisation/weights per band */
    HFOV_OPTIONS HFOVoption;        /**< see #HFOV_OPTIONS */
    ASPECT_RATIO_OPTIONS aspectRatioOption; /**< see #ASPECT_RATIO_OPTIONS */
    float covAvgCoeff;              /**< Covariance matrix averaging coefficient, [0..1] */
    float pmapAvgCoeff;             /**< Powermap averaging coefficient, [0..1] */
    int nSources;                   /**< Current number of sources (used for MUSIC) */
    POWERMAP_MODES pmap_mode;       /**< see #POWERMAP_MODES*/
    CH_ORDER chOrdering;            /**< Ambisonic channel order convention (see #CH_ORDER) */
    NORM_TYPES norm;                /**< Ambisonic normalisation convention (see #NORM_TYPES) */
    
} powermap_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/** Sets codec status (see #CODEC_STATUS enum) */
void powermap_setCodecStatus(void* const hPm, CODEC_STATUS newStatus);

/** Intialises the codec variables, based on current global/user parameters */
void powermap_initAna(void* const hPm);

/**
 * Initialise the filterbank used by powermap.
 *
 * @note Call this function before powermap_initAna()
 */
void powermap_initTFT(void* const hPm);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __POWERMAP_INTERNAL_H_INCLUDED__ */
