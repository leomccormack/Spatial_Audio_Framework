/*
 * Copyright 2019 Leo McCormack
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
 * @file dirass_internal.h
 * @brief A sound-field visualiser based on the directional re-assignment of
 *        beamformer energy based on local DoA estimates [1,2]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 * @see [2] McCormack, L., Politis, A., and Pulkki, V. (2019). "Sharpening of
 *          angular spectra based on a directional re-assignment approach for
 *          ambisonic sound-field visualisation". IEEE International Conference
 *          on Acoustics, Speech and Signal Processing (ICASSP).
 *
 * @author Leo McCormack
 * @date 21.02.2019
 * @license ISC
 */

#ifndef __DIRASS_INTERNAL_H_INCLUDED__
#define __DIRASS_INTERNAL_H_INCLUDED__

#include "dirass.h"        /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#if !defined(DIRASS_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define DIRASS_FRAME_SIZE ( FRAME_SIZE )   /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define DIRASS_FRAME_SIZE ( 1024 )         /**< Framesize, in time-domain samples */
# endif
#endif
#define MAX_DISPLAY_SH_ORDER ( 20 )          /**< Maximum display/upscaling SH order */
#define MAX_NUM_INPUT_SH_SIGNALS ( (MAX_SH_ORDER+1)*(MAX_SH_ORDER+1) )   /**< Maximum number of SH signals for the input */
#define MAX_NUM_DISPLAY_SH_SIGNALS ( (MAX_DISPLAY_SH_ORDER+1)*(MAX_DISPLAY_SH_ORDER+1) )  /**< Maximum number of SH signals for the display/upscaling SH output */
#define NUM_DISP_SLOTS ( 2 )                 /**< Number of display slots */

/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Contains variables for scanning grids, and sector beamforming
 */
typedef struct _dirass_codecPars
{
    /* scanning grid and intepolation table */
    float* grid_dirs_deg;     /**< scanning grid directions; FLAT: grid_nDirs x 2 */
    int grid_nDirs;           /**< number of grid directions */
    float* interp_dirs_deg;   /**< interpolation directions, in degrees; FLAT: interp_nDirs x 2 */
    float* interp_dirs_rad;   /**< interpolation directions, in radians; FLAT: interp_nDirs x 2 */
    float* interp_table;      /**< interpolation table (spherical->rectangular grid); FLAT: interp_nDirs x grid_nDirs */
    int interp_nDirs;         /**< number of interpolation directions */
    int interp_nTri;          /**< number of triangles in the spherical scanning grid mesh */
    float* ss;                /**< beamformer sector signals; FLAT: grid_nDirs x DIRASS_FRAME_SIZE */
    float* ssxyz;             /**< beamformer velocity signals; FLAT: 3 x DIRASS_FRAME_SIZE */
    int* est_dirs_idx;        /**< DoA indices, into the interpolation directions; grid_nDirs x 1 */
    float* prev_intensity;    /**< previous intensity vectors (for averaging); FLAT: grid_nDirs x 3 */
    float* prev_energy;       /**< previous energy (for averaging); FLAT: grid_nDirs x 1 */
    
    /* sector beamforming and upscaling */
    float* Cxyz;              /**< beamforming weights for velocity patterns; FLAT: nDirs x (order+1)^2 x 3 */
    float* Cw;                /**< beamforming weights; FLAT: nDirs x (order)^2 */
    float* Uw;                /**< beamforming weights; FLAT: nDirs x (upscaleOrder+1)^2 */
    float* Y_up;              /**< real SH weights for upscaling; FLAT: (upscaleOrder+1)^2 x grid_nDirs */
    float* est_dirs;          /**< estimated DoA per grid direction; grid_nDirs x 2 */
    
    /* regular beamforming */
    float* w;                 /**< beamforming weights; FLAT: nDirs x (order+1)^2 */
     
}dirass_codecPars;
    
/**
 * Main structure for dirass. Contains variables for audio buffers, filtering,
 * internal variables, flags, user parameters
 */
typedef struct _dirass
{
    /* FIFO buffers */
    int FIFO_idx;                           /**< FIFO buffer index */
    float inFIFO[MAX_NUM_INPUT_SH_SIGNALS][DIRASS_FRAME_SIZE]; /**< FIFO buffer */
    
    /* Buffers */
    float SHframeTD[MAX_NUM_INPUT_SH_SIGNALS][DIRASS_FRAME_SIZE];       /**< Input SH signals */
    float SHframe_upTD[MAX_NUM_DISPLAY_SH_SIGNALS][DIRASS_FRAME_SIZE];  /**< Upscaled SH signals */
    float fs;                               /**< host sampling rate */
    
    /* internal */ 
    int dispWidth;                          /**< number of interpolation points on the horizontal */
    float Wz12_hpf[MAX_NUM_INPUT_SH_SIGNALS][2]; /**< delayed elements used in the HPF */
    float Wz12_lpf[MAX_NUM_INPUT_SH_SIGNALS][2]; /**< delayed elements used in the LPF */
    int new_inputOrder;                     /**< New input/analysis order */
    int new_upscaleOrder;                   /**< New target upscale order */
    
    /* ana configuration */
    CODEC_STATUS codecStatus;               /**< see #CODEC_STATUS */
    PROC_STATUS procStatus;                 /**< see #PROC_STATUS */
    float progressBar0_1;                   /**< Current (re)initialisation progress, between [0..1] */
    char* progressBarText;                  /**< Current (re)initialisation step, string */
    dirass_codecPars* pars;                 /**< codec parameters */
    
    /* display */
    float* pmap;                            /**< grid_nDirs x 1 */
    float* pmap_grid[NUM_DISP_SLOTS];       /**< dirass interpolated to grid; interp_nDirs x 1 */
    int dispSlotIdx;                        /**< current display slot index */
    float pmap_grid_minVal;                 /**< minimum value in pmap */
    float pmap_grid_maxVal;                 /**< maximum value in pmap */
    int recalcPmap;                         /**< set this to 1 to generate a new image */
    int pmapReady;                          /**< 0: image generation not started yet, 1: image is ready for plotting*/
    
    /* User parameters */
    int inputOrder;                         /**< Current input/analysis order */
    STATIC_BEAM_TYPES beamType;             /**< beamformer type mode */
    DIRASS_REASS_MODES DirAssMode;          /**< see #DIRASS_REASS_MODES enum */
    int upscaleOrder;                       /**< Current target upscale order */
    DIRASS_GRID_OPTIONS gridOption;         /**< grid option */
    float pmapAvgCoeff;                     /**< averaging coefficient for the intensity vector per grid direction */
    float minFreq_hz;                       /**< minimum frequency to include in pmap generation, Hz */
    float maxFreq_hz;                       /**< maximum frequency to include in pmap generation, Hz */
    CH_ORDER chOrdering;                    /**< Ambisonic channel order convention (see #CH_ORDER) */
    NORM_TYPES norm;                        /**< Ambisonic normalisation convention (see #NORM_TYPES) */
    HFOV_OPTIONS HFOVoption;                /**< horizontal field-of-view option */
    ASPECT_RATIO_OPTIONS aspectRatioOption; /**< aspect ratio option */
    
} dirass_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/** Sets codec status (see #CODEC_STATUS enum) */
void dirass_setCodecStatus(void* const hDir, CODEC_STATUS newStatus);

/** Intialises the codec variables, based on current global/user parameters */
void dirass_initAna(void* const hDir);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __DIRASS_INTERNAL_H_INCLUDED__ */
