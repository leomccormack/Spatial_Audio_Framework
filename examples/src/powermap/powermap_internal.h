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
 *        input; note this code is a remnant from the work described in [1]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S. and Pulkki, V., 2017. Parametric
 *          acoustic camera for real-time sound capture, analysis and tracking.
 *          In Proceedings of the 20th International Conference on Digital Audio
 *          Effects (DAFx-17) (pp. 412-419)
 *
 * @author Leo McCormack
 * @date 26.04.2016
 */

#ifndef __POWERMAP_INTERNAL_H_INCLUDED__
#define __POWERMAP_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "powermap.h"

#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define FRAME_SIZE ( 1024 )
#define MAX_SH_ORDER ( 7 )
#define HOP_SIZE ( 128 )                   /* STFT hop size = nBands */
#define HYBRID_BANDS ( HOP_SIZE + 5 )      /* hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE ) /* Processing relies on fdHop = 16 */
#define NUM_DISP_SLOTS ( 2 )
#define MAX_COV_AVG_COEFF ( 0.45f )    /*  */
#ifndef M_PI
# define M_PI ( 3.14159265359f )
#endif
    
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Contains variables for scanning grids, and beamforming
 */
typedef struct _powermap_codecPars
{
    float* grid_dirs_deg; /* grid_nDirs x 2 */
    int grid_nDirs;
    float* interp_dirs_deg;
    float* interp_table;    /* interp_nDirs x grid_nDirs */
    int interp_nDirs;
    int interp_nTri;
    
    float* Y_grid[MAX_SH_ORDER];                 /* MAX_NUM_SH_SIGNALS x grid_nDirs */
    float_complex* Y_grid_cmplx[MAX_SH_ORDER];   /* MAX_NUM_SH_SIGNALS x grid_nDirs */
    
}powermap_codecPars;
    
/**
 * Main structure for powermap. Contains variables for audio buffers, internal
 * variables, flags, user parameters
 */
typedef struct _powermap
{
    /* FIFO buffers */
    int FIFO_idx;
    float inFIFO[MAX_NUM_SH_SIGNALS][FRAME_SIZE]; 

    /* TFT */
    float SHframeTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float_complex SHframeTF[HYBRID_BANDS][MAX_NUM_SH_SIGNALS][TIME_SLOTS];        
    void* hSTFT;
    complexVector* STFTInputFrameTF;
    float** tempHopFrameTD;
    float freqVector[HYBRID_BANDS];
    float fs;
    
    /* internal */
    float_complex Cx[HYBRID_BANDS][MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS];     /* cov matrices */
    int new_masterOrder;
    int dispWidth;
    
    /* ana configuration */
    CODEC_STATUS codecStatus;
    PROC_STATUS procStatus;
    float progressBar0_1;
    char* progressBarText;
    powermap_codecPars* pars;                                          /* codec parameters */
    
    /* display */
    float* pmap;                           /* grid_nDirs x 1 */
    float* prev_pmap;                      /* grid_nDirs x 1 */
    float* pmap_grid[NUM_DISP_SLOTS];      /* powermap interpolated to grid; interp_nDirs x 1 */
    int dispSlotIdx;
    float pmap_grid_minVal;
    float pmap_grid_maxVal;
    int recalcPmap;   /* set this to 1 to generate a new powermap */
    int pmapReady;    /* 0: powermap not started yet, 1: powermap is ready for plotting*/
    
    /* User parameters */
    int masterOrder;
    int analysisOrderPerBand[HYBRID_BANDS];
    float pmapEQ[HYBRID_BANDS]; 
    HFOV_OPTIONS HFOVoption;
    ASPECT_RATIO_OPTIONS aspectRatioOption;
    float covAvgCoeff;
    float pmapAvgCoeff;
    int nSources;
    POWERMAP_MODES pmap_mode;
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    
} powermap_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Sets codec status (see #_CODEC_STATUS enum)
 */
void powermap_setCodecStatus(void* const hPm, CODEC_STATUS newStatus);

/**
 * Intialises the codec variables, based on current global/user parameters
 */
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
