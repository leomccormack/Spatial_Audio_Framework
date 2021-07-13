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
 * @file sldoa_internal.h
 * @brief A spatially-localised active-intensity (SLAI) based direction-of-
 *        arrival estimator (SLDoA)
 *
 * VBAP gain patterns are imposed on the spherical harmonic signals, such that
 * the DoA can be estimated in a spatially-constrained region; thus mitigating
 * the effect of interferes and reflections arriving from other directions.
 * The DoA is estimated per sector for each frequency band.
 *
 * The algorithms within sldoa were developed in collaboration with Symeon
 * Delikaris-Manias and Angelo Farina, and are explained in more detail in [1,2]
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 * @see [2] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., "Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis," in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 *
 * @author Leo McCormack
 * @date 18.10.2017
 * @license ISC
 */

#ifndef __SLDOA_INTERNAL_H_INCLUDED__
#define __SLDOA_INTERNAL_H_INCLUDED__

#include "sldoa.h"         /* Include header for this example */
#include "sldoa_database.h"/* Database header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#if !defined(SLDOA_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define SLDOA_FRAME_SIZE ( FRAME_SIZE )  /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define SLDOA_FRAME_SIZE ( 512 )         /**< Framesize, in time-domain samples */
# endif
#endif
#define ORDER2NUMSECTORS(L) ( L*L )        /**< Macro to convert SH order to number of sectors */
#define HOP_SIZE ( 128 )                   /**< STFT hop size */
#define HYBRID_BANDS ( HOP_SIZE + 5 )      /**< hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( SLDOA_FRAME_SIZE / HOP_SIZE )          /**< Processing relies on fdHop = 16 */
#define MAX_NUM_SECTORS ( ORDER2NUMSECTORS(MAX_SH_ORDER) )  /**< maximum number of sectors */
#define NUM_DISP_SLOTS ( 2 )               /**< Number of display slots; needs to be at least 2. On slower systems that skip frames, consider adding more slots.  */

/* Checks: */
#if (SLDOA_FRAME_SIZE % HOP_SIZE != 0)
# error "SLDOA_FRAME_SIZE must be an integer multiple of HOP_SIZE"
#endif
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */
   
/** Main struct for sldoa */
typedef struct _sldoa
{
    /* FIFO buffers */
    int FIFO_idx;                    /**< FIFO buffer index */
    float inFIFO[MAX_NUM_SH_SIGNALS][SLDOA_FRAME_SIZE]; /**< FIFO buffer */

    /* TFT */
    float** SHframeTD;              /**< time-domain SH input frame; #MAX_NUM_SH_SIGNALS x #SLDOA_FRAME_SIZE */
    float_complex*** SHframeTF;     /**< time-frequency domain SH input frame; #HYBRID_BANDS x #MAX_NUM_SH_SIGNALS x #TIME_SLOTS */
    void* hSTFT;                    /**< afSTFT handle */
    float freqVector[HYBRID_BANDS]; /**< Frequency vector (filterbank centre frequencies) */
    float fs;                       /**< Host sampling rate, in Hz */
      
    /* ana configuration */
    CODEC_STATUS codecStatus;       /**< see #CODEC_STATUS */
    PROC_STATUS procStatus;         /**< see #PROC_STATUS */
    float progressBar0_1;           /**< Current (re)initialisation progress, between [0..1] */
    char* progressBarText;          /**< Current (re)initialisation step, string */
    
    /* internal */
    float grid_Y[64][NUM_GRID_DIRS];                 /**< SH basis */
    float grid_Y_dipoles_norm[3][NUM_GRID_DIRS];     /**< SH basis */
    float grid_dirs_deg[NUM_GRID_DIRS][2];           /**< Grid directions, in degrees */
    float_complex* secCoeffs[MAX_SH_ORDER-1];        /**< Sector beamforming weights/coefficients */
    float doa_rad[HYBRID_BANDS][MAX_NUM_SECTORS][2]; /**< Current DoA estimates per band and sector, in radians */
    float energy [HYBRID_BANDS][MAX_NUM_SECTORS];    /**< Current Sector energies */
    int nSectorsPerBand[HYBRID_BANDS];               /**< Number of sectors per band */
    int new_masterOrder;                             /**< New master/maximum analysis order (current value will be replaced by this after next re-init) */
    
    /* display */
    float* azi_deg[NUM_DISP_SLOTS];      /**< DoA azimuths, in degrees */
    float* elev_deg[NUM_DISP_SLOTS];     /**< DoA elevations, in degrees */
    float* colourScale[NUM_DISP_SLOTS];  /**< Values dictating each DoA marker colour */
    float* alphaScale[NUM_DISP_SLOTS];   /**< Values dictating each DoA marker transparency */
    int current_disp_idx;                /**< Current display slot */
    
    /* User parameters */
    int masterOrder;                     /**< Current master/maximum analysis order */
    int analysisOrderPerBand[HYBRID_BANDS]; /**< Analysis order MIN(anaPerBand, masterOrder) for each frequency band */
    float maxFreq;                       /**< Maximum display frequency, in Hz */
    float minFreq;                       /**< Minimum display frequency, in Hz */
    float avg_ms;                        /**< Temporal averaging, in ms */
    CH_ORDER chOrdering;                 /**< Ambisonic channel order convention (see #CH_ORDER) */
    NORM_TYPES norm;                     /**< Ambisonic normalisation convention (see #NORM_TYPES) */

} sldoa_data;
     

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/** Sets codec status (see #CODEC_STATUS enum) */
void sldoa_setCodecStatus(void* const hSld, CODEC_STATUS newStatus);

/**
 * Intialises the codec variables, based on current global/user parameters.
 *
 * The formulae for calculating the sector coefficients can be found in [1,2].
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., “Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis,” in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 * @see [2] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 */
void sldoa_initAna(void* const hSld);
    
/**
 * Initialise the filterbank used by sldoa.
 *
 * @note Call this function before sldoa_initAna()
 *
 * Input Arguments:
 *     hSld - sldoa handle
 */
void sldoa_initTFT(void* const hSld);
  
/**
 * Estimates the DoA using the active intensity vectors derived from spatially
 * localised sectors, as in [1,2].
 *
 * @note If anaOrder is 1, then the algorithm reverts to the standard active-
 *       intensity based DoA estimation.
 *
 * @param[in]  SHframeTF Input SH frame; MAX_NUM_SH_SIGNALS x TIME_SLOTS
 * @param[in]  anaOrder  Analysis order (1:AI, 2+: SLAI)
 * @param[in]  secCoeffs Sector coefficients for this order
 * @param[out] doa       Resulting DoA estimates per timeslot and sector
 * @param[out] energy    Resulting sector energies per time slot
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., “Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis,” in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 * @see [2] McCormack, L., Delikaris-Manias, S., Politis, A., Pavlidi, D.,
 *          Farina, A., Pinardi, D. and Pulkki, V., 2019. Applications of
 *          Spatially Localized Active-Intensity Vectors for Sound-Field
 *          Visualization. Journal of the Audio Engineering Society, 67(11),
 *          pp.840-854.
 */
void sldoa_estimateDoA(float_complex** SHframeTF,
                       int anaOrder,
                       float_complex* secCoeffs,
                       float doa[MAX_NUM_SECTORS][TIME_SLOTS][2],
                       float energy[MAX_NUM_SECTORS][TIME_SLOTS]);
    

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SLDOA_INTERNAL_H_INCLUDED__ */
