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

/*
 * Filename: sldoa_internal.h
 * --------------------------
 * A spatially-localised active-intensity based direction-of-arrival estimator
 * (SLDoA). VBAP gain patterns are imposed on the spherical harmonic signals,
 * such that the DoA can be estimated in a spatially-constrained region; thus
 * mitigating the effect of interferes and reflections arriving from other
 * directions. The DoA is estimated per sector for each frequency band.
 * The algorithms within sldoa were developed in collaboration with Symeon
 * Delikaris-Manias and Angelo Farina, and are explained in more detail in [1].
 * Dependencies:
 *     saf_utilities, afSTFTlib, saf_vbap, saf_sh
 * Author, date created:
 *     Leo McCormack, 18.10.2017
 *
 * [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki,
 *     V., “Real-time conversion of sensor array signals into spherical harmonic
 *     signals with applications to spatially localised sub-band sound-field
 *     analysis,” in Audio Engineering Society Convention 144, Audio Engineering
 *     Society, 2018.
 */

#ifndef __SLDOA_INTERNAL_H_INCLUDED__
#define __SLDOA_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "sldoa.h"
#include "sldoa_database.h"
#include "saf.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                               Internal Enums                               */
/* ========================================================================== */

/*
 * Enum: PROC_STATUS
 * -----------------
 * Current status of the processing loop.
 *
 * Options:
 *     PROC_STATUS_ONGOING     - Codec is processing input audio, and should not
 *                               be reinitialised at this time.
 *     PROC_STATUS_NOT_ONGOING - Codec is not processing input audio, and may
 *                               be reinitialised if needed.
 */
typedef enum _PROC_STATUS{
    PROC_STATUS_ONGOING = 0,
    PROC_STATUS_NOT_ONGOING
}PROC_STATUS;


/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define ORDER2NUMSIGS(L) ( (L+1) * (L+1) )
#define ELEV2INCL(E) ( (M_PI/2 - E) )
//#define ORDER2NUMSECTORS(L) ( 2*L )
#define ORDER2NUMSECTORS(L) ( L*L )

#define MAX_SH_ORDER ( 7 )
#define MAX_NUM_SH_SIGNALS ( (MAX_SH_ORDER + 1)*(MAX_SH_ORDER + 1)  )   /* (L+1)^2 */
#define HOP_SIZE ( 128 )                                    /* STFT hop size = nBands */
#define HYBRID_BANDS ( HOP_SIZE + 5 )                       /* hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )                /* Processing relies on fdHop = 16 */
#define MAX_NUM_SECTORS ( ORDER2NUMSECTORS(MAX_SH_ORDER) )      /* maximum number of sectors */
#define NUM_DISP_SLOTS ( 2 )                                /* needs to be at least 2. On slower systems that skip frames, consider more slots.  */
#ifndef M_PI
# define M_PI ( 3.14159265359f )
#endif
    
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */
    
typedef struct _sldoa
{
    /* TFT */
    float SHframeTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float_complex SHframeTF[HYBRID_BANDS][MAX_NUM_SH_SIGNALS][TIME_SLOTS];
    void* hSTFT;
    complexVector* STFTInputFrameTF;
    float** tempHopFrameTD;
    float freqVector[HYBRID_BANDS];
    float fs;
      
    /* ana configuration */
    CODEC_STATUS codecStatus;
    PROC_STATUS procStatus;
    float progressBar0_1;
    char* progressBarText;
    
    /* internal */
    float grid_Y[64][NUM_GRID_DIRS];
    float grid_Y_dipoles_norm[3][NUM_GRID_DIRS];
    float grid_dirs_deg[NUM_GRID_DIRS][2];
    float_complex* secCoeffs[MAX_SH_ORDER-1];
    float doa_rad[HYBRID_BANDS][MAX_NUM_SECTORS][2];
    float energy [HYBRID_BANDS][MAX_NUM_SECTORS];
    int nSectorsPerBand[HYBRID_BANDS];
    int new_masterOrder;
    
    /* display */
    float* azi_deg[NUM_DISP_SLOTS];
    float* elev_deg[NUM_DISP_SLOTS];
    float* colourScale[NUM_DISP_SLOTS];
    float* alphaScale[NUM_DISP_SLOTS]; 
    int current_disp_idx;
    
    /* User parameters */
    int masterOrder;
    int analysisOrderPerBand[HYBRID_BANDS];
    float maxFreq;
    float minFreq;
    float avg_ms;
    CH_ORDER chOrdering;
    NORM_TYPES norm;

} sldoa_data;
     

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/*
 * Function: sldoa_setCodecStatus
 * ------------------------------------
 * Sets codec status.
 *
 * Input Arguments:
 *     hSld      - sldoa handle
 *     newStatus - codec status (see 'CODEC_STATUS' enum)
 */
void sldoa_setCodecStatus(void* const hSld, CODEC_STATUS newStatus);

/*
 * sldoa_initAna
 * -------------
 * Intialises the codec variables, based on current global/user parameters.
 * The formulae for calculating the sector coefficients can be found in [1].
 *
 * Input Arguments:
 *     hSld - sldoa handle
 *
 * [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki,
 *     V., “Real-time conversion of sensor array signals into spherical harmonic
 *     signals with applications to spatially localised sub-band sound-field
 *     analysis,” in Audio Engineering Society Convention 144, Audio Engineering
 *     Society, 2018.
 */
void sldoa_initAna(void* const hSld);
    
/*
 * sldoa_initTFT
 * -------------
 * Initialise the filterbank used by sldoa.
 * Note: Call this function before sldoa_initAna
 *
 * Input Arguments:
 *     hSld - sldoa handle
 */
void sldoa_initTFT(void* const hSld);
  
/*
 * sldoa_estimateDoA
 * ----------------
 * Estimates the DoA using the active intensity vectors derived from spatially
 * localised sectors, as in [1].
 * Note: if anaOrder is 1, then the algorithm reverts to the standard active-
 * intensity based DoA estimation.
 *
 * Input Arguments:
 *     SHframeTF - input SH frame
 *     anaOrder  - analysis order (1:AI, 2+: SLAI)
 *     secCoeffs - sector coefficients for this order
 * Output Arguments:
 *     doa       - resulting DoA estimates per timeslot and sector
 *     energy    - resulting sector energies per time slot
 *
 * [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and Pulkki,
 *     V., “Real-time conversion of sensor array signals into spherical harmonic
 *     signals with applications to spatially localised sub-band sound-field
 *     analysis,” in Audio Engineering Society Convention 144, Audio Engineering
 *     Society, 2018.
 */
void sldoa_estimateDoA(float_complex SHframeTF[MAX_NUM_SH_SIGNALS][TIME_SLOTS],
                       int anaOrder,
                       float_complex* secCoeffs,
                       float doa[MAX_NUM_SECTORS][TIME_SLOTS][2],
                       float energy[MAX_NUM_SECTORS][TIME_SLOTS]);
    

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SLDOA_INTERNAL_H_INCLUDED__ */
