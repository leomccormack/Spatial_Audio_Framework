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
 * @file array2sh_internal.h
 * @brief Spatially encodes spherical or cylindrical sensor array signals into
 *        spherical harmonic signals utilising theoretical encoding filters.
 *
 * The algorithms within array2sh were pieced together and developed in
 * collaboration with Symeon Delikaris-Manias and Angelo Farina.
 * A detailed explanation of the algorithms within array2sh can be found in [1].
 * Also included, is a diffuse-field equalisation option for frequencies past
 * aliasing, developed in collaboration with Archontis Politis, 8.02.2019
 *
 * @note Since the algorithms are based on theory, only array designs where
 *       there are analytical solutions available are supported. i.e. only
 *       spherical or cylindrical arrays, which have phase-matched sensors.
 *
 * @see [1] McCormack, L., Delikaris-Manias, S., Farina, A., Pinardi, D., and
 *          Pulkki, V., "Real-time conversion of sensor array signals into
 *          spherical harmonic signals with applications to spatially localised
 *          sub-band sound-field analysis," in Audio Engineering Society
 *          Convention 144, Audio Engineering Society, 2018.
 *
 * @author Leo McCormack
 * @date 13.09.2017
 */

#ifndef __ARRAY2SH_INTERNAL_H_INCLUDED__
#define __ARRAY2SH_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "array2sh.h" 
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#ifndef FRAME_SIZE
# define FRAME_SIZE ( 128 ) 
#endif
#define HOP_SIZE ( 128 )                       /* STFT hop size = nBands */
#define HYBRID_BANDS ( HOP_SIZE + 5 )          /* hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )   /* 4/8/16 */
#define MAX_NUM_SENSORS ( ARRAY2SH_MAX_NUM_SENSORS ) /* Maximum permitted number of channels for the VST standard */
#define MAX_EVAL_FREQ_HZ ( 20e3f )             /* Up to which frequency should the evaluation be accurate */
#define MAX_NUM_SENSORS_IN_PRESET ( MAX_NUM_SENSORS )


/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Contains variables for describing the microphone/hydrophone array
 */
typedef struct _array2sh_arrayPars {
    int Q, newQ;                    /* number of sensors */
    float r;                        /* radius of sensors */
    float R;                        /* radius of scatterer (only for rigid arrays) */
    ARRAY2SH_ARRAY_TYPES arrayType;          /* array type, spherical/cylindrical */
    ARRAY2SH_WEIGHT_TYPES weightType;        /* open/rigid etc */
    float sensorCoords_rad[MAX_NUM_SENSORS][2];
    float sensorCoords_deg[MAX_NUM_SENSORS][2];
        
}array2sh_arrayPars;

/**
 * Main structure for array2sh. Contains variables for audio buffers, afSTFT,
 * encoding matrices, internal variables, flags, user parameters
 */
typedef struct _array2sh
{
    /* FIFO buffers */
    int FIFO_idx;
    float inFIFO[MAX_NUM_SENSORS][FRAME_SIZE];
    float outFIFO[MAX_NUM_SH_SIGNALS][FRAME_SIZE];

    /* audio buffers */
    float inputFrameTD[MAX_NUM_SENSORS][FRAME_SIZE];
    float SHframeTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float_complex inputframeTF[HYBRID_BANDS][MAX_NUM_SENSORS][TIME_SLOTS];
    float_complex SHframeTF[HYBRID_BANDS][MAX_NUM_SH_SIGNALS][TIME_SLOTS];
    complexVector* STFTInputFrameTF;
    complexVector* STFTOutputFrameTF;
    float** tempHopFrameTD_in;
    float** tempHopFrameTD_out;
    
    /* intermediates */
    double_complex bN_modal[HYBRID_BANDS][MAX_SH_ORDER + 1];
    double_complex* bN;
    double_complex bN_inv[HYBRID_BANDS][MAX_SH_ORDER + 1];
    double_complex bN_inv_R[HYBRID_BANDS][MAX_NUM_SH_SIGNALS]; 
    float_complex W[HYBRID_BANDS][MAX_NUM_SH_SIGNALS][MAX_NUM_SENSORS];
    float_complex W_diffEQ[HYBRID_BANDS][MAX_NUM_SH_SIGNALS][MAX_NUM_SENSORS];
    
    /* for displaying the bNs */
    float** bN_modal_dB;            /* modal responses / no regulaisation; HYBRID_BANDS x (MAX_SH_ORDER +1)  */
    float** bN_inv_dB;              /* modal responses / with regularisation; HYBRID_BANDS x (MAX_SH_ORDER +1)  */
    float* cSH;                     /* spatial correlation; HYBRID_BANDS x 1 */
    float* lSH;                     /* level difference; HYBRID_BANDS x 1 */ 
    
    /* time-frequency transform and array details */
    float freqVector[HYBRID_BANDS]; /* frequency vector */
    void* hSTFT;                    /* filterbank handle */
    void* arraySpecs;               /* array configuration */
    
    /* internal parameters */
    ARRAY2SH_EVAL_STATUS evalStatus;
    float progressBar0_1;
    char* progressBarText;
    int fs;                         /* sampling rate, hz */
    int new_order;                  /* new encoding order */
    
    /* flags */
    PROC_STATUS procStatus;
    int reinitSHTmatrixFLAG;        /* 0: do not reinit; 1: reinit; */
    int evalRequestedFLAG;          /* 0: do not reinit; 1: reinit; */
    
    /* additional user parameters that are not included in the array presets */
    int order;                      /* current encoding order */
    ARRAY2SH_MICROPHONE_ARRAY_PRESETS preset; /* currently selected MIC preset */
    ARRAY2SH_FILTER_TYPES filterType;  /* encoding filter approach */
    float regPar;                   /* regularisation upper gain limit, dB; */
    CH_ORDER chOrdering;   /* ACN */
    NORM_TYPES norm;       /* N3D/SN3D */
    float c;                        /* speed of sound, m/s */
    float gain_dB;                  /* post gain, dB */ 
    int enableDiffEQpastAliasing;   /* 0: disabled, 1: enabled */
    
} array2sh_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Initialise the filterbank used by array2sh.
 *
 * @note Call this function before array2sh_calculate_sht_matrix()
 */
void array2sh_initTFT(void* const hA2sh);

/**
 * Computes the spherical harmonic transform (SHT) matrix, to spatially encode
 * input microphone/hydrophone signals into spherical harmonic signals.
 */
void array2sh_calculate_sht_matrix(void* const hA2sh);

/**
 * Applies diffuse-field equalisation at frequencies above the spatial aliasing
 * limit.
 */
void array2sh_apply_diff_EQ(void* const hA2sh);

/**
 * Computes the magnitude responses of the equalisation filters; the
 * absolute values of the regularised inversed modal coefficients.
 */
void array2sh_calculate_mag_curves(void* const hA2sh);

/**
 * Evaluates the spherical harmonic transform performance with the currently
 * configured microphone/hydrophone array.
 *
 * @note This is based on an analytical model of the array, so may differ in
 *       practice (although, it is usually pretty close, and saves from having
 *       to measure the array)
 */
void array2sh_evaluateSHTfilters(void* hA2sh);

/**
 * Creates an instance of a struct, which contains the array configuration
 * data.
 *
 * @param[in] hPars (&) array configuration handle
 */
void array2sh_createArray(void ** const hPars);

/**
 * Destroys an instance of a struct, which contains the array configuration
 * data.
 *
 * @param[in] hPars (&) array configuration handle
 */
void array2sh_destroyArray(void ** const hPars);

/**
 * Intialises an instance of a struct based on a preset, which contains the
 * array configuration data
 *
 * @param[in] hPars         (&) array configuration handle
 * @param[in] preset        Array preset (see
 *                          #_ARRAY2SH_MICROPHONE_ARRAY_PRESETS enum)
 * @param[in] arrayOrder    (&) maximum encoding order of the current preset
 * @param[in] firstInitFLAG '1' this is the first time function is being called
 */
void array2sh_initArray(void * const hPars,
                        ARRAY2SH_MICROPHONE_ARRAY_PRESETS preset,
                        int* arrayOrder,
                        int firstInitFLAG);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __ARRAY2SH_INTERNAL_H_INCLUDED__ */
