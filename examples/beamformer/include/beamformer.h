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
 * @file beamformer.h
 * @brief Generates beamformers/virtual microphones in arbitrary directions
 *        with several different beam pattern to choose from
 * 
 * @author Leo McCormack
 * @date 17.05.2019
 */

#ifndef __BEAMFORMER_H_INCLUDED__
#define __BEAMFORMER_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/**
 * Available beamforming orders
 */
typedef enum _BEAMFORMER_BEAM_ORDERS{
    BEAM_ORDER_FIRST = 1, /**< First-order beamforming (4 channel input) */
    BEAM_ORDER_SECOND,    /**< Second-order beamforming (9 channel input) */
    BEAM_ORDER_THIRD,     /**< Third-order beamforming (16 channel input) */
    BEAM_ORDER_FOURTH,    /**< Fourth-order beamforming (25 channel input) */
    BEAM_ORDER_FIFTH,     /**< Fifth-order beamforming (36 channel input) */
    BEAM_ORDER_SIXTH,     /**< Sixth-order beamforming (49 channel input) */
    BEAM_ORDER_SEVENTH    /**< Seventh-order beamforming (64 channel input) */
    
}BEAMFORMER_BEAM_ORDERS;
    
/** Maximum supported Ambisonic order */
#define BEAMFORMER_MAX_SH_ORDER ( 7 )
    
/**
 * Available beamforming approaches
 */
typedef enum _BEAMFORMER_BEAM_TYPES {
    BEAM_TYPE_CARDIOID = 1,  /**< cardioid */
    BEAM_TYPE_HYPERCARDIOID, /**< hyper-cardioid */
    BEAM_TYPE_MAX_EV         /**< hyper-cardioid with max_rE weighting */
    
} BEAMFORMER_BEAM_TYPES;
    
/* Number of available beamformer types */
#define BEAMFORMER_NUM_BEAM_TYPES ( 3 )
 
/**
 * Available Ambisonic channel ordering conventions
 *
 * @note CH_FUMA only supported for 1st order input.
 */
typedef enum _BEAMFORMER_CH_ORDER {
    CH_ACN = 1, /**< Ambisonic Channel Numbering (ACN) */
    CH_FUMA     /**< (Obsolete) Furse-Malham/B-format (WXYZ) */
    
} BEAMFORMER_CH_ORDER;

/** Number of channel ordering options */
#define BEAMFORMER_NUM_CH_ORDERINGS ( 2 )

/**
 * Available Ambisonic normalisation conventions
 *
 * @note NORM_FUMA only supported for 1st order input and does NOT have the
 *       1/sqrt(2) scaling on the omni.
 */
typedef enum _BEAMFORMER_NORM_TYPES {
    NORM_N3D = 1, /**< orthonormalised (N3D) */
    NORM_SN3D,    /**< Schmidt semi-normalisation (SN3D) */
    NORM_FUMA     /**< (Obsolete) Same as NORM_SN3D for 1st order */
    
} BEAMFORMER_NORM_TYPES;

/** Number of normalisation options */
#define BEAMFORMER_NUM_NORM_TYPES ( 3 )
    
/** Maximum number of beams supported */
#define BEAMFORMER_MAX_NUM_BEAMS ( 64 )
    
    
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of beamformer
 *
 * @param[in] phBeam (&) address of beamformer handle
 */
void beamformer_create(void** const phBeam);

/**
 * Destroys an instance of beamformer
 *
 * @param[in] phBeam (&) address of beamformer handle
 */
void beamformer_destroy(void** const phBeam);

/**
 * Initialises an instance of beamformer with default settings
 *
 * @param[in] hBeam      beamformer handle
 * @param[in] samplerate Host samplerate.
 */
void beamformer_init(void* const hBeam,
                     int samplerate);

/**
 * Generates beamformers/virtual microphones in the specified directions
 *
 * @param[in] hBeam     beamformer handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs   Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nOutputs  Number of output channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 */
void beamformer_process(void* const hBeam,
                        float** const inputs,
                        float** const outputs,
                        int nInputs,
                        int nOutputs,
                        int nSamples);


/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1; re-initialising all settings/variables
 * as beamformer is currently configured, at next available opportunity.
 */
void beamformer_refreshSettings(void* const hBeam);
    
/**
 * Sets the beamforming order (see 'BEAMFORMER_BEAM_ORDERS' enum)
 *
 * If the beamforming order is higher than the
 * input signal order, the extra required channels are filled with zeros. If the
 * beamforming order is lower than the input signal order, the number input
 * signals is truncated accordingly.
 */
void beamformer_setBeamOrder(void* const hBeam,  int newValue);

/**
 * Sets a beamformer azimuth direction of a given index, in DEGREES
 */
void beamformer_setBeamAzi_deg(void* const hBeam, int index, float newAzi_deg);

/**
 * Sets a beamformer elevation direction for a given index, in DEGREES
 */
void beamformer_setBeamElev_deg(void* const hBeam, int index, float newElev_deg);

/**
 * Sets the number of beamformers to generate
 */
void beamformer_setNumBeams(void* const hBeam, int new_nBeams);

/**
 * Sets the Ambisonic channel ordering convention to decode with, in order to
 * match the convention employed by the input signals (see 'BEAMFORMER_CH_ORDER'
 * enum)
 */
void beamformer_setChOrder(void* const hBeam, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to decode with, in order to match
 * with the convention employed by the input signals (see 'BEAMFORMER_NORM_TYPE'
 * enum)
 */
void beamformer_setNormType(void* const hBeam, int newType);
    
/**
 * Sets the beamforming approach to employ (see 'BEAMFORMER_BEAM_TYPE' enum)
 */
void beamformer_setBeamType(void* const hBeam, int newID);

    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns tje beamforming order (see 'BEAMFORMER_BEAM_ORDERS' enum)
 */
int beamformer_getBeamOrder(void* const hBeam);
    
/**
 * Returns the number of frequency bands used by beamformer (only for adaptive
 * beamformer algorithms).
 */
int beamformer_getNumberOfBands(void);

/**
 * Returns the beamformer azimuth direction of a given index h, in DEGREES
 */
float beamformer_getBeamAzi_deg(void* const hBeam, int index);

/**
 * Returns the beamformer elevation direction of a given index, in DEGREES
 */
float beamformer_getBeamElev_deg(void* const hBeam, int index);

/**
 * Returns the number of beamformers being generated
 */
int beamformer_getNumBeams(void* const hBeam);
    
/**
 * Returns the maximum number of beamformers permitted
 */
int beamformer_getMaxNumBeams(void);
    
/**
 * Returns the number of spherical harmonic signals required by the currently
 * selected beamforming order: (current_order+1)^2
 */
int  beamformer_getNSHrequired(void* const hBeam);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * decode with, which should match the convention employed by the input signals
 * (see 'BEAMFORMER_CH_ORDER' enum)
 */
int beamformer_getChOrder(void* const hBeam);

/**
 * Returns the Ambisonic normalisation convention currently being usedto decode
 * with, which should match the convention employed by the input signals.
 * (see 'BEAMFORMER_NORM_TYPE' enum)
 */
int beamformer_getNormType(void* const hBeam);
    
/**
 * Returns the beamforming approach employed (see 'BEAMFORMER_BEAM_TYPE' enum)
 */
int beamformer_getBeamType(void* const hBeam);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int beamformer_getProcessingDelay(void);

    
    
#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __BEAMFORMER_H_INCLUDED__ */
