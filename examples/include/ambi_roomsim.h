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
 * @example ambi_roomsim.h
 * @brief A simple shoebox room Ambisonic encoder.
 *
 * ### Files
 * ambi_roomsim.h (include), ambi_roomsim_internal.h, ambi_roomsim.c,
 * ambi_roomsim_internal.c
 * ### Include Header
 */

/**
 * @file ambi_roomsim.h
 * @brief A simple shoebox room Ambisonic encoder.
 *
 * @author Leo McCormack
 * @date 10.08.2020
 * @license ISC
 */

#ifndef __AMBI_ROOMSIM_H_INCLUDED__
#define __AMBI_ROOMSIM_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/** Maximum supported number of receivers for the room sim example */
#define ROOM_SIM_MAX_NUM_RECEIVERS ( 16 )
/** Maximum supported number of sources for the room sim example */
#define ROOM_SIM_MAX_NUM_SOURCES ( 16 )

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of ambi_roomsim
 *
 * @param[in] phAmbi (&) address of ambi_roomsim handle
 */
void ambi_roomsim_create(void** const phAmbi);

/**
 * Destroys an instance of ambi_roomsim
 *
 * @param[in] phAmbi (&) address of ambi_roomsim handle
 */
void ambi_roomsim_destroy(void** const phAmbi);

/**
 * Initialises an instance of ambi_roomsim with default settings
 *
 * @param[in] hAmbi      ambi_roomsim handle
 * @param[in] samplerate Host samplerate.
 */
void ambi_roomsim_init(void* const hAmbi,
                     int samplerate);

/**
 * Processes audio
 *
 * @param[in] hAmbi    ambi_roomsim handle
 * @param[in] inputs   Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs  Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs  Number of input channels
 * @param[in] nOutputs Number of output channels
 * @param[in] nSamples Number of samples in 'inputs'/'output' matrices
 */
void ambi_roomsim_process(void* const hAmbi,
                          const float *const * inputs,
                          float** const outputs,
                          int nInputs,
                          int nOutputs,
                          int nSamples);

    
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/**
 * Sets all intialisation flags to 1; re-initialising all settings/variables
 * as ambi_roomsim is currently configured, at next available opportunity.
 */
void ambi_roomsim_refreshParams(void* const hAmbi);

/** Sets whether to include image sources (1) or not (0) */
void ambi_roomsim_setEnableIMSflag(void* const hAmbi, int newValue);

/** Sets the maximum reflection order */
void ambi_roomsim_setMaxReflectionOrder(void* const hAmbi, int newValue);
    
/** Sets the encoding order (see #SH_ORDERS enum) */
void ambi_roomsim_setOutputOrder(void* const hAmbi, int newValue);

/** Sets the number of input signals/sources to encode */
void ambi_roomsim_setNumSources(void* const hAmbi, int new_nSources);

/**
 * Sets the 'x' coordinate for a specific source index
 *
 * @param[in] hAmbi      ambi_roomsim handle
 * @param[in] index      Source index
 * @param[in] newValue   New 'x' coordinate, in metres
 */
void ambi_roomsim_setSourceX(void* const hAmbi, int index, float newValue);

/**
 * Sets the 'y' coordinate for a specific source index
 *
 * @param[in] hAmbi      ambi_roomsim handle
 * @param[in] index      Source index
 * @param[in] newValue   New 'y' coordinate, in metres
 */
void ambi_roomsim_setSourceY(void* const hAmbi, int index, float newValue);

/**
 * Sets the 'z' coordinate for a specific source index
 *
 * @param[in] hAmbi      ambi_roomsim handle
 * @param[in] index      Source index
 * @param[in] newValue   New 'z' coordinate, in metres
 */
void ambi_roomsim_setSourceZ(void* const hAmbi, int index, float newValue);

/** Sets the number of input SH receivers */
void ambi_roomsim_setNumReceivers(void* const hAmbi, int new_nReceievers);

/**
 * Sets the 'x' coordinate for a specific receiver index
 *
 * @param[in] hAmbi      ambi_roomsim handle
 * @param[in] index      Receiver index
 * @param[in] newValue   New 'x' coordinate, in metres
 */
void ambi_roomsim_setReceiverX(void* const hAmbi, int index, float newValue);

/**
 * Sets the 'y' coordinate for a specific receiver index
 *
 * @param[in] hAmbi      ambi_roomsim handle
 * @param[in] index      Receiver index
 * @param[in] newValue   New 'y' coordinate, in metres
 */
void ambi_roomsim_setReceiverY(void* const hAmbi, int index, float newValue);

/**
 * Sets the 'z' coordinate for a specific receiver index
 *
 * @param[in] hAmbi      ambi_roomsim handle
 * @param[in] index      Receiver index
 * @param[in] newValue   New 'z' coordinate, in metres
 */
void ambi_roomsim_setReceiverZ(void* const hAmbi, int index, float newValue);

/** Sets the room length along the x dimension */
void ambi_roomsim_setRoomDimX(void* const hAmbi, float newValue);

/** Sets the room length along the y dimension */
void ambi_roomsim_setRoomDimY(void* const hAmbi, float newValue);

/** Sets the room length along the z dimension */
void ambi_roomsim_setRoomDimZ(void* const hAmbi, float newValue);

/** Sets wall absorption coefficients */
void ambi_roomsim_setWallAbsCoeff(void* const hAmbi, int xyz_idx, int posNeg_idx, float new_value);

/**
 * Sets the Ambisonic channel ordering convention to encode with (see #CH_ORDER
 * enum)
 */
void ambi_roomsim_setChOrder(void* const hAmbi, int newOrder);

/**
 * Sets the Ambisonic normalisation convention to encode with (see #NORM_TYPES
 * enum)
 */
void ambi_roomsim_setNormType(void* const hAmbi, int newType);

    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int ambi_roomsim_getFrameSize(void);

/** Returns whether to include image sources (1) or not (0) */
int ambi_roomsim_getEnableIMSflag(void* const hAmbi);

/** Returns the maximum reflection order */
int ambi_roomsim_getMaxReflectionOrder(void* const hAmbi);

/**
 * Returns the decoding order (see #SH_ORDERS enum)
 *
 * If decoding order is higher than the input signal order, the extra required
 * channels are filled with zeros. If the decoding order is lower than the input
 * signal order, the number input signals is truncated accordingly.
 */
int ambi_roomsim_getOutputOrder(void* const hAmbi);

/** Returns the number of input signals/sources to encode. */
int ambi_roomsim_getNumSources(void* const hAmbi);

/** Returns the maximum number of input signals/sources that can be encoded. */
int ambi_roomsim_getMaxNumSources(void);

/**
 * Returns the number of spherical harmonic signals required by the current
 * decoding order: (current_order+1)^2
 */
int ambi_roomsim_getNSHrequired(void* const hAmbi);

/** Returns the 'x' coordinate for a specific source index */
float ambi_roomsim_getSourceX(void* const hAmbi, int index);

/** Returns the 'y' coordinate for a specific source index */
float ambi_roomsim_getSourceY(void* const hAmbi, int index);

/** Returns the 'z' coordinate for a specific source index */
float ambi_roomsim_getSourceZ(void* const hAmbi, int index);

/** Returns the number of SH receivers */
int ambi_roomsim_getNumReceivers(void* const hAmbi);

/** Returns the maximum number of receivers */
int ambi_roomsim_getMaxNumReceivers(void);

/** Returns the 'x' coordinate for a specific receiver index */
float ambi_roomsim_getReceiverX(void* const hAmbi, int index);

/** Returns the 'y' coordinate for a specific receiver index */
float ambi_roomsim_getReceiverY(void* const hAmbi, int index);

/** Returns the 'z' coordinate for a specific receiver index */
float ambi_roomsim_getReceiverZ(void* const hAmbi, int index);

/** Returns the room length along the x dimension */
float ambi_roomsim_getRoomDimX(void* const hAmbi);

/** Returns the room length along the y dimension */
float ambi_roomsim_getRoomDimY(void* const hAmbi);

/** Returns the room length along the z dimension */
float ambi_roomsim_getRoomDimZ(void* const hAmbi);

/** Returns wall absorption coefficients */
float ambi_roomsim_getWallAbsCoeff(void* const hAmbi,
                                   int xyz_idx,
                                   int posNeg_idx);

/**
 * Returns the Ambisonic channel ordering convention currently being used to
 * encode with (see #CH_ORDER enum)
 */
int ambi_roomsim_getChOrder(void* const hAmbi);
    
/**
 * Returns the Ambisonic normalisation convention currently being used to encode
 * with (see #NORM_TYPES enum)
 */
int ambi_roomsim_getNormType(void* const hAmbi);
 
/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int ambi_roomsim_getProcessingDelay(void);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __AMBI_ROOMSIM_H_INCLUDED__ */
