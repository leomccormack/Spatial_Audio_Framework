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
 * @file ambi_roomsim.h
 * @brief A simple shoebox room Ambisonic encoder.
 *
 * @author Leo McCormack
 * @date 10.08.2020
 */

#ifndef __AMBI_ROOMSIM_H_INCLUDED__
#define __AMBI_ROOMSIM_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

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
 * 
 *
 * @param[in] hAmbi    ambi_roomsim handle
 * @param[in] inputs   Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs  Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs  Number of input channels
 * @param[in] nOutputs Number of output channels
 * @param[in] nSamples Number of samples in 'inputs'/'output' matrices
 */
void ambi_roomsim_process(void* const hAmbi,
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
 * as ambi_roomsim is currently configured, at next available opportunity.
 */
void ambi_roomsim_refreshParams(void* const hAmbi);
    
/**
 * Sets the encoding order (see #SH_ORDERS enum)
 */
void ambi_roomsim_setOutputOrder(void* const hAmbi, int newValue);

/**
 * Sets the number of input signals/sources to encode.
 */
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


/**
 * Sets the input configuration preset (see #SOURCE_CONFIG_PRESETS enum)
 */
void ambi_roomsim_setInputConfigPreset(void* const hAmbi, int newPresetID);

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

/**
 * By default, ambi_roomsim will scale the output signals by the number of input
 * signals.
 */
void ambi_roomsim_setEnablePostScaling(void* const hAmbi, int newStatus);

    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int ambi_roomsim_getFrameSize(void);

/**
 * Returns the decoding order (see #SH_ORDERS enum)
 *
 * If decoding order is higher than the input signal order, the extra required
 * channels are filled with zeros. If the decoding order is lower than the input
 * signal order, the number input signals is truncated accordingly.
 */
int ambi_roomsim_getOutputOrder(void* const hAmbi);

/**
 * Returns the azimuth for a specific source, in DEGREES
 */
float ambi_roomsim_getSourceAzi_deg(void* const hAmbi, int index);

/**
 * Returns the elevation for a specific source, in DEGREES
 */
float ambi_roomsim_getSourceElev_deg(void* const hAmbi, int index);

/**
 * Returns the number of input signals/sources to encode.
 */
int ambi_roomsim_getNumSources(void* const hAmbi);

/**
 * Returns the maximum number of input signals/sources that can be encoded.
 */
int ambi_roomsim_getMaxNumSources(void);

/**
 * Returns the number of spherical harmonic signals required by the current
 * decoding order: (current_order+1)^2
 */
int ambi_roomsim_getNSHrequired(void* const hAmbi);

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
 * Returns 0: if post scaling is disabled, 1: if post scaling is enabled
 */
int ambi_roomsim_getEnablePostScaling(void* const hAmbi);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features)
 */
int ambi_roomsim_getProcessingDelay(void);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __AMBI_ROOMSIM_H_INCLUDED__ */
