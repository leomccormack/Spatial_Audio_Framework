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
 * @example panner.h
 * @brief A frequency-dependent 3D panner based on the Vector-base Amplitude
 *        Panning (VBAP) method, with an optional spread control
 * 
 * ### Files
 * panner.h (include), panner_internal.h, panner.c, panner_internal.c
 * ### Include Header
 */

/**
 * @file panner.h
 * @brief A frequency-dependent 3D panner based on the Vector-base Amplitude
 *        Panning (VBAP) method [1], with an optional spread control [2].
 *
 * Depending on the listening room, it may be beneficial to employ amplitude-
 * normalised gains for low frequencies, and energy-normalised gains for high
 * frequencies. Therefore, this VBAP implementation also uses the method
 * described in [3], to do just that.
 *
 * @see [1] Pulkki, V. (1997). Virtual sound source positioning using vector
 *          base amplitude panning. Journal of the audio engineering society,
 *          45(6), 456-466.
 * @see [2] Pulkki, V. (1999). Uniform spreading of amplitude panned virtual
 *          sources. In Proceedings of the 1999 IEEE Workshop on Applications of
 *          Signal Processing to Audio and Acoustics. WASPAA'99 (Cat. No.
 *          99TH8452) (pp. 187-190). IEEE.
 * @see [3] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V.
 *          (2014). Gain normalisation in amplitude panning as a function of
 *          frequency and room reverberance. 55th International Conference of
 *          the AES. Helsinki, Finland.
 *
 * @author Leo McCormack
 * @date 25.09.2017
 * @license ISC
 */

#ifndef __PANNER_H_INCLUDED__
#define __PANNER_H_INCLUDED__

#include "_common.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */

/** Minimum supported spread angle, degrees */
#define PANNER_SPREAD_MIN_VALUE ( 0.0f )

/** Maximum supported spread angle, degrees */
#define PANNER_SPREAD_MAX_VALUE ( 90.0f )
    

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * Creates an instance of the panner
 *
 * @param[in] phPan (&) address of panner handle
 */
void panner_create(void** const phPan);

/**
 * Destroys an instance of the panner
 *
 * @param[in] phPan (&) address of panner handle
 */
void panner_destroy(void** const phPan);

/**
 * Initialises an instance of panner with default settings
 *
 * @warning This should not be called while _process() is on-going!
 *
 * @param[in] hPan       panner handle
 * @param[in] samplerate Host samplerate.
 */
void panner_init(void* const hPan,
                 int samplerate);
    
/**
 * Intialises the codec variables, based on current global/user parameters
 *
 * @note This function is fully threadsafe. It can even be called periodically
 *       via a timer on one thread, while calling _process() on another thread.
 *       Since, if a set function is called (that warrants a re-init), then a
 *       flag is triggered internally and the next time this function is called,
 *       it will wait until the current process() function has completed before
 *       reinitialising the relevant parameters. If the _initCodec() takes
 *       longer than the time it takes for process() to be called again, then
 *       process() is simply bypassed until the codec is ready.
 * @note This function does nothing if no re-initialisations are required.
 *
 * @param[in] hPan panner handle
 */
void panner_initCodec(void* const hPan);

/**
 * Pans the input signals/sources to the loudspeaker channels using VBAP [1],
 * and optional spreading [2] and frequency-dependent normalisation as a
 * function of the room reverberation [3].
 *
 * @param[in] hPan      panner handle
 * @param[in] inputs    Input channel buffers; 2-D array: nInputs x nSamples
 * @param[in] outputs   Output channel buffers; 2-D array: nOutputs x nSamples
 * @param[in] nInputs   Number of input channels
 * @param[in] nOutputs  Number of output channels
 * @param[in] nSamples  Number of samples in 'inputs'/'output' matrices
 *
 * @see [1] Pulkki, V. (1997). Virtual sound source positioning using vector
 *          base amplitude panning. Journal of the audio engineering society,
 *          45(6), 456-466.
 * @see [2] Pulkki, V. (1999). Uniform spreading of amplitude panned virtual
 *          sources. In Proceedings of the 1999 IEEE Workshop on Applications of
 *          Signal Processing to Audio and Acoustics. WASPAA'99 (Cat. No.
 *          99TH8452) (pp. 187-190). IEEE.
 * @see [3] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V.
 *          (2014). Gain normalisation in amplitude panning as a function of
 *          frequency and room reverberance. 55th International Conference of
 *          the AES. Helsinki, Finland.
 */
void panner_process(void* const hPan,
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
 * as panner is currently configured, at next available opportunity.
 */
void panner_refreshSettings(void* const hPan);

/** Sets the azimuth of a specific input/source index, in DEGREES */
void panner_setSourceAzi_deg(void* const hPan, int index, float newAzi_deg);

/** Sets the elevation of a specific input/source index, in DEGREES */
void panner_setSourceElev_deg(void* const hPan, int index, float newElev_deg);

/** Sets the number of inputs/sources to pan */
void panner_setNumSources(void* const hPan, int new_nSources);

/** Sets the azimuth of a specific loudspeaker index, in DEGREES */
void panner_setLoudspeakerAzi_deg(void* const hPan, int index, float newAzi_deg);

/** Sets the elevation of a specific loudspeaker index, in DEGREES */
void panner_setLoudspeakerElev_deg(void* const hPan, int index, float newElev_deg);

/** Sets the number of loudspeakers to pan to */
void panner_setNumLoudspeakers(void* const hPan, int new_nLoudspeakers);

/**
 * Sets a preset for the output configuration (see #LOUDSPEAKER_ARRAY_PRESETS
 * enum)
 */
void panner_setOutputConfigPreset(void* const hPan, int newPresetID);

/**
 * Sets a preset for the input configuration (see #SOURCE_CONFIG_PRESETS enum)
 */
void panner_setInputConfigPreset(void* const hPan, int newPresetID);

/**
 * Sets the room coefficient value 0..1 [1]; 0: normal room, 0.5: dry listening
 * room, 1: anechoic
 *
 * @see [1] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V.
 *          (2014). Gain normalisation in amplitude panning as a function of
 *          frequency and room reverberance. 55th International Conference of
 *          the AES. Helsinki, Finland.
 */
void panner_setDTT(void* const hPan, float newValue);
    
/** Sets the degree of spread, in DEGREES */
void panner_setSpread(void* const hPan, float newValue);

/** Sets the 'yaw' rotation angle, in DEGREES */
void panner_setYaw(void* const hPan, float newYaw);

/** Sets the 'pitch' rotation angle, in DEGREES */
void panner_setPitch(void* const hPan, float newPitch);

/** Sets the 'roll' rotation angle, in DEGREES */
void panner_setRoll(void* const hPan, float newRoll);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'yaw' angle
 * (0: do not flip sign, 1: flip the sign)
 */
void panner_setFlipYaw(void* const hPan, int newState);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'pitch' angle
 * (0: do not flip sign, 1: flip the sign)
 */
void panner_setFlipPitch(void* const hPan, int newState);

/**
 * Sets a flag as to whether to "flip" the sign of the current 'roll' angle
 * (0: do not flip sign, 1: flip the sign)
 */
void panner_setFlipRoll(void* const hPan, int newState);
    
    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/**
 * Returns the processing framesize (i.e., number of samples processed with
 * every _process() call )
 */
int panner_getFrameSize(void);

/** Returns current codec status (see #CODEC_STATUS enum) */
CODEC_STATUS panner_getCodecStatus(void* const hPan);

/**
 * (Optional) Returns current intialisation/processing progress, between 0..1
 *  - 0: intialisation/processing has started
 *  - 1: intialisation/processing has ended
 */
float panner_getProgressBar0_1(void* const hPan);

/**
 * (Optional) Returns current intialisation/processing progress text
 *
 * @note "text" string should be (at least) of length:
 *       #PROGRESSBARTEXT_CHAR_LENGTH
 */
void panner_getProgressBarText(void* const hPan, char* text);
    
/** Returns the input/source azimuth for a given index, in DEGREES */
float panner_getSourceAzi_deg(void* const hPan, int index);

/** Returns the input/source elevation for a given index, in DEGREES */
float panner_getSourceElev_deg(void* const hPan, int index);

/** Returns the number of inputs/sources in the current layout */
int panner_getNumSources(void* const hPan);

/** Returns the maximum number of inputs/sources permitted by panner */
int panner_getMaxNumSources(void);

/** Returns the loudspeaker azimuth for a given index, in DEGREES */
float panner_getLoudspeakerAzi_deg(void* const hPan, int index);

/** Returns the loudspeaker elevation for a given index, in DEGREES */
float panner_getLoudspeakerElev_deg(void* const hPan, int index);

/** Returns the number of loudspeakers in the current layout */
int panner_getNumLoudspeakers(void* const hPan);

/** Returns the maximum number of loudspeakers permitted */
int panner_getMaxNumLoudspeakers(void);

/** Returns the DAW/Host sample rate */
int panner_getDAWsamplerate(void* const hPan);

/** Returns the room coefficient value 0..1 */
float panner_getDTT(void* const hPan);

/** Returns the spread value, in DEGREES */
float panner_getSpread(void* const hPan);

/** Returns the 'yaw' rotation angle, in DEGREES */
float panner_getYaw(void* const hPan);

/** Returns the 'pitch' rotation angle, in DEGREES */
float panner_getPitch(void* const hPan);

/** Returns the 'roll' rotation angle, in DEGREES */
float panner_getRoll(void* const hPan);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'yaw' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int panner_getFlipYaw(void* const hPan);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'pitch' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int panner_getFlipPitch(void* const hPan);

/**
 * Returns a flag as to whether to "flip" the sign of the current 'roll' angle
 * (0: do not flip sign, 1: flip the sign)
 */
int panner_getFlipRoll(void* const hPan);

/**
 * Returns the processing delay in samples (may be used for delay compensation
 * features) 
 */
int panner_getProcessingDelay(void);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __PANNER_H_INCLUDED__ */
