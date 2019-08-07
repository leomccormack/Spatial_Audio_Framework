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
 * Filename: panner.h (include header)
 * -----------------------------------
 * A frequency-dependent 3D panner, based on the Vector-base Amplitude Panning
 * (VBAP) method. Depending on the room, it may be beneficial to employ
 * amplitude-normalised gains for low frequencies, and energy-normalised gains
 * for high frequencies. Therefore, this VBAP implementation uses the method
 * described in [1], to do just that.
 *
 * Dependencies:
 *     saf_utilities, saf_vbap, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 25.09.2017
 *
 * [1] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V. (2014).
 *     Gain normalisation in amplitude panning as a function of frequency and
 *     room reverberance. 55th International Conference of the AES. Helsinki,
 *     Finland.
 */

#ifndef __PANNER_H_INCLUDED__
#define __PANNER_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                             Presets + Constants                            */
/* ========================================================================== */
     
#define ENABLE_MONO_PRESET
#define ENABLE_STEREO_PRESET
#define ENABLE_5PX_PRESET
#define ENABLE_7PX_PRESET
#define ENABLE_8PX_PRESET
#define ENABLE_9PX_PRESET
#define ENABLE_10PX_PRESET
#define ENABLE_11PX_PRESET
#define ENABLE_11PX_7_4_PRESET
#define ENABLE_13PX_PRESET
#define ENABLE_22PX_PRESET
#define ENABLE_AALTO_MCC_PRESET
#define ENABLE_AALTO_APAJA_PRESET
#define ENABLE_AALTO_APAJA2_PRESET
#define ENABLE_AALTO_LR_PRESET
#define ENABLE_DTU_AVIL_PRESET
#define ENABLE_T_DESIGN_4_PRESET
#define ENABLE_T_DESIGN_12_PRESET
#define ENABLE_T_DESIGN_24_PRESET
#define ENABLE_T_DESIGN_36_PRESET
#define ENABLE_T_DESIGN_48_PRESET
#define ENABLE_T_DESIGN_60_PRESET
    
typedef enum _PRESETS{
    PRESET_DEFAULT = 1
#ifdef ENABLE_MONO_PRESET
    , PRESET_MONO
#endif
#ifdef ENABLE_STEREO_PRESET
    , PRESET_STEREO
#endif
#ifdef ENABLE_5PX_PRESET
    , PRESET_5PX
#endif
#ifdef ENABLE_7PX_PRESET
    , PRESET_7PX
#endif
#ifdef ENABLE_8PX_PRESET
    , PRESET_8PX
#endif
#ifdef ENABLE_9PX_PRESET
    , PRESET_9PX
#endif
#ifdef ENABLE_10PX_PRESET
    , PRESET_10PX
#endif
#ifdef ENABLE_11PX_PRESET
    , PRESET_11PX
#endif
#ifdef ENABLE_11PX_7_4_PRESET
    , PRESET_11PX_7_4
#endif
#ifdef ENABLE_13PX_PRESET
    , PRESET_13PX
#endif
#ifdef ENABLE_22PX_PRESET
    , PRESET_22PX
#endif
#ifdef ENABLE_AALTO_MCC_PRESET
    , PRESET_AALTO_MCC
#endif
#ifdef ENABLE_AALTO_APAJA_PRESET
    , PRESET_AALTO_APAJA
#endif
#ifdef ENABLE_AALTO_APAJA2_PRESET
    , PRESET_AALTO_APAJA2
#endif
#ifdef ENABLE_AALTO_LR_PRESET
    , PRESET_AALTO_LR
#endif
#ifdef ENABLE_DTU_AVIL_PRESET
    , PRESET_DTU_AVIL
#endif
#ifdef ENABLE_T_DESIGN_4_PRESET
    , PRESET_T_DESIGN_4
#endif
#ifdef ENABLE_T_DESIGN_12_PRESET
    , PRESET_T_DESIGN_12
#endif
#ifdef ENABLE_T_DESIGN_24_PRESET
    , PRESET_T_DESIGN_24
#endif
#ifdef ENABLE_T_DESIGN_36_PRESET
    , PRESET_T_DESIGN_36
#endif
#ifdef ENABLE_T_DESIGN_48_PRESET
    , PRESET_T_DESIGN_48
#endif
#ifdef ENABLE_T_DESIGN_60_PRESET
    , PRESET_T_DESIGN_60
#endif
    
}PRESETS;
    
#define PANNER_MAX_NUM_INPUTS ( 64 )
#define PANNER_MAX_NUM_OUTPUTS ( 64 )
#define PANNER_SPREAD_MIN_VALUE ( 0.0f )
#define PANNER_SPREAD_MAX_VALUE ( 90.0f )
    

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/*
 * Function: panner_create
 * -----------------------
 * Creates an instance of the panner
 *
 * Input Arguments:
 *     phPan - & address of panner handle
 */
void panner_create(void** const phPan);

/*
 * Function: panner_destroy
 * ------------------------
 * Destroys an instance of the panner
 *
 * Input Arguments:
 *     phPan - & address of panner handle
 */
void panner_destroy(void** const phPan);

/*
 * Function: panner_init
 * ---------------------
 * Initialises an instance of panner with default settings
 *
 * Input Arguments:
 *     hPan       - panner handle
 *     samplerate - host samplerate.
 */
void panner_init(void* const hPan,
                 int samplerate);

/*
 * Function: panner_process
 * ------------------------
 * Pans the input signals/sources to the loudspeaker channels.
 *
 * Input Arguments:
 *     hPan      - panner handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     outputs   - Output channel buffers; 2-D array: nOutputs x nSamples
 *     nInputs   - number of input channels
 *     nOutputs  - number of output channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 *     isPlaying - flag to say if there is audio in the input buffers, 0: no
 *                 audio, reduced processing, 1: audio, full processing
 */
void panner_process(void* const hPan,
                    float** const inputs,
                    float** const outputs,
                    int nInputs,
                    int nOutputs,
                    int nSamples,
                    int isPlaying);
    
    
/* ========================================================================== */
/*                                Set Functions                               */
/* ========================================================================== */

/*
 * Function: panner_refreshSettings
 * --------------------------------
 * Sets all intialisation flags to 1. i.e. re-initialise all settings/variables
 * as panner is currently configured, at next available opportunity.
 *
 * Input Arguments:
 *     hPan - panner handle
 */
void panner_refreshSettings(void* const hPan);

/*
 * Function: panner_checkReInit
 * ----------------------------
 * Check if any reInit Flags are active, and reinitialise if they are.
 * Note: Only call when playback has stopped.
 *
 * Input Arguments:
 *     hPan - panner handle
 */
void panner_checkReInit(void* const hPan);

/*
 * Function: panner_setSourceAzi_deg
 * ---------------------------------
 * Sets the azimuth of a specific input/source
 *
 * Input Arguments:
 *     hPan       - panner handle
 *     index      - input/source index
 *     newAzi_deg - new azimuth in DEGREES
 */
void panner_setSourceAzi_deg(void* const hPan, int index, float newAzi_deg);

/*
 * Function: panner_setSourceElev_deg
 * ----------------------------------
 * Sets the elevation of a specific input/source
 *
 * Input Arguments:
 *     hPan        - panner handle
 *     index       - input/source index
 *     newElev_deg - new elevation in DEGREES
 */
void panner_setSourceElev_deg(void* const hPan, int index, float newElev_deg);

/*
 * Function: panner_setNumSources
 * ------------------------------
 * Sets the number of inputs/sources to pan.
 *
 * Input Arguments:
 *     hPan         - panner handle
 *     new_nSources - new number of inputs/sources
 */
void panner_setNumSources(void* const hPan, int new_nSources);

/*
 * Function: panner_setLoudspeakerAzi_deg
 * --------------------------------------
 * Sets the azimuth of a specific loudspeaker
 *
 * Input Arguments:
 *     hPan       - panner handle
 *     index      - loudspeaker index
 *     newAzi_deg - new azimuth in DEGREES
 */
void panner_setLoudspeakerAzi_deg(void* const hPan, int index, float newAzi_deg);

/*
 * Function: panner_setLoudspeakerElev_deg
 * ---------------------------------------
 * Sets the elevation of a specific loudspeaker
 *
 * Input Arguments:
 *     hPan        - panner handle
 *     index       - loudspeaker index
 *     newElev_deg - new elevation in DEGREES
 */
void panner_setLoudspeakerElev_deg(void* const hPan, int index, float newElev_deg);

/*
 * Function: panner_setNumLoudspeakers
 * -----------------------------------
 * Sets the number of loudspeakers to pan to.
 *
 * Input Arguments:
 *     hPan              - panner handle
 *     new_nLoudspeakers - new number of loudspeakers
 */
void panner_setNumLoudspeakers(void* const hPan, int new_nLoudspeakers);

/*
 * Function: panner_setOutputConfigPreset
 * --------------------------------------
 * Sets a preset for the output configuration. See "PRESET" enum
 *
 * Input Arguments:
 *     hPan        - panner handle
 *     newPresetID - See "PRESET" enum
 */
void panner_setOutputConfigPreset(void* const hPan, int newPresetID);

/*
 * Function: panner_setInputConfigPreset
 * -------------------------------------
 * Sets a preset for the input configuration. See "PRESET" enum
 *
 * Input Arguments:
 *     hPan        - panner handle
 *     newPresetID - See "PRESET" enum
 */
void panner_setInputConfigPreset(void* const hPan, int newPresetID);

/*
 * Function: panner_setDTT
 * -----------------------
 * Sets the room coefficient value [1].
 *
 * Input Arguments:
 *     hPan     - panner handle
 *     newValue - new room coefficient value; 0: normal room, 0.5: listening
 *                room, 1: anechoic
 *
 * [1] Laitinen, M., Vilkamo, J., Jussila, K., Politis, A., Pulkki, V. (2014).
 *     Gain normalisation in amplitude panning as a function of frequency and
 *     room reverberance. 55th International Conference of the AES. Helsinki,
 *     Finland.
 */
void panner_setDTT(void* const hPan, float newValue);
    
/*
 * Function: panner_setSpread
 * --------------------------
 * Sets the degree of spread
 *
 * Input Arguments:
 *     hPan     - panner handle
 *     newValue - spread, in DEGREES
 */
void panner_setSpread(void* const hPan, float newValue);

/*
 * Function: panner_setYaw
 * -----------------------
 * Sets the 'yaw' rotation angle
 *
 * Input Arguments:
 *     hPan   - panner handle
 *     newYaw - the 'yaw' rotation angle, in DEGREES
 */
void panner_setYaw(void* const hPan, float newYaw);

/*
 * Function: panner_setPitch
 * -------------------------
 * Sets the 'pitch' rotation angle
 *
 * Input Arguments:
 *     hPan     - panner handle
 *     newPitch - the 'pitch' rotation angle, in DEGREES
 */
void panner_setPitch(void* const hPan, float newPitch);

/*
 * Function: panner_setRoll
 * ------------------------
 * Sets the 'roll' rotation angle
 *
 * Input Arguments:
 *     hPan    - panner handle
 *     newRoll - the 'roll' rotation angle, in DEGREES
 */
void panner_setRoll(void* const hPan, float newRoll);

/*
 * Function: panner_setFlipYaw
 * ---------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'yaw' angle
 *
 * Input Arguments:
 *     hPan     - panner handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void panner_setFlipYaw(void* const hPan, int newState);

/*
 * Function: panner_setFlipPitch
 * -----------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'pitch' angle
 *
 * Input Arguments:
 *     hPan     - panner handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void panner_setFlipPitch(void* const hPan, int newState);

/*
 * Function: panner_setFlipRoll
 * ----------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'roll' angle
 *
 * Input Arguments:
 *     hPan     - panner handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void panner_setFlipRoll(void* const hPan, int newState);
    
    
/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/*
 * Function: panner_getSourceAzi_deg
 * ---------------------------------
 * Returns the input/source azimuth for a given index.
 *
 * Input Arguments:
 *     hPan - panner handle
 *     index - input/source index
 * Returns:
 *     input/source azimuth in DEGREES
 */
float panner_getSourceAzi_deg(void* const hPan, int index);

/*
 * Function: panner_getSourceElev_deg
 * ----------------------------------
 * Returns the input/source elevation for a given index.
 *
 * Input Arguments:
 *     hPan - panner handle
 *     index - input/source index
 * Returns:
 *     input/source elevation in DEGREES
 */
float panner_getSourceElev_deg(void* const hPan, int index);

/*
 * Function: panner_getNumSources
 * ------------------------------
 * Returns the number of inputs/sources in the current layout
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     number of inputs/sources
 */
int panner_getNumSources(void* const hPan);

/*
 * Function: panner_getMaxNumSources
 * ---------------------------------
 * Returns the maximum number of inputs/sources permitted
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     maximum number of inputs/sources
 */
int panner_getMaxNumSources(void);

/*
 * Function: panner_getLoudspeakerAzi_deg
 * --------------------------------------
 * Returns the loudspeaker azimuth for a given index.
 *
 * Input Arguments:
 *     hPan - panner handle
 *     index - loudspeaker index
 * Returns:
 *     loudspeaker azimuth in DEGREES
 */
float panner_getLoudspeakerAzi_deg(void* const hPan, int index);

/*
 * Function: panner_getLoudspeakerElev_deg
 * ---------------------------------------
 * Returns the loudspeaker elevation for a given index.
 *
 * Input Arguments:
 *     hPan - panner handle
 *     index - loudspeaker index
 * Returns:
 *     loudspeaker elevation in DEGREES
 */
float panner_getLoudspeakerElev_deg(void* const hPan, int index);

/*
 * Function: panner_getNumLoudspeakers
 * -----------------------------------
 * Returns the number of loudspeakers in the current layout
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     number of loudspeakers
 */
int panner_getNumLoudspeakers(void* const hPan);

/*
 * Function: panner_getNumLoudspeakers
 * -----------------------------------
 * Returns the maximum number of loudspeakers permitted
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     maximum number of loudspeakers
 */
int panner_getMaxNumLoudspeakers(void);

/*
 * Function: panner_getDAWsamplerate
 * ---------------------------------
 * Returns the DAW/Host sample rate
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     DAW/Host sampling rate
 */
int panner_getDAWsamplerate(void* const hPan);

/*
 * Function: panner_getDTT
 * -----------------------
 * Returns the room coefficient value
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     room coefficient value
 */
float panner_getDTT(void* const hPan);

/*
 * Function: panner_getSpread
 * --------------------------
 * Returns the spread value
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     spread, in DEGREES
 */
float panner_getSpread(void* const hPan);

/*
 * Function: panner_getYaw
 * -----------------------
 * Returns the 'yaw' rotation angle
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     the 'yaw' rotation angle, in DEGREES
 */
float panner_getYaw(void* const hPan);

/*
 * Function: panner_getPitch
 * -------------------------
 * Returns the 'pitch' rotation angle
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     the 'pitch' rotation angle, in DEGREES
 */
float panner_getPitch(void* const hPan);

/*
 * Function: panner_getRoll
 * ------------------------
 * Returns the 'roll' rotation angle
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     the 'roll' rotation angle, in DEGREES
 */
float panner_getRoll(void* const hPan);

/*
 * Function: panner_getFlipYaw
 * ---------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'yaw' angle
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int panner_getFlipYaw(void* const hPan);

/*
 * Function: panner_getFlipPitch
 * -----------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'pitch' angle
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int panner_getFlipPitch(void* const hPan);

/*
 * Function: panner_getFlipRoll
 * ----------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'roll' angle
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int panner_getFlipRoll(void* const hPan);

/*
 * Function: panner_getProcessingDelay
 * -----------------------------------
 * Returns the processing delay in samples. May be used for delay compensation
 * features
 *
 * Returns:
 *     processing delay in samples
 */
int panner_getProcessingDelay(void);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __PANNER_H_INCLUDED__ */
