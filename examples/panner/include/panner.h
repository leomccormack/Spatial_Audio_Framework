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
    
/*
 * Enum: PRESETS
 * -------------
 * Available source/loudspeaker configurations presets
 */
typedef enum _PRESETS{
    PRESET_DEFAULT = 1,
    PRESET_MONO,
    PRESET_STEREO,
    PRESET_5PX,
    PRESET_7PX,
    PRESET_8PX,
    PRESET_9PX,
    PRESET_10PX,
    PRESET_11PX,
    PRESET_11PX_7_4,
    PRESET_13PX,
    PRESET_22PX,
    PRESET_AALTO_MCC,
    PRESET_AALTO_APAJA,
    PRESET_AALTO_LR,
    PRESET_DTU_AVIL,
    PRESET_T_DESIGN_4,
    PRESET_T_DESIGN_12,
    PRESET_T_DESIGN_24,
    PRESET_T_DESIGN_36,
    PRESET_T_DESIGN_48,
    PRESET_T_DESIGN_60
    
}PRESETS;
    
/*
 * Enum: CODEC_STATUS
 * ------------------
 * Current status of the codec.
 *
 * Options:
 *     CODEC_STATUS_INITIALISED     - Codec is initialised and ready to process
 *                                    input audio.
 *     CODEC_STATUS_NOT_INITIALISED - Codec has not yet been initialised, or
 *                                    the codec configuration has changed. Input
 *                                    audio should not be processed.
 *     CODEC_STATUS_INITIALISING    - Codec is currently being initialised,
 *                                    input audio should not be processed.
 */
typedef enum _CODEC_STATUS{
    CODEC_STATUS_INITIALISED = 0,
    CODEC_STATUS_NOT_INITIALISED,
    CODEC_STATUS_INITIALISING
}CODEC_STATUS;
    
#define PANNER_MAX_NUM_INPUTS ( 64 )
#define PANNER_MAX_NUM_OUTPUTS ( 64 )
#define PANNER_SPREAD_MIN_VALUE ( 0.0f )
#define PANNER_SPREAD_MAX_VALUE ( 90.0f )
#define PANNER_PROGRESSBARTEXT_CHAR_LENGTH 256
    

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
 * Function: panner_initCodec
 * ----------------------------
 * Intialises the codec variables, based on current global/user parameters
 *
 * Input Arguments:
 *     hPan - panner handle
 */
void panner_initCodec(void* const hPan);

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
 */
void panner_process(void* const hPan,
                    float** const inputs,
                    float** const outputs,
                    int nInputs,
                    int nOutputs,
                    int nSamples);
    
    
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
 * Function: panner_getCodecStatus
 * -------------------------------
 * Returns current codec status.
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     codec status (see 'CODEC_STATUS' enum)
 */
CODEC_STATUS panner_getCodecStatus(void* const hPan);

/*
 * Function: panner_getProgressBar0_1
 * ----------------------------------
 * (Optional) Returns current intialisation/processing progress, between 0..1
 * 0: intialisation/processing has started
 * 1: intialisation/processing has ended
 *
 * Input Arguments:
 *     hPan - panner handle
 * Returns:
 *     current progress, 0..1
 */
float panner_getProgressBar0_1(void* const hPan);

/*
 * Function: panner_getProgressBarText
 * -----------------------------------
 * (Optional) Returns current intialisation/processing progress text
 * Note: "text" string should be (at least) of length:
 *     PANNER_PROGRESSBARTEXT_CHAR_LENGTH
 *
 * Input Arguments:
 *     hPan - panner handle
 * Output Arguments:
 *     text - process bar text; PANNER_PROGRESSBARTEXT_CHAR_LENGTH x 1
 */
void panner_getProgressBarText(void* const hPan, char* text);
    
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
