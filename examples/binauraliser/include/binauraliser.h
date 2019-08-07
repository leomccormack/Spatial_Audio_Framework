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
 * Filename: binauraliser.h (include header)
 * -----------------------------------------
 * Convolves input audio (up to 64 channels) with interpolated HRTFs in the
 * time-frequency domain. The HRTFs are interpolated by applying amplitude-
 * preserving VBAP gains to the HRTF magnitude responses and inter-aural time
 * differences (ITDs) individually, before being re-combined. The example also
 * allows the user to specify an external SOFA file for the convolution.
 *
 * Dependencies:
 *     saf_utilities, saf_hrir, saf_vbap, afSTFTlib
 * Author, date created:
 *     Leo McCormack, 25.09.2017
 */

#ifndef __BINAURALISER_H_INCLUDED__
#define __BINAURALISER_H_INCLUDED__

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
#define ENABLE_ZYLIA_LAB_PRESET
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
#ifdef ENABLE_ZYLIA_LAB_PRESET
    , PRESET_ZYLIA_LAB
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

typedef enum _INTERP_MODES{
    INTERP_TRI = 1  /* Triangular interpolation */
}INTERP_MODES;
    
#define BINAURALISER_MAX_NUM_INPUTS ( 64 )
    
    
/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */
    
/*
 * Function: binauraliser_create
 * -----------------------------
 * Creates an instance of the binauraliser
 *
 * Input Arguments:
 *     phBin - & address of binauraliser handle
 */
void binauraliser_create(void** const phBin);

/*
 * Function: binauraliser_destroy
 * -----------------------------
 * Destroys an instance of the binauraliser
 *
 * Input Arguments:
 *     phBin - & address of binauraliser handle
 */
void binauraliser_destroy(void** const phBin);

/*
 * Function: binauraliser_init
 * ---------------------------
 * Initialises an instance of binauraliser with default settings
 *
 * Input Arguments:
 *     hBin       - binauraliser handle
 *     samplerate - host samplerate.
 */
void binauraliser_init(void* const hBin,
                       int samplerate);

/*
 * Function: binauraliser_process
 * ------------------------------
 * Binauralises the input signals at the user specified directions
 *
 * Input Arguments:
 *     hBin      - binauraliser handle
 *     inputs    - input channel buffers; 2-D array: nInputs x nSamples
 *     outputs   - Output channel buffers; 2-D array: nOutputs x nSamples
 *     nInputs   - number of input channels
 *     nOutputs  - number of output channels
 *     nSamples  - number of samples in 'inputs'/'output' matrices
 *     isPlaying - flag to say if there is audio in the input buffers, 0: no
 *                 audio, reduced processing, 1: audio, full processing
 */
void binauraliser_process(void* const hBin,
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
 * Function: binauraliser_refreshSettings
 * --------------------------------------
 * Sets all intialisation flags to 1. i.e. re-initialise all settings/variables
 * as binauraliser is currently configured, at next available opportunity.
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 */
void binauraliser_refreshSettings(void* const hBin);

/*
 * Function: binauraliser_checkReInit
 * ----------------------------------
 * Check if any reInit Flags are active, and reinitialise if they are.
 * Note: Only call when playback has stopped.
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 */
void binauraliser_checkReInit(void* const hBin);
    
/*
 * Function: binauraliser_setSourceAzi_deg
 * ---------------------------------------
 * Sets the panning azimuth for a specific channel
 *
 * Input Arguments:
 *     hBin       - binauraliser handle
 *     index      - channel index
 *     newAzi_deg - new azimuth in DEGREES
 */
void binauraliser_setSourceAzi_deg(void* const hBin, int index, float newAzi_deg);

/*
 * Function: binauraliser_setSourceElev_deg
 * ----------------------------------------
 * Sets the panning elevation for a specific channel
 *
 * Input Arguments:
 *     hBin       - binauraliser handle
 *     index       - loudspeaker index
 *     newElev_deg - new elevation in DEGREES
 */
void binauraliser_setSourceElev_deg(void* const hBin, int index, float newElev_deg);

/*
 * Function: binauraliser_setNumSources
 * ------------------------------------
 * Sets the number of input channels/sources to binauralise.
 *
 * Input Arguments:
 *     hBin         - binauraliser handle
 *     new_nSources - new number of channels
 */
void binauraliser_setNumSources(void* const hBin, int new_nSources);

/*
 * Function: binauraliser_setUseDefaultHRIRsflag
 * ---------------------------------------------
 * Sets flag to dictate whether the default HRIRs in the Spatial_Audio_Framework
 * should be used, or a custom HRIR set loaded via a SOFA file.
 * Note: if the custom set failes to load correctly, binauraliser will revert to
 * the defualt set. Use 'binauraliser_getUseDefaultHRIRsflag' to check if
 * loading was successful.
 *
 * Input Arguments:
 *     hBin     - binauraliser handle
 *     newState - 0: use custom HRIR set, 1: use default HRIR set
 */
void binauraliser_setUseDefaultHRIRsflag(void* const hBin, int newState);

/*
 * Function: binauraliser_setSofaFilePath
 * --------------------------------------
 * Sets the file path for a .sofa file, in order to employ a custom HRIR set for
 * the decoding.
 * Note: if the custom set failes to load correctly, hcompass will revert to the
 * defualt set. Use 'binauraliser_getUseDefaultHRIRsflag' to check if loading was
 * successful.
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 *     path - file path to .sofa file (WITH file extension)
 */
void binauraliser_setSofaFilePath(void* const hBin, const char* path);

/*
 * Function: binauraliser_setInputConfigPreset
 * -------------------------------------------
 * Sets an input preset. e.g. 5.x, 9.x, 22.x etc.
 *
 * Input Arguments:
 *     hBin        - binauraliser handle
 *     newPresetID - new preset (see "PRESETS" enum)
 */
void binauraliser_setInputConfigPreset(void* const hBin, int newPresetID);

/*
 * Function: binauraliser_setEnableRotation
 * ----------------------------------------
 * Sets the flag to enable/disable rotation.
 *
 * Input Arguments:
 *     hBin     - binauraliser handle
 *     newState - 0: disable, 1: enable
 */
void binauraliser_setEnableRotation(void* const hBin, int newState);

/*
 * Function: binauraliser_setYaw
 * -----------------------------
 * Sets the 'yaw' rotation angle
 *
 * Input Arguments:
 *     hBin   - binauraliser handle
 *     newYaw - the 'yaw' rotation angle, in DEGREES
 */
void binauraliser_setYaw(void* const hBin, float newYaw);

/*
 * Function: binauraliser_setPitch
 * -------------------------------
 * Sets the 'pitch' rotation angle
 *
 * Input Arguments:
 *     hBin     - binauraliser handle
 *     newPitch - the 'pitch' rotation angle, in DEGREES
 */
void binauraliser_setPitch(void* const hBin, float newPitch);

/*
 * Function: binauraliser_setRoll
 * ------------------------------
 * Sets the 'roll' rotation angle
 *
 * Input Arguments:
 *     hBin    - binauraliser handle
 *     newRoll - the 'roll' rotation angle, in DEGREES
 */
void binauraliser_setRoll(void* const hBin, float newRoll);

/*
 * Function: binauraliser_setFlipYaw
 * ---------------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'yaw' angle
 *
 * Input Arguments:
 *     hBin     - binauraliser handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void binauraliser_setFlipYaw(void* const hBin, int newState);

/*
 * Function: binauraliser_setFlipPitch
 * -----------------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'pitch' angle
 *
 * Input Arguments:
 *     hBin     - binauraliser handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void binauraliser_setFlipPitch(void* const hBin, int newState);

/*
 * Function: binauraliser_setFlipRoll
 * ----------------------------------
 * Sets a flag as to whether to "flip" the sign of the current 'roll' angle
 *
 * Input Arguments:
 *     hBin     - binauraliser handle
 *     newState - 0: do not flip sign, 1: flip the sign
 */
void binauraliser_setFlipRoll(void* const hBin, int newState);

/*
 * Function: binauraliser_setRPYflag
 * ---------------------------------
 * Sets a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw" (1)
 * rotation order.
 *
 * Input Arguments:
 *     hBin     - binauraliser handle
 *     newState - 0: use "yaw-pitch-roll", 1: use "roll-pitch-yaw"
 */
void binauraliser_setRPYflag(void* const hBin, int newState);

    /* NOT IMPLEMENTED YET */
void binauraliser_setInterpMode(void* const hBin, int newMode);


/* ========================================================================== */
/*                                Get Functions                               */
/* ========================================================================== */

/*
 * Function: binauraliser_getSourceAzi_deg
 * ---------------------------------------
 * Returns the source azimuth for a given index.
 *
 * Input Arguments:
 *     hBin  - binauraliser handle
 *     index - source index
 * Returns:
 *     source azimuth in DEGREES
 */
float binauraliser_getSourceAzi_deg(void* const hBin, int index);

/*
 * Function: binauraliser_getSourceElev_deg
 * ----------------------------------------
 * Returns the source elevation for a given index.
 *
 * Input Arguments:
 *     hBin  - binauraliser handle
 *     index - source index
 * Returns:
 *     source elevation in DEGREES
 */
float binauraliser_getSourceElev_deg(void* const hBin, int index);

/*
 * Function: binauraliser_getNumSources
 * ------------------------------------
 * Returns the number of inputs/sources in the current layout
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     number of sources
 */
int binauraliser_getNumSources(void* const hBin);

/*
 * Function: binauraliser_getMaxNumSources
 * ---------------------------------------
 * Returns the maximum number of input sources supported by binauraliser
 *
 * Returns:
 *     maximum number of sources
 */
int binauraliser_getMaxNumSources(void);

/*
 * Function: binauraliser_getNumEars
 * ---------------------------------
 * Returns the number of ears possessed by the average homo sapien
 *
 * Returns:
 *     2
 */
int binauraliser_getNumEars(void);

/*
 * Function: binauraliser_getNDirs
 * -------------------------------
 * Returns the number of directions in the currently used HRIR set
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     number of HRIR directions
 */
int binauraliser_getNDirs(void* const hBin);

/*
 * Function: binauraliser_getNTriangles
 * ------------------------------------
 * Returns the number of triangular groupings (faces) returned by the Convex
 * Hull
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     number of faces returned by the Convex Hull of the HRIR grid
 */
int binauraliser_getNTriangles(void* const hBin);

/*
 * Function: binauraliser_getHRIRAzi_deg
 * -------------------------------------
 * Returns the HRIR/HRTF azimuth for a given index.
 *
 * Input Arguments:
 *     hBin  - binauraliser handle
 *     index - HRIR/HRTF index
 * Returns:
 *     HRIR/HRTF azimuth in DEGREES
 */
float binauraliser_getHRIRAzi_deg(void* const hBin, int index);

/*
 * Function: binauraliser_getHRIRElev_deg
 * --------------------------------------
 * Returns the HRIR/HRTF elevation for a given index.
 *
 * Input Arguments:
 *     hBin  - binauraliser handle
 *     index - HRIR/HRTF index
 * Returns:
 *     HRIR/HRTF elevation in DEGREES
 */
float binauraliser_getHRIRElev_deg(void* const hBin, int index);

/*
 * Function: binauraliser_getHRIRlength
 * ------------------------------------
 * Returns the length of HRIRs in time-domain samples
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     HRIR length in samples
 */
int binauraliser_getHRIRlength(void* const hBin);

/*
 * Function: binauraliser_getHRIRsamplerate
 * ----------------------------------------
 * Returns the HRIR sample rate
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     HRIR sampling rate
 */
int binauraliser_getHRIRsamplerate(void* const hBin);

/*
 * Function: binauraliser_getUseDefaultHRIRsflag
 * ---------------------------------------------
 * Returns the value of a flag used to dictate whether the default HRIRs in the
 * Spatial_Audio_Framework should be used, or a custom HRIR set loaded via a
 * SOFA file.
 * Note: if the custom set failes to load correctly, binauraliser will revert to
 * the defualt set.
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     0: use custom HRIR set, 1: use default HRIR set
 */
int binauraliser_getUseDefaultHRIRsflag(void* const hBin);

/*
 * Function: binauraliser_getSofaFilePath
 * --------------------------------------
 * Returns the file path for a .sofa file.
 * Note: if the custom set failes to load correctly, binauraliser will revert to
 * the defualt set. Use 'binauraliser_getUseDefaultHRIRsflag' to check if
 * loading was successful.
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *      file path to .sofa file (WITH file extension)
 */
char* binauraliser_getSofaFilePath(void* const hBin);

/*
 * Function: binauraliser_getDAWsamplerate
 * ---------------------------------------
 * Returns the DAW/Host sample rate
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     DAW/Host sampling rate
 */
int binauraliser_getDAWsamplerate(void* const hBin);

/*
 * Function: binauraliser_getEnableRotation
 * ----------------------------------------
 * Returns the flag value which dictates whether to enable/disable sound-field
 * rotation.
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     0: disabled, 1: enabled
 */
int binauraliser_getEnableRotation(void* const hBin);

/*
 * Function: binauraliser_getYaw
 * -----------------------------
 * Returns the 'yaw' rotation angle
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     the 'yaw' rotation angle, in DEGREES
 */
float binauraliser_getYaw(void* const hBin);

/*
 * Function: binauraliser_getPitch
 * -------------------------------
 * Returns the 'pitch' rotation angle
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     the 'pitch' rotation angle, in DEGREES
 */
float binauraliser_getPitch(void* const hBin);

/*
 * Function: binauraliser_getRoll
 * ------------------------------
 * Returns the 'roll' rotation angle
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     the 'roll' rotation angle, in DEGREES
 */
float binauraliser_getRoll(void* const hBin);

/*
 * Function: binauraliser_getFlipYaw
 * ---------------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'yaw' angle
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int binauraliser_getFlipYaw(void* const hBin);

/*
 * Function: binauraliser_getFlipPitch
 * -----------------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'pitch' angle
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int binauraliser_getFlipPitch(void* const hBin);

/*
 * Function: binauraliser_getFlipRoll
 * ----------------------------------
 * Returns a flag as to whether to "flip" the sign of the current 'roll' angle
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     0: do not flip sign, 1: flip the sign
 */
int binauraliser_getFlipRoll(void* const hBin);

/*
 * Function: binauraliser_getRPYflag
 * ---------------------------------
 * Returns a flag as to whether to use "yaw-pitch-roll" (0) or "roll-pitch-yaw"
 * (1) rotation order.
 *
 * Input Arguments:
 *     hBin - binauraliser handle
 * Returns:
 *     0: use "yaw-pitch-roll", 1: use "roll-pitch-yaw"
 */
int binauraliser_getRPYflag(void* const hBin);

    /* NOT IMPLEMENTED YET */
int binauraliser_getInterpMode(void* const hBin);

/*
 * Function: binauraliser_getProcessingDelay
 * -------------------------------------
 * Returns the processing delay in samples. May be used for delay compensation
 * features
 *
 * Returns:
 *     processing delay in samples
 */
int binauraliser_getProcessingDelay(void);
    

#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __BINAURALISER_H_INCLUDED__ */
