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
 * @file panner_internal.h
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

#ifndef __PANNER_INTERNAL_H_INCLUDED__
#define __PANNER_INTERNAL_H_INCLUDED__

#include "panner.h"        /* Include header for this example */
#include "saf.h"           /* Main include header for SAF */
#include "saf_externals.h" /* To also include SAF dependencies (cblas etc.) */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define FORCE_3D_LAYOUT /**< FLAG: Force 2D loudspeaker setups to also use 3D VBAP (i.e. with 2 virtual loudspeakers on the top/bottom) */
#if !defined(PANNER_FRAME_SIZE)
# if defined(FRAME_SIZE) /* Use the global framesize if it is specified: */
#  define PANNER_FRAME_SIZE ( FRAME_SIZE )          /**< Framesize, in time-domain samples */
# else /* Otherwise, the default framesize for this example is: */
#  define PANNER_FRAME_SIZE ( 128 )                 /**< Framesize, in time-domain samples */
# endif
#endif
#define HOP_SIZE ( 128 )                            /**< STFT hop size */
#define HYBRID_BANDS ( HOP_SIZE + 5 )               /**< Number of frequency bands */
#define TIME_SLOTS ( PANNER_FRAME_SIZE / HOP_SIZE ) /**< Number of STFT timeslots */

/* Checks: */
#if (PANNER_FRAME_SIZE % HOP_SIZE != 0)
# error "PANNER_FRAME_SIZE must be an integer multiple of HOP_SIZE"
#endif
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Main structure for panner. Contains variables for audio buffers, afSTFT,
 * internal variables, flags, user parameters
 */
typedef struct _panner
{
    /* audio buffers */
    float** inputFrameTD;           /**< Input signals, in the time-domain; #MAX_NUM_INPUTS x #PANNER_FRAME_SIZE */
    float** outputFrameTD;          /**< Output signals, in the time-domain; #MAX_NUM_OUTPUTS x #PANNER_FRAME_SIZE */
    float_complex*** inputframeTF;  /**< Input signals, in the time-frequency domain; #HYBRID_BANDS x #MAX_NUM_INPUTS x #TIME_SLOTS */
    float_complex*** outputframeTF; /**< Output signals, in the time-frequency domain; #HYBRID_BANDS x #MAX_NUM_OUTPUTS x #TIME_SLOTS */
    int fs;                         /**< Host sampling rate */
    
    /* time-frequency transform */
    float freqVector[HYBRID_BANDS]; /**< Frequency vector (centre frequencies) */
    void* hSTFT;                    /**< afSTFT handle */
    
    /* Internal */
    int vbapTableRes[2];            /**< [0] azimuth, and [1] elevation grid resolution, in degrees */
    float* vbap_gtable;             /**< Current VBAP gains; FLAT: N_hrtf_vbap_gtable x nLoudpkrs */
    int N_vbap_gtable;              /**< Number of directions in the VBAP gain table */
    float_complex G_src[HYBRID_BANDS][MAX_NUM_INPUTS][MAX_NUM_OUTPUTS];  /**< Current VBAP gains per source */
    
    /* flags */
    CODEC_STATUS codecStatus;       /**< see #CODEC_STATUS */
    PROC_STATUS procStatus;         /**< see #PROC_STATUS */
    float progressBar0_1;           /**< Current (re)initialisation progress, between [0..1] */
    char* progressBarText;          /**< Current (re)initialisation step, string */
    int recalc_gainsFLAG[MAX_NUM_INPUTS]; /**< 1: VBAP gains need to be recalculated for this source, 0: do not */
    int recalc_M_rotFLAG;           /**< 1: recalculate the rotation matrix, 0: do not */
    int reInitGainTables;           /**< 1: reinitialise the VBAP gain table, 0: do not */
    
    /* misc. */
    float src_dirs_rot_deg[MAX_NUM_INPUTS][2]; /**< Intermediate rotated source directions, in degrees */
    float src_dirs_rot_xyz[MAX_NUM_INPUTS][3]; /**< Intermediate rotated source directions, as unit-length Cartesian coordinates */
    float src_dirs_xyz[MAX_NUM_INPUTS][3];     /**< Intermediate source directions, as unit-length Cartesian coordinates */
    int nTriangles;                 /**< Number of loudspeaker triangles */
    int output_nDims;               /**< Dimensionality of the loudspeaker array, 2: 2-D, 3: 3-D */
    int new_nLoudpkrs;              /**< New number of loudspeakers in the array */
    int new_nSources;               /**< New number of inputs/sources */
    
    /* pValue */
    float pValue[HYBRID_BANDS];     /**< Used for the frequency-dependent panning normalisation */
    
    /* user parameters */
    int nSources;                   /**< Current number of inputs/sources */
    float src_dirs_deg[MAX_NUM_INPUTS][2]; /**< Current source directions */
    float DTT;                      /**< Room coefficient [3] */
    float spread_deg;               /**< Source spread/MDAP [2] */
    int nLoudpkrs;                  /**< Current number of loudspeakers in the array */
    float loudpkrs_dirs_deg[MAX_NUM_OUTPUTS][2]; /**< Current loudspeaker directions */
    float yaw;                      /**< yaw (Euler) rotation angle, in degrees */
    float roll;                     /**< roll (Euler) rotation angle, in degrees */
    float pitch;                    /**< pitch (Euler) rotation angle, in degrees */
    int bFlipYaw;                   /**< flag to flip the sign of the yaw rotation angle */
    int bFlipPitch;                 /**< flag to flip the sign of the pitch rotation angle */
    int bFlipRoll;                  /**< flag to flip the sign of the roll rotation angle */
    
} panner_data;
     

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/** Sets codec status (see #CODEC_STATUS enum) */
void panner_setCodecStatus(void* const hPan, CODEC_STATUS newStatus);
    
/**
 * Intialises the VBAP gain table used for panning.
 *
 * @note Call ambi_dec_initTFT() (if needed) before calling this function
 */
void panner_initGainTables(void* const hPan);
    
/**
 * Initialise the filterbank used by panner.
 *
 * @note Call this function before panner_initGainTables()
 */
void panner_initTFT(void* const hPan);
    
/**
 * Loads source directions from preset
 *
 * @param[in]  preset   See #SOURCE_CONFIG_PRESETS enum
 * @param[out] dirs_deg Source/loudspeaker directions
 * @param[out] newNCH   (&) new number of channels
 * @param[out] nDims    (&) estimate of the number of dimensions (2 or 3)
 */
void panner_loadSourcePreset(SOURCE_CONFIG_PRESETS preset,
                             float dirs_deg[MAX_NUM_INPUTS][2],
                             int* newNCH,
                             int* nDims);

/**
 * Loads source/loudspeaker directions from preset
 *
 * @param[in]  preset   See #LOUDSPEAKER_ARRAY_PRESETS enum
 * @param[out] dirs_deg Source/loudspeaker directions
 * @param[out] newNCH   (&) new number of channels
 * @param[out] nDims    (&) estimate of the number of dimensions (2 or 3)
 */
void panner_loadLoudspeakerPreset(LOUDSPEAKER_ARRAY_PRESETS preset,
                                  float dirs_deg[MAX_NUM_INPUTS][2],
                                  int* newNCH,
                                  int* nDims);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __PANNER_INTERNAL_H_INCLUDED__ */
