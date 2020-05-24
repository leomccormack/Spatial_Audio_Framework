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
 */

#ifndef __PANNER_INTERNAL_H_INCLUDED__
#define __PANNER_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "panner.h"
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
    
/* ========================================================================== */
/*                               Internal Enums                               */
/* ========================================================================== */

/**
 * Current status of the processing loop.
 */
typedef enum _PANNER_PROC_STATUS{
    PROC_STATUS_ONGOING = 0, /**< Codec is processing input audio, and should
                              *   not be reinitialised at this time. */
    PROC_STATUS_NOT_ONGOING  /**< Codec is not processing input audio, and may
                              *   be reinitialised if needed. */
}PANNER_PROC_STATUS;
    

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define FORCE_3D_LAYOUT /* Even 2D loudspeaker setups will use 3D VBAP, with 2 virtual loudspeakers on the top/bottom */

#define FRAME_SIZE ( 128 )
#define HOP_SIZE ( 128 )                            /* STFT hop size = nBands */
#define HYBRID_BANDS ( HOP_SIZE + 5 )               /* hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )        /* 4/8/16 */
#define MAX_NUM_INPUTS ( PANNER_MAX_NUM_INPUTS )    /* Maximum permited channels for the VST standard */
#define MAX_NUM_OUTPUTS ( PANNER_MAX_NUM_OUTPUTS )  /* Maximum permited channels for the VST standard */
#ifndef DEG2RAD
# define DEG2RAD(x) (x * SAF_PI / 180.0f)
#endif
#ifndef RAD2DEG
# define RAD2DEG(x) (x * 180.0f / SAF_PI)
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
    /* FIFO buffers */
    int FIFO_idx;
    float inFIFO[MAX_NUM_INPUTS][FRAME_SIZE];
    float outFIFO[MAX_NUM_OUTPUTS][FRAME_SIZE];

    /* audio buffers */
    float inputFrameTD[MAX_NUM_INPUTS][FRAME_SIZE];
    float_complex inputframeTF[HYBRID_BANDS][MAX_NUM_INPUTS][TIME_SLOTS];
    float_complex outputframeTF[HYBRID_BANDS][MAX_NUM_OUTPUTS][TIME_SLOTS];
    complexVector* STFTInputFrameTF;
    complexVector* STFTOutputFrameTF;
    float** tempHopFrameTD;
    int fs;
    
    /* time-frequency transform */
    float freqVector[HYBRID_BANDS];
    void* hSTFT;
    
    /* Internal */
    int vbapTableRes[2];
    float* vbap_gtable; /**< N_hrtf_vbap_gtable x nLoudpkrs */
    int N_vbap_gtable;
    float_complex G_src[HYBRID_BANDS][MAX_NUM_INPUTS][MAX_NUM_OUTPUTS];
    
    /* flags */
    PANNER_CODEC_STATUS codecStatus;
    PANNER_PROC_STATUS procStatus;
    float progressBar0_1;
    char* progressBarText;
    int recalc_gainsFLAG[MAX_NUM_INPUTS];
    int recalc_M_rotFLAG;
    int reInitGainTables;
    
    /* misc. */
    float src_dirs_rot_deg[MAX_NUM_INPUTS][2];
    float src_dirs_rot_xyz[MAX_NUM_INPUTS][3];
    float src_dirs_xyz[MAX_NUM_INPUTS][3]; 
    int nTriangles;
    int output_nDims; /**< 2: 2-D, 3: 3-D */
    
    /* pValue */
    float pValue[HYBRID_BANDS];
    
    /* user parameters */
    int nSources, new_nSources;
    float src_dirs_deg[MAX_NUM_INPUTS][2];
    float DTT, spread_deg;
    int nLoudpkrs, new_nLoudpkrs;
    float loudpkrs_dirs_deg[MAX_NUM_OUTPUTS][2];
    float yaw, roll, pitch;                  /**< rotation angles in degrees */
    int bFlipYaw, bFlipPitch, bFlipRoll;     /**< flag to flip the sign of the individual rotation angles */
    
} panner_data;
     

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Sets codec status (see 'PANNER_CODEC_STATUS' enum)
 */
void panner_setCodecStatus(void* const hPan, PANNER_CODEC_STATUS newStatus);
    
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
 * Loads source/loudspeaker directions from preset
 *
 * @param[in]  preset   See PANNER_PRESET enum
 * @param[out] dirs_deg Source/loudspeaker directions
 * @param[out] newNCH   (&) new number of channels
 * @param[out] nDims    (&) estimate of the number of dimensions (2 or 3)
 */
void panner_loadPreset(PANNER_PRESETS preset,
                       float dirs_deg[MAX_NUM_INPUTS][2],
                       int* newNCH,
                       int* nDims);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __PANNER_INTERNAL_H_INCLUDED__ */
