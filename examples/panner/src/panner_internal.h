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
 * Filename: panner_internal.h
 * ---------------------------
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

#ifndef __PANNER_INTERNAL_H_INCLUDED__
#define __PANNER_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "panner.h"
#define SAF_ENABLE_AFSTFT /* for time-frequency transform */
#define SAF_ENABLE_VBAP   /* for VBAP gains */
#define SAF_ENABLE_SH     /* for rotations */
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#define FORCE_3D_LAYOUT /* Even 2D loudspeaker setups will use 3D VBAP, with 2 virtual loudspeakers on the top/bottom */
    
#define HOP_SIZE ( 128 )                            /* STFT hop size = nBands */
#define HYBRID_BANDS ( HOP_SIZE + 5 )               /* hybrid mode incurs an additional 5 bands  */
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )        /* 4/8/16 */
#define MAX_NUM_INPUTS ( PANNER_MAX_NUM_INPUTS )    /* Maximum permited channels for the VST standard */
#define MAX_NUM_OUTPUTS ( PANNER_MAX_NUM_OUTPUTS )  /* Maximum permited channels for the VST standard */
#ifndef DEG2RAD
# define DEG2RAD(x) (x * PI / 180.0f)
#endif
#ifndef RAD2DEG
# define RAD2DEG(x) (x * 180.0f / PI)
#endif
    
    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/*
 * Struct: panner_data
 * -------------------
 * Main structure for panner. Contains variables for audio buffers, afSTFT,
 * internal variables, flags, user parameters
 */
typedef struct _panner
{
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
    float* vbap_gtable; /* N_hrtf_vbap_gtable x nLoudpkrs */
    int N_vbap_gtable;
    float_complex G_src[HYBRID_BANDS][MAX_NUM_INPUTS][MAX_NUM_OUTPUTS];
    
    /* flags */
    int recalc_gainsFLAG[MAX_NUM_INPUTS];
    int reInitGainTables;
    int reInitTFT;
    int recalc_M_rotFLAG;
    
    /* misc. */
    float src_dirs_rot_deg[MAX_NUM_INPUTS][2];
    float src_dirs_rot_xyz[MAX_NUM_INPUTS][3];
    float src_dirs_xyz[MAX_NUM_INPUTS][3]; 
    int nTriangles;
    int output_nDims; /* 2: 2-D, 3: 3-D */
    
    /* pValue */
    float pValue[HYBRID_BANDS];
    
    /* user parameters */
    int nSources, new_nSources;
    float src_dirs_deg[MAX_NUM_INPUTS][2];
    float DTT, spread_deg;
    int nLoudpkrs, new_nLoudpkrs;
    float loudpkrs_dirs_deg[MAX_NUM_OUTPUTS][2];
    float yaw, roll, pitch;                  /* rotation angles in degrees */
    int bFlipYaw, bFlipPitch, bFlipRoll;     /* flag to flip the sign of the individual rotation angles */
    
} panner_data;
     

/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/*
 * panner_initGainTables
 * ---------------------
 * Intialises the VBAP gain table used for panning.
 * Note: call "ambi_dec_initTFT" (if needed) before calling this function
 *
 * Input Arguments:
 *     hPan - panner handle
 */
void panner_initGainTables(void* const hPan);
    
/*
 * panner_initTFT
 * --------------
 * Initialise the filterbank used by panner.
 * Note: Call this function before "panner_initGainTables"
 *
 * Input Arguments:
 *     hPan - panner handle
 */
void panner_initTFT(void* const hPan);
    
/*
 * panner_loadPreset
 * -----------------
 * Loads source/loudspeaker directions from preset
 *
 * Input Arguments:
 *     preset   - see PRESET enum
 * Output Arguments:
 *     dirs_deg - source/loudspeaker directions
 *     newNCH   - & new number of channels
 *     nDims    - & estimate of the number of dimensions (2 or 3)
 */
void panner_loadPreset(PRESETS preset,
                       float dirs_deg[MAX_NUM_INPUTS][2],
                       int* newNCH,
                       int* nDims);


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __PANNER_INTERNAL_H_INCLUDED__ */
