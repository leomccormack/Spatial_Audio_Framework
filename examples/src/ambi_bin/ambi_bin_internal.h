/*
 * Copyright 2018 Leo McCormack
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
 * @file ambi_bin_internal.h
 * @brief A binaural Ambisonic decoder for reproducing ambisonic signals over
 *        headphones
 *
 * The decoder includes many historic and current state-of-the-art decoding
 * approaches. It also supports sound-field rotation for head-tracking and may
 * also accomodate custom HRIR sets via the SOFA standard.
 *
 * @author Leo McCormack
 * @date 14.04.2018
 */

#ifndef __AMBI_BIN_INTERNAL_H_INCLUDED__
#define __AMBI_BIN_INTERNAL_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ambi_bin.h"   
#include "saf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                            Internal Parameters                             */
/* ========================================================================== */

#ifndef FRAME_SIZE
# define FRAME_SIZE ( 128 )
#endif
#define HOP_SIZE ( 128 ) /* STFT hop size */
#define HYBRID_BANDS ( 133 )
#define TIME_SLOTS ( FRAME_SIZE / HOP_SIZE )  
#define POST_GAIN ( -9.0f )   /* dB */
#ifndef DEG2RAD
# define DEG2RAD(x) (x * M_PI / 180.0f)
#endif
#ifndef RAD2DEG
# define RAD2DEG(x) (x * 180.0f / M_PI)
#endif

    
/* ========================================================================== */
/*                                 Structures                                 */
/* ========================================================================== */

/**
 * Contains variables for sofa file loading, HRIRs, and the binaural decoder.
 */
typedef struct _ambi_bin_codecPars
{
    /* Decoder */
    float_complex M_dec[HYBRID_BANDS][NUM_EARS][MAX_NUM_SH_SIGNALS];
    
    /* sofa file info */
    char* sofa_filepath;    /**< absolute/relevative file path for a sofa file */
    float* hrirs;           /**< time domain HRIRs; FLAT: N_hrir_dirs x 2 x hrir_len */
    float* hrir_dirs_deg;   /**< directions of the HRIRs in degrees [azi elev]; FLAT: N_hrir_dirs x 2 */
    int N_hrir_dirs;        /**< number of HRIR directions in the current sofa file */
    int hrir_len;           /**< length of the HRIRs, this can be truncated, see "saf_sofa_reader.h" */
    int hrir_fs;            /**< sampling rate of the HRIRs, should ideally match the host sampling rate, although not required */
   
    /* hrtf filterbank coefficients */
    float* itds_s;          /**< interaural-time differences for each HRIR (in seconds); N_hrirs x 1 */
    float_complex* hrtf_fb; /**< HRTF filterbank coeffs; FLAT: nBands x nCH x N_hrirs */
    
}ambi_bin_codecPars;
    
/**
 * Main structure for ambi_bin. Contains variables for audio buffers, afSTFT,
 * rotation matrices, internal variables, flags, user parameters
 */
typedef struct _ambi_bin
{
    /* audio buffers + afSTFT time-frequency transform handle */
    int fs;                         /**< host sampling rate */ 
    float SHFrameTD[MAX_NUM_SH_SIGNALS][FRAME_SIZE];
    float_complex SHframeTF[HYBRID_BANDS][MAX_NUM_SH_SIGNALS][TIME_SLOTS];
    float_complex SHframeTF_rot[HYBRID_BANDS][MAX_NUM_SH_SIGNALS][TIME_SLOTS];
    float_complex binframeTF[HYBRID_BANDS][NUM_EARS][TIME_SLOTS];
    complexVector* STFTInputFrameTF;
    complexVector* STFTOutputFrameTF;
    void* hSTFT;                    /**< afSTFT handle */
    int afSTFTdelay;                /**< for host delay compensation */
    float** tempHopFrameTD;         /**< temporary multi-channel time-domain buffer of size "HOP_SIZE". */
    float freqVector[HYBRID_BANDS]; /**< frequency vector for time-frequency transform, in Hz */
     
    /* our codec configuration */
    CODEC_STATUS codecStatus;
    float progressBar0_1;
    char* progressBarText;
    ambi_bin_codecPars* pars;
    
    /* internal variables */
    PROC_STATUS procStatus;
    float_complex M_rot[MAX_NUM_SH_SIGNALS][MAX_NUM_SH_SIGNALS]; 
    int new_order;                  /**< new decoding order */
    int nSH;                        /**< number of spherical harmonic signals */
    
    /* flags */ 
    int recalc_M_rotFLAG;           /**< 0: no init required, 1: init required */
    int reinit_hrtfsFLAG;           /**< 0: no init required, 1: init required */
    
    /* user parameters */
    int order;                      /**< current decoding order */
    int enableMaxRE;                /**< 0: disabled, 1: enabled */
    int enableDiffuseMatching;      /**< 0: disabled, 1: enabled */
    int enablePhaseWarping;         /**< 0: disabled, 1: enabled */
    AMBI_BIN_DECODING_METHODS method; /* current decoding method */
    float EQ[HYBRID_BANDS];         /**< EQ curve */
    int useDefaultHRIRsFLAG;        /**< 1: use default HRIRs in database, 0: use those from SOFA file */
    CH_ORDER chOrdering;
    NORM_TYPES norm;
    int enableRotation;
    float yaw, roll, pitch;         /**< rotation angles in degrees */
    int bFlipYaw, bFlipPitch, bFlipRoll; /**< flag to flip the sign of the individual rotation angles */
    int useRollPitchYawFlag;        /**< rotation order flag, 1: r-p-y, 0: y-p-r */
    
} ambi_bin_data;


/* ========================================================================== */
/*                             Internal Functions                             */
/* ========================================================================== */

/**
 * Sets codec status. 
 */
void ambi_bin_setCodecStatus(void* const hAmbi,
                             CODEC_STATUS newStatus);


#ifdef __cplusplus
} /* extern "C" { */
#endif /* __cplusplus */

#endif /* __AMBI_BIN_INTERNAL_H_INCLUDED__ */
