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
 *@addtogroup SOFA_Reader
 *@{
 * @file saf_sofa_reader.h 
 * @brief Main header for the sofa reader module (#SAF_SOFA_READER_MODULE)
 *
 * @warning This (optional) SOFA reader, requires netcdf to be linked to your
 *          project!
 *
 * @author Leo McCormack
 * @date 21.11.2017
 */

#ifndef __SAF_SOFA_READER_H_INCLUDED__
#define __SAF_SOFA_READER_H_INCLUDED__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                           Public Structures/Enums                          */
/* ========================================================================== */

/** SOFA container struct, as laid down in the SOFA 1.0 FIR standard:
 *      https://www.sofaconventions.org/mediawiki/index.php/GeneralFIR */
typedef struct _saf_sofa_container{
    /* Possible SOFA/NetCDF variables */
    int nSources;                 /**< Number of source/measurement positions */
    int nReceivers;               /**< Number of ears/number of mics etc. */
    int DataLengthIR;             /**< Length of the IRs, in samples */
    float* DataIR;                /**< The impulse response (IR) Data;
                                   *   FLAT:nSources x nReceivers x DataLengthIR
                                   */
    float DataSamplingRate;       /**< Sampling rate used to measure the IRs */
    int* DataDelay;               /**< Delay in samples; nReceivers x 1 */
    float* SourcePosition;        /**< Source positions [azi,elev,radius]
                                   *   (degrees), or [x,y,z] (meters);
                                   *   FLAT: nSources x 3 */
    float* ReceiverPosition;      /**< Receiver positions [azi,elev,radius]
                                   *   (degrees), or [x,y,z] (meters);
                                   *   FLAT: nReceivers x 3 */
    int nListeners;               /**< Number of listener positions (Always 1)*/
    int nEmitters;                /**< Number of emitter positions */
    float* ListenerPosition;      /**< Listener position (The object
                                   *   incorporating all the receivers),
                                   *   [azi,elev,radius] (degrees), or [x,y,z]
                                   *   (meters); 3 x 1  */
    float* ListenerUp;            /**< Vector pointing upwards from the listener
                                   *   position; 3 x 1  */
    float* ListenerView;          /**< Vector pointing forwards from the listner
                                   *   position; 3 x 1 */
    float* EmitterPosition;       /**< Positions of acoustic excitation used for
                                   *   the measurement; numEmitterPosition x 3
                                   */

    /* Possible SOFA/NetCDF Attributes */
    char* ListenerPositionType;   /**< (default=NULL) */
    char* ListenerPositionUnits;  /**< (default=NULL) */
    char* ReceiverPositionType;   /**< (default=NULL) */
    char* ReceiverPositionUnits;  /**< (default=NULL) */
    char* SourcePositionType;     /**< (default=NULL) */
    char* SourcePositionUnits;    /**< (default=NULL) */
    char* EmitterPositionType;    /**< (default=NULL) */
    char* EmitterPositionUnits;   /**< (default=NULL) */
    char* DataSamplingRateUnits;  /**< (default=NULL) */

    /* Possible SOFA/NetCDF global Attributes */
    char* Conventions;            /**< (default=NULL) */
    char* Version;                /**< (default=NULL) */
    char* SOFAConventions;        /**< (default=NULL) */
    char* SOFAConventionsVersion; /**< (default=NULL) */
    char* APIName;                /**< (default=NULL) */
    char* APIVersion;             /**< (default=NULL) */
    char* ApplicationName;        /**< (default=NULL) */
    char* ApplicationVersion;     /**< (default=NULL) */
    char* AuthorContact;          /**< (default=NULL) */
    char* Comment;                /**< (default=NULL) */
    char* DataType;               /**< (default=NULL) */
    char* History;                /**< (default=NULL) */
    char* License;                /**< (default=NULL) */
    char* Organisation;           /**< (default=NULL) */
    char* References;             /**< (default=NULL) */
    char* RoomType;               /**< (default=NULL) */
    char* Origin;                 /**< (default=NULL) */
    char* DateCreated;            /**< (default=NULL) */
    char* DateModified;           /**< (default=NULL) */
    char* Title;                  /**< (default=NULL) */
    char* DatabaseName;           /**< (default=NULL) */
    char* ListenerShortName;      /**< (default=NULL) */

}saf_sofa_container;

/** SOFA loader error codes */
typedef enum{
    /** None of the error checks failed */
    SAF_SOFA_OK,
    /** Not a SOFA file, or no such file was found in the specified location */
    SAF_SOFA_ERROR_INVALID_FILE_OR_FILE_PATH,
    /** Dimensions of the SOFA data were not as expected */
    SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED,
    /** The data-type of the SOFA data was not as expected */
    SAF_SOFA_ERROR_FORMAT_UNEXPECTED

} SAF_SOFA_ERROR_CODES;


/* ========================================================================== */
/*                              Main Functions                                */
/* ========================================================================== */

/**
 * Fills a sofa_container with the specified SOFA file data
 *
 * @warning Not all SOFA files conform to the proposed standard [1,2]! Certain
 *          accomodatations for off-standard SOFA files (the ones that are
 *          wide-spread) have been built into this SOFA loader. However, if you
 *          encounter a SOFA file that this SOFA loader cannot load, (or it
 *          misses some informtion) then please send it to the developers :-)
 * @warning This loader currently supports only FIR data (not TF data)!
 *
 * @param[in] hSOFA         The sofa_container
 * @param[in] sofa_filepath SOFA file path (including .sofa extension)
 *
 * @see [1] Majdak, P., Iwaya, Y., Carpentier, T., Nicol, R., Parmentier, M.,
 *          Roginska, A., Suzuki, Y., Watanabe, K., Wierstorf, H., Ziegelwanger,
 *          H. and Noisternig, M., 2013, May. Spatially oriented format for
 *          acoustics: A data exchange format representing head-related transfer
 *          functions. In Audio Engineering Society Convention 134. Audio
 *          Engineering Society.
 * @see [2] https://www.sofaconventions.org/mediawiki/index.php/GeneralFIR
 */
SAF_SOFA_ERROR_CODES saf_sofa_open(saf_sofa_container* hSOFA,
                                   char* sofa_filepath);

/**
 * Frees all SOFA data in a sofa_container
 *
 * @param[in] hSOFA The sofa_container
 */
void saf_sofa_close(saf_sofa_container* hSOFA);


/* ========================================================================== */
/*                            Deprecated Functions                            */
/* ========================================================================== */

/**
 * A bare-bones SOFA file reader
 *
 * Allocates memory and copies the values of the essential data contained in a
 * SOFA file to the output arguments.
 *
 * @warning This function is deprecated, use saf_sofa_open().
 * @warning This function assumes the SOFA file comprises HRIR data! (i.e.
 *          not general IR measurement data).
 * @note The hrirs are returned as NULL if the file does not exist.
 *
 * @param[in]  sofa_filepath Directory/file_name of the SOFA file you wish to
 *                           load. Optionally, you may set this as NULL, and the
 *                           function will return the default HRIR data.
 * @param[out] hrirs         (&) the HRIR data;
 *                           FLAT: N_hrir_dirs x #NUM_EARS x hrir_len
 * @param[out] hrir_dirs_deg (&) the HRIR positions; FLAT: N_hrir_dirs x 2
 * @param[out] N_hrir_dirs   (&) number of HRIR positions
 * @param[out] hrir_len      (&) length of the HRIRs, in samples
 * @param[out] hrir_fs       (&) sampling rate of the HRIRs
 */ 
void loadSofaFile(/* Input Arguments */
                  char* sofa_filepath,
                  /* Output Arguments */
                  float** hrirs,
                  float** hrir_dirs_deg,
                  int* N_hrir_dirs,
                  int* hrir_len,
                  int* hrir_fs );


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_SOFA_READER_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup SOFA_Reader */
