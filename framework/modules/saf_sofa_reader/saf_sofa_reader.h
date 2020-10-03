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
 * @note This (optional) SOFA reader, which returns only the bare minimum,
 *       requires netcdf to be linked to the project.
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


typedef struct _saf_sofa_container{
    /* Main data: */
    int nSources;                 /**< Number of measurement positions */
    int nReceivers;               /**< Number of ears/number of mics etc. */
    int DataLengthIR;             /**< Length of IR in samples */
    float* DataIR;                /**< IRs; FLAT: nSources x nReceivers x DataLengthIR */
    float DataSamplingRate;       /**< Sampling rate used to measure IRs */
    float* DataDelay;             /**< Delay in seconds; nReceivers x 1 */
    float* SourcePosition;        /**< */
    float* ReceiverPosition;      /**< */

    /* Also parsing this */ 
    int numListenerPosition;      /**< */
    int numListenerUp;            /**< */
    int numListenerView;          /**< */
    int numEmitterPosition;       /**< */
    float* ListenerPosition;      /**< Listener position [azi,elev,radius] (degrees) */
    float* ListenerUp;            /**< Vector pointing upwards from the listener position */
    float* ListenerView;          /**< Vector pointing forwards from the listner position */
    float* EmitterPosition;       /**< No idea what this could be */

    /* SOFA file Attributes (only if "pullAttributesFLAG" is set to 1) */
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
    char* Organization;           /**< (default=NULL) */
    char* References;             /**< (default=NULL) */
    char* RoomType;               /**< (default=NULL) */
    char* Origin;                 /**< (default=NULL) */
    char* DateCreated;            /**< (default=NULL) */
    char* DateModified;           /**< (default=NULL) */
    char* Title;                  /**< (default=NULL) */
    char* DatabaseName;           /**< (default=NULL) */
    char* ListenerShortName;      /**< (default=NULL) */

}saf_sofa_container;

typedef enum{
    SAF_SOFA_OK,
    SAF_SOFA_ERROR_INVALID_FILE_OR_FILE_PATH,
    SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED,
    SAF_SOFA_ERROR_FORMAT_UNEXPECTED
} SAF_SOFA_ERROR_CODES;


void saf_SOFAcontainer_create(saf_sofa_container** phCon);

SAF_SOFA_ERROR_CODES saf_SOFAcontainer_load(saf_sofa_container* hCon,
                                            char* sofa_filepath,
                                            int pullAttributesFLAG);

void saf_SOFAcontainer_destroy(saf_sofa_container** phCon);

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * A bare-bones SOFA file reader
 *
 * Allocates memory and copies the values of the essential data contained in a
 * SOFA file to the output arguments.
 *
 * @warning This function assumes the SOFA file comprises HRIR data!
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
