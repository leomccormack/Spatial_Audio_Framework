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
 *          project! Refer to docs/SOFA_READER_MODULE_DEPENDENCIES.md for
 *          more information.
 *
 * @author Leo McCormack
 * @date 21.11.2017
 */

#ifndef __SAF_SOFA_READER_H_INCLUDED__
#define __SAF_SOFA_READER_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef SAF_ENABLE_SOFA_READER_MODULE

/* ========================================================================== */
/*                          Public Structures/Enums                           */
/* ========================================================================== */

/**
 * SOFA container struct comprising all possible data that can be extracted
 * from SOFA 1.0 files; as laid down in the GeneralFIR and SimpleFreeFieldHRIR
 * specifications:
 *    https://www.sofaconventions.org/mediawiki/index.php/GeneralFIR
 *    https://www.sofaconventions.org/mediawiki/index.php/SimpleFreeFieldHRIR
 */
typedef struct _saf_sofa_container{
    /* All possible SOFA variables (defaults={-1|NULL}) */
    int nSources;                 /**< Number of source/measurement positions */
    int nReceivers;               /**< Number of ears/number of mics etc. */
    int DataLengthIR;             /**< Length of the IRs, in samples */
    float* DataIR;                /**< The impulse response (IR) Data;
                                   * FLAT:nSources x nReceivers x DataLengthIR*/
    float DataSamplingRate;       /**< Sampling rate used to measure the IRs */
    int* DataDelay;               /**< Delay in samples; nReceivers x 1 */
    float* SourcePosition;        /**< Source positions (refer to
                                   *   SourcePositionType & SourcePositionUnits
                                   *   for the convention and units);
                                   *   FLAT: nSources x 3 */
    float* ReceiverPosition;      /**< Receiver positions (refer to
                                   *   ReceiverPositionType &
                                   *   ReceiverPositionUnits for the convention
                                   *   and units);
                                   *   FLAT: nReceivers x 3 */
    int nListeners;               /**< Number of listener positions */
    int nEmitters;                /**< Number of emitter positions */
    float* ListenerPosition;      /**< Listener position (The object
                                   *   incorporating all receivers; refer to
                                   *   ListenerPositionType &
                                   *   ListenerPositionUnits for the convention
                                   *   and units); FLAT: nListeners x 3  */
    float* ListenerUp;            /**< Vector pointing upwards from the listener
                                   *   position (Cartesian); 1 x 3 or FLAT: nListeners x 3  */
    float* ListenerView;          /**< Vector pointing forwards from the
                                   *   listener position (Cartesian); 3 x 1 */
    float* EmitterPosition;       /**< Positions of acoustic excitation used for
                                   *   the measurement (refer to
                                   *   EmitterPositionType &
                                   *   EmitterPositionUnits for the convention
                                   *   and units); FLAT: nEmitters x 3 */

    /* All possible SOFA variable attributes (defaults=NULL) */
    char* ListenerPositionType;   /**< {'cartesian'|'spherical'} */
    char* ListenerPositionUnits;  /**< {'degree, degree, metre'|'metre'} */
    char* ListenerViewType;       /**< {'cartesian'|'spherical'} */
    char* ListenerViewUnits;      /**< {'degree, degree, metre'|'metre'} */
    char* ReceiverPositionType;   /**< {'cartesian'|'spherical'} */
    char* ReceiverPositionUnits;  /**< {'degree, degree, metre'|'metre'} */
    char* SourcePositionType;     /**< {'cartesian'|'spherical'} */
    char* SourcePositionUnits;    /**< {'degree, degree, metre'|'metre'} */
    char* EmitterPositionType;    /**< {'cartesian'|'spherical'} */
    char* EmitterPositionUnits;   /**< {'degree, degree, metre'|'metre'} */
    char* DataSamplingRateUnits;  /**< {'hertz'} */

    /* All possible SOFA global attributes (defaults=NULL) */
    char* Conventions;            /**< {'SOFA'} */
    char* Version;                /**< Version number */
    char* SOFAConventions;        /**< {'GeneralFIR'|'GeneralTF'|
                                   *   'SimpleFreeFieldHRIR'} */
    char* SOFAConventionsVersion; /**< SOFA convention number */
    char* APIName;                /**< API name */
    char* APIVersion;             /**< API version */
    char* ApplicationName;        /**< Name of Application that created file */
    char* ApplicationVersion;     /**< Ver. of Application that created file */
    char* AuthorContact;          /**< Contact information */
    char* Comment;                /**< File comments */
    char* DataType;               /**< {'FIR'|'TF'} */
    char* History;                /**< History information */
    char* License;                /**< License under which file is provided */
    char* Organisation;           /**< Organisation reponsible for the file */
    char* References;             /**< References */
    char* RoomType;               /**< Room type (free field etc.) */
    char* Origin;                 /**< Where this file came from */
    char* DateCreated;            /**< Date file was created */
    char* DateModified;           /**< Date file was modified */
    char* Title;                  /**< Title of file */
    char* DatabaseName;           /**< Name of database this file belongs to */
    char* ListenerShortName;      /**< Name of the listener/dummyhead/mic etc.*/

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
    SAF_SOFA_ERROR_FORMAT_UNEXPECTED,
    /** NetCDF is not thread safe! */
    SAF_SOFA_ERROR_NETCDF_IN_USE

} SAF_SOFA_ERROR_CODES;


/* ========================================================================== */
/*                              Main Functions                                */
/* ========================================================================== */

/**
 * Fills a 'sofa_container' with data found in a SOFA file (GeneralFIR or
 * SimpleFreeFieldHRIR), as detailed in the SOFA 1.0 standard [1,2,3]
 *
 * @warning This loader currently does not support TF SOFA files!
 * @note If you encounter a SOFA file that this SOFA loader cannot load, (or it
 *       misses some of the data) then please send it to the developers :-)
 *
 * @param[in] hSOFA         The sofa_container
 * @param[in] sofa_filepath SOFA file path (including .sofa extension)
 * @returns An error code (see #SAF_SOFA_ERROR_CODES)
 *
 * @see [1] Majdak, P., Iwaya, Y., Carpentier, T., Nicol, R., Parmentier, M.,
 *          Roginska, A., Suzuki, Y., Watanabe, K., Wierstorf, H., Ziegelwanger,
 *          H. and Noisternig, M., 2013, May. Spatially oriented format for
 *          acoustics: A data exchange format representing head-related transfer
 *          functions. In Audio Engineering Society Convention 134. Audio
 *          Engineering Society.
 * @see [2] https://www.sofaconventions.org/mediawiki/index.php/GeneralFIR
 * @see [3] https://www.sofaconventions.org/mediawiki/index.php/SimpleFreeFieldHRIR
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

#endif /* SAF_ENABLE_SOFA_READER_MODULE */


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_SOFA_READER_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup SOFA_Reader */
