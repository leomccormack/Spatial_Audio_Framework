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
 * @note This SOFA reader may optionally use netcdf if "SAF_ENABLE_NETCDF" is
 *       defined. Otherwise, the reader will use libmysofa [1] in conjunction
 *       with zlib, which is included in framework/resources/zlib; i.e. no
 *       external libraries need to be linked by default.
 *
 * @see [1] https://github.com/hoene/libmysofa (BSD-3-Clause license)
 *
 * @author Leo McCormack
 * @date 21.11.2017
 * @license ISC
 */

#ifndef __SAF_SOFA_READER_H_INCLUDED__
#define __SAF_SOFA_READER_H_INCLUDED__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef SAF_ENABLE_SOFA_READER_MODULE

#include "libmysofa/mysofa.h"

/** SOFA file reader options */
typedef enum{
    /** The default option is #SAF_SOFA_READER_OPTION_LIBMYSOFA */
    SAF_SOFA_READER_OPTION_DEFAULT,

    /** This option uses the libmysofa library to load SOFA files, which is
     *  adopted from: https://github.com/hoene/libmysofa (BSD-3-Clause license)
     *
     *  The benefits of this option is that it only depends on zlib, which is
     *  included in SAF. While the downsides of this option, is that zlib has
     *  file size limits for each chunk (<4GB) and it is quite slow at
     *  decompressing large files. */
    SAF_SOFA_READER_OPTION_LIBMYSOFA,

    /** If SAF_ENABLE_NETCDF is defined, then an alternative SOFA reader may be
     *  used. This version requires netcdf to be linked to SAF, along with its
     *  dependencies. The netcdf loader gets around the file size limits of
     *  the libmysofa loader and is also approximately 3 times faster.
     *  Therefore, if you intend to load many large SOFA files
     *  (especially microphone arrays or Ambisonic IRs), then this alternative
     *  SOFA reader is either required (to get around the file size limit) or
     *  may be preferred due to the shorter loading times. The downsides of
     *  using the netcdf option is that it is NOT thread-safe! and requires
     *  these additional external libraries to be linked to SAF. */
    SAF_SOFA_READER_OPTION_NETCDF

} SAF_SOFA_READER_OPTIONS;


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
    float* DataDelay;             /**< Delay in samples; nReceivers x 1 */
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

    /* libmysofa handle, which is used if SAF_ENABLE_NETCDF is not defined */
    void* hLMSOFA;                /**< libmysofa handle */

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
 * @test test__saf_sofa_open(), test__mysofa_load(), test__sofa_comparison()
 *
 * @param[in] hSOFA         The sofa_container
 * @param[in] sofa_filepath SOFA file path (including .sofa extension)
 * @param[in] option        See #SAF_SOFA_READER_OPTIONS
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
                                   char* sofa_filepath,
                                   SAF_SOFA_READER_OPTIONS option);

/**
 * Frees all SOFA data in a sofa_container
 *
 * @param[in] hSOFA The sofa_container
 */
void saf_sofa_close(saf_sofa_container* hSOFA);

#endif /* SAF_ENABLE_SOFA_READER_MODULE */


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __SAF_SOFA_READER_H_INCLUDED__ */

/**@} */ /* doxygen addtogroup SOFA_Reader */
