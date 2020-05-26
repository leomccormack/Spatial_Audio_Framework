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
 * @file saf_sofa_reader.h 
 * @brief A simple sofa reader, which returns only the bare minimum
 *
 * @note This (optional) SOFA reader, requires netcdf to be linked to the
 *       project.
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
#ifdef SAF_ENABLE_SOFA_READER
/* If your compiler stopped at this point, then please add the path to the
 * netcdf include file to your project include header paths.
 * Instructions for linking the required "netcdf" library may also be found
 * here: https://github.com/leomccormack/Spatial_Audio_Framework */
# include <netcdf.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* ========================================================================== */
/*                               Main Functions                               */
/* ========================================================================== */

/**
 * A bare-bones SOFA file reader
 *
 * Allocates memory and copies the values of the essential data contained in a
 * SOFA file to the output arguments.
 *
 * @note The hrirs are returned as NULL if the file does not exist.
 *
 * @param[in]  sofa_filepath Directory/file_name of the SOFA file you wish to
 *                           load. Optionally, you may set this as NULL, and the
 *                           function will return the default HRIR data.
 * @param[out] hrirs         (&) the HRIR data;
 *                           FLAT: N_hrir_dirs x 2 x hrir_len
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
