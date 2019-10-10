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
 * Filename: saf_sofa_reader.h (include header)
 * --------------------------------------------
 * A simple sofa reader, which returns only the bare minimum.
 *
 * Dependencies:
 *     netcdf
 * Author, date created:
 *     Leo McCormack, 21.11.2017
 */

#ifdef SAF_ENABLE_SOFA_READER
 
#include "saf_sofa_reader.h"
#include "saf_hrir.h"

void loadSofaFile
(
    char* sofa_filepath,
    float** hrirs,
    float** hrir_dirs_deg,
    int* N_hrir_dirs,
    int* hrir_len,
    int* hrir_fs
)
{
    int i, j, k, retval,dimid[6], *dimids, ndimsp,ncid, varid, is0_360;
    size_t dimlength[6], IR_dims[3], SourcePosition_dims[2];
    size_t hrir_dims[3], hrir_pos_dims[2];
    char dimname[6];
    const char* errorMessage;
    double* IR, *SourcePosition, IR_fs;
    
    /* free any existing memory */
    free1d((void**)&(*hrirs));
    free1d((void**)&(*hrir_dirs_deg));
    
    /* open sofa file */
    /* (returns error value if sofa_filepath==NULL (intentional), or if the file
     * path/name was not found (unintentional). */
    if ((retval = nc_open(sofa_filepath, NC_NOWRITE, &ncid)))
        errorMessage = nc_strerror(retval);
    
    /* if error: */
    if(retval!=NC_NOERR){
        is0_360 = 0;
        /* return default HRIR data */
        (*N_hrir_dirs) = __default_N_hrir_dirs;
        (*hrir_len) = __default_hrir_len;
        (*hrir_fs) = __default_hrir_fs;
        (*hrirs) = malloc1d((*N_hrir_dirs) * 2 * (*hrir_len)*sizeof(float));
        for(i=0; i<(*N_hrir_dirs); i++)
            for(j=0; j<2; j++)
                for(k=0; k< (*hrir_len); k++)
                    (*hrirs)[i*2*(*hrir_len) + j*(*hrir_len) + k] = (float)__default_hrirs[i][j][k];
        (*hrir_dirs_deg) = malloc1d((*N_hrir_dirs) * 2 * sizeof(float));
        for(i=0; i<(*N_hrir_dirs); i++){
            for(j=0; j<2; j++)
                (*hrir_dirs_deg)[i*2+j] = (float)__default_hrir_dirs_deg[i][j];
            if((*hrir_dirs_deg)[2*i+0]>=181.0f)
                is0_360 = 1;
        }
        /* convert to -180..180, if azi is 0..360 */
        if(is0_360)
            for(i=0; i<(*N_hrir_dirs); i++)
                (*hrir_dirs_deg)[2*i+0] = (*hrir_dirs_deg)[2*i+0]>180.0f ? (*hrir_dirs_deg)[2*i+0] -360.0f : (*hrir_dirs_deg)[2*i+0];
        
#ifndef NDEBUG
        /* also output warning message, if encountering this error value was
         * unintentional (i.e. sofa_filepath!=NULL) */
        if(sofa_filepath!=NULL)
            saf_error_print(SAF_WARNING__SOFA_FILE_NOT_FOUND);
#endif
        return;
    }
 
    /* Determine dimension IDs and lengths */
    /* Note: there are 6 possible dimension lengths in the sofa standard */
    for (i=0; i<6; i++){
        retval = nc_inq_dim(ncid, i, &dimname[i], &dimlength[i]);
        errorMessage = nc_strerror(retval);
        retval = nc_inq_dimid(ncid, &dimname[i], &dimid[i]);
        errorMessage = nc_strerror(retval);
    }
    
    /* Extract IR data */
    if ((retval = nc_inq_varid(ncid, "Data.IR", &varid)))
        errorMessage = nc_strerror(retval);
    if ((retval =  nc_inq_varndims (ncid, varid, &ndimsp)))
        errorMessage = nc_strerror(retval);
    dimids = malloc1d(ndimsp*sizeof(int));
    if ((retval = nc_inq_vardimid(ncid, varid, dimids)))
        errorMessage = nc_strerror(retval);
    for(i=0; i<3; i++)
        IR_dims[i] = dimlength[dimid[dimids[i]]];
    free(dimids);
    IR = malloc1d(IR_dims[0]*IR_dims[1]*IR_dims[2]*sizeof(double));
    if ((retval = nc_get_var(ncid, varid, IR)))
        errorMessage = nc_strerror(retval);
    if ((retval = nc_inq_varid(ncid, "Data.SamplingRate", &varid)))
        errorMessage = nc_strerror(retval);
    if ((retval = nc_get_var(ncid, varid, &IR_fs)))
        errorMessage = nc_strerror(retval);
    
    /* Extract positional data */
    if ((retval = nc_inq_varid(ncid, "SourcePosition", &varid)))
        errorMessage = nc_strerror(retval);
    if ((retval =  nc_inq_varndims (ncid, varid, &ndimsp)))
        errorMessage = nc_strerror(retval);
    dimids = malloc1d(ndimsp*sizeof(int));
    if ((retval = nc_inq_vardimid(ncid, varid, dimids)))
        errorMessage = nc_strerror(retval);
    for(i=0; i<2; i++)
        SourcePosition_dims[i] = dimlength[dimid[dimids[i]]];
    free(dimids);
    SourcePosition = malloc1d(SourcePosition_dims[0]*SourcePosition_dims[1]*sizeof(double));
    if ((retval = nc_get_var(ncid, varid, SourcePosition)))
        errorMessage = nc_strerror(retval);
    
    /* Close the file, freeing all resources. */
    if ((retval = nc_close(ncid)))
        errorMessage = nc_strerror(retval);
    
    /* Allocate sufficient memory */
    (*hrir_len) = MIN((int)IR_dims[2], MAX_HRIR_LENGTH); /* truncate the HRIR length (1024 should be plenty) */
    (*hrirs) = malloc1d(IR_dims[0]*IR_dims[1]*(*hrir_len) *sizeof(float));
    (*hrir_dirs_deg) = malloc1d(SourcePosition_dims[0]*2*sizeof(float));
    
    /* Store relevent info in handle */
    (*hrir_fs) = (int)(IR_fs+0.5);
    hrir_dims[0] = IR_dims[0];
    hrir_dims[1] = IR_dims[1];
    hrir_dims[2] = (size_t)(*hrir_len);
    (*N_hrir_dirs) = (int)hrir_dims[0];
    for(i=0; i<2; i++)
        hrir_pos_dims[i] = SourcePosition_dims[i];
    for(i=0; i<hrir_dims[0]; i++){
        for(j=0; j<hrir_dims[1]; j++){
            for(k=0; k<hrir_dims[2]; k++){
                /* trunctate IRs and store in floating point precision */
                (*hrirs)[i*hrir_dims[1]*hrir_dims[2] + j*hrir_dims[2] + k] = (float)IR[i*IR_dims[1]*IR_dims[2] + j*IR_dims[2] + k];
            }
        }
    }
    
    /* store in floating point precision */
    is0_360 = 0;
    for(i=0; i<SourcePosition_dims[0]; i++){
        (*hrir_dirs_deg)[2*i+0] = (float)SourcePosition[i*SourcePosition_dims[1]+0];
        (*hrir_dirs_deg)[2*i+1] = (float)SourcePosition[i*SourcePosition_dims[1]+1];
        if((*hrir_dirs_deg)[2*i+0]>=181.0f)
            is0_360 = 1;
    }
    
    /* convert to -180..180, if azi is 0..360 */
    if(is0_360)
        for(i=0; i<SourcePosition_dims[0]; i++)
            (*hrir_dirs_deg)[2*i+0] = (*hrir_dirs_deg)[2*i+0]>180.0f ? (*hrir_dirs_deg)[2*i+0] -360.0f : (*hrir_dirs_deg)[2*i+0];
            
    free(IR);
    free(SourcePosition);
}

#else
/* To conform to the ISO standard, a '.c' file should not be empty.
 * This typedef is a lie: */
typedef int saf_sofa_reader_is_empty;

#endif /* SAF_ENABLE_SOFA_READER */
