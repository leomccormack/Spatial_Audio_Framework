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
 * @file saf_sofa_reader.c
 * @ingroup SOFA_Reader
 * @brief Public source for the sofa reader module (#SAF_SOFA_READER_MODULE)
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

#include "saf_sofa_reader.h"
#include "../saf_utilities/saf_utilities.h"
#include "saf_externals.h"

#ifdef SAF_ENABLE_SOFA_READER_MODULE

/* ========================================================================== */
/*                              Main Functions                                */
/* ========================================================================== */

SAF_SOFA_ERROR_CODES saf_sofa_open
(
    saf_sofa_container* h,
    char* sofa_filepath,
    SAF_SOFA_READER_OPTIONS option
)
{
    /* libmysofa reader vars */
    int err;
    MYSOFA_HRTF *hrtf;
    MYSOFA_ATTRIBUTE* tmp_a;
#ifdef SAF_ENABLE_NETCDF
    /* NetCDF reader vars */
    int varid, attnum, i, j, ncid, retval, ndimsp, nvarsp, nattsp, unlimdimidp, varnattsp;
    size_t tmp_size, lenp;
    char varname[NC_MAX_NAME+1], attname[NC_MAX_NAME+1];
    char* dimname;
    double* tmp_data;
    int* dimids, *dimid;
    nc_type typep;
    size_t* dimlength;
#endif /* SAF_ENABLE_NETCDF */

    /* Default variables */
    h->nSources = h->nReceivers = h->DataLengthIR = -1;
    h->DataSamplingRate = 0.0f;
    h->nEmitters = h->nListeners = -1;
    h->DataIR = h->SourcePosition = h->ReceiverPosition = h->ListenerPosition =
    h->ListenerUp = h->ListenerView = h->EmitterPosition = NULL;
    h->DataDelay = NULL;

    /* Default variable attributes */
    h->ListenerPositionType = h->ListenerPositionUnits = h->ReceiverPositionType
    = h->ReceiverPositionUnits = h->SourcePositionType = h->SourcePositionUnits
    = h->EmitterPositionType = h->EmitterPositionUnits = h->DataSamplingRateUnits
    = h->ListenerViewType = h->ListenerViewUnits = NULL;

    /* Default global attributes */
    h->Conventions = h->Version = h->SOFAConventions = h->SOFAConventionsVersion
    = h->APIName = h->APIVersion = h->ApplicationName = h->ApplicationVersion
    = h->AuthorContact = h->Comment = h->DataType = h->History = h->License
    = h->Organisation = h->References = h->RoomType = h->Origin = h->DateCreated
    = h->DateModified = h->Title = h->DatabaseName = h->ListenerShortName = NULL;

    /* Read the SOFA file */
    switch(option){
        case SAF_SOFA_READER_OPTION_DEFAULT: /* fall through */
        case SAF_SOFA_READER_OPTION_LIBMYSOFA:
            /* Load SOFA file using the libmysofa library: */
            hrtf = mysofa_load(sofa_filepath, &err);
            h->hLMSOFA = (void*)hrtf;
            switch(err){
                case MYSOFA_OK:
                    /* Copy variables and pointers to data: */
                    h->nSources = hrtf->M;
                    h->nReceivers = hrtf->R;
                    h->DataLengthIR = hrtf->N;
                    h->DataSamplingRate = hrtf->DataSamplingRate.values[0];
                    h->nEmitters = hrtf->E;
                    h->nListeners = hrtf->M; // changed to M for multiple listeners
                    h->DataIR = hrtf->DataIR.values;
                    h->DataDelay = hrtf->DataDelay.values;
                    h->SourcePosition = hrtf->SourcePosition.values;
                    h->ReceiverPosition = hrtf->ReceiverPosition.values;
                    h->ListenerPosition = hrtf->ListenerPosition.values;
                    h->ListenerUp = hrtf->ListenerUp.values;
                    h->ListenerView = hrtf->ListenerView.values;
                    h->EmitterPosition = hrtf->EmitterPosition.values;

                    /* Variable attributes */
                    tmp_a = hrtf->ListenerPosition.attributes;
                    while (tmp_a) {
                        if (!strcmp((char*)tmp_a->name,"Type"))
                            h->ListenerPositionType = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"Units"))
                            h->ListenerPositionUnits = tmp_a->value;
                        tmp_a = tmp_a->next;
                    }
                    tmp_a = hrtf->ReceiverPosition.attributes;
                    while (tmp_a) {
                        if (!strcmp((char*)tmp_a->name,"Type"))
                            h->ReceiverPositionType = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"Units"))
                            h->ReceiverPositionUnits = tmp_a->value;
                        tmp_a = tmp_a->next;
                    }
                    tmp_a = hrtf->SourcePosition.attributes;
                    while (tmp_a) {
                        if (!strcmp((char*)tmp_a->name,"Type"))
                            h->SourcePositionType = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"Units"))
                            h->SourcePositionUnits = tmp_a->value;
                        tmp_a = tmp_a->next;
                    }
                    tmp_a = hrtf->EmitterPosition.attributes;
                    while (tmp_a) {
                        if (!strcmp((char*)tmp_a->name,"Type"))
                            h->EmitterPositionType = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"Units"))
                            h->EmitterPositionUnits = tmp_a->value;
                        tmp_a = tmp_a->next;
                    }
                    tmp_a = hrtf->ListenerView.attributes;
                    while (tmp_a) {
                        if (!strcmp((char*)tmp_a->name,"Type"))
                            h->ListenerViewType = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"Units"))
                            h->ListenerViewUnits = tmp_a->value;
                        tmp_a = tmp_a->next;
                    }
                    tmp_a = hrtf->DataSamplingRate.attributes;
                    while (tmp_a) {
                        if (!strcmp((char*)tmp_a->name,"Units"))
                            h->DataSamplingRateUnits = tmp_a->value;
                        tmp_a = tmp_a->next;
                    }

                    /* Global attributes */
                    tmp_a = hrtf->attributes;
                    while (tmp_a) {
                        if (!strcmp((char*)tmp_a->name,"Conventions"))
                            h->Conventions = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"Version"))
                            h->Version = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"SOFAConventions"))
                            h->SOFAConventions = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"SOFAConventionsVersion"))
                            h->SOFAConventionsVersion = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"APIName"))
                            h->APIName = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"APIVersion"))
                            h->APIVersion = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"ApplicationName"))
                            h->ApplicationName = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"ApplicationVersion"))
                            h->ApplicationVersion = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"AuthorContact"))
                            h->AuthorContact = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"Comment"))
                            h->Comment = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"DataType"))
                            h->DataType = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"History"))
                            h->History = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"License"))
                            h->License = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"Organization"))
                            h->Organisation = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"References"))
                            h->References = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"RoomType"))
                            h->RoomType = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"Origin"))
                            h->Origin = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"DateCreated"))
                            h->DateCreated = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"DateModified"))
                            h->DateModified = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"Title"))
                            h->Title = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"DatabaseName"))
                            h->DatabaseName = tmp_a->value;
                        else if (!strcmp((char*)tmp_a->name,"ListenerShortName"))
                            h->ListenerShortName = tmp_a->value;
                        tmp_a = tmp_a->next;
                    }
                    break;
                case MYSOFA_READ_ERROR:
                    return SAF_SOFA_ERROR_INVALID_FILE_OR_FILE_PATH;
                case MYSOFA_INVALID_DIMENSIONS:
                    return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED;
                default:
                    return SAF_SOFA_ERROR_FORMAT_UNEXPECTED;
            }
            break;

        case SAF_SOFA_READER_OPTION_NETCDF:
#ifdef SAF_ENABLE_NETCDF
            h->hLMSOFA = NULL; /* Not used */

            /* Open NetCDF file */
            if ((retval = nc_open(sofa_filepath, NC_NOWRITE, &ncid)))
                nc_strerror(retval);
            if(retval!=NC_NOERR)/* if error: */
                return SAF_SOFA_ERROR_INVALID_FILE_OR_FILE_PATH;
            nc_inq(ncid, &ndimsp, &nvarsp, &nattsp, &unlimdimidp); /* find number of possible dimensions, variables, and attributes */

            /* Find dimension IDs and lengths */
            dimid = malloc1d(ndimsp*sizeof(int));
            dimlength = malloc1d(ndimsp*sizeof(size_t));
            dimname = malloc1d(ndimsp*(NC_MAX_NAME+1)*sizeof(char)); /* +1 as NC_MAX_NAME does not include null-termination */
            for (i=0; i<ndimsp; i++){
                nc_inq_dim(ncid, i, &dimname[i*(NC_MAX_NAME+1)], &dimlength[i]);
                nc_inq_dimid(ncid, &dimname[i*(NC_MAX_NAME+1)], &dimid[i]);
            }

            /* Loop over the variables and pull the data accordingly */
            dimids = NULL;
            tmp_data = NULL;
            for(varid=0; varid<nvarsp; varid++){
                nc_inq_var(ncid, varid, (char*)varname, NULL, NULL, NULL, NULL); /* Variable name */
                nc_inq_varndims(ncid, varid, &ndimsp);                           /* Variable dimensionality */
                dimids = realloc1d(dimids, ndimsp*sizeof(int));
                nc_inq_vardimid(ncid, varid, dimids);                            /* Variable dimension IDs */
                nc_inq_vartype(ncid, varid, &typep);                             /* Variable data type */
                nc_inq_varnatts(ncid, varid, &varnattsp);                        /* Variable number of associated attributes */

                if (!strcmp((char*)varname,"Data.IR")){
                    /* Checks */
                    if(ndimsp!=3) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if(h->nReceivers!=-1 && (int)dimlength[dimid[dimids[1]]] != h->nReceivers) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if(typep!=NC_DOUBLE) { return SAF_SOFA_ERROR_FORMAT_UNEXPECTED; }

                    /* Pull data */
                    h->nSources =     (int)dimlength[dimid[dimids[0]]];
                    h->nReceivers =   (int)dimlength[dimid[dimids[1]]];
                    h->DataLengthIR = (int)dimlength[dimid[dimids[2]]];
                    tmp_size = dimlength[dimid[dimids[0]]] * dimlength[dimid[dimids[1]]] * dimlength[dimid[dimids[2]]];
                    tmp_data = realloc1d(tmp_data, tmp_size*sizeof(double));
                    nc_get_var(ncid, varid, tmp_data);
                    h->DataIR = malloc1d(tmp_size*sizeof(float));
                    for(i=0; i<(int)tmp_size; i++)
                        h->DataIR[i] = (float)tmp_data[i];
                }
                else if(!strcmp((char*)varname,"Data.SamplingRate")){
                    /* Checks */
                    if(!(ndimsp==1 && dimlength[dimid[dimids[0]]] == 1)) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if(typep!=NC_DOUBLE) { return SAF_SOFA_ERROR_FORMAT_UNEXPECTED; }

                    /* Pull data */
                    tmp_data = realloc1d(tmp_data, 1*sizeof(double));
                    nc_get_var(ncid, varid, tmp_data);
                    h->DataSamplingRate = (float)tmp_data[0];

                    /* Pull associated attributes */
                    for(attnum=0; attnum<varnattsp; attnum++){
                        nc_inq_attname(ncid, varid, attnum, attname);
                        nc_inq_attlen(ncid, varid, attname, &lenp);
                        if (!strcmp((char*)attname,"Units")){
                            h->DataSamplingRateUnits = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                            nc_get_att(ncid, varid, attname, h->DataSamplingRateUnits);
                        }
                    }
                }
                else if (!strcmp((char*)varname,"Data.Delay")){
                    /* Checks */
                    if(!(ndimsp==2 || ndimsp==3)) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if(h->nReceivers!=-1 && !((int)dimlength[dimid[dimids[1]]] == h->nReceivers || (int)dimlength[dimid[dimids[0]]] == h->nReceivers)) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if((int)dimlength[dimid[dimids[0]]] != 1 && (int)dimlength[dimid[dimids[1]]] != 1) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if(typep!=NC_DOUBLE) { return SAF_SOFA_ERROR_FORMAT_UNEXPECTED; }

                    /* Pull data */
                    tmp_size = dimlength[dimid[dimids[0]]] * dimlength[dimid[dimids[1]]];
                    tmp_data = realloc1d(tmp_data, tmp_size*sizeof(double));
                    nc_get_var(ncid, varid, tmp_data);
                    h->DataDelay = malloc1d(tmp_size*sizeof(int));
                    for(i=0; i<(int)tmp_size; i++)
                        h->DataDelay[i] = (float)tmp_data[i];
                }
                else if (!strcmp((char*)varname,"SourcePosition")){
                    /* Checks */
                    if(ndimsp!=2) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if(h->nSources!=-1 && (int)dimlength[dimid[dimids[0]]] != h->nSources) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if((int)dimlength[dimid[dimids[1]]] != 3) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if(typep!=NC_DOUBLE) { return SAF_SOFA_ERROR_FORMAT_UNEXPECTED; }

                    /* Pull data */
                    tmp_size = dimlength[dimid[dimids[0]]] * dimlength[dimid[dimids[1]]];
                    h->nSources = (int)dimlength[dimid[dimids[0]]];
                    tmp_data = realloc1d(tmp_data, tmp_size*sizeof(double));
                    nc_get_var(ncid, varid, tmp_data);
                    h->SourcePosition = malloc1d(tmp_size*sizeof(float));
                    for(i=0; i<(int)tmp_size; i++)
                        h->SourcePosition[i] = (float)tmp_data[i];

                    /* Pull associated attributes */
                    for(attnum=0; attnum<varnattsp; attnum++){
                        nc_inq_attname(ncid, varid, attnum, attname);
                        nc_inq_attlen(ncid, varid, attname, &lenp);
                        if (!strcmp((char*)attname,"Type")){
                            h->SourcePositionType = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                            nc_get_att(ncid, varid, attname, h->SourcePositionType);
                        }
                        else if (!strcmp((char*)attname,"Units")){
                            h->SourcePositionUnits = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                            nc_get_att(ncid, varid, attname, h->SourcePositionUnits);
                        }
                    }
                }
                else if (!strcmp((char*)varname,"ReceiverPosition")){
                    switch(ndimsp){
                        /* Many SOFA files have the "ReceiverPosition" variable with the following dimensions: nReceivers x 3  */
                        case 2:
                            /* Checks */
                            if(h->nReceivers!=-1 && (int)dimlength[dimid[dimids[0]]] != h->nReceivers) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                            if((int)dimlength[dimid[dimids[1]]] != 3) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                            if(typep!=NC_DOUBLE) { return SAF_SOFA_ERROR_FORMAT_UNEXPECTED; }

                            /* Pull data */
                            h->nReceivers = (int)dimlength[dimid[dimids[0]]];
                            tmp_data = realloc1d(tmp_data, h->nReceivers*3*sizeof(double));
                            nc_get_var(ncid, varid, tmp_data);
                            h->ReceiverPosition = malloc1d(h->nReceivers*3*sizeof(float));
                            for(i=0; i<h->nReceivers*3; i++)
                                    h->ReceiverPosition[i] = (float)tmp_data[i];
                            break;

                        /* Some SOFA files have the "ReceiverPosition" variable with the following dimensions: 1 x nReceivers x 3
                         * This is the reason for the switch case found here {2,3}, as it is not fully understood if this '1' is
                         * for the number of emmiters or listeners? - Until this is made clear, the
                         * following code will just pull the first one (i.e. nReceivers x 3). Therefore, if you have an example
                         * of a file that has "ReceiverPosition" with the dimensions: N x nReceivers x 3  (where N>1)
                         * then please send it to the developers :-) */
                        case 3:
                            /* Checks */
                            if(h->nReceivers!=-1 && (int)dimlength[dimid[dimids[0]]] != h->nReceivers) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                            if((int)dimlength[dimid[dimids[1]]] != 3) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                            if(typep!=NC_DOUBLE) { return SAF_SOFA_ERROR_FORMAT_UNEXPECTED; }

                            /* Pull data */
                            h->nReceivers = (int)dimlength[dimid[dimids[0]]];
                            tmp_data = realloc1d(tmp_data, h->nReceivers*3*sizeof(double));
                            nc_get_var(ncid, varid, tmp_data);
                            h->ReceiverPosition = malloc1d(h->nReceivers*3*sizeof(float));
                            for(i=0; i<h->nReceivers*3; i++)
                                h->ReceiverPosition[i] = (float)tmp_data[i];
                            break;
                        default:
                            return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED;
                    }

                    /* Pull associated attributes */
                    for(attnum=0; attnum<varnattsp; attnum++){
                        nc_inq_attname(ncid, varid, attnum, attname);
                        nc_inq_attlen(ncid, varid, attname, &lenp);
                        if (!strcmp((char*)attname,"Type")){
                            h->ReceiverPositionType = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                            nc_get_att(ncid, varid, attname, h->ReceiverPositionType);
                        }
                        else if (!strcmp((char*)attname,"Units")){
                            h->ReceiverPositionUnits = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                            nc_get_att(ncid, varid, attname, h->ReceiverPositionUnits);
                        }
                    }
                }
                else if (!strcmp((char*)varname,"ListenerPosition")){
                    /* Checks */
                    if(ndimsp!=2) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if((int)dimlength[dimid[dimids[1]]] != 3 && (int)dimlength[dimid[dimids[0]]] != 3) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    //if((int)dimlength[dimid[dimids[1]]] != 1 && (int)dimlength[dimid[dimids[0]]] != 1) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if(typep!=NC_DOUBLE) { return SAF_SOFA_ERROR_FORMAT_UNEXPECTED; }

                    /* Pull data */
                    tmp_size = dimlength[dimid[dimids[0]]] * dimlength[dimid[dimids[1]]];
                    h->nListeners = (int)dimlength[dimid[dimids[0]]];
                    tmp_data = realloc1d(tmp_data, tmp_size*sizeof(double));
                    nc_get_var(ncid, varid, tmp_data);
                    h->ListenerPosition = malloc1d(tmp_size*sizeof(float));
                    for(j=0; j<(int)tmp_size; j++)
                        h->ListenerPosition[j] = (float)tmp_data[j];

                    /* Pull associated attributes */
                    for(attnum=0; attnum<varnattsp; attnum++){
                        nc_inq_attname(ncid, varid, attnum, attname);
                        nc_inq_attlen(ncid, varid, attname, &lenp);
                        if (!strcmp((char*)attname,"Type")){
                            h->ListenerPositionType = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                            nc_get_att(ncid, varid, attname, h->ListenerPositionType);
                        }
                        else if (!strcmp((char*)attname,"Units")){
                            h->ListenerPositionUnits = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                            nc_get_att(ncid, varid, attname, h->ListenerPositionUnits);
                        }
                    }
                }
                else if (!strcmp((char*)varname,"ListenerUp")){
                    /* Checks */
                    if(ndimsp!=2) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if((int)dimlength[dimid[dimids[1]]] != 3 && (int)dimlength[dimid[dimids[0]]] != 3) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if((int)dimlength[dimid[dimids[1]]] != 1 && (int)dimlength[dimid[dimids[0]]] != 1) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if(typep!=NC_DOUBLE) { return SAF_SOFA_ERROR_FORMAT_UNEXPECTED; }

                    /*  Pull data */
                    tmp_size = dimlength[dimid[dimids[0]]] * dimlength[dimid[dimids[1]]];
                    tmp_data = realloc1d(tmp_data, tmp_size*sizeof(double));
                    nc_get_var(ncid, varid, tmp_data);
                    h->ListenerUp = malloc1d(tmp_size*sizeof(float));
                    for(j=0; j<(int)tmp_size; j++)
                        h->ListenerUp[j] = (float)tmp_data[j];
                }
                else if (!strcmp((char*)varname,"ListenerView")){
                    /* Checks */
                    if(ndimsp!=2) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if((int)dimlength[dimid[dimids[1]]] != 3 && (int)dimlength[dimid[dimids[0]]] != 3) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if((int)dimlength[dimid[dimids[1]]] != 1 && (int)dimlength[dimid[dimids[0]]] != 1) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                    if(typep!=NC_DOUBLE) { return SAF_SOFA_ERROR_FORMAT_UNEXPECTED; }

                    /* Pull data */
                    tmp_size = dimlength[dimid[dimids[0]]] * dimlength[dimid[dimids[1]]];
                    tmp_data = realloc1d(tmp_data, tmp_size*sizeof(double));
                    nc_get_var(ncid, varid, tmp_data);
                    h->ListenerView = malloc1d(tmp_size*sizeof(float));
                    for(j=0; j<(int)tmp_size; j++)
                        h->ListenerView[j] = (float)tmp_data[j];

                    /* Pull associated attributes */
                    for(attnum=0; attnum<varnattsp; attnum++){
                        nc_inq_attname(ncid, varid, attnum, attname);
                        nc_inq_attlen(ncid, varid, attname, &lenp);
                        if (!strcmp((char*)attname,"Type")){
                            h->ListenerViewType = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                            nc_get_att(ncid, varid, attname, h->ListenerViewType);
                        }
                        else if (!strcmp((char*)attname,"Units")){
                            h->ListenerViewUnits = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                            nc_get_att(ncid, varid, attname, h->ListenerViewUnits);
                        }
                    }
                }
                else if (!strcmp((char*)varname,"EmitterPosition")){
                    switch(ndimsp){
                        /* Many SOFA files have the "EmitterPosition" variable with the following dimensions: nEmitters x 3  */
                        case 2:
                            /* Checks */
                            if((int)dimlength[dimid[dimids[1]]] != 3 && (int)dimlength[dimid[dimids[0]]] != 3) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                            if(typep!=NC_DOUBLE) { return SAF_SOFA_ERROR_FORMAT_UNEXPECTED; }

                            /* Pull data */
                            tmp_size = dimlength[dimid[dimids[0]]] * dimlength[dimid[dimids[1]]];
                            tmp_size *= ndimsp==3 ? dimlength[dimid[dimids[2]]] : 1;
                            h->nEmitters = dimlength[dimid[dimids[1]]] == 3 ? (int)dimlength[dimid[dimids[0]]] : (int)dimlength[dimid[dimids[1]]];
                            tmp_data = realloc1d(tmp_data, tmp_size*sizeof(double));
                            nc_get_var(ncid, varid, tmp_data);
                            h->EmitterPosition = malloc1d(tmp_size*sizeof(float));
                            for(i=0; i<(int)tmp_size; i++)
                                h->EmitterPosition[i] = (float)tmp_data[i];
                            break;

                        /* Some SOFA files have the "EmitterPosition" variable with the following dimensions: 1 x nEmitters x 3
                         * This is the reason for the switch case found here {2,3}, as it is not fully understood if this '1' is
                         * for the number of listeners? - Until this is made clear, the
                         * following code will just pull the first one (i.e. nEmitters x 3). Therefore, if you have an example
                         * of a file that has "EmitterPosition" with the dimensions: N x nEmitters x 3  (where N>1)
                         * then please send it to the developers :-) */
                        case 3:
                            /* Checks */
                            if((int)dimlength[dimid[dimids[1]]] != 3 && (int)dimlength[dimid[dimids[0]]] != 3) { return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED; }
                            if(typep!=NC_DOUBLE) { return SAF_SOFA_ERROR_FORMAT_UNEXPECTED; }

                            /* Pull data */
                            tmp_size = dimlength[dimid[dimids[0]]] * dimlength[dimid[dimids[1]]];
                            tmp_size *= ndimsp==3 ? dimlength[dimid[dimids[2]]] : 1;
                            h->nEmitters = dimlength[dimid[dimids[1]]] == 3 ? (int)dimlength[dimid[dimids[0]]] : (int)dimlength[dimid[dimids[1]]];
                            tmp_data = realloc1d(tmp_data, tmp_size*sizeof(double));
                            nc_get_var(ncid, varid, tmp_data);
                            h->EmitterPosition = malloc1d(tmp_size*sizeof(float));
                            for(i=0; i<(int)tmp_size; i++)
                                h->EmitterPosition[i] = (float)tmp_data[i];
                            break;
                        default:
                            return SAF_SOFA_ERROR_DIMENSIONS_UNEXPECTED;
                    }

                    /* Pull associated attributes */
                    for(attnum=0; attnum<varnattsp; attnum++){
                        nc_inq_attname(ncid, varid, attnum, attname);
                        nc_inq_attlen(ncid, varid, attname, &lenp);
                        if (!strcmp((char*)attname,"Type")){
                            h->EmitterPositionType = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                            nc_get_att(ncid, varid, attname, h->EmitterPositionType);
                        }
                        else if (!strcmp((char*)attname,"Units")){
                            h->EmitterPositionUnits = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                            nc_get_att(ncid, varid, attname, h->EmitterPositionUnits);
                        }
                    }
                }
            }

            /* Loop over the attributes and pull the info accordingly */
            for(attnum=0; attnum<nattsp; attnum++){
                nc_inq_attname(ncid, -1, attnum, attname);
                nc_inq_attlen(ncid, -1, attname, &lenp);

                if (!strcmp((char*)attname,"DataType")){
                    h->DataType = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, -1, attname, h->DataType);
                }
                else if (!strcmp((char*)attname,"Conventions")){
                    h->Conventions = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, -1, attname, h->Conventions);
                }
                else if (!strcmp((char*)attname,"Version")){
                    h->Version = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->Version);
                }
                else if (!strcmp((char*)attname,"SOFAConventions")){
                    h->SOFAConventions = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->SOFAConventions);
                }
                else if (!strcmp((char*)attname,"SOFAConventionsVersion")){
                    h->SOFAConventionsVersion = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->SOFAConventionsVersion);
                }
                else if (!strcmp((char*)attname,"APIName")){
                    h->APIName = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->APIName);
                }
                else if (!strcmp((char*)attname,"APIVersion")){
                    h->APIVersion = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->APIVersion);
                }
                else if (!strcmp((char*)attname,"ApplicationName")){
                    h->ApplicationName = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->ApplicationName);
                }
                else if (!strcmp((char*)attname,"ApplicationVersion")){
                    h->ApplicationVersion = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->ApplicationVersion);
                }
                else if (!strcmp((char*)attname,"AuthorContact")){
                    h->AuthorContact = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->AuthorContact);
                }
                else if (!strcmp((char*)attname,"Comment")){
                    h->Comment = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->Comment);
                }
                else if (!strcmp((char*)attname,"History")){
                    h->History = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->History);
                }
                else if (!strcmp((char*)attname,"License")){
                    h->License = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->License);
                }
                else if (!strcmp((char*)attname,"Organization")||!strcmp((char*)attname,"Organisation")){
                    h->Organisation = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->Organisation);
                }
                else if (!strcmp((char*)attname,"References")){
                    h->References = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->References);
                }
                else if (!strcmp((char*)attname,"RoomType")){
                    h->RoomType = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->RoomType);
                }
                else if (!strcmp((char*)attname,"Origin")){
                    h->Origin = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->Origin);
                }
                else if (!strcmp((char*)attname,"DateCreated")){
                    h->DateCreated = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->DateCreated);
                }
                else if (!strcmp((char*)attname,"DateModified")){
                    h->DateModified = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->DateModified);
                }
                else if (!strcmp((char*)attname,"Title")){
                    h->Title = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->Title);
                }
                else if (!strcmp((char*)attname,"DatabaseName")){
                    h->DatabaseName = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->DatabaseName);
                }
                else if (!strcmp((char*)attname,"ListenerShortName")){
                    h->ListenerShortName = calloc1d(SAF_MAX((NC_MAX_NAME+1),lenp), sizeof(char));
                    nc_get_att(ncid, NC_GLOBAL, attname, h->ListenerShortName);
                }
            }

            /* Close the file and clean-up */
            nc_close(ncid);
            free(dimid);
            free(dimlength);
            free(dimname);
            free(dimids);
            free(tmp_data);
#else
            saf_print_error("SAF_ENABLE_NETCDF must be defined to use this SOFA reader!");
#endif /* SAF_ENABLE_NETCDF */
            break;
    }

    return SAF_SOFA_OK;
}

void saf_sofa_close
(
    saf_sofa_container* c
)
{
    /* If NetCDF was used: */
    if (c->hLMSOFA == NULL){
        /* Vars */
        free(c->DataIR);
        free(c->SourcePosition);
        free(c->ReceiverPosition);
        free(c->DataDelay);
        free(c->ListenerPosition);
        free(c->ListenerView);
        free(c->ListenerUp);
        free(c->EmitterPosition);

        /* Var Atts */
        free(c->ListenerPositionType);
        free(c->ListenerPositionUnits);
        free(c->ListenerViewType);
        free(c->ListenerViewUnits);
        free(c->ReceiverPositionType);
        free(c->ReceiverPositionUnits);
        free(c->SourcePositionType);
        free(c->SourcePositionUnits);
        free(c->EmitterPositionType);
        free(c->EmitterPositionUnits);
        free(c->DataSamplingRateUnits);

        /* Global Atts */
        free(c->Conventions);
        free(c->Version);
        free(c->SOFAConventions);
        free(c->SOFAConventionsVersion);
        free(c->APIName);
        free(c->APIVersion);
        free(c->ApplicationName);
        free(c->ApplicationVersion);
        free(c->AuthorContact);
        free(c->Comment);
        free(c->DataType);
        free(c->History);
        free(c->License);
        free(c->Organisation);
        free(c->References);
        free(c->RoomType);
        free(c->Origin);
        free(c->DateCreated);
        free(c->DateModified);
        free(c->Title);
        free(c->DatabaseName);
        free(c->ListenerShortName);
    }
    else
        mysofa_free((MYSOFA_HRTF*)c->hLMSOFA);
}

#endif /* SAF_ENABLE_SOFA_READER_MODULE */
