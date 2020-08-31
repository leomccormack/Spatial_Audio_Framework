/*
 * Copyright 2020 Leo McCormack
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
 * @file safmex_faf_IIRFilterbank.c
 * @brief MEX wrapper for faf_IIRFilterbank (see the .m file of the same name
 *        for documentation)
 * @author Leo McCormack
 * @date 06.08.2020
 */

#include "safmex.h"

/* ===================================================================== */
/*                                Config                                 */
/* ===================================================================== */

#define NUM_INPUT_ARGS_CREATE  ( 4 )
const MEX_DATA_TYPES inputDataTypes_create[NUM_INPUT_ARGS_CREATE] = {
    SM_INT32,                 
    SM_DOUBLE_REAL_1D,  
    SM_INT32,
    SM_DOUBLE_REAL
}; 
#define NUM_INPUT_ARGS_APPLY ( 1 )
#define NUM_OUTPUT_ARGS_APPLY ( 1 ) 
const MEX_DATA_TYPES inputDataTypes_apply[NUM_INPUT_ARGS_APPLY] = {
    SM_DOUBLE_REAL_1D      
};   
const MEX_DATA_TYPES outputDataTypes_apply[NUM_OUTPUT_ARGS_APPLY] = {
    SM_DOUBLE_REAL_2D 
}; 


/* ===================================================================== */
/*                                 Vars                                  */
/* ===================================================================== */

/* user arguments */
int order;          /* filter order, 1 or 3 */
float* fc = NULL;   /* filter cutoff frequencies */
int lSig;           /* Number of samples to process at a time */
float fs;           /* sampling rate */

/* internal parameters */
void* hFaF = NULL;               /* faf handle */
float* data_in = NULL;
float** data_out = NULL; 
int nCutoffFreqs;


/* ===================================================================== */
/*                              MEX Wrapper                              */
/* ===================================================================== */

void mexFunction
(
    int nlhs,             /* Number of output argments */
    mxArray *plhs[],      /* Pointers for output arguments */
    int nrhs,             /* Number of input argments */
    const mxArray *prhs[] /* Pointers for input arguments */
)
{  
    /* mex variables */
    int i, nDims;
    int *pDims = NULL;
     
    /* DESTROY */
    if(nrhs == 0){
        if(hFaF!=NULL){
            mexPrintf("Destroying FaF filterbank.\n");
            faf_IIRFilterbank_destroy(&hFaF);
            free(fc);       fc = NULL;
            free(data_in);  data_in = NULL;
            free(data_out); data_out = NULL;
            hFaF = NULL;
        } 
        else
            mexPrintf("FaF filterbank is already dead!\n"); 
    }
    
    /* CREATE */
    else if( nrhs==NUM_INPUT_ARGS_CREATE && nlhs==0 ){
        if(hFaF!=NULL)
            mexErrMsgIdAndTxt("MyToolbox:inputError","safmex_faf_IIRFilterbank is already initialised! First destroy it if you want to change its configuration.");
        
        /* Check input argument datatypes are as expected */ 
        checkArgDataTypes((mxArray**)prhs, (MEX_DATA_TYPES*)inputDataTypes_create, NUM_INPUT_ARGS_CREATE); 
        
        /* Copy user arguments */
        order = (int)mxGetScalar(prhs[0]); 
        lSig = (int)mxGetScalar(prhs[2]); 
        fs = (float)mxGetScalar(prhs[3]); 
        MEXdouble2SAFsingle(prhs[1], &fc, &nDims, &pDims);
        
        /* Extra checks */
        nCutoffFreqs = pDims[0];
        if( !((order==1) || (order==3)) )
            mexErrMsgIdAndTxt("MyToolbox:inputError","'order' must be either 1 or 3."); 
        if( nCutoffFreqs<=1 )
            mexErrMsgIdAndTxt("MyToolbox:inputError","cut-off frequency vector must be longer than 1 element."); 
        
        /* Create an instance of the hFaF filterbank */ 
        faf_IIRFilterbank_create(&hFaF, order, fc, nCutoffFreqs, fs, lSig);  
        
        /* Allocate buffers */
        data_in = malloc1d(lSig * sizeof(float)); 
        data_out = (float**)malloc2d(nCutoffFreqs+1, lSig, sizeof(float)); 
         
        /* Mainly just for debugging... */
        mexPrintf("Creating FaF filterbank:");
        snprintf(message, MSG_STR_LENGTH, " filter order = %d,", order); mexPrintf(message);
        snprintf(message, MSG_STR_LENGTH, " signal length = %d,", lSig); mexPrintf(message); 
        mexPrintf(" filter cut-off frequencies = ["); 
        for(i=0; i<nCutoffFreqs; i++){
            snprintf(message, MSG_STR_LENGTH, " %.2f ", fc[i]); 
            mexPrintf(message);
        }  
        mexPrintf("]\n");
    }
    
    /* TRANSFORM */
    else if(nrhs == 1 && nlhs == 1){
        if(hFaF==NULL)
            mexErrMsgIdAndTxt("MyToolbox:inputError","safmex_faf_IIRFilterbank is uninitialised!");
          
        /* Check input argument datatypes are as expected */ 
        checkArgDataTypes((mxArray**)prhs, (MEX_DATA_TYPES*)inputDataTypes_apply, NUM_INPUT_ARGS_APPLY); 
 
        /* Find dimensionality of input */
        mwSize nDims_mx;
        const mwSize *pDims_mx;
        nDims_mx = mxGetNumberOfDimensions(prhs[0]);
        pDims_mx = mxGetDimensions(prhs[0]); 
       
        /* extra checks */
        if( pDims_mx[0] != (mwSize)lSig ){
            snprintf(message, MSG_STR_LENGTH, "Was expecting %d samples.", lSig);
            mexErrMsgIdAndTxt("MyToolbox:inputError", message);
        }
        if( nDims_mx>1 )
            if(pDims_mx[1] != 1 || nDims_mx>2 )
                mexErrMsgIdAndTxt("MyToolbox:inputError", "Was expecting just one input channel.");
             
        /* Apply filterbank */
        MEXdouble2SAFsingle(prhs[0], &data_in, &nDims, &pDims);  
        faf_IIRFilterbank_apply(hFaF, data_in, data_out, lSig);
            
        /* output */
        nDims = 2;
        pDims = realloc1d(pDims, nDims*sizeof(int));
        pDims[0] = nCutoffFreqs+1; pDims[1] = lSig; 
        SAFsingle2MEXdouble(FLATTEN2D(data_out), nDims, pDims, &plhs[0]);

        /* Check output argument datatypes are as expected */ 
        checkArgDataTypes((mxArray**)plhs, (MEX_DATA_TYPES*)outputDataTypes_apply, NUM_OUTPUT_ARGS_APPLY); 
    }
 
    /* ERROR */
    else 
        mexErrMsgIdAndTxt("MyToolbox:inputError","Unrecognised input/output configuration, refer to help instructions.");
}
