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
 * @file safmex_generateVBAPgainTable3D.c
 * @brief MEX wrapper for generateVBAPgainTable3D (see the .m file of the same
 *        name for documentation)
 * @author Leo McCormack
 * @date 06.08.2020
 */
 
#include "safmex.h"

/* ===================================================================== */
/*                                Config                                 */
/* ===================================================================== */

#define NUM_INPUT_ARGS  ( 6 )
#define NUM_OUTPUT_ARGS ( 1 )
const MEX_DATA_TYPES inputDataTypes[NUM_INPUT_ARGS] = {
    SM_DOUBLE_REAL_2D,
    SM_INT32,
    SM_INT32,
    SM_INT32,
    SM_INT32,
    SM_DOUBLE_REAL
};
const MEX_DATA_TYPES outputDataTypes[NUM_OUTPUT_ARGS] = { 
    SM_DOUBLE_REAL_2D
};


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
    /* check for proper number of arguments and input argument datatypes */
    checkNumInOutArgs(nrhs, nlhs, NUM_INPUT_ARGS, NUM_OUTPUT_ARGS); 
    checkArgDataTypes((mxArray**)prhs, (MEX_DATA_TYPES*)inputDataTypes, NUM_INPUT_ARGS); 
      
    /* mex variables */
    int nDims;
    int *pDims = NULL;

    /* saf variables */
    float* ls_dirs_deg = NULL;
    float* gtable = NULL;  
    int L, az_res_deg, el_res_deg, omitLargeTriangles, enableDummies;
    int N_gtable, nTriangles;
    float spread;
    
    /* prep */ 
    MEXdouble2SAFsingle(prhs[0], &ls_dirs_deg, &nDims, &pDims);
    L = pDims[0];
    az_res_deg = (int)mxGetScalar(prhs[1]);
    el_res_deg = (int)mxGetScalar(prhs[2]);
    omitLargeTriangles = (int)mxGetScalar(prhs[3]);
    enableDummies = (int)mxGetScalar(prhs[4]);
    spread = (float)mxGetScalar(prhs[5]); 
     
    /* any extra checks */  
    if(pDims[1]!=2)
        mexErrMsgIdAndTxt("MyToolbox:inputError","the second dimension of the first argument should be of size: 2");
    
    /* call SAF function */
//     snprintf(message, MSG_STR_LENGTH, "L = %d,\n", L); mexPrintf(message);
//     snprintf(message, MSG_STR_LENGTH, "az_res_deg = %d,\n", az_res_deg); mexPrintf(message);
//     snprintf(message, MSG_STR_LENGTH, "el_res_deg = %d,\n", el_res_deg); mexPrintf(message);
//     snprintf(message, MSG_STR_LENGTH, "omitLargeTriangles = %d,\n", omitLargeTriangles); mexPrintf(message);
//     snprintf(message, MSG_STR_LENGTH, "enableDummies = %d,\n", enableDummies); mexPrintf(message);
//     snprintf(message, MSG_STR_LENGTH, "spread = %.5f,\n", spread); mexPrintf(message);
     generateVBAPgainTable3D(ls_dirs_deg, L, az_res_deg, el_res_deg, omitLargeTriangles,
                             enableDummies, spread, &gtable, &N_gtable, &nTriangles);
    
    /* output */
    nDims = 2;
    pDims = realloc1d(pDims, nDims*sizeof(int));
    pDims[0] = N_gtable;
    pDims[1] = L;
    SAFsingle2MEXdouble(gtable, nDims, pDims, &plhs[0]);
            
    /* clean-up */
    free(pDims);
    free(ls_dirs_deg); 
    free(gtable); 
    
    /* check output argument datatypes */
    checkArgDataTypes((mxArray**)plhs, (MEX_DATA_TYPES*)outputDataTypes, NUM_OUTPUT_ARGS); 
}
