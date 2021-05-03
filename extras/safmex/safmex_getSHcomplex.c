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
 * @file safmex_getSHcomplex.c
 * @brief MEX wrapper for getSHcomplex (see the .m file of the same name for
 *        documentation)
 * @author Leo McCormack
 * @date 06.08.2020
 */
 
#include "safmex.h" 

/* ===================================================================== */
/*                                Config                                 */
/* ===================================================================== */

#define NUM_INPUT_ARGS  ( 2 )
#define NUM_OUTPUT_ARGS ( 1 )
const MEX_DATA_TYPES inputDataTypes[NUM_INPUT_ARGS] = {
    SM_INT32,
    SM_DOUBLE_REAL_1D_OR_2D
};
const MEX_DATA_TYPES outputDataTypes[NUM_OUTPUT_ARGS] = { 
    SM_DOUBLE_COMPLEX_1D_OR_2D
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
    float* dirs_rad = NULL;
    float_complex* Y = NULL;  
    int order, nSH, nDirs;
    
    /* prep */
    order = (int)mxGetScalar(prhs[0]);
    nSH = ORDER2NSH(order); 
    MEXdouble2SAFsingle(prhs[1], &dirs_rad, &nDims, &pDims);
    nDirs = pDims[0];
    Y = malloc1d(nSH*nDirs*sizeof(float_complex));
     
    /* any extra checks */  
    if(pDims[1]!=2)
        mexErrMsgIdAndTxt("MyToolbox:inputError","the second dimension of the second argument should be of size: 2");
    
    /* call SAF function */
    getSHcomplex(order, dirs_rad, nDirs, Y);
    
    /* output */
    nDims = 2;
    pDims = realloc1d(pDims, 2*sizeof(int));
    pDims[0] = nSH;
    pDims[1] = nDirs;
    SAFsingle2MEXdouble_complex(Y, nDims, pDims, &plhs[0]);
 
    /* clean-up */
    free(dirs_rad);
    free(pDims);
    free(Y); 
    
    /* check output argument datatypes */
    checkArgDataTypes((mxArray**)plhs, (MEX_DATA_TYPES*)outputDataTypes, NUM_OUTPUT_ARGS); 
}
