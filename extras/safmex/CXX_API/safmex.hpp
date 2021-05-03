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
 * @file safmex.hpp
 * @brief Main header for the SAFMEX C++ API
 * @author Leo McCormack
 * @date 22.10.2020
 */

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "saf.h" 

using namespace matlab::data;
using namespace matlab::mex;

/** Can be passed to assert_isMatrix() to indicate that the matrix can have any dimension length */
#define SAFMEX_ANY_LENGTH ( -1 )

/** Available permute options for SAF2MEX_array() and MEX2SAF_array() */
typedef enum{
    SAFMEX_NO_PERMUTE,  /**< No permutation */
    SAFMEX_INV_PERMUTE  /**< Intended for row<->column major reordering, i.e.: ^T for 2-D arrays, permute(A, [3 2 1]) for 3-D arrays etc. */
}SAFMEX_PERMUTE_OPTIONS;

/** Available SAFMEX message options */
typedef enum{
    SAFMEX_MESSAGE,  /**< Message printed unchanged */
    SAFMEX_WARNING,  /**< Adds the "SAF WARNING: " prefix to message */
    SAFMEX_ERROR     /**< Adds the "SAF ERROR: " prefix to message and also triggers a MATLAB error */
}SAFMEX_MESSAGE_TYPE;

/** SAFMEX base class for interfacing with SAF within C++MEX wrappers; it comprises a bunch of helpful functions for: 
 *  argument parsing, argument checking, and printing messages to the command window etc. */ 
class safmex : public Function {
public:  
    /** Default constructor */
    safmex(){ 
        matlabPtr = getEngine();
        stream.str(""); /* Clear */
    }
    
    /** Converts a 2-D MEX array into a 2-D SAF array (real double->real float)
     *  @note Always check that the dimensions of 'arg' are what you expect; e.g. by first calling: assert_realMatrix()
     *  @note Always initialise (*pOutM) as NULL for the first call */
    void MEX2SAF_array(Array& arg,  SAFMEX_PERMUTE_OPTIONS tOpt, float*** pOutM){
        ArrayDimensions dims = arg.getDimensions(); 
        int nDims = dims.size();   
        switch(tOpt){
            case SAFMEX_INV_PERMUTE: {
                (*pOutM) = (float**)realloc2d((void**)(*pOutM), dims[1], dims[0], sizeof(float));
                for(int i=0; i<dims[0]; i++)
                    for(int j=0; j<dims[1]; j++)
                        (*pOutM)[j][i] = (float)arg[i][j];
                break;
            }
            case SAFMEX_NO_PERMUTE: {
                (*pOutM) = (float**)realloc2d((void**)(*pOutM), dims[0], dims[1], sizeof(float));
                for(int i=0; i<dims[0]; i++)
                    for(int j=0; j<dims[1]; j++)
                        (*pOutM)[i][j] = (float)arg[i][j];
                break;
            }
        }
    }
    
    /** Converts a 3-D MEX array into a 3-D SAF array (complex double->complex float)
     *  @note Always check that the dimensions of 'arg' are what you expect; e.g. by first calling: assert_realMatrix()
     *  @note Always initialise (*pOutM) as NULL for the first call */
    void MEX2SAF_array(Array& arg,  SAFMEX_PERMUTE_OPTIONS tOpt, float_complex**** pOutM){ 
        ArrayDimensions dims = arg.getDimensions();  
        TypedArray<std::complex<double>> ar(std::move(arg));
        std::complex<double>* cArg = ar.release().get();
        
        int nDims = dims.size();   
        switch(tOpt){
            case SAFMEX_INV_PERMUTE: {
                (*pOutM) = (float_complex***)realloc3d((void***)(*pOutM), dims[2], dims[1], dims[0], sizeof(float_complex));
                for(int i=0; i<dims[0]; i++)
                    for(int j=0; j<dims[1]; j++)
                        for(int k=0; k<dims[2]; k++)
                            (*pOutM)[k][j][i] = cmplxf((float)std::real(cArg[i*dims[1]*dims[2]+j*dims[2]+k]), (float)std::imag(cArg[i*dims[1]*dims[2]+j*dims[2]+k]));
                break;
            }
            case SAFMEX_NO_PERMUTE: {
                (*pOutM) = (float_complex***)realloc3d((void***)(*pOutM), dims[0], dims[1], dims[2], sizeof(float_complex));
                for(int i=0; i<dims[0]; i++)
                    for(int j=0; j<dims[1]; j++)
                        for(int k=0; k<dims[1]; k++)
                            (*pOutM)[i][j][k] = cmplxf((float)std::real(cArg[i*dims[1]*dims[2]+j*dims[2]+k]), (float)std::imag(cArg[i*dims[1]*dims[2]+j*dims[2]+k]));
                break;
            }
        }
    }
    
    /** Converts a 2-D SAF array into a 2-D MEX array (real float->real double) */
    void SAF2MEX_array(float** inM, int dim1, int dim2, SAFMEX_PERMUTE_OPTIONS tOpt, Array* arg){ 
        switch(tOpt){
            case SAFMEX_INV_PERMUTE: {
                TypedArray<double> outM_inv = factory.createArray<double>({ static_cast<unsigned long>(dim2), static_cast<unsigned long>(dim1) }); 
                for(int i=0; i<dim1; i++)
                    for(int j=0; j<dim2; j++)
                        outM_inv[j][i] = (double)inM[i][j]; 
                (*arg) = outM_inv;
                break;
            } 
            case SAFMEX_NO_PERMUTE: {
                TypedArray<double> outM = factory.createArray<double>({ static_cast<unsigned long>(dim1), static_cast<unsigned long>(dim2) }); 
                for(int i=0; i<dim1; i++)
                    for(int j=0; j<dim2; j++)
                        outM[i][j] = (double)inM[i][j]; 
                (*arg) = outM;
                break;
            }
        }    
    }
    
    /** Converts a FLAT 2-D SAF array into a 2-D MEX array (real float->real double) */
    void SAF2MEX_array(float* inM, int dim1, int dim2, SAFMEX_PERMUTE_OPTIONS tOpt, Array* arg){
        switch(tOpt){ 
            case SAFMEX_INV_PERMUTE: {
                TypedArray<double> outM_inv = factory.createArray<double>({ static_cast<unsigned long>(dim2), static_cast<unsigned long>(dim1) }); 
                for(int i=0; i<dim1; i++)
                    for(int j=0; j<dim2; j++)
                        outM_inv[j][i] = (double)inM[i*dim2+j]; 
                (*arg) = outM_inv;
                break;
            }
            case SAFMEX_NO_PERMUTE: {
                TypedArray<double> outM = factory.createArray<double>({ static_cast<unsigned long>(dim1), static_cast<unsigned long>(dim2) }); 
                for(int i=0; i<dim1; i++)
                    for(int j=0; j<dim2; j++)
                        outM[i][j] = (double)inM[i*dim2+j]; 
                (*arg) = outM;
                break;
            }
        }     
    }
    
    /** Converts a FLAT 2-D SAF array into a 2-D MEX array (int->double) */
    void SAF2MEX_array(int* inM, int dim1, int dim2, SAFMEX_PERMUTE_OPTIONS tOpt, Array* arg){
        switch(tOpt){
            case SAFMEX_INV_PERMUTE: {
                TypedArray<double> outM_inv = factory.createArray<double>({ static_cast<unsigned long>(dim2), static_cast<unsigned long>(dim1) }); 
                for(int i=0; i<dim1; i++)
                    for(int j=0; j<dim2; j++)
                        outM_inv[j][i] = (double)inM[i*dim2+j]; 
                (*arg) = outM_inv;
                break;
            }
            case SAFMEX_NO_PERMUTE: {
                TypedArray<double> outM = factory.createArray<double>({ static_cast<unsigned long>(dim1), static_cast<unsigned long>(dim2) }); 
                for(int i=0; i<dim1; i++)
                    for(int j=0; j<dim2; j++)
                        outM[i][j] = (double)inM[i*dim2+j]; 
                (*arg) = outM;
                break;
            }
        }   
    }
    
    /** Converts a 3-D SAF array into a 3-D MEX array (complex float->complex double) */
    void SAF2MEX_array(float_complex*** inM, int dim1, int dim2, int dim3, SAFMEX_PERMUTE_OPTIONS tOpt, Array* arg){
        switch(tOpt){
            case SAFMEX_INV_PERMUTE: {
                TypedArray<std::complex<double>> outM_inv = factory.createArray<std::complex<double>>({ static_cast<unsigned long>(dim3), static_cast<unsigned long>(dim2), static_cast<unsigned long>(dim1) }); 
                for(int i=0; i<dim1; i++)
                    for(int j=0; j<dim2; j++)
                        for(int k=0; k<dim3; k++)
                            outM_inv[k][j][i] = std::complex<double>(std::real(inM[i][j][k]), std::imag(inM[i][j][k])); 
                (*arg) = outM_inv;
                break;
            }
            case SAFMEX_NO_PERMUTE: {
                TypedArray<std::complex<double>> outM = factory.createArray<std::complex<double>>({ static_cast<unsigned long>(dim1), static_cast<unsigned long>(dim2), static_cast<unsigned long>(dim3) }); 
                for(int i=0; i<dim1; i++)
                    for(int j=0; j<dim2; j++)
                        for(int k=0; k<dim3; k++)
                            outM[i][j][k] = std::complex<double>(std::real(inM[i][j][k]), std::imag(inM[i][j][k])); 
                (*arg) = outM;
                break;
            }
        }    
    }
    
    /** Asserts that "args[argInd]" is a scalar and within the specified min and max values */
    void assert_isScalar(ArgumentList& args, int argInd, double minVal, double maxVal, ArrayType format){
        if (args[argInd].getNumberOfElements() != 1) {
            stream << "Argument " << argInd+1 << " must be a scalar!";
            print2commandWindow(SAFMEX_ERROR);
        } 
        if (args[argInd].getType() != format){
            stream << "Argument: " << argInd+1 << " is of an unsupported data type!";
            print2commandWindow(SAFMEX_ERROR);
        }
        TypedArray<double> val = args[argInd];
        if (val[0]>maxVal || val[0]<minVal){
            stream << "Argument: " << argInd+1 << " must be in the range ["<< minVal << ".." << maxVal << "]!";
            print2commandWindow(SAFMEX_ERROR);
        }
    }
    
    /** Asserts that "args[argInd]" is a 2-D matrix of dimensions: dim1 x dim2 */
    void assert_isMatrix(ArgumentList& args, int argInd, int dim1, int dim2, ArrayType format){
        ArrayDimensions dims = args[argInd].getDimensions(); 
        int nDims = dims.size(); 
        if (args[argInd].getType() != format){
            stream << "Argument: " << argInd+1 << " is of an unsupported data type!";
            print2commandWindow(SAFMEX_ERROR);
        }
        if(nDims!=2){
            stream << "Argument: " << argInd+1 << " must be a 2-D matrix!";
            print2commandWindow(SAFMEX_ERROR);
        }
        if((dim1!=SAFMEX_ANY_LENGTH && dims[0]!=dim1) || (dim2!=SAFMEX_ANY_LENGTH && dims[1]!=dim2)){
            stream << "Argument: " << argInd+1 << " must be a 2-D matrix with dimensions " << (dim1==SAFMEX_ANY_LENGTH ? "?" : std::to_string(dim1)) << " x " <<  (dim2==SAFMEX_ANY_LENGTH ? "?" : std::to_string(dim2)) << " !";
            print2commandWindow(SAFMEX_ERROR);
        }
        MemoryLayout memLayout = args[argInd].getMemoryLayout(); 
        if(memLayout != MemoryLayout::COLUMN_MAJOR){
            stream << "Argument: " << argInd+1 << " must be a 2-D matrix in column-major memory layout!";
            print2commandWindow(SAFMEX_ERROR);
        } 
    }
    
    /** Asserts that "args[argInd]" is a 3-D matrix of dimensions: dim1 x dim2 */
    void assert_isMatrix(ArgumentList& args, int argInd, int dim1, int dim2, int dim3, ArrayType format){
        ArrayDimensions dims = args[argInd].getDimensions(); 
        int nDims = dims.size(); 
        if (args[argInd].getType() != format){
            stream << "Argument: " << argInd+1 << " is of an unsupported data type!";
            print2commandWindow(SAFMEX_ERROR);
        }
        if(nDims!=3){
            stream << "Argument: " << argInd+1 << " must be a 3-D matrix!";
            print2commandWindow(SAFMEX_ERROR);
        }
        if((dim1!=SAFMEX_ANY_LENGTH && dims[0]!=dim1) || (dim2!=SAFMEX_ANY_LENGTH && dims[1]!=dim2) || (dim3!=SAFMEX_ANY_LENGTH && dims[2]!=dim3)){
            stream << "Argument: " << argInd+1 << " must be a 3-D matrix with dimensions " << (dim1==SAFMEX_ANY_LENGTH ? "?" : std::to_string(dim1)) << " x " << (dim2==SAFMEX_ANY_LENGTH ? "?" : std::to_string(dim2)) << " x " << (dim3==SAFMEX_ANY_LENGTH ? "?" : std::to_string(dim3)) << " !";
            print2commandWindow(SAFMEX_ERROR);
        }
        MemoryLayout memLayout = args[argInd].getMemoryLayout(); 
        if(memLayout != MemoryLayout::COLUMN_MAJOR){
            stream << "Argument: " << argInd+1 << " must be a 3-D matrix in column-major memory layout!";
            print2commandWindow(SAFMEX_ERROR);
        } 
    }
         
    /** Prints messages, warnings, and errors to the MATLAB command window */
    void print2commandWindow(SAFMEX_MESSAGE_TYPE message) { 
        stream << std::endl; /* New line */ 
        switch(message){
            case SAFMEX_MESSAGE: matlabPtr->feval(u"fprintf", 0, std::vector<Array>({ factory.createScalar(stream.str()) }));  break;
            case SAFMEX_WARNING: matlabPtr->feval(u"fprintf", 0, std::vector<Array>({ factory.createScalar("SAFMEX Warning: " + stream.str() ) })); break;
            case SAFMEX_ERROR:   matlabPtr->feval(u"error",   0, std::vector<Array>({ factory.createScalar("SAFMEX Error: " + stream.str()) })); break;
        }  
        stream.str(""); /* Clear */
    }
     
    /* Variables common to all SAFMEX wrappers */ 
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr; /**< Pointer to MATLAB engine */
    ArrayFactory factory;      /**< For handling arrays */
    std::ostringstream stream; /**< For printing messages */
};