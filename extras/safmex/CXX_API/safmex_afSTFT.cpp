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
 * @file safmex_afSTFT.cpp
 * @brief SAFMEX C++ wrapper for afSTFT
 * @author Leo McCormack
 * @date 22.10.2020
 */

#include "safmex.hpp"  

/** C++MEX Wrapper */
class MexFunction : public safmex {
private: 
    /* Local copy of user parameters */
    int nCHin;                   /**< Number of input channels */
    int nCHout;                  /**< Number of output channels */
    int hopsize;                 /**< Hop size, in samples */ 
    int hybridmode;              /**< 0: disabled, 1: hybrid-filtering enabled */
    AFSTFT_FDDATA_FORMAT format; /**< Frequency-domain frame format, 0:nBands x nChannels x nTimeHops 1: nTimeHops x nChannels x nBands */
    float fs;                    /**< Samplerate */
      
    /* Internals */
    void* hSTFT = NULL;          /**< afSTFT handle */ 
    float* freqVector;           /**< frequency vector; nBands x 1 */
    int nBands;                  /**< Number of frequency bands */
    int procDelay;               /**< Processing delay in samples */
    int timeSlots;               /**< Number of down-sampled time frames */
    float** dataTD_in;           /**< Input time-domain Buffer; nCHin x blocksize */
    float** dataTD_out;          /**< Output time-domain Buffer; nCHout x blocksize */
    float_complex*** dataFD_in;  /**< Input frequency-domain Buffer; nBands x nCHin x timeslots */
    float_complex*** dataFD_out; /**< Output frequency-domain Buffer; nBands x nCHout x timeslots  */
    
public:
    /** Default constructor */
    MexFunction(){
        hSTFT = NULL;
        freqVector = NULL;
        dataTD_in = NULL;
        dataTD_out = NULL;
        dataFD_in = NULL;
        dataFD_out = NULL;
    }
    
    /** Point of entry */
    void operator()(ArgumentList outputs, ArgumentList inputs) {   
        
        /* CREATE */
        if(inputs.size() == 6){ 
            /* Checks */  
            if(hSTFT!=NULL){
                stream << "Object has already been initialised! First destroy it in order to intialise another...";
                print2commandWindow(SAFMEX_ERROR);
            } 
            assert_isScalar(inputs, 0, 1, 64, ArrayType::DOUBLE);
            assert_isScalar(inputs, 1, 1, 64, ArrayType::DOUBLE);
            assert_isScalar(inputs, 2, 16, 1024, ArrayType::DOUBLE);
            assert_isScalar(inputs, 3, 0, 1, ArrayType::DOUBLE); 
            assert_isScalar(inputs, 4, 0, 1, ArrayType::DOUBLE); 
            assert_isScalar(inputs, 5, 100, 1e6, ArrayType::DOUBLE);
            
            /* Parse/Convert */
            nCHin = (int)inputs[0][0];
            nCHout = (int)inputs[1][0];
            hopsize = (int)inputs[2][0]; 
            hybridmode = (int)inputs[3][0];
            format = (AFSTFT_FDDATA_FORMAT)((int)inputs[4][0]); 
            fs  = inputs[5][0];
 
            /* Create */
            stream << "Creating and initialising an instance of afSTFT: nCHin=" << nCHin << ", nCHout=" << nCHout << ", hopsize=" << hopsize << ", hybridmode=" << hybridmode << ", format=" << format << ", fs=" << fs;
            print2commandWindow(SAFMEX_MESSAGE);
            afSTFT_create(&hSTFT, nCHin, nCHout, hopsize, 0, hybridmode, format);
            
            /* Store these variables locally */
            nBands = afSTFT_getNBands(hSTFT); 
            
            /* Optional return parameters */
            if(outputs.size() > 0){
                freqVector = (float*)realloc1d((void*)freqVector, nBands*sizeof(float));
                afSTFT_getCentreFreqs(hSTFT, fs, nBands, freqVector);  
                SAF2MEX_array(freqVector, nBands, 1, SAFMEX_NO_PERMUTE, &outputs[0]);  
            }
            if(outputs.size() > 1){
                int procDelay =  afSTFT_getProcDelay(hSTFT);
                SAF2MEX_array(&procDelay, 1, 1, SAFMEX_NO_PERMUTE, &outputs[1]); 
            }
        }
        
        /* PROCESS */
        else if(inputs.size() == 1 && outputs.size() == 1){  
            /* Checks */  
            if(hSTFT==NULL){
                stream << "Object has not yet been created and initialised!";
                print2commandWindow(SAFMEX_ERROR);
            } 
            
            /* FORWARD */
            if(inputs[0].getType() == ArrayType::DOUBLE){ // COMPLEX_DOUBLE 
                /* Checks */
                assert_isMatrix(inputs, 0, SAFMEX_ANY_LENGTH, nCHin, ArrayType::DOUBLE);
                int framesize = inputs[0].getDimensions()[0];
                if(framesize % hopsize != 0){
                    stream << "Input blocksize '" << framesize  << "' is not divisable by hopsize '" << hopsize << "'";
                    print2commandWindow(SAFMEX_ERROR);
                }  
                
                /* Parse/Convert */
                MEX2SAF_array(inputs[0], SAFMEX_INV_PERMUTE, &dataTD_in); 
                int nHops = framesize/hopsize;
            
                /* Process */
                switch(format){
                    case AFSTFT_BANDS_CH_TIME: dataFD_in = (float_complex***)realloc3d((void***)dataFD_in, nBands, nCHin, nHops, sizeof(float_complex)); break;
                    case AFSTFT_TIME_CH_BANDS: dataFD_in = (float_complex***)realloc3d((void***)dataFD_in, nHops, nCHin, nBands, sizeof(float_complex)); break;
                }
                afSTFT_forward(hSTFT, dataTD_in, framesize, dataFD_in);
                
                /* Output */
                switch(format){
                    case AFSTFT_BANDS_CH_TIME: SAF2MEX_array(dataFD_in, nBands, nCHin, nHops, SAFMEX_NO_PERMUTE, &outputs[0]);  ; break;
                    case AFSTFT_TIME_CH_BANDS: SAF2MEX_array(dataFD_in, nHops, nCHin, nBands, SAFMEX_NO_PERMUTE, &outputs[0]);  ; break;
                }
            }
            /* BACKWARD */
            else if(inputs[0].getType() == ArrayType::COMPLEX_DOUBLE){ 
                /* Checks */
                int framesize;
                switch(format){
                    case AFSTFT_BANDS_CH_TIME: 
                        assert_isMatrix(inputs, 0, nBands, nCHout, SAFMEX_ANY_LENGTH, ArrayType::COMPLEX_DOUBLE); 
                        framesize = inputs[0].getDimensions()[2]*hopsize;
                        break;
                    case AFSTFT_TIME_CH_BANDS: 
                        assert_isMatrix(inputs, 0, SAFMEX_ANY_LENGTH, nCHout, nBands, ArrayType::COMPLEX_DOUBLE);  
                        framesize = inputs[0].getDimensions()[0]*hopsize;
                        break;
                } 
                
                /* Parse/Convert */
                MEX2SAF_array(inputs[0], SAFMEX_NO_PERMUTE, &dataFD_out); 
                  
                /* Process */
                dataTD_out = (float**)realloc2d((void**)dataTD_out, framesize, nCHout, sizeof(float));
                afSTFT_backward(hSTFT, dataFD_out, framesize, dataTD_out);
               
                /* Output */
                SAF2MEX_array(dataTD_out, nCHout, framesize, SAFMEX_INV_PERMUTE, &outputs[0]);  
            }
            
            /* ERROR */
            else{
                stream << "Input/output argument configuration was unexpected.";
                print2commandWindow(SAFMEX_ERROR);
            } 
        } 
        
        /* DESTROY */
        else if(inputs.size() == 0 && outputs.size() == 0){
            stream << "Destroying afSTFT instance.";
            print2commandWindow(SAFMEX_MESSAGE);
            afSTFT_destroy(&hSTFT); hSTFT = NULL; 
            free(freqVector); freqVector = NULL;
            free(dataTD_in);  dataTD_in = NULL;
            free(dataTD_out); dataTD_out = NULL;
            free(dataFD_in);  dataFD_in = NULL;
            free(dataFD_out); dataFD_out = NULL;
        }
        
        /* ERROR */
        else{
            stream << "Input/output argument configuration was unexpected.";
            print2commandWindow(SAFMEX_ERROR);
        }
    } 
};