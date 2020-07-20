


#define SAMPLERATE_PERMIT_SAF_INVASION

#include "samplerate_interface.h"
#include "samplerate.h"

#include <assert.h>






void sampleratelib_resample
(
    float* insig,
    int length_insig,
    int length_outsig,
    int input_fs,
    int output_fs,
    int nChannels,
    RESAMPLE_QUALITY_OPTIONS quality,
    float* outsig
)
{
    SRC_DATA data;
    int err;


    data.data_in = insig;
    data.data_out = outsig;
    data.input_frames = length_insig; /* frames refers to number of samples */
    data.output_frames = length_outsig;
    data.src_ratio = (double)output_fs/(double)input_fs;
    err = src_simple(&data, (int)quality, nChannels); 


    const char* error = src_strerror (err) ;

    assert(err==0);

}
