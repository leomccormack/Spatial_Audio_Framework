![GitHub release (latest by date)](https://img.shields.io/github/v/release/leomccormack/Spatial_Audio_Framework) ![GitHub Release Date](https://img.shields.io/github/release-date/leomccormack/Spatial_Audio_Framework) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4707945.svg)](https://doi.org/10.5281/zenodo.4118286)  

<img src="saf.svg"> 

* git: [https://github.com/leomccormack/Spatial_Audio_Framework](https://github.com/leomccormack/Spatial_Audio_Framework)
* doxygen: [https://leomccormack.github.io/Spatial_Audio_Framework/](https://leomccormack.github.io/Spatial_Audio_Framework/)

# About

The Spatial_Audio_Framework (SAF) is an open-source and cross-platform framework for developing spatial audio related algorithms and software in C/C++. While originally intended as a resource for researchers in the field, the framework has gradually grown into a rather large codebase comprising a number of distinct [**modules**](framework/modules); with each module targeting a specific sub-field of spatial audio (e.g. Ambisonics encoding/decoding, spherical array processing, amplitude-panning, HRIR processing, room simulation, etc.). The framework also makes use of highly optimised linear algebra libraries (such as Intel MKL, Apple Accelerate, OpenBLAS) as well as SIMD intrinsics (SSE, AVX, AVX-512). Several [**examples**](examples/include) are also included in the repository, which serve to demonstrate the functionality of the framework and may also act as a starting point for new projects.

Owing to its modular design, expanding the framework is relatively straightforward, and contributions from researchers and developers of spatial audio technologies is actively encouraged! :-) 

# Prerequisites

The framework requires the following external libraries:
* Any library (or libraries) conforming to the [CBLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms#Implementations) and [LAPACK](https://en.wikipedia.org/wiki/LAPACK) standards
* (**Optional**) Intel's [Integrated Performance Primitives (IPP)](https://software.intel.com/content/www/us/en/develop/tools/integrated-performance-primitives.html) for the FFT and/or resampler
* (**Optional**) [FFTW](https://www.fftw.org/) for the FFT
* (**Optional**) a SSE, AVX, AVX-512 supporting CPU

In order to inform SAF which CBLAS/LAPACK supporting library/libraries you have linked to your project, simply add **one** of the following global pre-processor definitions:
```
SAF_USE_INTEL_MKL_LP64        # great option, but only for x86 architectures (using the LP64 config [int32])
SAF_USE_INTEL_MKL_ILP64       # great option, but only for x86 architectures (using the ILP64 config [int64])
SAF_USE_OPEN_BLAS_AND_LAPACKE # good option, works on everything
SAF_USE_APPLE_ACCELERATE      # good option (x86 and ARM), faster than OpenBLAS, but MacOSX only & slower than MKL
SAF_USE_ATLAS                 # bad option (x86 and ARM), many LAPACK functions are missing
SAF_USE...                    # please get in touch if you use something else! :-)
```
[Detailed instructions regarding how to build and link these libraries can be found here](docs/PERFORMANCE_LIBRARY_INSTRUCTIONS.md).

## Framework structure

The [framework](docs/FRAMEWORK_STRUCTURE.md) comprises the following core modules (**ISC** License):
* **saf_hoa** - a collection of higher-order Ambisonics binaural and loudspeaker decoders.
* **saf_sh** - spherical harmonic and spherical array processing related functions.
* **saf_vbap** - Vector-base Amplitude Panning (VBAP) functions.
* **saf_cdf4sap** - Covariance Domain Framework for Spatial Audio Processing (CDF4SAP).
* **saf_hrir** - HRIR/HRTF related functions (estimating ITDs, HRTF interpolation, diffuse-field EQ etc.).
* **saf_reverb** - a collection of reverbs and room simulation algorithms.
* **saf_utilities** - a collection of useful utility functions and cross-platform wrappers.

The framework also includes the following optional modules:
* **saf_sofa_reader** - a simple SOFA file reader (**ISC** License).
* **saf_tracker** - a particle-filtering based tracker (**GPLv2** License).
* **saf_hades** - for binaural rendering of Hearing-Assistive/Augmented-reality Devices (HADES)  (**GPLv2** License).

To enable optional framework modules, simply add the relevant pre-processor definition:
```
SAF_ENABLE_SOFA_READER_MODULE  # to enable saf_sofa_reader
SAF_ENABLE_TRACKER_MODULE      # to enable saf_tracker
SAF_ENABLE_HADES_MODULE        # to enable saf_hades
```

### Additional options

The framework can be configured further, with the following options:
```
SAF_USE_INTEL_IPP # To use Intel IPP for performing the DFT/FFT and resampling
SAF_USE_FFTW      # To use the FFTW library for performing the DFT/FFT 
SAF_ENABLE_SIMD   # To enable SIMD (SSE, AVX, AVX512) intrinsics for certain vector operations
```

# Using the framework

Once a CBLAS/LAPACK flag is defined (and the correct libraries are linked to your project), add the files found in the **framework** folder to your project and add the following directory to your project's header search paths:

```
Spatial_Audio_Framework/framework/include  
```

The framework's master include header is then:

```c
#include "saf.h"
#include "saf_externals.h"  /* (Optional) To also carry over CBLAS/LAPACK routines and other external functions. */
```

## Building with CMake 

![CMake Build](https://github.com/leomccormack/Spatial_Audio_Framework/workflows/CMake%20Build/badge.svg?branch=master)

The framework may also be included within an existing CMake workflow with simply:
```
add_subdirectory(Spatial_Audio_Framework)
target_link_libraries(${PROJECT_NAME} PRIVATE saf)
```

The available SAF-related CMake options (and their default values) are:
```
-DSAF_PERFORMANCE_LIB=SAF_USE_INTEL_MKL_LP64 # performance library to employ
-DSAF_ENABLE_SOFA_READER_MODULE=0            # enable/disable the saf_sofa_reader module 
-DSAF_ENABLE_TRACKER_MODULE=0                # enable/disable the saf_tracker module 
-DSAF_ENABLE_HADES_MODULE=0                  # enable/disable the saf_hades module 
-DSAF_BUILD_EXAMPLES=1                       # build saf examples
-DSAF_BUILD_EXTRAS=0                         # build safmex etc.
-DSAF_BUILD_TESTS=1                          # build unit testing program
-DSAF_USE_INTEL_IPP=0                        # link and use Intel IPP for the FFT, resampler, etc.
-DSAF_ENABLE_SIMD=0                          # enable/disable SSE, AVX, and/or AVX-512 support
-DSAF_ENABLE_NETCDF=0                        # enable the use of NetCDF (requires external libs)
-DSAF_ENABLE_FAST_MATH_FLAG=1                # enable the -ffast-math compiler flag on clang/gcc
```

If using e.g. **SAF_USE_INTEL_MKL_LP64** as the performance library, note that the default header and library search paths may be overridden [according to your setup](docs/PERFORMANCE_LIBRARY_INSTRUCTIONS.md) with:
``` 
-DINTEL_MKL_HEADER_PATH="path/to/mkl/headers"
-DINTEL_MKL_LIB="path/to/mkl/libs/mkl_rt(.so/.dylib/.lib)"   # OR:
-DINTEL_MKL_LIB="path/to/custom/mkl/lib/saf_mkl_custom_lp64(.so/.dylib/.lib)"
```

For Linux/MacOS users: the framework, examples, and unit testing program may be built as follows:
```
cmake -S . -B build 
# Or to also enable e.g. SSE3 and AVX2 intrinsics (for both C and C++ code)
cmake -S . -B build -DSAF_ENABLE_SIMD=1 -DCMAKE_CXX_FLAGS="-msse3 -mavx2" -DCMAKE_C_FLAGS="-msse3 -mavx2"
cd build
make
test/saf_test # To run the unit testing program
```

Or for Visual Studio users (using x64 Native Tools Command Prompt for VS):
```
# e.g. for VS2019:
cmake -S . -B build -G "Visual Studio 16" -A x64  
# e.g. for VS2017:
cmake -S . -B build -G "Visual Studio 15 Win64"
cd build
msbuild ALL_BUILD.vcxproj /p:Configuration=Release /m
cd test/Release
saf_test.exe  # To run the unit testing program
```

## Documentation 

![Doxygen Generate](https://github.com/leomccormack/Spatial_Audio_Framework/workflows/Doxygen%20Generate/badge.svg?branch=master)

[Doxygen](http://www.doxygen.nl/index.html)-based documentation is generated via a GitHub Action everytime a commit is pushed to the master branch. The documentation is hosted [here](https://leomccormack.github.io/Spatial_Audio_Framework/).

Alternatively, you may generate the documentation yourself (e.g. for the other branches) with the following commands:
```
cd docs/doxygen
doxygen doxygen_config
# (optional) to also build the pdf version:
cd latex
make
```

## Examples

Several **examples** have also been included in the repository, which may serve as a starting point for learning how to use the framework:

* **ambi_bin** - a binaural Ambisonic decoder with built-in rotator. It supports the following decoding options: least-squares (LS), spatial re-sampling (SPR), Time-alignment (TA), Magnitude Least-Squares (MagLS).
* **ambi_dec** - a frequency-dependent Ambisonic decoder. It supports the following decoding options: sampling Ambisonic decoder (SAD), AllRAD, Energy-Preserving decoder (EPAD), Mode-Matching decoder (MMD).
* **ambi_drc** - a frequency-dependent dynamic range compressor (DRC) for Ambisonic signals.
* **ambi_enc** - a basic Ambisonic encoder.
* **ambi_roomsim** - an Ambisonic encoder that also includes reflections and source distance based on an image-source model of a shoebox room. Multiple sources and Ambisonic receivers are supported.
* **array2sh** - converts microphone array signals into spherical harmonic signals (aka Ambisonic signals), based on theoretical descriptions of the array configuration and construction.
* **beamformer** - a beamformer/virtual microphone generator for Ambisonic signals, with several different beam pattern options.
* **binauraliser** - convolves input audio with interpolated HRTFs, which can be optionally loaded from a SOFA file.
* **binauraliser_nf** - binauraliser, with the addition of proximity filtering for near field sound sources.
* **decorrelator** - a basic multi-channel signal decorrelator.
* **dirass** - a sound-field visualiser based on re-assigning the energy of beamformers. This re-assignment is based on the DoA estimates extracted from spatially-localised active-intensity vectors, which are biased towards each beamformer direction.
* **matrixconv** - a basic matrix convolver with an optional partitioned convolution mode. 
* **multiconv** - a basic multi-channel convolver with an optional partitioned convolution mode. 
* **tvconv** - a time-varying partitioned convolution multi-channel convolver for SOFA files containing RIRs with multiple listener positions.
* **panner** - a frequency-dependent VBAP panner, which accommodates a source loudness compensation (as a function of the room) option.
* **pitch_shifter** - a basic multi-channel pitch shifter, based on the phase vocoder approach.
* **powermap** - sound-field visualiser based on beamformer (PWD, MVDR) energy or subspace methods (MUSIC).
* **rotator** - rotates spherical harmonic signals (aka Ambisonic signals) given yaw-pitch-roll Euler rotation angles.
* **sldoa** - a sound-field visualiser based on directly depicting the DoA estimates extracted from multiple spatially-localised active-intensity vectors for multiple frequencies. 
* **spreader** - an arbitrary array panner (HRIRs, microphone array IRs, etc.) with coherent and incoherent spreading modes.

Many of these examples have also been released as VST audio plug-ins under the [SPARTA](https://github.com/leomccormack/SPARTA) banner. The following open-source projects also employ the framework: [Super-Hearing](https://github.com/leomccormack/Super-Hearing/), [HO-SIRR-GUI](https://github.com/leomccormack/HO-SIRR-GUI), and [CroPaC-Binaural](https://github.com/leomccormack/CroPaC-Binaural).

## Extras

The repository also includes the following **extras**:

* [matlab](extras/matlab/SAF_MATLAB_CODE.md) - a bunch of MATLAB scripts/functions to accompany the framework (a script to generate saf_default_hrirs.c, MATLAB versions of certain SAF functions, etc.).
* [safmex](extras/safmex/SAFMEX.md) - a bunch of MATLAB MEX wrappers, which allow certain SAF functions to be used within MATLAB.
* [safpy](extras/safpy/SAFPY.md) - a work-in-progress initiative to bring SAF functionality to Python.
* [safwwise](extras/safwwise/SAFWWISE.md) - a proof of concept regarding how one might integrate SAF into Wwise.

## Contributing

Suggestions and contributions to the code are both welcomed and encouraged. It should be highlighted that the framework has been designed to be highly modular with plenty of room for expansion. Therefore:
* if you are researcher who has developed a spatial-audio related method and want to integrate it into the framework... or
* if you notice that an existing piece of code can be rewritten to make it clearer, faster, or to fix a bug...

then please feel free to do so and submit a pull request. We may also be able to help with the implementation if needed, just get in touch :- )

## Contributors

* **Leo McCormack** - C programming and algorithm design (contact: leo.mccormack(at)aalto.fi)
* **Symeon Delikaris-Manias** - algorithm design
* **Archontis Politis** - algorithm design
* **Ville Pulkki** - algorithm design
* **Juhani Paasonen** - C programming
* **Chris Hold** - C programming and algorithm design 
* **Janani Fernandez** - C programming and algorithm design 

# License

This software is dual-licensed. By default, this software is provided permissively under the terms of the [ISC License](https://choosealicense.com/licenses/isc/); since all of the core (non-optional) modules are licensed as such. However, including and enabling certain optional modules, which are instead provided under the copy-left [GNU GPLv2 License](https://choosealicense.com/licenses/gpl-2.0/), will mean that the use of this software is instead governed by the GNU GPLv2 licencing terms.

For full licensing terms see [LICENSE.md](LICENSE.md).

Note that while we do not impose any copyleft licensing philosophies for the ISC licensed modules, we still appreciate it when improvements and/or bug fixes are also merged into this public repository where possible :-)
