# Spatial_Audio_Framework

A cross-platform framework for developing spatial audio related applications in C/C++.

![](saf.png)

* git: [https://github.com/leomccormack/Spatial_Audio_Framework](https://github.com/leomccormack/Spatial_Audio_Framework)
* doxygen: [http://research.spa.aalto.fi/projects/spatial_audio_framework](http://research.spa.aalto.fi/projects/spatial_audio_framework)

## Prerequisites

The Spatial_Audio_Framework (SAF) requires the following libraries:
* A library (or libraries) conforming to the [CBLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms#Implementations) and [LAPACK](https://en.wikipedia.org/wiki/LAPACK) standards
* (Optional) [netCDF](https://www.unidata.ucar.edu/software/netcdf/) for reading [SOFA](https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)) files

Add one of the following pre-processor definitions to inform SAF as to which performance library/libraries you have linked to your project:
```
SAF_USE_INTEL_MKL             # great option, but only for x86 architectures    
SAF_USE_OPEN_BLAS_AND_LAPACKE # good option, works on everything
SAF_USE_APPLE_ACCELERATE      # solid option (x86 and ARM), but MacOS only and slower than MKL
SAF_USE_ATLAS                 # bad option (x86 and ARM), many LAPACK functions are missing
``` 
Detailed instructions regarding how to build and link these performance libraries can be found [here](dependencies/PERFORMANCE_LIBRARY_INSTRUCTIONS.md).

## Framework structure

The framework comprises the following core modules:
* **saf_hoa** - a collection of higher-order Ambisonics binaural and loudspeaker decoders.
* **saf_sh** - spherical harmonic and spherical array processing related functions.
* **saf_vbap** - Vector-base Amplitude Panning (VBAP) functions.
* **saf_cdf4sap** - Covariance Domain Framework for Spatial Audio Processing (CDF4SAP).
* **saf_hrir** - HRIR/HRTF related functions (estimating ITDs, HRTF interpolation, diffuse-field EQ etc.).
* **saf_reverb** - a collection of reverbs and room simulation algorithms.
* **saf_utilities** - a collection of useful utility functions and cross-platform wrappers.

The framework also includes the following optional modules:
* **saf_sofa_reader** - a simple SOFA file reader.

To enable optional framework modules, simply add the relevant pre-processor definitions:
```
SAF_ENABLE_SOFA_READER_MODULE  # to enable saf_sofa_reader
```
Note that the **saf_sofa_reader** module also requires [netCDF](https://www.unidata.ucar.edu/software/netcdf/) to be linked to your project. Instructions on how to install/link this dependency can be found [here](dependencies/SOFA_READER_MODULE_DEPENDENCIES.md). 


## Using the framework

Once a CBLAS/LAPACK flag is defined (and the correct libraries are linked to your project), add the files found in the **framework** folder to your project and add the following directory to your project's header search paths:

```
Spatial_Audio_Framework/framework/include  
```

The framework's master include header is:

```c
#include "saf.h"
```

## Building with CMake 

The framework may be included within an existing CMake workflow with simply:
```
add_subdirectory(Spatial_Audio_Framework)
target_link_libraries(${PROJECT_NAME} PRIVATE saf)
```

The available SAF-related build options (and their default values) are:
```
-DSAF_PERFORMANCE_LIB=SAF_USE_INTEL_MKL    # performance library to employ
-DSAF_ENABLE_SOFA_READER_MODULE=0          # enable/disable the saf_sofa_reader module 
-DSAF_BUILD_EXAMPLES=1                     # build saf examples
-DSAF_BUILD_TESTS=1                        # build unit testing program
```

If using **SAF_USE_INTEL_MKL** as the performance library, note that the default header and library search paths may be overridden with:
``` 
-DINTEL_MKL_HEADER_PATH="path/to/mkl/headers"
-DINTEL_MKL_LIB="path/to/mkl/libs/mkl_rt(.so/.dylib/.lib)"   # OR:
-DINTEL_MKL_LIB="path/to/custom/mkl/lib/saf_mkl_custom(.so/.dylib/.lib)"
```

If the **saf_sofa_reader** is enabled, CMake will use the statically built dependencies found in **dependencies** for MacOSX and MSVC users by default. Linux and MSYS2 users may instead install a shared [netcdf library](dependencies/SOFA_READER_MODULE_DEPENDENCIES.md) and inform CMake of its location via:
```
# e.g. Linux users:
-DNETCDF_LIBRARY="/usr/lib/x86_64-linux-gnu/libnetcdf.so"
# e.g. MSYS2 users
-DNETCDF_LIBRARY="/c/msys64/mingw64/lib/libnetcdf.dll.a"
```

The framework, examples, and unit testing program may be built with:
```
mkdir build 
cmake -S . -B build 
cd build
make
test/saf_test # To run the unit testing program
```

## Documentation

Documentation is periodically generated using [Doxygen](http://www.doxygen.nl/index.html) and hosted [here](http://research.spa.aalto.fi/projects/spatial_audio_framework/index.html).

Alternatively, you may compile the most recent documentation (HTML) yourself with the following commands:
```
cd doxygen
doxygen doxygen_config
# (optional) to also build the pdf version:
cd latex
make
```

## Examples

Several examples have also been included in the repository, which may serve as a starting point for learning the framework:

* **ambi_bin** - a binaural Ambisonic decoder with built-in rotator. It supports the following decoding options: least-squares (LS), spatial re-sampling (SPR), Time-alignment (TA) [1], Magnitude Least-Squares (MagLS) [2].
* **ambi_dec** - a frequency-dependent Ambisonic decoder. It supports the following decoding options: sampling Ambisonic decoder (SAD), AllRAD [3], Energy-Preserving decoder (EPAD) [4], Mode-Matching decoder (MMD).
* **ambi_drc** - a frequency-dependent dynamic range compressor (DRC) for Ambisonic signals, based on the design proposed in [5].
* **ambi_enc** - a basic Ambisonic encoder.
* **array2sh** - converts microphone array signals into spherical harmonic signals (aka Ambisonic signals), based on theoretical descriptions [6,7]; more details found in [8].
* **beamformer** - a beamformer example with several different beamforming options.
* **binauraliser** - convolves input audio with interpolated HRTFs, which can be optionally loaded from a SOFA file.
* **decorrelator** - a basic multi-channel signal decorrelator.
* **dirass** - a sound-field visualiser based on re-assigning the energy of beamformers. This re-assignment is based on the DoA estimates extracted from spatially-localised active-intensity vectors, which are biased towards each beamformer direction [9].
* **matrixconv** - a basic matrix convolver with an optional partitioned convolution mode. 
* **multiconv** - a basic multi-channel convolver with an optional partitioned convolution mode. 
* **panner** - a frequency-dependent VBAP panner [10], which permits source loudness compensation as a function of the room [11].
* **pitch_shifter** - a basic multi-channel pitch shifter, based on the phase vocoder approach.
* **powermap** - sound-field visualiser based on beamformer (PWD, MVDR) energy or sub-space-based methods (MUSIC).
* **rotator** - rotates spherical harmonic signals (aka Ambisonic signals) given yaw-pitch-roll angles [12].
* **sldoa** - a sound-field visualiser based on directly depicting the DoA estimates extracted from multiple spatially-localised active-intensity vectors, at multiple frequencies; as proposed in [8]. 

Many of these examples have also been released as VST audio plug-ins under the [SPARTA](https://github.com/leomccormack/SPARTA) banner.

The following open-source projects also employ the framework: [HO-SIRR-GUI](https://github.com/leomccormack/HO-SIRR-GUI), and [CroPaC-Binaural](https://github.com/leomccormack/CroPaC-Binaural).

## Contributing

Suggestions and contributions to the code are both welcomed and encouraged. It should be highlighted that the framework has been designed to be highly modular with plenty of room for expansion. Therefore:
* if you are researcher who has developed a spatial-audio related method and want to integrate it into the framework... or
* if you notice that an existing piece of code can be rewritten to make it clearer, faster, or to fix a bug...

then please feel free to do so and submit a pull request. Note, however, that if the changes/additions are major, then maybe consider first discussing it via a github "issue" or by contacting the developers directly via email. We may also be able to help with the implementation if needed :- )

## Contributors

* **Leo McCormack** - C programmer and algorithm design (contact: leo.mccormack(at)aalto.fi)
* **Symeon Delikaris-Manias** - algorithm design
* **Archontis Politis** - algorithm design
* **Ville Pulkki** - algorithm design
* **Juhani Paasonen** - C programmer

## License

This framework is provided under the permissive [ISC license](https://choosealicense.com/licenses/isc/). It also includes the following resources: a modified version of the [alias-free STFT](https://github.com/jvilkamo/afSTFT) (MIT license); [kissFFT](https://github.com/mborgerding/kissfft) (BSD 3-clause license); [convhull_3d](https://github.com/leomccormack/convhull_3d) (MIT license); and [md_malloc](https://github.com/leomccormack/md_malloc) (MIT license).

Please be aware of and respect the license conditions found at the top of each file. Furthermore, while we do not impose any copyleft licensing philosophies, we would still appreciate it if any improvements and/or bug fixes are also merged into this public repository where possible.

## References

[1] Zaunschirm M, Scho"rkhuber C, Ho"ldrich R. **Binaural rendering of Ambisonic signals by head-related impulse response time alignment and a diffuseness constraint**.
The Journal of the Acoustical Society of America. 2018 Jun 19;143(6):3616-27.

[2] Scho"rkhuber C, Zaunschirm M, Ho"ldrich R. **Binaural Rendering of Ambisonic Signals via Magnitude Least Squares**.
InProceedings of the DAGA 2018 (Vol. 44, pp. 339-342).

[3] Zotter F, Frank M. **All-round ambisonic panning and decoding**.
Journal of the audio engineering society. 2012 Nov 26;60(10):807-20.

[4] Zotter F, Pomberger H, Noisternig M. **Energy-preserving ambisonic decoding**.
Acta Acustica united with Acustica. 2012 Jan 1;98(1):37-47.

[5] McCormack L, Va"lima"ki V. **FFT-based Dynamic Range Compression**.
In Proceedings of the 14th Sound and Music Computing Conference, July 5--8, Espoo, Finland, At Espoo, Finland 2017

[6] Williams EG. **Fourier acoustics: sound radiation and nearfield acoustical holography**.
Elsevier; 1999 Jun 10.

[7] Rafaely B. **Fundamentals of spherical array processing**.
Berlin: Springer; 2015 Feb 18.

[8] McCormack L, Delikaris-Manias S, Farina A, Pinardi D, Pulkki V. **Real-Time Conversion of Sensor Array Signals into Spherical Harmonic Signals with Applications to Spatially Localized Sub-Band Sound-Field Analysis**.
In Audio Engineering Society Convention 144 2018 May 14. Audio Engineering Society.

[9] McCormack L, Politis A, Pulkki V. **Sharpening of Angular Spectra Based on a Directional Re-assignment Approach for Ambisonic Sound-field Visualisation**.
In ICASSP 2019-2019 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP) 2019 Apr 17 (pp. 576-580). IEEE.

[10] Pulkki V. **Virtual sound source positioning using vector base amplitude panning**.
Journal of the audio engineering society. 1997 Jun 1;45(6):456-66.

[11] Laitinen MV, Vilkamo J, Jussila K, Politis A, Pulkki V. **Gain normalization in amplitude panning as a function of frequency and room reverberance**.
In Audio Engineering Society Conference: 55th International Conference: Spatial Audio 2014 Aug 26. Audio Engineering Society.

[12] Ivanic J, Ruedenberg K. **Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion**.
The Journal of Physical Chemistry A. 1998 Nov 5;102(45):9099-100.
