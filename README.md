# Spatial_Audio_Framework

A cross-platform Spatial Audio Framework (SAF) written in C. The framework includes functions for performing Vector-Base Amplitude Panning (VBAP), Spherical Array Processing, and Higher-order Ambisonics (HOA); to name a few.

![](saf.png)


## Prerequisites

The framework requires the following libraries:
* Any Library/Libraries conforming to the [CBLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms#Implementations) and [LAPACK](https://en.wikipedia.org/wiki/LAPACK) standards
* (Optional) [netCDF](https://www.unidata.ucar.edu/software/netcdf/) for reading [SOFA](https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)) files

The rationale for the former requirement is that the framework employs the use of CBLAS/LAPACK routines for tackling all of the linear algebra operations, which are used quite prolifically throughout the code. Therefore, a performance library, which conforms to the CBLAS/LAPACK standards, is required by most of the framework modules. In principle, any such library (or combination of libraries) may be employed for this task, and if you've worked with such libraries before, then you probably already know what to do. Although, using a [custom Intel MKL library](CUSTOM_INTEL_MKL_INTRUCTIONS.md) for this requirement is still generally recommended, as this is the approach used by the developers.

## Available CBLAS/LAPACK options

Define one of the following preprocessor definitions, which will dictate which library/libraries you should link to your project:

```
SAF_USE_INTEL_MKL
SAF_USE_OPEN_BLAS_AND_LAPACKE
SAF_USE_ATLAS
```

* SAF_USE_INTEL_MKL - to use [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries), or a [**custom Intel MKL library**](CUSTOM_INTEL_MKL_INTRUCTIONS.md)  (recommended for x86_64/amd64).
* SAF_USE_OPEN_BLAS_AND_LAPACKE - to use [OpenBLAS](https://github.com/xianyi/OpenBLAS) and the LAPACKE interface (recommended for ARM)
* SAF_USE_ATLAS - to use [ALTAS](http://math-atlas.sourceforge.net/) which is not recommended, as some LAPACK functions are missing. However, if you don't mind loosing some SAF functionality then it may still be a good choice for your particular project.

**MacOSX users only**: If you do not define one of the above flags, then SAF will use [Apple Accelerate](https://developer.apple.com/documentation/accelerate) for CBLAS/LAPACK and also vDSP for the FFT. However, note that Intel MKL is still the more recommended option, as it is usually faster than Accelerate.

## Enable SOFA support (Optional)

In order to use the built-in [SOFA](https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)) reader (framework/modules/saf_hrir/saf_sofa_reader.h), your project must also link against the [netCDF](https://www.unidata.ucar.edu/software/netcdf/) library (including its dependencies). For those already familar with building and linking this particular library, you know what to do. However, for convenience, suggested platform specfic instructions have been provided below.

Also, the following preprocessor definition is required:
```
SAF_ENABLE_SOFA_READER
```

### Windows (64-bit) users

For convenience, the following statically built libraries are included in "dependencies/Win64/"; simply link your project against them:

```
libszip.lib; libzlib.lib; libhdf5.lib; libhdf5_hl.lib; netcdf.lib;
```

And also add the following two directories to your project's header and library search paths, respectively:

``` 
Spatial_Audio_Framework/dependencies/Win64/include
Spatial_Audio_Framework/dependencies/Win64/lib
```

### MacOSX users 

For convenience, the following statically built libraries are included in "dependencies/MacOSX/"; simply link your project against them:

```
netcdf; hdf5; hdf5_hl; z; 
```

And also add the following two directories to your project's header and library search paths, respectively:

``` 
Spatial_Audio_Framework/dependencies/MacOSX/include
Spatial_Audio_Framework/dependencies/MacOSX/lib
```

###  Linux (amd64) and Raspberry Pi (ARM) users

For ubuntu based distros, you may install [netCDF](https://www.unidata.ucar.edu/software/netcdf/) and its dependencies with these terminal commands:

```
sudo apt-get install libhdf5-dev
sudo apt-get install libnetcdf-dev libnetcdff-dev
```

Then simply add the following directory to the header search path:

```
/usr/include  
```

And add this linker flag to your project (or wherever it was installed):

```
-L/lib/x86_64-linux-gnu -lnetcdf
```

## Using the framework

Add all of the files found in the "framework" folder to your project (you may simply drag the folder into your IDE if you wish), and add the following directory to your header search paths:

```
Spatial_Audio_Framework/framework/include 
```

The framework's master include header is:

```c
#include "saf.h"
```

Detailed instructions regarding how to use the functions offered by each of the framework's module, is provided in the main header file for the respective module (e.g. "/modules/saf_**sh**/saf_**sh**.h", or  "/modules/saf_**vbap**/saf_**vbap**.h").

### Examples

Many examples have been included in the repository, which may also serve as a starting point for learning the framework:

* **ambi_bin** - a binaural Ambisonic decoder with built-in rotator. It includes the following decoding approaches: least-squares (LS), spatial re-sampling (SPR), Time-alignment (TA) [1], Magnitude Least-Squares (MagLS) [2].
* **ambi_dec** - a frequency-dependent Ambisonic decoder. Including the following decoding approaches: sampling ambisonic decoder (SAD), AllRAD [3], Energy-Preserving decoder (EPAD) [4], Mode-Matching decoder (MMD).
* **ambi_drc** - a frequency-dependent dynamic range compressor (DRC) for Ambisonic signals, based on the design proposed in [5].
* **ambi_enc** - a simple Ambisonic encoder.
* **array2sh** - converts microphone array signals into spherical harmonic signals (aka Ambisonic signals), based on theoretical descriptions [6,7]. More details in [8].
* **beamformer** - a beamforming example with several different beamforming options.
* **binauraliser** - convolves input audio with interpolated HRTFs, which can be optionally loaded from a SOFA file.
* **dirass** - a sound-field visualiser based on re-assigning the energy of beamformers. This re-assignment is based on the DoA estimates extracted from spatially-localised active-intensity vectors, which are biased towards each beamformer direction [9].
* **panner** - a frequency-dependent VBAP panner [10], which allows for more consistent source loudness as a function of the room [11].
* **matrixconv** - a basic matrix convolver with an optional partitioned convolution mode. 
* **multiconv** - a basic multi-channel convolver with an optional partitioned convolution mode. 
* **powermap** - sound-field visualiser using beamformers (PWD, MVDR) or sub-space methods (MUSIC).
* **rotator** - rotates spherical harmonic signals (aka Ambisonic signals) given yaw-pitch-roll angles [12].
* **sldoa** - a sound-field visualiser based on directly depicting the DoA estimates extracted from multiple spatially-localised active-intensity vectors; as proposed in [8].
* **upmix** - a (soon to be) collection of upmixing algorithms (currently only stereo to 5.x upmixing).

Many of these examples have also been integrated into VST audio plug-ins using the JUCE framework and can be found [here](http://research.spa.aalto.fi/projects/sparta_vsts/).

## Developers

* **Leo McCormack** - C programmer and algorithm design (contact: leo.mccormack@aalto.fi)
* **Symeon Delikaris-Manias** - algorithm design
* **Archontis Politis** - algorithm design

## License

This framework is provided under the [ISC license](https://choosealicense.com/licenses/isc/). However, it also includes a modified version of the ['alias-free STFT'](https://github.com/jvilkamo/afSTFT) implementation by Juha Vilkamo (MIT license); [kissFFT](https://github.com/mborgerding/kissfft) (BSD 3-clause license) by Mark Borgerding; and the ['convhull_3d'](https://github.com/leomccormack/convhull_3d) 3-D Convex Hull implementation by Leo McCormack (MIT license).

## References

[1] Zaunschirm M, Schörkhuber C, Höldrich R. **Binaural rendering of Ambisonic signals by head-related impulse response time alignment and a diffuseness constraint**.  
The Journal of the Acoustical Society of America. 2018 Jun 19;143(6):3616-27.

[2] Schörkhuber C, Zaunschirm M, Höldrich R. **Binaural Rendering of Ambisonic Signals via Magnitude Least Squares**.   
InProceedings of the DAGA 2018 (Vol. 44, pp. 339-342).

[3] Zotter F, Frank M. **All-round ambisonic panning and decoding**.   
Journal of the audio engineering society. 2012 Nov 26;60(10):807-20.

[4] Zotter F, Pomberger H, Noisternig M. **Energy-preserving ambisonic decoding**.   
Acta Acustica united with Acustica. 2012 Jan 1;98(1):37-47.

[5] McCormack L, Välimäki V. **FFT-based Dynamic Range Compression**.   
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
