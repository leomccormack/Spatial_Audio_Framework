# Spatial_Audio_Framework

A cross-platform Spatial Audio Framework (SAF) written in C. The framework includes functions for performing Vector-Base Amplitude Panning (VBAP), Spherical Array Processing, and Higher-order Ambisonics (HOA); to name a few.

![](saf.png)

## Getting Started

### Prerequisites

The framework requires the following libraries:
* Any Library/Libraries conforming to the [BLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms#Implementations) and [LAPACK](https://en.wikipedia.org/wiki/LAPACK) standards
* (Optional) [netCDF](https://www.unidata.ucar.edu/software/netcdf/) for reading [SOFA](https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)) files

The rationale for the former requirement is that the framework employs the use of BLAS/LAPACK routines for tackling all of the linear algebra operations, which are used quite prolifically throughout the code. Therefore, a performance library, which conforms to the BLAS/LAPACK standards, is required by most of the framework modules. In principle, any such library (or combination of libraries) may be employed for this task, and if you've worked with such libraries before, then you probably already know what to do. However, it is still generally recommended to use a custom [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries) library, as this is the approach used by the developers; and the framework will also employ Intel MKL's FFT, which is pretty damn fast.

In order to make this a painless endevour, detailed instructions regarding the acquisition and linking of this custom Intel MKL library have been tailored for specific operating systems and provided below. 

## Acquiring and linking a custom [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries) library tailored for SAF (Recommended)

### Windows (64-bit) users

Note: use the "x64 Developer Command Prompt for VS.exe" (open as administrator) to run the following commands.

1. Install [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries). 

2. The required custom library may be obtained by first copying the included "**dependencies/saf_mkl_list**" file into the MKL "builder" folder:

```
xcopy saf_mkl_list C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl\tools\builder /R
```

3. EITHER (The blue pill): to generate and copy the "**saf_mkl_custom.dll**" library to a system path folder, and the "**saf_mkl_custom.lib**" for you to use locally, enter the following commands:

```
cd /Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/tools/builder
nmake intel64 interface=lp64 threading=sequential name=saf_mkl_custom export=saf_mkl_list
xcopy saf_mkl_custom.dll C:\Windows\System32 /R
xcopy saf_mkl_custom.lib C:\Users\[YOUR WORKING DIRECTORY]\Spatial_Audio_Framework\dependencies\Win64\lib /R
```

3. OR (The red pill): you may instead build a threaded version of the library (which involves some additional steps):

```
cd /Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/tools/builder
nmake intel64 interface=lp64 threading=parallel name=saf_mkl_custom export=saf_mkl_list
xcopy saf_mkl_custom.dll C:\Windows\System32 /R
xcopy saf_mkl_custom.lib C:\Users\[YOUR WORKING DIRECTORY]\Spatial_Audio_Framework\dependencies\Win64\lib /R
cd C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/compiler/lib/intel64/ 
xcopy libiomp5md.lib C:\Users\[YOUR WORKING DIRECTORY]\Spatial_Audio_Framework\dependencies\Win64\lib /R
cd C:\Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/redist/intel64/compiler
xcopy libiomp5md.dll C:\Windows\System32 /R
```

4. Add the following header search path to your project:

```
C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include
```

5. Add the following library search path to your project:

```
C:/Users/[YOUR WORKING DIRECTORY]/SDKs/Spatial_Audio_Framework/dependencies/Win64/lib
```

6. link your project against the following libraries (note that the second library is only needed if you built the threaded version):
```
saf_mkl_custom.lib
libiomp5md.lib
```

7. Add "SAF_USE_INTEL_MKL" to your pre-processor definitions.

Note: If you built the threaded version of the library, then there are some more pre-processors defintions you can use. See below.


### MacOSX users 

 By default, the framework will use Apple's Accelerate library for the BLAS/LAPACK routines and FFT, so you may ignore all of these steps if you wish. However, Mac users may still elect to use Intel MKL, as is it often faster than Accelerate. 

1. Install [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries).  

1. 1. Optionally, you may want to add the MKL global environment variables using this command:

```
source /opt/intel/compilers_and_libraries/mac/mkl/bin/mklvars.sh intel64
```

2. The required custom library may be obtained by first copying the included "**dependencies/saf_mkl_list**" file into the MKL "builder" folder:

```
sudo cp saf_mkl_list /opt/intel/compilers_and_libraries/mac/mkl/tools/builder
```

3. EITHER (The blue pill): to generate and copy the "**saf_mkl_custom.dylib**" library to a system path folder, ready for you to use, enter the following commands:

```
cd /opt/intel/compilers_and_libraries/mac/mkl/tools/builder
sudo make intel64 interface=lp64 threading=sequential name=libsaf_mkl_custom export=saf_mkl_list
sudo cp saf_mkl_custom.dylib /usr/local/lib
```

3. OR (The red pill): you may instead build a threaded version of the library (which involves some additional steps):

```
cd /opt/intel/compilers_and_libraries/mac/mkl/tools/builder
sudo make intel64 interface=lp64 threading=parallel name=libsaf_mkl_custom export=saf_mkl_list
sudo cp saf_mkl_custom.dylib /usr/local/lib
sudo cp /opt/intel/compilers_and_libraries/mac/compiler/lib/libiomp5.dylib /usr/local/lib
```

4. Add the following header search path to your project (where the X's are the version numbers):

```
/opt/intel/compilers_and_libraries_20XX.X.XXX/mac/mkl/include
```

5. Then add the following linker flags to your project (note that the second library is only needed if you built the threaded version):

```
-L/usr/local/lib -lsaf_mkl_custom
-L/usr/local/lib -liomp5    
```

6. Finally, add "SAF_USE_INTEL_MKL" to your project's pre-processor definitions.

Note: If you built the threaded version of the library, then there are some more pre-processors defintions you can use. See below.


### Linux (x86_64) users

0. Add the following directories to the header search paths:

```
Spatial_Audio_Framework/framework/include
Spatial_Audio_Framework/dependencies/Linux/include
```

1. Add the following directory to the library search paths:

```
Spatial_Audio_Framework/dependencies/Linux/lib
```

2. Install [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries). 

3. The required "**saf_mkl_custom.so**" file may be generated using Intel's custom dll builder. 

4. First copy the included "dependencies/saf_mkl_list" file to this folder:

```
/opt/intel/compilers_and_libraries/linux/mkl/tools/builder
```

5. To generate and copy this library to a system path folder, use the following commands (root permissions required):

```
cd /opt/intel/compilers_and_libraries/linux/mkl/tools/builder
sudo make intel64 interface=lp64 threading=sequential name=libsaf_mkl_custom export=saf_mkl_list
sudo cp saf_mkl_custom.so /usr/lib

```

6. Then add the following linker flag to your project:

```
-L/usr/lib -lsaf_mkl_custom
```

7. Add "SAF_USE_INTEL_MKL" to your pre-processor definitions.


### Intel MKL threading options (Optional)

If you built a threaded version of the custom library, then there are some additional options you may specify. More information can be found [here](https://software.intel.com/en-us/articles/recommended-settings-for-calling-intel-mkl-routines-from-multi-threaded-applications).

Note by default: MKL_NUM_THREADS = [number of CPU cores], MKL_DYNAMIC = 1. 

You may also change these at run time using the following functions:

```
/* for example: */
MKL_Set_Num_Threads(2);
MKL_Set_Dynamic(1);
```

## Enable SOFA support (Optional)

In order to use the built-in [SOFA](https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)) reader (framework/modules/saf_hrir/saf_sofa_reader.h), your project must also link against the [netCDF](https://www.unidata.ucar.edu/software/netcdf/) library (including its dependencies). For those already familar with building and linking this particular library, you know what to do. However, for convenience, suggested platform specfic instructions have been provided below.

### Windows (64-bit) users

For convenience, the following statically built libraries are included in "dependencies/Win64/"; simply link your project against them:

```
libszip.lib; libzlib.lib; libhdf5.lib; libhdf5_hl.lib; netcdf.lib;
```

Also add the following two directories to your project's header and library search paths, respectively:

``` 
Spatial_Audio_Framework/dependencies/Win64/include
Spatial_Audio_Framework/dependencies/Win64/lib
```

### MacOSX users 

For convenience, the following statically built libraries are included in "dependencies/MacOSX/"; simply link your project against them:

```
netcdf; hdf5; hdf5_hl; z; 
```

Also add the following two directories to your project's header and library search paths, respectively:

``` 
Spatial_Audio_Framework/dependencies/MacOSX/include
Spatial_Audio_Framework/dependencies/MacOSX/lib
```

###  Linux (x86_64) users

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

Add this directory to your header search paths:

```
Spatial_Audio_Framework/framework/include 
```

Note the framework is divided into individual modules and instructions on how to enable these modules is provided in the main header include file (saf.h). However, the general idea is that you enable modules by defining specific preprocessor flags:

```c
SAF_ENABLE_AFSTFT       - to enable use of the alias-free STFT library 
SAF_ENABLE_CDF4SAP      - to enable use of the covariance-domain framework module 
SAF_ENABLE_HOA          - to enable use of the higher-order Ambisonics module 
SAF_ENABLE_SH           - to enable use of the spherical harmonic domain related stuff module 
SAF_ENABLE_HRIR         - to enable use of the HRIR related stuff module 
SAF_ENABLE_VBAP         - to enable use of the VBAP module 
SAF_ENABLE_SOFA_READER  - to enable the SOFA reader (requires netcdf) 

#include "saf.h"
```

Detailed instructions regarding how to use the functions offered by each module is provided in the main header file for the respective module (e.g. "/modules/saf_sh/saf_sh.h", or  "/modules/saf_vbap/saf_vbap.h").


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

This framework is provided under the [ISC license](https://choosealicense.com/licenses/isc/). However, it also includes a modified version of the ['alias-free STFT'](https://github.com/jvilkamo/afSTFT) implementation by Juha Vilkamo (MIT license); and the ['convhull_3d'](https://github.com/leomccormack/convhull_3d) header only 3-D Convex Hull implementation by Leo McCormack (MIT license).

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
In Proceedings of the 14th Sound and Music Computing Conference, July 5--8, Espoo, Finland, At Espoo, Finland 2017 Jul 5.

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
