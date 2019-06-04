# Spatial_Audio_Framework

A cross-platform Spatial Audio Framework (SAF) written in C. The framework includes functions for performing Vector-Base Amplitude Panning (VBAP), Spherical Harmonic Transforms (SHT), Beamforming, and Higher-order Ambisonics (HOA); to name a few.

![](saf.png)

## Getting Started

To use this framework in your project, first add the following directories to the header search paths:

```
Spatial_Audio_Framework/framework/include
Spatial_Audio_Framework/dependencies/.../include
```
Then add the following directory to the library search paths:

```
Spatial_Audio_Framework/dependencies/.../lib
```

The remaining instructions are tailored for individual operating systems and described below.

### For Windows users

Windows users must link against a custom [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries) library. The "dependencies/Win64/lib/**saf_mkl_custom.lib**" and "**saf_mkl_custom.dll**" files may be generated using Intel's custom dll builder. 

First copy the included "dependencies/saf_mkl_list" file to this folder:

```
C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/tools/builder
```

Then enter the following commands into "x64 Developer Command Prompt for VS.exe" (admin rights required):

```
cd /Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/tools/builder
nmake intel64 interface=lp64 threading=sequential name=saf_mkl_custom export=saf_mkl_list
```

The generated "saf_mkl_custom.dll" file should be placed in a suitable system PATH folder, for example:

```
C:/Windows/System32
```

(Optional) To enable the [SOFA](https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)) reading feature, your project must also link against the [netCDF](https://www.unidata.ucar.edu/software/netcdf/) library (including its dependencies). For convenience, the following statically built libraries are included in "dependencies/Win64/":

```
libszip.lib; libzlib.lib; libhdf5.lib; libhdf5_hl.lib; netcdf.lib;
```

### For Mac OSX users (Optional)

By default, the framework uses Apple's Accelerate library for the linear algebra speed-ups. However, Mac users may choose to instead use [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries) via the global pre-processor definition: "SAF_USE_INTEL_MKL" (often faster than Accelerate). 

The required "**saf_mkl_custom.dylib**" may be obtained by first placing the "dependencies/saf_mkl_list" file in:

```
/opt⁩/intel⁩/compilers_and_libraries/mac⁩/mkl⁩/tools⁩/builder⁩ 
```

To generate and copy the library to a system path folder, use the following commands (requires root permissions):

```
cd /opt⁩/intel⁩/compilers_and_libraries/mac⁩/mkl⁩/tools⁩/builder⁩
sudo make intel64 interface=lp64 threading=sequential name=libsaf_mkl_custom export=saf_mkl_list
sudo cp saf_mkl_custom.dylib /usr/local/lib

```

Then add the following linker flag to your project:

```
-L/usr/local/lib -lsaf_mkl_custom
```

(Optional) To enable the [SOFA](https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)) reading feature, your project must also link against the [netCDF](https://www.unidata.ucar.edu/software/netcdf/) library (including its dependencies). For convenience, the following statically built libraries are included in "dependencies/MacOSX/":

```
netcdf; hdf5; hdf5_hl; z; 
```

### For Linux users

Linux users must also link against a custom [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries) library, and your project must include the global pre-processor definition: "SAF_USE_INTEL_MKL". 

The required "**saf_mkl_custom.so**" may be obtained by first placing the "dependencies/saf_mkl_list" file in:

```
/opt⁩/intel⁩/compilers_and_libraries/linux/mkl⁩/tools⁩/builder⁩ 
```

To generate and copy this library to a system path folder, use the following commands (requires root permissions):

```
cd /opt⁩/intel⁩/compilers_and_libraries/linux/mkl⁩/tools⁩/builder⁩
sudo make intel64 interface=lp64 threading=sequential name=libsaf_mkl_custom export=saf_mkl_list
sudo cp saf_mkl_custom.so /usr/lib

```

Then add the following linker flag to your project:

```
-L/usr/lib -lsaf_mkl_custom
```

(Optional) To enable the SOFA loading feature, you must install netcdf and hdf5 on your system. For ubuntu based distros, this is simply:

```
sudo apt-get install libhdf5-dev
sudo apt-get install libnetcdf-dev libnetcdff-dev
```

Then add the following to the header search path:

```
/usr/include  
```

And finally, add this linker flag:

```
-L/usr/lib -lnetcdf
```

## Examples

Many examples are also included in the repository:
* **ambi_bin** - a binaural Ambisonic decoder with built-in rotator. It includes the following decoding approaches: least-squares (LS), spatial re-sampling (SPR), Time-alignment (TA) [1], Magnitude Least-Squares (MagLS) [2].
* **ambi_dec** - a frequency-dependent Ambisonic decoder. Including the following decoding approaches: sampling ambisonic decoder (SAD), AllRAD [3], Energy-Preserving decoder (EPAD) [4], Mode-Matching decoder (MMD).
* **ambi_drc** - a frequency-dependent dynamic range compressor (DRC) for Ambisonic signals, based on the design proposed in [5].
* **ambi_enc** - a simple Ambisonic encoder.
* **array2sh** - converts microphone array signals into spherical harmonic signals (aka Ambisonic signals), based on theoretical descriptions [6,7]. More details in [8].
* **beamformer** - a beamforming example with several different beamforming options.
* **binauraliser** - convolves input audio with interpolated HRTFs, which can be optionally loaded from a SOFA file.
* **dirass** - a sound-field visualiser based on re-assigning the energy of beamformers. This re-assignment is based on the DoA estimates extracted from spatially-localised active-intensity vectors, which are biased towards each beamformer direction [9].
* **panner** - a frequency-dependent VBAP panner [10], which allows for more consistent source loudness as a function of the room [11].
* **powermap** - sound-field visualiser using beamformers (PWD, MVDR) or sub-space methods (MUSIC).
* **rotator** - rotates spherical harmonic signals (aka Ambisonic signals) given yaw-pitch-roll angles [12].
* **sldoa** - a sound-field visualiser based on directly depicting the DoA estimates extracted from multiple spatially-localised active-intensity vectors; as proposed in [8].
* **upmix** - a (soon to be) collection of upmixing algorithms (currently only stereo to 5.x upmixing).

### GUI implementations

Many of these examples have been integrated into VST audio plug-ins using the JUCE framework and can be found [here](http://research.spa.aalto.fi/projects/sparta_vsts/).

## Authors

* **Leo McCormack** - C programmer and DSP researcher (contact: leo.mccormack@aalto.fi)
* **Symeon Delikaris-Manias** - DSP researcher
* **Archontis Politis** - DSP researcher

## License

This framework is provided under the [ISC license](https://choosealicense.com/licenses/isc/). However, it also includes a modified version of the ['alias-free STFT'](https://github.com/jvilkamo/afSTFT) implementation by Juha Vilkamo (MIT license); and the ['convhull_3d'](https://github.com/leomccormack/convhull_3d) header only 3-D Convex Hull implementation by Leo McCormack (MIT license).

## References

[1] Zaunschirm M, Schörkhuber C, Höldrich R. Binaural rendering of Ambisonic signals by head-related impulse response time alignment and a diffuseness constraint. The Journal of the Acoustical Society of America. 2018 Jun 19;143(6):3616-27.

[2] Schörkhuber C, Zaunschirm M, Höldrich R. Binaural Rendering of Ambisonic Signals via Magnitude Least Squares. InProceedings of the DAGA 2018 (Vol. 44, pp. 339-342).

[3] Zotter F, Frank M. All-round ambisonic panning and decoding. Journal of the audio engineering society. 2012 Nov 26;60(10):807-20.

[4] Zotter F, Pomberger H, Noisternig M. Energy-preserving ambisonic decoding. Acta Acustica united with Acustica. 2012 Jan 1;98(1):37-47.

[5] McCormack L, Välimäki V. FFT-based Dynamic Range Compression. InProceedings of the 14th Sound and Music Computing Conference, July 5--8, Espoo, Finland, At Espoo, Finland 2017 Jul 5.

[6] Williams EG. Fourier acoustics: sound radiation and nearfield acoustical holography. Elsevier; 1999 Jun 10.

[7] Rafaely B. Fundamentals of spherical array processing. Berlin: Springer; 2015 Feb 18.

[8] McCormack L, Delikaris-Manias S, Farina A, Pinardi D, Pulkki V. Real-Time Conversion of Sensor Array Signals into Spherical Harmonic Signals with Applications to Spatially Localized Sub-Band Sound-Field Analysis. InAudio Engineering Society Convention 144 2018 May 14. Audio Engineering Society.

[9] McCormack L, Politis A, Pulkki V. Sharpening of Angular Spectra Based on a Directional Re-assignment Approach for Ambisonic Sound-field Visualisation. InICASSP 2019-2019 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP) 2019 Apr 17 (pp. 576-580). IEEE.

[10] Pulkki V. Virtual sound source positioning using vector base amplitude panning. Journal of the audio engineering society. 1997 Jun 1;45(6):456-66.

[11] Laitinen MV, Vilkamo J, Jussila K, Politis A, Pulkki V. Gain normalization in amplitude panning as a function of frequency and room reverberance. InAudio Engineering Society Conference: 55th International Conference: Spatial Audio 2014 Aug 26. Audio Engineering Society.

[12] Ivanic J, Ruedenberg K. Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion. The Journal of Physical Chemistry A. 1998 Nov 5;102(45):9099-100.




