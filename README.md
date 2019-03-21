# Spatial_Audio_Framework

A Spatial Audio Framework (SAF) written in C. The framework includes functions for performing Vector-Base Amplitude Panning (VBAP), Spherical Harmonic Transforms (SHT), and Higher-order Ambisonics (HOA); among others.

![](saf.png)

## Getting Started

To use this framework, add the following directories to the relevant header and library search paths:

```
Spatial_Audio_Framework/framework/include
Spatial_Audio_Framework/dependencies/.../include
Spatial_Audio_Framework/dependencies/.../lib
```

To enable the SOFA loading feature, your project must also link against the following included static libraries:

```
netcdf; hdf5; hdf5_hl; z; (for MacOSX users)
libszip.lib; libzlib.lib; libhdf5.lib; libhdf5_hl.lib; netcdf.lib; (for Windows users)
```

### For Windows users only

Windows users must also link against a custom [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries) library. The "dependencies/Win64/lib/saf_mkl_custom.lib" and "saf_mkl_custom.dll" files may be generated using Intel's custom dll builder, by first copying the included "dependencies/saf_mkl_list" file into this folder:

```
C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/tools/builder
```

and enter the following commands into "x64 Developer Command Prompt for VS.exe" (admin rights required):

```
cd /Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/tools/builder
nmake intel64 interface=lp64 threading=sequential name=saf_mkl_custom export=saf_mkl_list
```

The generated "saf_mkl_custom.dll" file should be placed in a suitable system PATH folder, for example:

```
C:/Windows/System32
```

Alternatively, this .dll is installed in a system PATH folder alongside the software found [here](http://research.spa.aalto.fi/projects/sparta_vsts/download/); however, this may not always be up-to-date. 

### For Mac OSX users only (Optional)

By default, the framework uses Apple's Accelerate library for the linear algebra speed-ups. However, Mac users may choose to instead use Intel's MKL via the global pre-processor definition: "SAF_USE_INTEL_MKL" (often around 20-40% faster than Accelerate). 

The required "saf_mkl_custom.dylib" may be obtained by installing the software found [here](http://research.spa.aalto.fi/projects/sparta_vsts/download/), or by first placing the "dependencies/saf_mkl_list" file in:

```
/opt⁩/intel⁩/compilers_and_libraries/mac⁩/mkl⁩/tools⁩/builder⁩ 
```

and generating it using the following commands via the terminal (admin rights required):

```
cd /opt⁩/intel⁩/compilers_and_libraries/mac⁩/mkl⁩/tools⁩/builder⁩
make intel64 interface=lp64 threading=sequential name=libsaf_mkl_custom export=saf_mkl_list
```

The generated "saf_mkl_custom.dylib" file should be placed in a suitable system path folder, for example:

```
/usr/local/lib
```

## Examples

Several examples are also included in the repository:
* **ambi_bin** - a binaural Ambisonic decoder with built-in rotator
* **ambi_dec** - a frequency-dependent Ambisonic decoder (AllRAD, EPAD, MMD etc)
* **ambi_drc** - a frequency-dependent dynamic range compressor (DRC) for spherical harmonic signals (aka Ambisonic signals)
* **ambi_enc** - a simple Ambisonic encoder
* **array2sh** - converts microphone array signals into spherical harmonic signals (aka Ambisonic signals)
* **binauraliser** - convolves input audio with interpolated HRTFs, which can be optionally loaded from a SOFA file
* **dirass** - a sound-field visualiser based on re-assigning the energy of beamformers. This re-assignment is based on the DoA estimates extracted from spatially-localised active-intensity vectors, which are biased towards each beamformer direction.
* **panner** - a frequency-dependent VBAP panner
* **powermap** - sound-field visualiser using beamformers (PWD, MVDR) or sub-space methods (MUSIC)
* **rotator** - rotates spherical harmonic signals (aka Ambisonic signals) given yaw-pitch-roll angles
* **sldoa** - a sound-field visualiser based on directly depicting the DoA estimates extracted from multiple spatially-localised active-intensity vectors
* **upmix** - a (soon to be) collection of upmixing algorithms (currently only stereo to 5.x upmixing)

### GUI implementations

Many of these examples have been integrated into VST audio plug-ins using the JUCE framework and can be found [here](http://research.spa.aalto.fi/projects/sparta_vsts/).

## Authors

* **Leo McCormack** - C programmer and DSP researcher (contact: leo.mccormack@aalto.fi)
* **Symeon Delikaris-Manias** - DSP researcher
* **Archontis Politis** - DSP researcher

## License

This framework is provided under the [ISC license](https://choosealicense.com/licenses/isc/). However, it also includes a modified version of the ['alias-free STFT'](https://github.com/jvilkamo/afSTFT) implementation by Juha Vilkamo (MIT license); and the ['convhull_3d'](https://github.com/leomccormack/convhull_3d) header only 3-D Convex Hull implementation by Leo McCormack (MIT license).

