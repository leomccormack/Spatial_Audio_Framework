# Spatial_Audio_Framework

A Spatial Audio Framework (SAF) written in C. The framework includes functions for performing Vector-Base Amplitude Panning (VBAP), Spherical Harmonic Transforms (SHT), and Higher-order Ambisonics (HOA); among others.

![](saf.png)

## Getting Started

To use this framework in your project, first add the following folder and all of its contained files:

```
Spatial_Audio_Framework/framework
```

Add the following directories to the relevant header and library search paths:

```
Spatial_Audio_Framework/framework/include
Spatial_Audio_Framework/dependencies/.../include
Spatial_Audio_Framework/dependencies/.../lib
```

To enable the SOFA loading feature, your project must link against the following included static libraries:

```
netcdf; hdf5; hdf5_hl; z; (for MacOSX users)
libszip.lib; libzlib.lib; libhdf5.lib; libhdf5_hl.lib; netcdf.lib; (for Windows users)
```

### For Windows users only

Windows users must also link against an included custom [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries) library:

```
saf_mkl_custom.lib; (for Windows users)
```

The "saf_mkl_custom.dll" may be generated with Intel's custom dll builder and the included "saf_mkl_list", with the following command:

```c
nmake dllintel64 interface=lp64 threading=sequential name=saf_mkl_custom export=saf_mkl_list
```
Alternatively, this .dll is also installed alongside the software found [here](http://research.spa.aalto.fi/projects/sparta_vsts/download/).

### For Mac OSX users only (Optional)

By default, the framework uses Apple's Accelerate library for the linear algebra speed-ups. However, Mac users may also choose to use Intel's MKL instead, with the global pre-processor definition: "SAF_USE_INTEL_MKL". 

The "saf_mkl_custom.dylib" may be obtained in similar manner as in the above, by either installing the software found [here](http://research.spa.aalto.fi/projects/sparta_vsts/download/), or by generating it with the following command:

```c
make intel64 interface=lp64 threading=sequential name=saf_mkl_custom export=saf_mkl_list
```

## Examples

Several examples have also been included:
* **ambi_bin** - a binaural Ambisonic decoder with built-in rotator
* **ambi_dec** - a frequency-dependent Ambisonic decoder (AllRAD, EPAD, MMD etc)
* **ambi_drc** - a frequency-dependent dynamic range compressor (DRC) for spherical harmonic signals (aka Ambisonic signals)
* **ambi_enc** - a simple Ambisonic encoder/panner
* **array2sh** - converts microphone array signals into spherical harmonic signals (aka Ambisonic signals)
* **binauraliser** - convolves input audio with interpolated HRTFs, which can be optionally loaded from a SOFA file
* **panner** - a frequency-dependent VBAP panner
* **powermap** - sound-field visualiser using beamformers (PWD, MVDR) or sub-space methods (MUSIC)
* **rotator** - rotates spherical harmonic signals (aka Ambisonic signals) given yaw-pitch-roll angles
* **sldoa** - spatially-localised direction of arrival estimator
* **upmix** - a (soon to be) collection of upmixing algorithms (currently only stereo to 5.x upmixing)

### GUI implementations

Many of these examples have been integrated into VST audio plug-ins using the JUCE framework and can be found [here](http://research.spa.aalto.fi/projects/sparta_vsts/).

## Authors

* **Leo McCormack** - C programmer and DSP researcher (contact: leo.mccormack@aalto.fi)
* **Symeon Delikaris-Manias** - DSP researcher
* **Archontis Politis** - DSP researcher

## License

This framework is licensed under the [ISC license](https://choosealicense.com/licenses/isc/). However, it also includes the 'alias-free STFT' implementation by Juha Vilkamo (MIT license), which can be found [here](https://github.com/jvilkamo/afSTFT); and the 'convhull_3d' header only Convex Hull implementation by Leo McCormack (MIT license), which can be found [here](https://github.com/leomccormack/convhull_3d).

