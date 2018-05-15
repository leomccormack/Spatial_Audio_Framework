# Spatial_Audio_Framework

A Spatial Audio Framework (SAF) written in C. The framework includes functions for Vector-Base Amplitude Panning (VBAP), Spherical Harmonic Transforms (SHT), Higher-order Ambisonics (HOA) etc.

## Getting Started

To include this framework in a project, simply add the following code

```
Spatial_Audio_Framework/framework
```

And add

```
Spatial_Audio_Framework/framework/include
```

To the header search paths.

Windows users must also install Intel's MKL, which can be freely acquired from
* [Intel MKL](https://software.intel.com/en-us/articles/free-ipsxe-tools-and-libraries)

## Examples

Several examples have also been included

```
Spatial_Audio_Framework/examples/ambi_dec
Spatial_Audio_Framework/examples/ambi_enc
Spatial_Audio_Framework/examples/array2sh
Spatial_Audio_Framework/examples/binauraliser
Spatial_Audio_Framework/examples/panner
Spatial_Audio_Framework/examples/powermap
Spatial_Audio_Framework/examples/rotator
Spatial_Audio_Framework/examples/sldoa
Spatial_Audio_Framework/examples/upmix
```

* **ambi_dec** - a frequency-dependent Ambisonic decoder (AllRAD, EPAD, MMD etc)
* **ambi_enc** - a simple Ambisonic encoder/panner
* **array2sh** - converts microphone array signals into spherical harmonic signals (aka HOA signals)
* **binauraliser** - convolves input audio with interpolated HRTFs, which can be optionally loaded from a SOFA file
* **panner** - a frequency-dependent VBAP panner
* **powermap** - sound-field visualiser using beamformers (PWD, MVDR) or sub-space methods (MUSIC)
* **rotator** - rotates spherical harmonic signals (HOA signals) given yaw-pitch-roll angles
* **sldoa** - spatially-localised direction of arrival estimator
* **upmix** - a (soon to be) collection of upmixing algorithms (currently only stereo to 5.x upmixing)

### GUI implementations

Many of these examples have been intergrated into VST audio plug-ins using the JUCE framework and can be found here:
http://research.spa.aalto.fi/projects/sparta_vsts/

## Authors

* **Leo McCormack** - C programmer and DSP researcher (contact: leo.mccormack@aalto.fi)
* **Symeon Delikaris-Manias** - DSP researcher
* **Archontis Politis** - DSP researcher

## License

This framework is licensed under the ISC License. However, it also includes a slightly modified version of the alias-free STFT implementation by Juha Vilkamo (MIT license).

