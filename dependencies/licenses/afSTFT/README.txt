Juha Vilkamo
Software developed at Aalto Univesity, present e-mail: juha.vilkamo@nokia.com 


Alias-free short-time Fourier transform – a robust time-frequency transform for audio processing

- use compile_afSTFT_mex.m to compile a Matlab MEX file, or use the precompiled MEX files.
- find Matlab usage example in afSTFT_usage_instruction.m


*** What it is ***

The software does processing similar to the short-time Fourier transform (STFT), but applies the longer and more overlapping analysis/synthesis window to make the time-frequency processing fundamentally robust for aliasing. The difference is relevant for very active time-frequency processing techniques, such as DirAC or optimized adaptive beamformers.  

There are also use cases in which the aliasing effects are not an issue, such as smooth equalizers, spectral analyzers and such, in which this library should perform similarly to a traditional STFT.

The design principles of this filter bank are not fundamentally different than those of the complex-modulated QMF bank found in many MPEG systems. The main difference is that the MPEG filter bank employs a 0.5 band shift to preserve compatibility to real-valued transforms. 

In this afSTFT implementation, the band center frequencies are according to the familiar STFT processing (unless the “hybrid” mode is applied). For most use cases, you can directly replace any existing STFT processing with this afSTFT library.

The software can be applied in C-programs and also in Matlab.


*** The package contains ***

The package contains:

-- afSTFT.mexmaci64 & afSTFT.mexw64: The pre-compiled Matlab MEX files for Mac and Windows

-- compile_afSTFT_mex.m: A matlab file that compiles and tests the MEX-file (compiler needs to be installed) —> Although not tested, this should work on linux too

-- afSTFT_mex.c: The mex-C-file that is compiled by the above.

-- afSTFTlib.c & .h, vecTools.c & .h, afSTFT_protofilter.h, fft4g.c: The main files for the library, that can be compiled as part of any project.

-- afSTFT_usage_instruction.m: Usage instruction for the MEX file in Matlab. Not maybe super clear presently, but should be understandable for a person in the field.

-- afAnalyze.m, afSynthesize.m: For quick single-shot usage of the mex file. Continuous-time programs apply afSTFT differently, check the afSTFT_usage_instruction.m


*** Documentation ***

A detailed documentation of the software is tbd.


