# safmex

**extras/safmex** is a collection of MEX wrappers which allow SAF functions to be called within MATLAB. 

Even for functions that have already been implemented in MATLAB, swapping them out for these wrappers has proved beneficial simply for the faster computations afforded by them. These wrappers also serve to bring SAF functionality to MATLAB, which is not already found in existing MATLAB libraries. 

## Getting started

The easiest way to build the MEX wrappers is by enabling the **SAF_BUILD_EXTRAS** flag and building them via CMake.

Unix users:
```
cmake -S . -B build -DSAF_BUILD_EXTRAS=1 -DSAF_ENABLE_TRACKER_MODULE=1
cd build && make -j6
```

Windows users:
```
cmake -S . -B build -G "Visual Studio 15 Win64" -DSAF_BUILD_EXTRAS=1 -DSAF_ENABLE_TRACKER_MODULE=1
cd build
msbuild ALL_BUILD.vcxproj /p:Configuration=Release /m
```

You can install/copy the mex files into this folder with
```
make install
```

Alternatively, the **BUILD_SAFMEX.m** script may be used to build the wrappers. However, note that since the SAF framework must be built via CMake anyway (in order to run this script), you might as well just use CMake...

## Folder structure

**extras/safmex/** comprises the following:
* A script, **BUILD_SAFMEX.m**, to build all of the safmex wrappers via MATLAB.
* **CMakeLists.txt** for building the safmex wrappers via CMake.
* A folder, **safmex_tests**, containing unit tests for the safmex wrappers. 
* '.c' files which are the actual MEX wrappers written in C.
* '.m' files (with the same name as their '.c' counterparts) which serve as documentation for the wrappers. In the MATLAB command window, simply call: "help safmex_xxx" to view its help information.


## Dependencies

You will need MATLAB version 2017b or higher installed on your system. The remaining dependencies are then the same as the requirements when building the SAF framework.

Note that a few of the unit tests involve comparing the output of the safmex wrapper with that of an existing MATLAB implementation of the same function. Therefore, the following MATLAB libraries are required for running these unit tests:
* [Spherical-Harmonic-Transform](https://github.com/polarch/Spherical-Harmonic-Transform)
* [Higher-Order-Ambisonics](https://github.com/polarch/Higher-Order-Ambisonics)
* [Vector-Base-Amplitude-Panning](https://github.com/polarch/Vector-Base-Amplitude-Panning)

## Contributing

Contributions are very much welcomed and encouraged. Please feel free to add more wrappers!
