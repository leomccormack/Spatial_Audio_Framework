# Framework Structure

<img src="framework.svg"> 

## Include

The framework's master include header is:

```c
#include "saf.h"
```

Optionally, if you would like to carry over the CBLAS/LAPACK routines (and the functions of any other SAF external library) to use in your own project, then also include the following header:
```c
#include "saf_externals.h"  
```

## Modules

The framework comprises the following core modules (**ISC**):
* **saf_hoa** - a collection of higher-order Ambisonics binaural and loudspeaker decoders.
* **saf_sh** - spherical harmonic and spherical array processing related functions.
* **saf_vbap** - Vector-base Amplitude Panning (VBAP) functions.
* **saf_cdf4sap** - Covariance Domain Framework for Spatial Audio Processing (CDF4SAP).
* **saf_hrir** - HRIR/HRTF related functions (estimating ITDs, HRTF interpolation, diffuse-field EQ etc.).
* **saf_reverb** - a collection of reverbs and room simulation algorithms.
* **saf_utilities** - a collection of useful utility functions and cross-platform wrappers.

The framework also includes the following optional modules:
* **saf_sofa_reader** - a simple SOFA file reader (**ISC**).
* **saf_tracker** - a particle-filtering based tracker (**GPLv2**).

To enable optional framework modules, simply add the relevant pre-processor definition:
```
SAF_ENABLE_SOFA_READER_MODULE  # to enable saf_sofa_reader
SAF_ENABLE_TRACKER_MODULE      # to enable saf_tracker
```

## Resources

The **saf_utilities** module also inherits and offers use of the following third-party **resources**:
* [**afSTFT**](https://github.com/jvilkamo/afSTFT) - a slightly modified version of the alias-free Short-Time Fourier Transform (afSTFT) filterbank (**MIT**).
* [**md_malloc**](https://github.com/leomccormack/md_malloc) - a utility for allocating contiguous multi-dimensional arrays (**MIT**).
* [**convhull_3d**](https://github.com/leomccormack/convhull_3d) - for building 3-D convex hulls (**MIT**).
* [**kissFFT**](https://github.com/mborgerding/kissfft) - the default FFT implementation for when no optimised alternative is linked (**BSD-3-Clause**).

## Dependencies

The only hard dependency for SAF is a library (or a combination of libraries) which supports the [CBLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms#Implementations) and [LAPACK](https://en.wikipedia.org/wiki/LAPACK) standards. 

For more details, refer to: [**PERFORMANCE_LIBRARY_INSTRUCTIONS.md**](PERFORMANCE_LIBRARY_INSTRUCTIONS.md)

### Optional externals

In order to use the optional built-in [SOFA](https://www.sofaconventions.org/mediawiki/index.php/SOFA_(Spatially_Oriented_Format_for_Acoustics)) reader module (**saf_sofa_reader**), your project must also link against the [netCDF](https://www.unidata.ucar.edu/software/netcdf/) library (including its dependencies). 

For more details, refer to: [**SOFA_READER_MODULE_DEPENDENCIES.md**](SOFA_READER_MODULE_DEPENDENCIES.md)

Optionally, Intel's [Integrated Performance Primitives (IPP)](https://software.intel.com/content/www/us/en/develop/tools/integrated-performance-primitives.html) may be employed for the FFT with the following definition:
```
SAF_USE_INTEL_IPP            
```

## Contributing

Contributions are very much welcomed and encouraged. Please feel free to make suggestions, pull-requests, or get in touch (via leo.mccormack(at)aalto.fi or github "issues") if you would like to discuss additions to the framework. These additions can come in many forms; including:
* bug fixes or optimisations to existing code
* adding new functionality (which falls under the scope of an existing module) to an existing module
* adding a new utility to the **saf_utilities** module or an entirely new module (see below)
* using the framework to create new **examples**, **extras**, or unit tests
* documentation improvements

## Module structure

All modules should reside in the **framework/modules** folder, named in lower-case and with the **saf_** prefix. The main include header for the module should also have the same name as this folder; e.g.:

```
framework/modules/saf_new_module/saf_new_module.h
framework/modules/saf_sh/saf_sh.h
```

The rest of the module can be organised as you wish. However, one suggestion is to use the following structure:

```
saf_new_module.h           # Main header comprising functions/structs/enums that the user can interface with
saf_new_module_internal.h  # Internal header comprising functions/structs/enums that only the module uses
saf_new_module.c           # Source code for the functions declared in the main header
saf_new_module_internal.c  # Source code for the functions declared in the internal header
```

### Adding a module to the framework

All core modules are included by the **framework/include/saf.h** master include header using its relative path, and prefaced with a definition of the same name (except upper case). Doxygen documentation should then be added for this definition, which: describes the module, its license, and any dependencies. For example:
```c
/* *  
 * SAF Module: HOA
 *
 * A collection of higher-order Ambisonics related functions; many of which are
 * derived from the Matlab library found in [1] (BSD-3-Clause License).
 *
 * ## Dependencies
 *   saf_utilities.h, saf_vbap.h, saf_sh.h
 *
 * @see [1] https://github.com/polarch/Higher-Order-Ambisonics
 *
 * License: ISC
 */
#define SAF_HOA_MODULE
#include "../modules/saf_hoa/saf_hoa.h"
```

Optional modules are also included by **framework/include/saf.h** in a similar manner as with the core modules, except they should be guarded by a **SAF_ENABLE_X_MODULE** definition and instructions should be provided regarding how to enable it, along with a link to instructions on how to build/link any dependencies. For example:
```c
/* *
 * SAF Module: SOFA_Reader
 *
 * A simple SOFA file reader that returns only the bare minimum needed to
 * load HRIR data.
 *
 * ## Enable instructions
 *   Add this pre-processor definition to your project:
 *       SAF_ENABLE_SOFA_READER_MODULE
 *   and ensure that the netcdf library is also linked to your project. More
 *   information can be found in:
 *       docs/SOFA_READER_MODULE_DEPENDENCIES.md
 * ## Dependencies
 *   saf_utilities.h, saf_hrir.h, netcdf
 *
 * License: ISC
 */
#define SAF_SOFA_READER_MODULE
#ifdef  SAF_ENABLE_SOFA_READER_MODULE
# include "../modules/saf_sofa_reader/saf_sofa_reader.h"
#endif /* SAF_ENABLE_SOFA_READER_MODULE */
```

### Coding style

All code should be written in C and compile under C99 compliant compilers and the ancient MSVC compiler. Note that complex number wrappers (**saf_utilities/saf_utility_complex.h**) have been provided to make this latter requirement a bit easier, but we can also help port C99-only code to be MSVC compliant if needed. 

As for coding-style, there are no strict requirements, only suggestions. In general, taking an existing SAF module and using it as a basis for developing a new module is probably the easiest way to go. However, suggested coding styles are also listed below:
* Indentations should be 4 spaces in length - no tabs please!
* Comments are preferred with the /* comment */ syntax, rather than //
* Comments and function declarations in header files should respect an 80 character margin
* Doxygen syntax should be used to document all functions/structures open to the user (i.e., the main header file of the module). 
* Functions/structures open to the user should be particularly verbose and make reference to the appropriate publications.
* Using Doxygen syntax to also document internal functions is preferred, but we'll take any documentation over no documentation.
* Breaking up source code to describe what is going on with comments will make it easier to find and fix bugs in the future.
* Open curly brackets should come after "for" loops and "if" statements on the same line.
* Open curly brackets should be omitted for "for" loops and "if" statements that comprise only one line of execution.
* Structure handles are preferred to be of type "void*" if the structure is only visible internally to the module (i.e., defined in the module's internal header).
* Structure handles should be of their declaration type if the structure is visible to the user (i.e., defined in the main module header).
* Using the SAF-native "float_complex" and "double_complex" structures for handling complex numbers is favoured for all user interface functions. 
* Writing appropriate unit tests (which also demonstrate the functionality of your added code) to accompany your module is highly recommended. See the **test** folder for more details.

### Dependencies

If the new module requires third-party dependencies, (which are not already included in the framework), then there are a few options to add them:
* If the dependency is permissively licensed and is a relatively small library (comprising 1-3 source files), then these may be added to the module itself or, if the dependency may be useful to other modules, it can be added to the **framework/resources** folder and inherited by the **saf_utilities** module.  
* If the dependency is **not** permissively licensed and is a relatively small library (comprising 1-3 source files), then these may be added to the module itself and included only internally by the module.
* If the dependency is large, then it should be included via the **framework/include/saf_externals.h** header, and detailed instructions regarding how to build/link these libraries should be added to the **docs** folder. 

Note: your module can only be added as a **core** module if it is free from large dependencies and comprises only permissively licensed code (see below).

### Licensing

All core modules must be released under the [ISC License](https://choosealicense.com/licenses/isc/) and only use third-party code that is provided under similar permissive licensing terms (ISC, MIT, BSD etc.). 

Optional modules may be released under alternative licenses. If you would prefer (or are required) to release your module under a copyleft license, then you may do so. For example, the **saf_tracker** module is provided under the terms of the copyleft [GNU GPLv2 License](https://choosealicense.com/licenses/gpl-2.0/). However, please discuss with us if you wish to use a license other than ISC or GPLv2.
