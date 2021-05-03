#!/bin/bash
set -e

# Help message, if no input arguments were given
if [[ $# -eq 0 ]]; then
    cat <<EOT
A script to build and install a custom Intel MKL library tailored to SAF.

Usage:"

  sudo ./$(basename $0) build_type

  build_type must be "sequential" or "threaded".

  Examples:
      sudo ./$(basename $0) sequential
      sudo ./$(basename $0) threaded
EOT
    exit 1
fi

# Check build_type argument is valid
if [ "$1" == "sequential" ]; then
    build_type=${1}
elif [ "$1" == "threaded" ]; then
    build_type=${1}
else
    echo "Error: the build_type argument must be \"sequential\" or \"threaded\"."
    exit 1
fi
shift

# Define output and MKL build directories
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    output_dir="/usr/lib/"
    mkl_builder_dir="~/intel/compilers_and_libraries/linux/mkl/tools/builder"
    if ! [ -d ${mkl_builder_dir} ]; then
        mkl_builder_dir="/opt/intel/compilers_and_libraries/linux/mkl/tools/builder"
    fi
    if ! [ -d ${mkl_builder_dir} ]; then
        echo "Error: Intel MKL not installed"
        exit 1
    fi
    iomp5_dir="~/intel/compilers_and_libraries/linux/lib/intel64"
    if ! [ -d ${iomp5_dir} ]; then
        iomp5_dir="/opt/intel/compilers_and_libraries/linux/lib/intel64"
    fi

elif [[ "$OSTYPE" == "darwin"* ]]; then
    output_dir="/usr/local/lib/"
    mkl_builder_dir="~/intel/compilers_and_libraries/mac/mkl/tools/builder"
    if ! [ -d ${mkl_builder_dir} ]; then
        mkl_builder_dir="/opt/intel/compilers_and_libraries/mac/mkl/tools/builder"
    fi
    if ! [ -d ${mkl_builder_dir} ]; then
        echo "Error: Intel MKL not installed"
        exit 1
    fi
    iomp5_dir="~/intel/compilers_and_libraries/mac/lib"
    if ! [ -d ${iomp5_dir} ]; then
        iomp5_dir="/opt/intel/compilers_and_libraries/mac/lib"
    fi
else
    echo "Error: unknown OS"
    exit 1
fi

# Current path
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# copy saf_mkl_list
cp ${parent_path}/saf_mkl_list ${mkl_builder_dir}

#echo ${build_type}
#echo ${parent_path}
#echo ${output_dir}
#echo ${mkl_builder_dir}

# build custom library
if [[ ${build_type} == "sequential"* ]]; then
    (cd ${mkl_builder_dir} && make intel64 interface=lp64 threading=sequential name=libsaf_mkl_custom export=saf_mkl_list)
elif [[ ${build_type} == "threaded"* ]]; then
    (cd ${mkl_builder_dir} && make intel64 interface=lp64 threading=parallel name=libsaf_mkl_custom export=saf_mkl_list)
fi

# copy library
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    (cd ${mkl_builder_dir} && cp libsaf_mkl_custom.so ${output_dir})
    if [[ ${build_type} == "threaded"* ]]; then
        (cd ${iomp5_dir} && cp libiomp5.so ${output_dir})
    fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
    (cd ${mkl_builder_dir} && cp libsaf_mkl_custom.dylib ${output_dir})
    if [[ ${build_type} == "threaded"* ]]; then
        (cd ${iomp5_dir} && cp libiomp5.dylib ${output_dir})
    fi
fi

echo "Installed libsaf_mkl_custom into ${output_dir}"

set +e
