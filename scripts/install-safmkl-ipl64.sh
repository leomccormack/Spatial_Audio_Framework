#!/bin/bash
set -e

# Help message, if no input arguments were given
if [[ $# -eq 0 ]]; then
    cat <<EOT
A script to build and install a custom Intel MKL library tailored to SAF.

Usage:"

  sudo ./$(basename $0) build_type [mkl_builder_dir]

  build_type must be "sequential" or "threaded", with optional "mkl_builder_dir".

  Examples:
      sudo ./$(basename $0) sequential
      sudo ./$(basename $0) threaded
      sudo ./$(basename $0) sequential /opt/intel/oneapi/mkl/latest/tools/builder
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
# Check if MKL build directory is provided
if [ ! -z "$2" ]; then
    mkl_builder_dir=${2}
else
    echo "Using default MKL builder path"
    mkl_builder_dir="/opt/intel/oneapi/mkl/latest/tools/builder"
fi

# Define output dir
output_dir="/usr/local/lib/"

# MKL Interface (lp64 / ilp64)
mkl_interface="ilp64"

# Define output and MKL build directories
if [[ "$OSTYPE" == "linux"* ]]; then
    if ! [ -d ${mkl_builder_dir} ]; then
        echo "Error: Intel MKL not installed"
        exit 1
    fi

elif [[ "$OSTYPE" == "darwin"* ]]; then
    if ! [ -d ${mkl_builder_dir} ]; then
        echo "Error: Intel MKL not installed"
        exit 1
    fi
else
    echo "Error: unknown OS"
    exit 1
fi

# Current path
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# copy saf_mkl_list
cp ${parent_path}/saf_mkl_list ${mkl_builder_dir}

echo "Configuration: Builder ${mkl_builder_dir}, ${build_type} ${mkl_interface}"
#echo ${parent_path}
#echo ${output_dir}

# build custom library
if [[ ${build_type} == "sequential" ]]; then
    (cd ${mkl_builder_dir} && make libintel64 interface=${mkl_interface} threading=sequential name=libsaf_mkl_custom export=saf_mkl_list)
elif [[ ${build_type} == "threaded" ]]; then
    (cd ${mkl_builder_dir} && make libintel64 interface=${mkl_interface} threading=parallel name=libsaf_mkl_custom export=saf_mkl_list)
fi

# copy library
if [[ "$OSTYPE" == "linux"* ]]; then
    (cd ${mkl_builder_dir} && cp libsaf_mkl_custom.so ${output_dir})
elif [[ "$OSTYPE" == "darwin"* ]]; then
    (cd ${mkl_builder_dir} && cp libsaf_mkl_custom.dylib ${output_dir})
fi

echo "Installed libsaf_mkl_custom into ${output_dir}"

set +e
