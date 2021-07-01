#!/bin/bash
set -e

# Help message, if no input arguments were given
if [[ $# -eq 0 ]]; then
    cat <<EOT
A script to build and install a custom Intel MKL library tailored to SAF.

Usage:"

  sudo ./$(basename $0) build_type [mkl_interface] [mkl_builder_dir]

  build_type must be either "sequential" or "threaded",
  (optional) mkl_interface must be either "lp64" (default) or "ilp64"
  (optional) mkl_builder_dir the path for the MKL "builder" folder

  Examples:
      sudo ./$(basename $0) sequential
      sudo ./$(basename $0) threaded
      sudo ./$(basename $0) sequential ilp64 /opt/intel/oneapi/mkl/latest/tools/builder
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

# Check if MKL interface is specified (and if so, is valid)
if [ ! -z "$2" ]; then
    if [ "$2" == "lp64" ]; then
        mkl_interface=${2}
    elif [ "$2" == "ilp64" ]; then
        mkl_interface=${2}
    else
        echo "Error: the mkl_interface argument must be \"lp64\" or \"ilp64\"."
        exit 1
    fi
else
    echo "Using default MKL interface configuration (lp64)"
    mkl_interface="lp64"
fi

# Check if MKL build directory is provided
if [ ! -z "$3" ]; then
    mkl_builder_dir=${3}
else
    mkl_builder_dir="/opt/intel/oneapi/mkl/latest/tools/builder"
    echo "Using default MKL builder path (${mkl_builder_dir})"
fi

# Define output dir
output_dir="/usr/local/lib/"

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
    (cd ${mkl_builder_dir} && make libintel64 interface=${mkl_interface} threading=sequential name="libsaf_mkl_custom_${mkl_interface}" export=saf_mkl_list)
elif [[ ${build_type} == "threaded" ]]; then
    (cd ${mkl_builder_dir} && make libintel64 interface=${mkl_interface} threading=parallel name="libsaf_mkl_custom_${mkl_interface}" export=saf_mkl_list)
fi

# copy library
if [[ "$OSTYPE" == "linux"* ]]; then
    (cd ${mkl_builder_dir} && mv "libsaf_mkl_custom_${mkl_interface}.so" ${output_dir})
elif [[ "$OSTYPE" == "darwin"* ]]; then
    (cd ${mkl_builder_dir} && mv "libsaf_mkl_custom_${mkl_interface}.dylib" ${output_dir})
fi

echo "Installed libsaf_mkl_custom_${mkl_interface} into ${output_dir}"

set +e
