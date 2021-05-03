#!/bin/bash
set -e

# Help message, if no input arguments were given
if [[ $# -eq 0 ]]; then
    cat <<EOT
A script to build and install a custom Intel IPP library tailored to SAF.

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
    ipp_builder_dir="~/intel/compilers_and_libraries/linux/ipp/tools/custom_library_tool_python"
    if ! [ -d ${ipp_builder_dir} ]; then
        ipp_builder_dir="/opt/intel/compilers_and_libraries/linux/ipp/tools/custom_library_tool_python"
    fi
    if ! [ -d ${ipp_builder_dir} ]; then
        echo "Error: Intel IPP not installed"
        exit 1
    fi

elif [[ "$OSTYPE" == "darwin"* ]]; then
    output_dir="/usr/local/lib/"
    ipp_builder_dir="~/intel/compilers_and_libraries/mac/ipp/tools/custom_library_tool_python"
    if ! [ -d ${ipp_builder_dir} ]; then
        ipp_builder_dir="/opt/intel/compilers_and_libraries/mac/ipp/tools/custom_library_tool_python"
    fi
    if ! [ -d ${ipp_builder_dir} ]; then
        echo "Error: Intel IPP not installed"
        exit 1
    fi
else
    echo "Error: unknown OS"
    exit 1
fi

# Current path
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# copy saf_mkl_list
cp ${parent_path}/saf_ipp_list ${ipp_builder_dir}

# build and copy custom library into output dir
if [[ ${build_type} == "sequential"* ]]; then
    (cd ${ipp_builder_dir} && python main.py -c -g
    –n saf_ipp_custom
    -p “${output_dir}”
    -ff “saf_ipp_list”
    -d avx512bw
    -arch=intel64)
elif [[ ${build_type} == "threaded"* ]]; then
    (cd ${ipp_builder_dir} && python main.py -c -g
    –n saf_ipp_custom
    -p “${output_dir}”
    -ff “saf_ipp_list”
    -d avx512bw
    -arch=intel64 -mt )
fi

echo "Installed libsaf_ipp_custom into ${output_dir}"

set +e
