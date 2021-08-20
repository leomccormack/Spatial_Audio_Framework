#!/bin/bash
set -e

# Help message, if no input arguments were given
if [[ $# -eq 0 ]]; then
    cat <<EOT
A script to build and install a custom Intel IPP library tailored to SAF.

Usage:"

  sudo ./$(basename $0) build_type [ipp_builder_dir]

  build_type must be "sequential" or "threaded"
  (optional) ipp_builder_dir the path for the IPP "custom_library_tool_python" folder

  Examples:
      sudo ./$(basename $0) sequential
      sudo ./$(basename $0) threaded
      sudo ./$(basename $0) sequential /opt/intel/oneapi/ipp/latest/tools/custom_library_tool_python
      sudo ./$(basename $0) sequential /opt/intel/compilers_and_libraries/linux/ipp/tools/custom_library_tool_python
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

# Check if IPP build directory is provided
if [ ! -z "$2" ]; then
    ipp_builder_dir=${2}
else
    ipp_builder_dir="/opt/intel/oneapi/ipp/latest/tools/custom_library_tool_python"
    echo "Using default IPP builder path (${ipp_builder_dir})"
fi

# Define output and MKL build directories
if [[ "$OSTYPE" == "linux"* ]]; then
    if ! [ -d ${ipp_builder_dir} ]; then
        echo "Error: Intel IPP not installed"
        exit 1
    fi

elif [[ "$OSTYPE" == "darwin"* ]]; then
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

# copy saf_ipp_list
#cp "${parent_path}/saf_ipp_list.txt" ${ipp_builder_dir}

echo "Configuration: ${build_type}"
echo "IPP list: ${parent_path}/saf_ipp_list.txt"
#source "/opt/intel/oneapi/setvars.sh" -arch intel64 --force

# Create the build script
if [[ ${build_type} == "sequential"* ]]; then
    (cd ${ipp_builder_dir} && python3 main.py -c -g -n saf_ipp_custom -p "${parent_path}" -ff "${parent_path}/saf_ipp_list.txt" -arch=intel64)
elif [[ ${build_type} == "threaded"* ]]; then
    (cd ${ipp_builder_dir} && python3 main.py -c -g -n saf_ipp_custom -p "${parent_path}" -ff "${parent_path}/saf_ipp_list.txt" -arch=intel64 -mt)
fi

# Run the build script
./build_saf_ipp_custom_intel64.sh

# Define output dir
output_dir="/usr/local/lib/"

# copy library
if [[ "$OSTYPE" == "linux"* ]]; then
    (mv "libsaf_ipp_custom.so" ${output_dir})
elif [[ "$OSTYPE" == "darwin"* ]]; then
    (mv "libsaf_ipp_custom.dylib" ${output_dir})
fi

echo "Installed libsaf_ipp_custom into ${output_dir}"

# clean-up
if [[ "$OSTYPE" == "linux"* ]]; then
    rm export.def

elif [[ "$OSTYPE" == "darwin"* ]]; then
    rm "export.lib-export"
fi
rm main.c
rm main.obj

set +e
