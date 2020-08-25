%
% Copyright 2020 Leo McCormack 
%
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
% REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
% INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
% LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
% OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
% PERFORMANCE OF THIS SOFTWARE.
%
clear all, close all %#ok

% 1. compile SAF via CMake prior to running this script
% 2. configure user options below as needed
% If you encounter any problems during compilation, please refer to the FAQ 
% at the bottom of this script before submitting a bug report.

% NOTE: This stuff is very much WIP. You may need to change the
% configurations below. If you would like to contribute, then please go 
% ahead and do so, or get in touch :-)


%% User Options 
saf_perf_lib = []; % leave empty to scan cache, 'SAF_USE_INTEL_MKL', 'SAF_USE_APPLE_ACCELERATE', 'SAF_USE_OPEN_BLAS_AND_LAPACKE'
run_tests = 1;     % 0: disabled, 1: enabled
cmake_build_folder = '../../build/'; % path to CMake build folder
outputFolder = './';


%% Pull settings from CMakeCache.txt if not defined in user options 
if isempty(saf_perf_lib)
    text = fileread([cmake_build_folder 'CMakeCache.txt']);
    txt = textscan(text,'%s','delimiter','\n'); 
    cac = [txt{:}];

    idx =[];
    if isempty(idx)
        idx = find(strcmp( cac, 'SAF_PERFORMANCE_LIB:STRING=SAF_USE_INTEL_MKL' ));
        if isempty(idx), idx = find(strcmp( cac, 'SAF_PERFORMANCE_LIB:UNINITIALIZED=SAF_USE_INTEL_MKL' )); end
        if ~isempty(idx), saf_perf_lib = 'SAF_USE_INTEL_MKL'; end 
    end
    if isempty(idx)
        idx = find(strcmp( cac, 'SAF_PERFORMANCE_LIB:STRING=SAF_USE_APPLE_ACCELERATE' ));
        if isempty(idx), idx = find(strcmp( cac, 'SAF_PERFORMANCE_LIB:UNINITIALIZED=SAF_USE_APPLE_ACCELERATE' )); end
        if ~isempty(idx), saf_perf_lib = 'SAF_USE_APPLE_ACCELERATE'; end 
    end
    if isempty(idx)
        idx = find(strcmp( cac, 'SAF_PERFORMANCE_LIB:UNINITIALIZED=SAF_USE_OPEN_BLAS_AND_LAPACKE' ));   %SAF_PERFORMANCE_LIB:UNINITIALIZED=SAF_USE_OPEN_BLAS_AND_LAPACKE
        if isempty(idx), idx = find(strcmp( cac, 'SAF_PERFORMANCE_LIB:UNINITIALIZED=SAF_USE_OPEN_BLAS_AND_LAPACKE' )); end
        if ~isempty(idx), saf_perf_lib = 'SAF_USE_OPEN_BLAS_AND_LAPACKE'; end 
    end  
end


%% Configure
header_paths = {'../../framework/include'};

libs = {};
defines = {saf_perf_lib};
defines{end+1} = 'SAF_ENABLE_TRACKER_MODULE';

% Platform specific configurations
if ismac 
    lib_paths = {[ cmake_build_folder 'framework' ]};
    libs{end+1} = 'saf';
    
    % Error checking
    if ~exist([cmake_build_folder 'framework/libsaf.a'], 'file')
        if ~exist([cmake_build_folder 'framework/libsaf.dylib'], 'file'), error('libsaf.a/dylib is missing, run CMake and build saf first'), end
    end
    
    % extra header+lib paths
    switch saf_perf_lib
        case 'SAF_USE_APPLE_ACCELERATE' 
            % none needed
        case 'SAF_USE_INTEL_MKL'
            header_paths{end+1} = '/opt/intel/compilers_and_libraries/mac/mkl/include';
            lib_paths{end+1} = '/usr/local/lib';
            libs{end+1} = 'saf_mkl_custom';
            saf_mkl = [lib_paths{end} '/lib' libs{end} '.dylib' ];
            if ~exist(saf_mkl, 'file'), error([ saf_mkl ' is missing']), end
        otherwise
            assert(0); % currently unsupported, please contribute!
    end 
    
elseif isunix
    lib_paths = {[ cmake_build_folder 'framework' ]};
    libs{end+1} = 'saf';
    
    % Error checking
    if ~exist([cmake_build_folder 'framework/libsaf.a'], 'file')
        if ~exist([cmake_build_folder 'framework/libsaf.so'], 'file'), error('libsaf.a/so is missing, run CMake and build saf first'), end
    end
    
    % extra header+lib paths
    switch saf_perf_lib
        case 'SAF_USE_INTEL_MKL'
            assert(0); % At least for me, MKL ver.2020.1.217 is fucked on ubuntu 20.04 LTS and crashes matlab, no idea why. Needs investigating...
            
            % If you remove this assert and proceed anyway, please let me
            % know how it goes.
            
            %header_paths{end+1} = '/opt/intel/compilers_and_libraries/linux/mkl/include';
            header_paths{end+1} = '~/intel/compilers_and_libraries/linux/mkl/include';
            lib_paths{end+1} = '/usr/lib';
            libs{end+1} = 'saf_mkl_custom';
            saf_mkl = [lib_paths{end} '/lib' libs{end} '.so' ];
            if ~exist(saf_mkl, 'file'), error([ saf_mkl ' is missing']), end
        case 'SAF_USE_OPEN_BLAS_AND_LAPACKE' 
            header_paths{end+1} = '/usr/include/x86_64-linux-gnu';
            lib_paths{end+1} = '/usr/lib/x86_64-linux-gnu';
            libs{end+1} = 'openblas'; 
            openblas = [lib_paths{end} '/lib' libs{end} '.so' ];
            if ~exist(openblas, 'file'), error([ openblas ' is missing']), end
            libs{end+1} = 'lapacke'; 
            lapacke = [lib_paths{end} '/lib' libs{end} '.so' ];
            if ~exist(lapacke, 'file'), error([ lapacke ' is missing']), end
        otherwise
            assert(0); % currently unsupported, please contribute!
    end 
elseif ispc
    lib_paths = {[ cmake_build_folder 'framework/Release' ]};
    libs{end+1} = 'saf';
    
    switch saf_perf_lib
        case 'SAF_USE_INTEL_MKL'
            header_paths{end+1} = 'C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include';
            lib_paths{end+1} = '../../dependencies/Win64/lib';
            libs{end+1} = 'saf_mkl_custom';
            saf_mkl = [lib_paths{end} '/' libs{end} '.lib' ];
            if ~exist(saf_mkl, 'file'), error([ saf_mkl ' is missing']), end
        otherwise
            assert(0); % currently unsupported, please contribute!
    end
end 

% prepends and stack
header_paths = cellfun(@(c)['-I' c],header_paths,'uni',false);
lib_paths = cellfun(@(c)['-L' c],lib_paths,'uni',false);
libs = cellfun(@(c)['-l' c],libs,'uni',false);
defines = cellfun(@(c)['-D' c],defines,'uni',false);
compilerFlags = {header_paths{:} lib_paths{:} libs{:} defines{:}}; %#ok
 
% output folder
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end


%% Build MEX 
mex('-v', compilerFlags{:}, './safmex_tracker3d.c', '-outdir', outputFolder);
mex('-v', compilerFlags{:}, './safmex_latticeDecorrelator.c', '-outdir', outputFolder);
mex('-v', compilerFlags{:}, './safmex_faf_IIRFilterbank.c', '-outdir', outputFolder);
mex('-v', compilerFlags{:}, './safmex_afSTFT.c', '-outdir', outputFolder);
mex('-v', compilerFlags{:}, './safmex_qmf.c', '-outdir', outputFolder);
mex('-v', compilerFlags{:}, './safmex_generateVBAPgainTable3D.c', '-outdir', outputFolder);
mex('-v', compilerFlags{:}, './safmex_getSHreal.c', '-outdir', outputFolder);
mex('-v', compilerFlags{:}, './safmex_getSHcomplex.c', '-outdir', outputFolder);


%% TESTS
if ~run_tests, return; end
addpath('safmex_tests/')
SAFMEX_TESTS


%% FAQ
% If you get a bug like: 
% "xcodebuild: error: SDK "macosx10.15.4" cannot be located."
% 
% Enter the following commands:
%     compilerCfg = mex.getCompilerConfigurations;
%     open(compilerCfg(1).MexOpt)
%     open(compilerCfg(2).MexOpt)
% 
% Then replace:
%     <SDKVER>
%         <cmdReturns name="xcrun -sdk macosx --show-sdk-version"/>
%     </SDKVER>
%
% With:
%     <SDKVER>
%         <cmdReturns name="xcrun -sdk macosx --show-sdk-version | cut -c1-5"/>
%     </SDKVER>
% 
% And run:
%     mex -setup C++
%     mex -setup C
