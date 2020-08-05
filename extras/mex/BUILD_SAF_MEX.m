%
% Copyright 2020 Leo McCormack
%
% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
%
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
% REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
% INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
% LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
% OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
% PERFORMANCE OF THIS SOFTWARE.
%
clear all %#ok

% 1. compile SAF via CMake prior to running this script
% 2. configure user options below as needed
% If you encounter any problems during compilation, please refer to the FAQ 
% at the bottom of this script before submitting a bug report.


%% User Options 
saf_perf_lib = []; % leave empty to scan cache, 'SAF_USE_INTEL_MKL', 'SAF_USE_APPLE_ACCELERATE', 'SAF_USE_OPEN_BLAS_AND_LAPACKE'
run_tests = 1;     % 0: disabled, 1: enabled
cmake_build_folder = '../../build/'; % path to CMake build folder
outputFolder = 'safmex/';


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
    elseif isempty(idx)
        idx = find(strcmp( cac, 'SAF_PERFORMANCE_LIB:STRING=SAF_USE_APPLE_ACCELERATE' ));
        if isempty(idx), idx = find(strcmp( cac, 'SAF_PERFORMANCE_LIB:UNINITIALIZED=SAF_USE_APPLE_ACCELERATE' )); end
        if ~isempty(idx), saf_perf_lib = 'SAF_USE_APPLE_ACCELERATE'; end 
    elseif isempty(idx)
        idx = find(strcmp( cac, 'SAF_PERFORMANCE_LIB:STRING=SAF_USE_OPEN_BLAS_AND_LAPACKE' ));
        if isempty(idx), idx = find(strcmp( cac, 'SAF_PERFORMANCE_LIB:UNINITIALIZED=SAF_USE_OPEN_BLAS_AND_LAPACKE' )); end
        if ~isempty(idx), saf_perf_lib = 'SAF_USE_OPEN_BLAS_AND_LAPACKE'; end 
    end  
end


%% Configure
header_paths = {'../../framework/include'};
lib_paths = {[ cmake_build_folder 'framework' ]};
libs = {};
defines = {saf_perf_lib};

% Platform specific configurations
if ismac 
    libs{end+1} = 'saf';
    
    % Error checking
    if ~exist([cmake_build_folder 'framework/libsaf.a'], 'file'), error('libsaf.a is missing, run CMake and build saf first'), end
    
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
    end 
    
elseif isunix
    assert(0); %currently unsupported 
elseif ispc
    assert(0); %currently unsupported 
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
mex('-v', compilerFlags{:}, './src/safmex_afSTFT.c', '-outdir', outputFolder);
%mex('-v', compilerFlags{:}, './src/safmex_qmf.c', '-outdir', outputFolder);
%mex('-v', compilerFlags{:}, './src/safmex_generateVBAPgainTable3D.c', '-outdir', outputFolder);
%mex('-v', compilerFlags{:}, './src/safmex_getSHreal.c', '-outdir', outputFolder);
%mex('-v', compilerFlags{:}, './src/safmex_getSHcomplex.c', '-outdir', outputFolder);


%% TESTS
if ~run_tests, return; end
addpath(outputFolder)
nTests = 0; nPass = 0; nFail = 0;
tol = 1e-5; % FLT_EPSILON

% safmex_afSTFT
nCHin = 6;
nCHout = 13;
hopsize = 128;
blocksize = 2048*24;
[freqVector, procDelay] = safmex_afSTFT(nCHin, nCHout, hopsize, blocksize, 1, 1, 48e3); 
dataTD_ref = randn(blocksize, nCHin); 
dataFD = safmex_afSTFT(dataTD_ref.');
dataFD = repmat(dataFD(:,1,:), [1 nCHout 1]); % copy 1st input channel to all output channels
dataTD = safmex_afSTFT(dataFD);
dataTD_ref = dataTD_ref(1:end-procDelay,1);
dataTD = dataTD(1,procDelay+1:end).';
nTests = nTests+1;
if max(abs(dataTD(:,1)-dataTD_ref(:,1)))<0.01, nPass=nPass+1; else, nFail=nFail+1; end 
safmex_afSTFT();

% safmex_qmf
nCHin = 10;
nCHout = 4;
hopsize = 128;
blocksize = 2048*20;
[freqVector2, procDelay] = safmex_qmf(nCHin, nCHout, hopsize, blocksize, 1, 0, 48e3); 
dataTD_ref = randn(blocksize, nCHin); 
dataFD = safmex_qmf(dataTD_ref.');
dataFD = repmat(dataFD(:,1,:), [1 nCHout 1]); % copy 1st input channel to all output channels
dataTD = safmex_qmf(dataFD);
dataTD_ref = dataTD_ref(1:end-procDelay,1);
dataTD = dataTD(1,procDelay+1:end).';
nTests = nTests+1;
if max(abs(dataTD(:,1)-dataTD_ref(:,1)))<0.01, nPass=nPass+1; else, nFail=nFail+1; end 
safmex_qmf();

% safmex_generateVBAPgainTable3D
[~,ls_dirs] = getTdesign(10);
ls_dirs = ls_dirs*180/pi;
aziElevRes = [5 5];
gtable = safmex_generateVBAPgainTable3D(ls_dirs, aziElevRes(1), aziElevRes(2), 1, 0, 0);
gtable_ref = getGainTable(ls_dirs,aziElevRes,0,'vbap');
nTests = nTests+1;
if max(abs(gtable_ref(:)-gtable(:)))<tol, nPass=nPass+1; else, nFail=nFail+1; end 

% safmex_getSHreal
order = 5;
dirs = randn(1000,2); 
Y = safmex_getSHreal(order,dirs); 
Y_ref = getSH(order, dirs, 'real').'; 
nTests = nTests+1;
if max(abs(Y(:)-Y_ref(:)))<tol, nPass=nPass+1; else, nFail=nFail+1; end 

% safmex_getSHcomplex
order = 4;
dirs = randn(400,2); 
Y = safmex_getSHcomplex(order,dirs);
Y_ref = getSH(order, dirs, 'complex').';
nTests = nTests+1;
if max(abs(Y(:)-Y_ref(:)))<tol, nPass=nPass+1; else, nFail=nFail+1; end 

% Output
display(['Number of tests: ' num2str(nTests)]);
display(['Number of failures: ' num2str(nFail)]);


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
