clear all %#ok


%% User Options 
saf_perf_lib = 'SAF_USE_INTEL_MKL'; % 'SAF_USE_INTEL_MKL', 'SAF_USE_APPLE_ACCELERATE', 'SAF_USE_OPEN_BLAS_AND_LAPACKE'
run_tests = 1; % 0:disable, 1:enable
cmake_build_folder = '../../build'; % absolute/relative path to CMake build folder (compile saf prior to running script)


%% Configure
header_paths = {'../../framework/include'};
lib_paths = {[ cmake_build_folder '/framework' ]};
libs = {};
defines = {saf_perf_lib};

% Platform specific configurations
if ismac 
    libs{end+1} = 'saf';
    
    % Error checking
    if ~exist([cmake_build_folder '/framework/libsaf.a'], 'file'), error('libsaf.a is missing, run CMake and build saf first'), end
    
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

% prepends
header_paths = cellfun(@(c)['-I' c],header_paths,'uni',false);
lib_paths = cellfun(@(c)['-L' c],lib_paths,'uni',false);
libs = cellfun(@(c)['-l' c],libs,'uni',false);
defines = cellfun(@(c)['-D' c],defines,'uni',false);


%% Build MEX
mex('-v', header_paths{:}, lib_paths{:}, libs{:}, defines{:}, 'saf_getSHreal.c');


%% TESTS
A = ones(10,1).';
s = 2;
B = saf_getSHreal(s,A) 


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
