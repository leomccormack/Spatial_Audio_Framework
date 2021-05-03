function mat2c( matfiles, filename, dataPrepend, headerdef, headerNotice, anyExtraStuff )
%mat2c - converts matlab structs to a C header and source file
% note: does not currently support cells within structs, or 5-dimensional+
% arrays; structs inside structs are unwrapped into one struct and any data
% elements that have the same name may be distinguished by the '_' suffix;
%
% Syntax: mat2c( matfiles, filename, dataPrepend, headerdef,
%         headerNotice )
%
% Inputs:
%    matfiles - cell array of matlab structs
%    filename - file name string ('filename.h' & 'filename.c')
%    dataPrepend - prepend for the names of all of the data elements
%    headerdef - header definition
%    headerNotice - e.g. copyright notice and/or description of contents
%    anyExtraStuff - Literally any extra stuff [string]...
%
% Outputs:
%    The header and source files are created in the current file directory
%
% Example:
%    matfiles = {load('example.mat'), load('example2.mat')}; % etc.
%    filename = 'example_database';
%    dataPrepend = 'ex_';
%    headerdef = '__EXAMPLE_DATABASE_INCLUDED__';
%    headerNotice = {
%    '    Leo McCormack\n'
%    '    An example C database using 'mat2c function
%    ['    Generated using MatLab on: ' date '\n']
%    '    Copyright X 2016\n'};
%
%    mat2c( matfiles, filename, dataPrepend, headerdef, headerNotice, '' );
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Leo McCormack, M.Sc, Acoustics and Audio Technology
% Aalto University, Dept. of Acoustics and Signal Processing
% email address: leo.mccormack@aalto.fi
% October 2016; Last revision: 13-Oct-2016
% License: ISC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1; error('not enough input arguments'); end
% defaults:
if nargin < 2; filename = 'database'; end
if nargin < 3; dataPrepend = ''; end
if nargin < 4; headerdef = '__DATABASE_H_INCLUDED__'; end
if nargin < 5; headerNotice = 'This is a database'; end


%% Unwrap any structs within structs
for m = 1: size(matfiles,2)
    matfile = matfiles(m);
    matfile = cell2mat(matfile);
    foundStruct = 1;
    while foundStruct == 1
        % recursively pull the contents from substructs to the main struct
        [foundStruct, matfile] = unwrapstructs(matfile);
    end
    matfiles{m} = matfile;
end


%% Generate Header File
fd=fopen([filename '.h'],'wt');

fprintf(fd, '/*\n');
for i=1:size(headerNotice,1)
    str = cell2mat(headerNotice(i));
    fprintf(fd, str);
end
fprintf(fd, '*/\n');

fprintf(fd, '\n');
fprintf(fd, ['#ifndef ' headerdef '\n']);
fprintf(fd, ['#define ' headerdef '\n']);
fprintf(fd, '\n');

fprintf(fd,anyExtraStuff);

for m = 1: size(matfiles,2)
    matfile = matfiles(m);
    matfile = cell2mat(matfile);
    names = fieldnames(matfile);

    for i = 1:numel(names)
        subfield = matfile.(names{i});
        if (isreal(subfield))
            fprintf(fd, ['extern const ' dataType(subfield) ' ' dataPrepend cell2mat(names(i))]);
            for dim = 1:ndims(subfield)
                size_dim = size(subfield, dim);
                if size_dim ~=1
                    fprintf(fd, '[%d]', size_dim);
                end
            end
            fprintf(fd, ';\n');
        else
            fprintf(fd, ['extern const ' dataType(subfield) ' ' dataPrepend cell2mat(names(i)) '_re']);
            for dim = 1:ndims(subfield)
                size_dim = size(subfield, dim);
                if size_dim ~=1
                    fprintf(fd, '[%d]', size_dim);
                end

            end
            fprintf(fd, ';\n');

            fprintf(fd, ['extern const ' dataType(subfield) ' ' dataPrepend cell2mat(names(i)) '_im']);
            for dim = 1:ndims(subfield)
                size_dim = size(subfield, dim);
                if size_dim ~=1
                    fprintf(fd, '[%d]', size_dim);
                end
            end
            fprintf(fd, ';\n');
        end
    end
end

fprintf(fd, '\n');
fprintf(fd, ['#endif /* ' headerdef ' */']);
fclose(fd);


%% Generate Source File
fd=fopen([filename '.c'],'wt');

fprintf(fd,['/* Generated using MatLab on: ' date ', using a script by Leo McCormack */\n']);
fprintf(fd,'/* (leo.t.mccormack@gmail.com) */');
fprintf(fd,'\n');


for m = 1: size(matfiles,2)
    matfile = matfiles(m);
    matfile = cell2mat(matfile);
    names = fieldnames(matfile);

    for i = 1:numel(names)
        subfield = matfile.(names{i});
        if (isreal(subfield))
            fprintf(fd, [ 'const ' dataType(subfield) ' ' dataPrepend cell2mat(names(i))]);
            dim_idx = ndims(subfield);
            for dim = 1:ndims(subfield)
                size_dim = size(subfield, dim);
                if size_dim ~= 1
                    fprintf(fd, '[%d]', size_dim);
                else
                    dim_idx = dim_idx - 1;
                end
            end

            fprintf(fd, ' = ');
            dim_idx = 1:dim_idx;
            if size(dim_idx,2) ~= 0
                fprintf(fd, '\n{    ');
            end
            switch size(dim_idx,2)
                case 0
                    print0d(fd, subfield);
                case 1
                    print1d(fd, dim_idx, subfield);
                case 2
                    print2d(fd, dim_idx, subfield);
                case 3
                    print3d(fd, dim_idx, subfield);
                case 4
                    print4d(fd, dim_idx, subfield);
                % currently only supports <= 4-dimensional arrays
            end
            if size(dim_idx,2) ~= 0
                fprintf(fd, '    }');
            end
            fprintf(fd, ';\n\n');
        else
            %real part
            fprintf(fd, [ 'const ' dataType(subfield) ' ' dataPrepend cell2mat(names(i)) '_re']);
            dim_idx = ndims(subfield);
            for dim = 1:ndims(subfield)
                size_dim = size(subfield, dim);
                if size_dim ~= 1
                    fprintf(fd, '[%d]', size_dim);
                else
                    dim_idx = dim_idx - 1;
                end
            end

            fprintf(fd, ' = ');
            dim_idx = 1:dim_idx;
            if size(dim_idx,2) ~= 0
                fprintf(fd, '\n{    ');
            end
            switch size(dim_idx,2)
                case 0
                    print0d(fd, real(subfield));
                case 1
                    print1d(fd, dim_idx, real(subfield));
                case 2
                    print2d(fd, dim_idx, real(subfield));
                case 3
                    print3d(fd, dim_idx, real(subfield));
                case 4
                    print4d(fd, dim_idx, real(subfield));
                % currently only supports <= 4-dimensional arrays
            end
            if size(dim_idx,2) ~= 0
                fprintf(fd, '    }');
            end
            fprintf(fd, ';\n\n');

            %imaginary part
            fprintf(fd, [ 'const ' dataType(subfield) ' ' dataPrepend cell2mat(names(i)) '_im']);
            dim_idx = ndims(subfield);
            for dim = 1:ndims(subfield)
                size_dim = size(subfield, dim);
                if size_dim ~= 1
                    fprintf(fd, '[%d]', size_dim);
                else
                    dim_idx = dim_idx - 1;
                end
            end

            fprintf(fd, ' = ');
            dim_idx = 1:dim_idx;
            if size(dim_idx,2) ~= 0
                fprintf(fd, '\n{    ');
            end
            switch size(dim_idx,2)
                case 0
                    print0d(fd, real(subfield));
                case 1
                    print1d(fd, dim_idx, imag(subfield));
                case 2
                    print2d(fd, dim_idx, imag(subfield));
                case 3
                    print3d(fd, dim_idx, imag(subfield));
                case 4
                    print4d(fd, dim_idx, imag(subfield));
                % currently only supports <= 4-dimensional arrays
            end
            if size(dim_idx,2) ~= 0
                fprintf(fd, '    }');
            end
            fprintf(fd, ';\n\n');
        end
    end
end


fclose(fd);


end % database_generator


%% the guts
function string = dataType(data)

if strcmp(class(data), 'int32')
    string = 'int';
elseif strcmp(class(data), 'single')
    string = 'float';
else
    string = class(data);
end

end

function print2f(fd, data)

if strcmp(class(data), 'int32') %#ok
    fprintf(fd,'%d', round(data));
elseif strcmp(class(data), 'single')  %#ok
    fprintf(fd,'%.9ff', data);
elseif strcmp(class(data), 'double')  %#ok
    fprintf(fd,'%.13f', data);
end

end


function [foundStruct, matfile] = unwrapstructs(matfile)
n = fieldnames(matfile);
foundStruct = 0;
for i = 1:numel(n)
    subfield = matfile.(n{i});

    if isstruct(subfield)
        nn = fieldnames(subfield);
        for j = 1:numel(nn)
            if isfield(matfile, nn{j})
                matfile.([nn{j} '_']) = subfield.(nn{j}); % add '_' suffix, if field name already exists
            else
                matfile.(nn{j}) = subfield.(nn{j});
            end
        end
        foundStruct = 1;
        matfile = rmfield(matfile, n{i});
    end
end

end

function print0d(fd, subfield)
print2f(fd, subfield);
end

function print1d(fd, dim_idx, subfield)
if ischar(subfield)
    fprintf(fd,' "%s" ', subfield);
else
    if isrow(subfield)
        subfield = subfield.'; %force column
    end

    for j = 1: size(subfield, max(dim_idx))
        print2f(fd, subfield(j,1));
        if j~=size(subfield, max(dim_idx));  fprintf(fd,', ');  end
    end
end
end

function print2d(fd, dim_idx, subfield)
for j=1:size(subfield, dim_idx(1))
    fprintf(fd,'\n    { ');
    for k=1:size(subfield, dim_idx(2))
        print2f(fd, subfield(j,k));
        if k~=size(subfield, dim_idx(2));  fprintf(fd,', ');  end
    end
    fprintf(fd,'} ');
    if j~=size(subfield, dim_idx(1));  fprintf(fd,', ');  end
end
end

function print3d(fd, dim_idx, subfield)
for j=1:size(subfield, dim_idx(1))
    fprintf(fd,'\n    { ');
    for k=1:size(subfield, dim_idx(2))
        fprintf(fd,'{ ');
        for l=1:size(subfield, dim_idx(3))
            print2f(fd, subfield(j,k,l));
            if l~=size(subfield, dim_idx(3));  fprintf(fd,', ');  end
        end
        fprintf(fd,'} ');
        if k~=size(subfield, dim_idx(2));  fprintf(fd,', ');  end
    end
    fprintf(fd,'} ');
    if j~=size(subfield, dim_idx(1));  fprintf(fd,', ');  end
end
end

function print4d(fd, dim_idx, subfield)
for j=1:size(subfield, dim_idx(1))
    fprintf(fd,'\n    { ');
    for k=1:size(subfield, dim_idx(2))
        fprintf(fd,'{ ');
        for l=1:size(subfield, dim_idx(3))
            fprintf(fd,'{ ');
            for m=1:size(subfield, dim_idx(4))
                print2f(fd, subfield(j,k,l,m));
                if m~=size(subfield, dim_idx(4));  fprintf(fd,', ');  end
            end
            fprintf(fd,'} ');
            if l~=size(subfield, dim_idx(3));  fprintf(fd,', ');  end
        end
        fprintf(fd,'} ');
        if k~=size(subfield, dim_idx(2));  fprintf(fd,', ');  end
    end
    fprintf(fd,'} ');
    if j~=size(subfield, dim_idx(1));  fprintf(fd,', ');  end
end
end

% currently only supports <= 4-dimensional arrays

