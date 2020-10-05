% A script to load HRIR data from a sofa file, and generate the .h/.c
% files to be used as the default HRIR data within SAF.
%
% Author: Leo McCormack
% Date: 04.10.2020
% License: ISC

close all, clear all %#ok
addpath('../utils');                  % for mat2c.m
addpath('../saf_sofa_reader_module'); % for saf_sofa_open.m

% Load SOFA file
sofa = saf_sofa_open('/Users/mccorml1/Documents/HRIRs_SOFA/KemarAuralID.sofa');

% Permute as needed
hrir_dirs_deg = sofa.SourcePosition(1:2,:).';
hrirs = permute(sofa.Data_IR, [3 2 1]);

% Copy data into a struct
structToPass.default_hrirs = single(hrirs);  % nDirs x 2 x len
structToPass.default_hrir_dirs_deg = single(hrir_dirs_deg);
structToPass.default_N_hrir_dirs =  int32(size(hrir_dirs_deg,1));
structToPass.default_hrir_len = int32(size(hrirs,3));
structToPass.default_hrir_fs = int32(sofa.Data_SamplingRate);

% Specify .h/.c file export configurations
filename = 'saf_default_hrirs';
dataPrepend = '__';
headerdef = '__SAF_DEFAULT_HRIRS_INCLUDED__';
headerNotice = { ...
 ' * Copyright 2020 Leo McCormack                                                  \n'
 ' *                                                                               \n'
 ' * Permission to use, copy, modify, and/or distribute this software for any      \n'
 ' * purpose with or without fee is hereby granted, provided that the above        \n'
 ' * copyright notice and this permission notice appear in all copies.             \n'
 ' *                                                                               \n'
 ' * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH \n'
 ' * REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY   \n'
 ' * AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,  \n'
 ' * INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM   \n'
 ' * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR \n'
 ' * OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR        \n'
 ' * PERFORMANCE OF THIS SOFTWARE.                                                 \n'
 ' */                                                                              \n'
 '                                                                                 \n'
 '/**                                                                              \n'
 ' * @file saf_default_hrirs.h                                                     \n'
 ' * @ingroup HRIR                                                                 \n'
 ' * @brief Default HRIR data                                                      \n'
 ' *                                                                               \n'
 ' * The default HRIR set is a Genelec Aural ID of a KEMAR Dummy Head (@48kHz).    \n'
 ' * Kindly provided by Aki Mäkivirta and Jaan Johansson                           \n'
 ' *                                                                               \n'
 ' * @author Leo McCormack                                                         \n'
 ' * @date 17.04.2020                                                              \n'};

% Generate the .h and .c file
passAsCell = {structToPass};
mat2c(passAsCell, filename, dataPrepend, headerdef, headerNotice, '' );




