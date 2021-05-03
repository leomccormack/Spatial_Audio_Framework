%SAFMEX_QMF Quadrature-Mirror Filterbank (QMF)
%   [FREQ_VECTOR, PROC_DELAY] = SAFMEX_QMF(NCH_IN, NCH_OUT, HOPSIZE, 
%   BLOCKSIZE, HYBRIDMODE, FORMAT, FS) creates an the QMF filterbank, with
%   NCH_IN and NCH_OUT number of input and output channels, respectively.
%
%   The filterbank employs a hopsize of HOPSIZE, and will process BLOCKSIZE
%   number of samples at a time. The flag HYBRIDMODE dictates whether the
%   optional hybrid-filtering is enabled, which sub-divides the lowest 3 
%   QMF bands into 10 hybrid-QMF bands. The FORMAT flag dictates whether  
%   the time-frequency data is given in nBands x nCH x time, or time x nCH  
%   x nBands.
%
%   SAFMEX_QMF() destroys the QMF filterbank.
%
%   Y = SAFMEX_QMF(X) performs the forward transform (QMF analysis) if X is
%   a real-valued time-domain (2D) input - resulting in complex-valued 
%   time-frequency (3D) output Y. If Y is complex-valued time-frequency
%   (3D) input, then the backward transform (QMF synthesis) is performed -
%   resulting in real-valued time-domain (2D) output Y.
%
% INPUT ARGUMENTS 
%   NCH_IN:      Number of input channels
%   NCH_OUT:     Number of output channels
%   HOPSIZE:     Hop size, in samples
%   BLOCKSIZE:   Block size, in samples
%   HYBRIDMODE:  '0' disabled, '1' hybrid-mode is enabled 
%   FORMAT       '0' nBands x nCH x time, '1' time x nCH x nBands
%   FS           Samplerate, in Hz
% OUTPUT ARGUMENTS 
%   FREQ_VECTOR: Frequency vector/centre frequencies; nBands x 1
%   PROC_DELAY:  Processing delay, in samples
%
% EXAMPLE
%
%   % User params
%   nCHin = 10;
%   nCHout = nCHin;
%   hopsize = 128;
%   blocksize = 2048*10;
%   hybridmode = 1;
%   format = 0;
%   fs = 48e3;
%
%   % Create
%   [freqVector, procDelay] = safmex_qmf(nCHin, nCHout, hopsize, blocksize, hybridmode, 0, fs); 
%
%   % Forward
%   dataTD_ref = randn(blocksize, nCHin); 
%   dataFD = safmex_qmf(dataTD_ref.'); 
%
%   % Backward
%   dataTD = safmex_qmf(dataFD);
%
%   % Destroy
%   safmex_qmf()
%
%   % dataTD_ref will approx. (-50dB) equal dataTD, shifted by procDelay
%   % samples
%

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
