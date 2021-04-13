%SAFMEX_LATTICEDECORRELATOR Lattice all-pass filter based decorrelator
%   SAFMEX_LATTICEDECORRELATOR(FS, HOPSIZE, NCH, ORDERS, FREQCUTOFFS, 
%   MAXDELAY, FREQVECTOR, TIMESLOTS) creates the latticeDecorrelator for 
%   NCH number of channels. 
%
%   The decorrelator groups the input signals based on the FREQCUTOFFS
%   divisions, and applies a fixed delay to each grouping defined by vector
%   FIXEDDELAYS and lattice all-pass filters of order ORDERS.
%
%   SAFMEX_LATTICEDECORRELATOR() destroys the latticeDecorrelator.
%
%   Y = SAFMEX_LATTICEDECORRELATOR(X) performs the decorrelation on the
%   complex-valued time-frequency (3D) input X. Resulting in the
%   decorrelated complex-valued time-frequency (3D) output Y. The
%   dimensions of both X and Y, are nBands x nCH x timeslots; where nBands
%   is length(FREQVECTOR).
%
% INPUT ARGUMENTS 
%   FS:          Sampling rate
%   HOPSIZE:     Hop size in samples
%   NCH:         Number of channels
%   ORDERS       Lattice all-pass filter orders (2,3,4,6,8,10,12,14,15,16,18
%                , or 20) per band grouping; nCutoffs x 1
%   FREQCUTOFFS  Frequency cut-offs defining the band groupings; 
%                nCutoffs x 1
%   MAXDELAY:    Maximum delay in hops (TIMESLOTS)
%   FREQVECTOR   Frequency vector; nBands x 1 */ 
%   TIMESLOTS    Number of TF frames to process at a time  
%
% EXAMPLE
%
%   % User params
%   fs = 48e3;
%   nCH = 6; 
%   hopsize = 128;
%   blocksize = 2048*24; 
%   orders = [20 15 6 3].'; 
%   freqCutoffs = [ 900 2.4e3 8e3, 16e3].';
%   timeslots = blocksize/hopsize;
%
%   % Create
%   [freqVector, procDelay] = safmex_afSTFT(nCH, nCH, hopsize, blocksize, 1, 0, 48e3); 
%   safmex_latticeDecorrelator(fs, hopsize, nCH, orders, freqCutoffs, maxDelay, freqVector, timeslots);
%
%   % Forward TF transform
%   dataTD_ref = randn(blocksize, nCH); 
%   dataFD = safmex_afSTFT(dataTD_ref.'); 
%
%   % Decorrelate
%   dataFD = safmex_latticeDecorrelator(dataFD);
%
%   % Inverse TF transform
%   dataTD = safmex_afSTFT(dataFD);
%   dataTD_ref = dataTD_ref(1:end-procDelay,1);
%   dataTD = dataTD(1,procDelay+1:end).';
%
%   % Destroy
%   safmex_afSTFT();
%   safmex_latticeDecorrelator();
%
%   % the interchannel coherence (ICC) between dataTD_ref and dataTD be
%   around ~0
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
