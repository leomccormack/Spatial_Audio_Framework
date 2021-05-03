%SAFMEX_FAF_IIRFILTERBANK Power-complementary IIR filterbank
%   SAFMEX_FAF_IIRFILTERBANK(ORDER, CUTOFFS, SIGLEN, FS) creates the FaF 
%   filterbank.
%
%   This filterbank employs power-complementary butterworth IIR filters of 
%   order ORDER to divide a mono input signal of lenth SIGLEN. The result 
%   is frequency band signals (also of length SIGLEN), which lay within the
%   cut-off frequencies CUTOFFS.
%
%   SAFMEX_FAF_IIRFILTERBANK() destroys the FaF filterbank.
%
%   Y = SAFMEX_FAF_IIRFILTERBANK(X) performs the filtering on mono signal X
%   to obtain frequency band signals Y.
%
% INPUT ARGUMENTS 
%   ORDER:    IIR filter order (must be 1 or 3)
%   CUTOFFS:  Cut-off frequencies, in Hz, as column vector
%   SIGLEN:   Signal length, in samples 
%   FS        Samplerate, in Hz 
%
% EXAMPLE
%
%   % User params
%   order = 1;
%   lSig = 2048*10;
%   fs = 48e3;
%   cutoffFreqs = [125 250 500 1000 2000 4000]*2/sqrt(2);
%
%   % Create
%   safmex_faf_IIRFilterbank(order, cutoffFreqs.', lSig, fs);  
%
%   % Apply
%   Impulse = zeros(lSig, 1);  Impulse(1,1) = 1;
%   IRs = safmex_faf_IIRFilterbank(Impulse).';
%
%   % The magnitude responses of IRs, should sum to ~<0.1 dB at all
%   % frequencies
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
