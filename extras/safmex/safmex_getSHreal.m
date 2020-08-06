%SAFMEX_GETSHREAL compute real spherical harmonics
%   Y = SAFMEX_GETSHREAL(ORDER, DIRS_RAD) computes spherical harmonics Y of
%   order ORDER for all DIRS_RAD directions.
%
%   The spherical harmonics are computed WITH the 1/sqrt(4*pi) term.
%
% INPUT ARGUMENTS 
%   ORDER:    Spherical harmonic order
%   DIRS_RAD: Directions on the sphere [azi, INCLINATION] convention, in 
%             RADIANS; nDirs x 2
% OUTPUT ARGUMENTS 
%   Y: Spherical harmonics; (order+1)^2 x nDirs
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
