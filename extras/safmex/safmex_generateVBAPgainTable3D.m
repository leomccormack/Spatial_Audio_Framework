%SAFMEX_GENERATEVBAPGAINTABLE3D computes a VBAP gain table
%   GTABLE = SAFMEX_GENERATEVBAPGAINTABLE3D(DIRS_DEG, AZI_RES, ELEV_RES, 
%   OMIT_LARGE_TRIS, ENABLE_DUMMIES, SPREAD) computes a VBAP gain table
%   GTABLE for loudspeaker directions DIRS_DEG, for every azimuth and 
%   elevation step AZI_RES and ELEV_RES, repectively. 
%
%   The OMIT_LARGE_TRIS flag will remove any triangles which wrap around 
%   the sphere; ENABLE_DUMMIES will place dummy loudspeakers at +90 or -90
%   degrees elevation if no loudspeakers exist above +60 or -60 degrees 
%   elevation, respectively; and SPREAD (in degrees) introduces source 
%   spreading.
%
% INPUT ARGUMENTS 
%   DIRS_DEG:        Loudspeaker directions in degrees; L x 2 
%   AZI_RES:         Azimuthal resolution in degrees
%   ELEV_RES:        Elevation resolution in degrees
%   OMIT_LARGE_TRIS: '0' normal triangulation, '1' remove large triangles
%   ENABLE_DUMMIES:  '0' disabled, '1' enabled. Dummies are placed at +/-90
%                    elevation if required
%   SPREAD           Spreading factor in degrees, 0: VBAP, >0: MDAP
% OUTPUT ARGUMENTS 
%   GTABLE: Gain Table; nSources x L
%
% EXAMPLE
%
%   N_azi = 360 / AZI_RES + 1;
%   aziIndex = round(mod(azi+180,360)/AZI_RES);
%   elevIndex = round((elev+90)/ELEV_RES);
%   idx3D = elevIndex*N_azi+aziIndex+1;
%   gains3D = GTABLE(idx3D,:);
%
%   % Will return the gains required to pan a source to [azi, elev].

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
