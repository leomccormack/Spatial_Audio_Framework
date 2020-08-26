%SAFMEX_TRACKER3D Particle-filtering based tracker
%   SAFMEX_TRACKER3D(TPARS) creates the tracker
% 
%   SAFMEX_TRACKER3D() destroys the tracker
%
%   [TARGET_STATES, TARGET_IDS] = SAFMEX_TRACKER3D(MEASUREMENTS) takes 
%   measurements/observations MEASUREMENTS (nObs x 3) and feeds them to the
%   tracker; and returns the TARGET_STATES and TARGET_IDS per target. 
%
%   Note that TARGET_STATES and TARGET_IDS will be empty ([]) if no targets
%   are being tracked. SAFMEX_TRACKER3D(MEASUREMENTS) must also be called
%   for every step in time (TPARS.dt). If there are no new measurements for
%   a specific time step, then still call the function but pass 
%   MEASUREMENTS=[].
%
% INPUT ARGUMENTS 
%   TPARS.Np:                Number of Monte Carlo samples/particles 1..100
%   TPARS.maxNactiveTargets: Maximum number of targets to track 1..100
%   TPARS.noiseLikelihood:   Likelihood of clutter/noise [0..1]
%   TPARS.measNoiseSD:       Measurement noise
%   TPARS.noiseSpecDen       Noise spectral density
%   TPARS.ALLOW_MULTI_DEATH: FLAG - whether to allow for multiple target 
%                            deaths in the same tracker prediction step 
%   TPARS.init_birth:        Prior probability of birth [0..1]
%   TPARS.alpha_death = 1:   Prior probability of death >= 1 
%   TPARS.beta_death = 1:    Prior probability of death >= 1
%   TPARS.dt:                Elapsed time (in seconds) between observations 
%   TPARS.W_avg_coeff:       Averaging coeff for importance weights [0..1]
%   TPARS.FORCE_KILL_TARGETS Kill targets that are too close to eachother
%   TPARS.forceKillDistance: Distance at which to start killing...
%   TPARS.M0:                Mean position and velocity priors; 6 x 1
%   TPARS.P0:                Position and velocity variance priors; 6 x 6
%   TPARS.cd:                Prior probability of noise, [0..1]
%
% EXAMPLE
% 
%   % Configure the tracker 
%   tpars.Np = 20;
%   tpars.maxNactiveTargets = 4; 
%   tpars.noiseLikelihood = 0.2; 
%   tpars.measNoiseSD = 0.2340; 
%   tpars.noiseSpecDen = 1.5231e-06;
%   tpars.ALLOW_MULTI_DEATH = 1;
%   tpars.init_birth = 0.1;  
%   tpars.alpha_death = 1; 
%   tpars.beta_death = 1;  
%   tpars.dt = 0.02;  
%   tpars.W_avg_coeff = 0.0;
%   tpars.FORCE_KILL_TARGETS = 1;
%   tpars.forceKillDistance = 0.25;
%   tpars.M0(1,1) = 1; tpars.M0(2,1) = 0; tpars.M0(3,1) = 0; 
%   tpars.M0(4,1) = 0; tpars.M0(5,1) = 0; tpars.M0(6,1) = 0;
%   tpars.P0 = zeros(6); 
%   tpars.P0(1,1) = 4; tpars.P0(2,2) = 4; tpars.P0(3,3) = 4;
%   tpars.P0(4,4) = 0.0014;
%   tpars.P0(5,5) = 0.0014;
%   tpars.P0(6,6) = 0.0007;
%   tpars.cd = 1/(4*pi);
%
%   % Create
%   safmex_tracker3d(tpars);
%
%   % Track
%   [target_xyz, target_IDs] = safmex_tracker3d(measurements_xyz_time_0); 
%   [target_xyz, target_IDs] = safmex_tracker3d(measurements_xyz_time_1); 
%   ...
%   [target_xyz, target_IDs] = safmex_tracker3d(measurements_xyz_time_N); 
% 
%   % Destroy
%   safmex_tracker3d();
%  

%
% This file is part of the saf_tracker module.
% Copyright (c) 2020 - Leo McCormack
%
% The saf_tracker module is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% The saf_tracker module is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
%
% See <http://www.gnu.org/licenses/> for a copy of the GNU General Public
% License.
%
