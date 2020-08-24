
fs = 48e3;
hopsize = 128;


%% Configure the tracker 
% Number of Monte Carlo samples/particles. The more complex the
% distribution is, the more particles required (but also, the more
% computationally expensive the tracker becomes). 
tpars.Np = 20;
tpars.maxNactiveTargets = 4; % about 2 higher than expected is good 
% Likelihood of an estimate being noise/clutter 
tpars.noiseLikelihood = 0.2; % between [0..1] 
% Measurement noise - e.g. to assume that estimates within the range +/-20
% degrees belong to the same target, set SDmnoise_deg = 20 
measNoiseSD_deg = 20;
tpars.measNoiseSD = 1-cos(measNoiseSD_deg*pi/180);
% Noise spectral density - not fully understood. But it influences the
% smoothness of the target tracks 
noiseSpecDen_deg = 1;
tpars.noiseSpecDen = 1-cos(noiseSpecDen_deg*pi/180);
% FLAG - whether to allow for multiple target deaths in the same tracker
% prediction step 
tpars.ALLOW_MULTI_DEATH = 1;
% Probability of birth and death 
tpars.init_birth = 0.5; % value between [0 1] - Prior probability of birth 
tpars.alpha_death = 1.0; % always >= 1; 1 is good 
tpars.beta_death = 20.0; % always >= 1; 1 is good 
% Elapsed time (in seconds) between observations 
tpars.dt = 1.0/(fs/hopsize); % Hop length of frames 
% Whether or not to allow multiple active sources for each update  
% Real-time tracking is based on the particle with highest weight. A
% one-pole averaging filter is used to smooth the weights over time. 
tpars.W_avg_coeff = 0.0;
% Force kill targets that are close to another target. In these cases, the
% target that has been 'alive' for the least amount of time, is killed 
tpars.FORCE_KILL_TARGETS = 1;
tpars.forceKillDistance = 0.2;
% Mean position priors x,y,z (assuming directly in-front) 
tpars.M0(1,1) = 1.0; tpars.M0(2,1) = 0.0; tpars.M0(3,1) = 0.0;
% Mean Velocity priors x,y,z (assuming stationary) 
tpars.M0(4,1) = 0.0; tpars.M0(5,1) = 0.0; tpars.M0(6,1) = 0.0;
% Target velocity - e.g. to assume that a target can move 20 degrees in
% two seconds along the horizontal, set V_azi = 20/2 
Vazi_deg = 3.0;  % Velocity of target on azimuthal plane 
Vele_deg = 3.0;  % Velocity of target on median plane 
tpars.P0 = zeros(6); 
% Variance PRIORs of estimates along the x,y,z axes, respectively. Assuming
% coordinates will lay on the unit sphere +/- x,y,z, so a range of 2, and
% hence a variance of 2^2: 
tpars.P0(1,1) = 4.0; tpars.P0(2,2) = 4.0; tpars.P0(3,3) = 4.0;
% Velocity PRIORs of estimates the x,y,z axes 
tpars.P0(4,4) = 1.0-cos(Vazi_deg*pi/180.0); % x 
tpars.P0(5,5) = tpars.P0(4,4);              % y 
tpars.P0(6,6) = 1.0-cos(Vele_deg*pi/180.0); % z 
% PRIOR probabilities of noise. (Assuming the noise is uniformly
% distributed in the entire spatial grid). 
tpars.cd = 1.0/(4.0*pi);


%%
safmex_tracker3d(tpars);

safmex_tracker3d();


