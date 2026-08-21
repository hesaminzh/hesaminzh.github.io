function [values] = para_init()
% --- System & Simulation Parameters ---
params.M              = 32;       % # BS antennas
params.Mx             = 8;
params.My             = 4;
params.N              = 16;       % # RIS elements
params.Nx             = 4;
params.Ny             = 4;
params.K              = 2;        % # users
params.L              = 2;        % NLOS path nums
params.Ng             = 4;        % group size for GC
params.P_tx           = 1;        % total BS power (normalized)
params.noise          = 3.16e-12; % noise variance
params.tolerance      = 1e-3;
params.weights        = ones(params.K,1);
%params.tau         = 0.5;      % outage threshold in bps/Hz
nMonte                = 500;      % # Monte Carlo trials
%p_list                = 0:0.1:1;  % blockage probabilities
