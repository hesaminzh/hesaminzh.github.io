%% main_simulate_blockage_effects.m
% Generates Fig.7 (WSR vs. blockage) and Fig.8 (Outage vs. blockage)

clear; close all; clc;
rng(3,'twister');

% --- System & Simulation Parameters ---
params.M       = 32;       % # BS antennas
params.N       = 16;       % # RIS elements
params.K       = 2;        % # users
params.Ng      = 4;        % group size for GC
params.P_tx    = 1;        % total BS power (normalized)
params.noise   = 3.16e-12; % noise variance
params.weights = ones(params.K,1);
%params.tau     = 0.5;      % outage threshold in bps/Hz
nMonte         = 500;      % # Monte Carlo trials
p_list         = 0:0.1:1;  % blockage probabilities

schemes = { ...
  'FC' ...               % fully-connected RIS (AO)
  %'GC', ...              % group-connected RIS (AO, Ng=4)
  %'SC', ...              % single-connected RIS (optimized AO)
  %'RandPhase', ...       % random-phase RIS (SC)  
  %'NoRIS', ...           % direct path only
};



%%%%%%%%%%%%
%% find the correct value of tau so that outage vs blockage shows the S-curve in the figure.
%
%
%% 1) Pilot simulation at p=0.5 to estimate the WSR distribution
p_pilot  = 0.5;
nPilot   = 200;          % enough to get a decent median
wsr_pilot = zeros(nPilot,1);

for it=1:nPilot
  % generate channels & blockages exactly as in simulateBlockagePoint,
  % but with p_blk = 0.5
  [G,Hr,Hd] = generateChannels(params.M,params.N,params.K);
  I_d = rand(params.K,1)>p_pilot;
  I_r = rand(params.K,1)>p_pilot;
  for k=1:params.K
    if ~I_d(k), Hd(:,k)=0; end
    if ~I_r(k), Hr(:,k)=0; end
  end
  % run your AO‐WMMSE + Manopt (or whatever scheme you want
  % to pick τ for—often you choose the *best* scheme here)
  Theta0 = diag(exp(1i*2*pi*rand(params.N,1)));  % e.g. random SC
  W0     = update_W_WMMSE_closed_form(G,Hr,Hd,Theta0,...
               params.weights,params.P_tx,params.noise);
  wsr_pilot(it) = computeWSR(Hd,Hr,G,Theta0,W0,...
                   params.noise,params.weights);
end

%% 2) Choose tau as the median of that distribution
tau = median(wsr_pilot);
fprintf('Chosen outage threshold \\tau = %.3f bps/Hz\n', tau);

% now plug this tau into your main simulate_blockage_effects.m
params.tau = tau;


%%%%%%%%%%%%



% Pre-allocate
WSR_avg   = NaN(numel(p_list),numel(schemes));
Outage    = NaN(size(WSR_avg));

% Main loop over blockage probability
for ip = 1:numel(p_list)
  p_blk = p_list(ip);
  fprintf('Simulating p=%.1f …\n', p_blk);
  for is = 1:numel(schemes)
    scheme = schemes{is};
    [rWSR, rOut] = simulateBlockagePoint(p_blk,scheme,params,nMonte);
    WSR_avg(ip,is) = rWSR;
    Outage(ip,is)  = rOut;
  end
end

% --- Plot Figure 7: Average WSR vs. Blockage Probability ---
figure; hold on; box on;
for is = 1:numel(schemes)
  plot(p_list, WSR_avg(:,is),'LineWidth',1.5);
end
xlabel('Link blockage probability, p');
ylabel('Average WSR (bit/s/Hz)');
legend(schemes,'Location','southwest');
title('Fig.7: WSR vs. Blockage');
ylim([0, max(WSR_avg(:))*1.1]);

% --- Plot Figure 8: Outage vs. Blockage Probability ---
figure; hold on; box on;
for is = 1:numel(schemes)
  plot(p_list, Outage(:,is),'LineWidth',1.5);
end
xlabel('Link blockage probability, p');
ylabel('Outage probability (R < 0.5 bps/Hz)');
legend(schemes,'Location','northwest');
title('Fig.8: Outage vs. Blockage');
ylim([0,1]);
