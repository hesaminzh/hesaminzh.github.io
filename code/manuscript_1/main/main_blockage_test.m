%% main_simulate_blockage_effects.m
% Generates Fig.7 (WSR vs. blockage) and Fig.8 (Outage vs. blockage)

clear; close all; clc;
rng(3,'twister');

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
p_list                = 0:0.1:1;  % blockage probabilities

schemes = { ...
  'FC' ...               % fully-connected RIS (AO)
  %'GC', ...              % group-connected RIS (AO, Ng=4)
  %'SC', ...              % single-connected RIS (optimized AO)
  %'RandPhase', ...       % random-phase RIS (SC)  
  %'NoRIS', ...           % direct path only
};

%% distances
if exist('distances.mat', 'file')
    load('distances.mat', 'd_BR', 'd_RU', 'd_BU');
else
    d_BR = 50;
    cord = [(2*rand(params.K,1)-1)+4, 2*rand(params.K,1)-1]*10;
    for ii=1:params.K
        c_IU(ii,:) = cord(ii,:)-[40, 30];
    end

    d_RU = sqrt(sum(c_IU.^2, 2));
    d_BU = sqrt(sum(cord.^2, 2));
    save('distances.mat', 'd_BR', 'd_RU', 'd_BU');
end
dist = struct();
dist.d_BR = d_BR;
dist.d_RU = d_RU;
dist.d_BU = d_BU;
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the correct value of tau so that outage vs blockage shows the S-curve in the figure.
nPilot   = 200;          % enough to get a decent median
%nPilot   = 1;          % enough to get a decent median

%% Create channels here then load wherever needed
H_BS_RIS_all   = cell(nPilot, 1);
H_RIS_user_all = cell(nPilot, 1);
H_BS_user_all  = cell(nPilot, 1);

for iPilot = 1:nPilot
    H_BS_RIS_all{iPilot,1}   = BS_IRS(params.Nx, params.Ny(1), params.Mx, params.My, params.L, dist.d_BR);
    H_RIS_user_all{iPilot,1} = FF_User(params.Nx, params.Ny(1), params.K, dist.d_RU);
    H_BS_user_all{iPilot,1}  = BS_UE_v2(params.Mx, params.My, params.K, params.L, dist.d_BU);
end

%
%
%% 1) Pilot simulation at p=0.5 to estimate the WSR distribution
p_pilot  = 0.5;
%nPilot   = 200;          % enough to get a decent median
parpool('local');
wsr_pilot = NaN(nPilot,1);

parfor it=1:nPilot
    % generate channels & blockages exactly as in simulateBlockagePoint,
  % but with p_blk = 0.5
    %[H_Bs_RIS,H_RIS_user,H_BS_user] = generateChannels(params.M,params.N,params.K);
    %% BS-IRS channel matrix
    H_BS_RIS = H_BS_RIS_all{it,1};

    %% IRS-user channel matrix
    H_RIS_user = H_RIS_user_all{it,1};

    %% BS-user channel matrix
    H_BS_user = H_BS_user_all{it,1};
    I_d = rand(1,params.K)>p_pilot;
    I_r = rand(1,params.K)>p_pilot;
    for k=1:params.K
        if ~I_d(k), H_BS_user(:,k)  = 0; end
        if ~I_r(k), H_RIS_user(:,k) = 0; end
    end
    tmpWSR = NaN;
    % run your AO‐WMMSE + Manopt (or whatever scheme you want
  % to pick τ for—often you choose the *best* scheme here)
    %Theta0 = diag(exp(1i*2*pi*rand(params.N,1)));  % e.g. random SC
    %W0     = update_W_WMMSE_closed_form(H_Bs_RIS,H_RIS_user,H_BS_user,Theta0,...
    %params.weights,params.P_tx,params.noise);
    [tmpWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(H_BS_user,H_RIS_user,H_BS_RIS,params.tolerance,params.P_tx, params.noise, params.weights);
    wsr_pilot(it) = tmpWSR;
    %wsr_pilot(it) = computeWSR(H_BS_user,H_RIS_user,H_Bs_RIS,Theta0,W0,...
    %params.noise,params.weights);
end

%% 2) Choose tau as the median of that distribution
tau = median(wsr_pilot, 'omitnan');
fprintf('Chosen outage threshold \\tau = %.3f bps/Hz\n', tau);

% now plug this tau into your main simulate_blockage_effects.m
params.tau = tau;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create channels here and load them inside the loop
if exist('channels.mat', 'file')
    load('channels.mat', 'H_BS_RIS_all', 'H_RIS_user_all', 'H_BS_user_all');
else
    H_BS_RIS_all   = cell(nMonte, 1);
    H_RIS_user_all = cell(nMonte, 1);
    H_BS_user_all  = cell(nMonte, 1);

    for iMonte = 1:nMonte
        H_BS_RIS_all{iMonte,1}   = BS_IRS(params.Nx, params.Ny(1), params.Mx, params.My, params.L, dist.d_BR);
        H_RIS_user_all{iMonte,1} = FF_User(params.Nx, params.Ny(1), params.K, dist.d_RU);
        H_BS_user_all{iMonte,1}  = BS_UE_v2(params.Mx, params.My, params.K, params.L, dist.d_BU);
    end
end
channels.H_BS_RIS_all   = H_BS_RIS_all;
channels.H_RIS_user_all = H_RIS_user_all;
channels.H_BS_user_all  = H_BS_user_all;


% Pre-allocate
WSR_avg   = NaN(numel(p_list),numel(schemes));
Outage    = NaN(size(WSR_avg));

% Main loop over blockage probability
if isempty(gcp('nocreate'))
    parpool('local');
end

for ip = 1:numel(p_list)
  p_blk = p_list(ip);
  fprintf('Simulating p=%.1f …\n', p_blk);
  for is = 1:numel(schemes)
    scheme = schemes{is};
    [rWSR, rOut] = simulateBlockagePoint(p_blk,scheme,params,channels,nMonte);
    WSR_avg(ip,is) = rWSR;
    Outage(ip,is)  = rOut;
  end
end

%% save the result
save('plot_Outage_FC__N_16', 'WSR_avg', 'Outage', 'p_list', 'schemes');
%load('plot_Outage_FC__N_16', 'WSR_avg', 'Outage', 'p_list', 'schemes');

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
