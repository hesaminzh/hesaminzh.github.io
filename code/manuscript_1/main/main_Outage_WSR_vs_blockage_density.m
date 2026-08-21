%% main_simulate_blockage_effects.m
% Generates Fig.7 (WSR vs. blockage) and Fig.8 (Outage vs. blockage)

clear; close all; clc;
rng(3,'twister');

params = para_init();
% --- System & Simulation Parameters ---
M              = 32;       % # BS antennas
Mx             = 8;
My             = 4;
N              = 16;       % # RIS elements
Nx             = 4;
Ny             = 4;
K              = 2;        % # users
L              = 2;        % NLOS path nums
Ng             = 4;        % group size for GC
P_tx           = 1;        % total BS power (normalized)
noise          = 3.16e-12; % noise variance
tolerance      = 1e-3;
weights        = ones(K,1);
%tau         = 0.5;      % outage threshold in bps/Hz
nMonte         = 500;      % # Monte Carlo trials
%p_list         = 0:0.1:1;  % blockage probabilities

%schemes = { ...
%  'FC' ...               % fully-connected RIS (AO)
%  'GC', ...              % group-connected RIS (AO, Ng=4)
%  'SC', ...              % single-connected RIS (optimized AO)
%  'RandPhase', ...       % random-phase RIS (SC)  
%  'NoRIS', ...           % direct path only
%};

schemes = { ...
  'FC' ...               % fully-connected RIS (AO)
  'GC', ...              % group-connected RIS (AO, Ng=4)
  'SC', ...              % single-connected RIS (optimized AO)
  'NoRIS', ...           % direct path only
};



%%% distances
%if exist('distances.mat', 'file')
%    load('distances.mat', 'd_BR', 'd_RU', 'd_BU');
%else
%    d_BR = 50;
%    cord = [(2*rand(K,1)-1)+4, 2*rand(K,1)-1]*10;
%    for ii=1:K
%        c_IU(ii,:) = cord(ii,:)-[40, 30];
%    end
%
%    d_RU = sqrt(sum(c_IU.^2, 2));
%    d_BU = sqrt(sum(cord.^2, 2));
%    save('distances.mat', 'd_BR', 'd_RU', 'd_BU');
%end

para.BS_loc = [0,0]; % location of BS
para.RIS_loc = [50,0]; % location of RIS in m
para.user_center = [60, 0];
para.user_range = [0, 10]; % range of user location from RIS in m



%% distances
if exist('distances.mat', 'file')
    load('distances.mat', 'd_BR', 'd_RU', 'd_BU');
else
    [user_loc, user_angle, d_BR, d_RU, d_BU] = generate_user_location(para,K);
    save('distances.mat', 'd_BR', 'd_RU', 'd_BU');
end
%
%

%%% Blockage density
%% 1) Desired window of RIS–user blockage probabilities
p_min = 0.1;   % 10 % blocked
p_max = 0.9;   % 90 % blocked

%% 2) Your RIS–user distances
%d_RU = [10.1161,  2.9741];   % in metres

%% 3) Invert p = 1−exp(−η⋅d) to get η_low and η_high
d_RU_max = max(d_RU);  % the longest RIS–user path
d_RU_min = min(d_RU);  % the shortest RIS–user path

eta_low  = -log(1 - p_min) / d_RU_max;  % ensures at η_low, p_RU≈0.1
eta_high = -log(1 - p_max) / d_RU_min;  % ensures at η_high, p_RU≈0.9

%% 4) Build your η grid
numPoints = 9;                                   % however many you like
eta_list  = linspace(eta_low, eta_high, numPoints);

%% 5) Check the actual RIS–user blockage sweep
p_RU = 1 - exp(-eta_list .* d_RU);
p_BU = 1 - exp(-eta_list .* d_BU);



%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the correct value of tau so that outage vs blockage shows the S-curve in the figure.
nPilot   = 200;          % enough to get a decent median
%nPilot   = 1;          % enough to get a decent median

%% Create channels here then load wherever needed
H_BS_RIS_all   = cell(nPilot, 1);
H_RIS_user_all = cell(nPilot, 1);
H_BS_user_all  = cell(nPilot, 1);

for iPilot = 1:nPilot
    H_BS_RIS_all{iPilot,1}   = BS_IRS(Nx, Ny(1), Mx, My, L, d_BR);
    H_RIS_user_all{iPilot,1} = FF_User(Nx, Ny(1), K, d_RU);
    H_BS_user_all{iPilot,1}  = BS_UE_v2(Mx, My, K, L, d_BU);
end

%
%
%% 1) Pilot simulation at p=0.5 to estimate the WSR distribution
%p_pilot  = 0.5;
%nPilot   = 200;          % enough to get a decent median
eta_pilot = log(2)/mean(d_RU);
p_d = 1 - exp(-eta_pilot*d_BU);
p_r = 1 - exp(-eta_pilot*d_RU);

if isempty(gcp('nocreate'))
    parpool('local');
end
wsr_pilot = NaN(nPilot,1);

parfor it=1:nPilot
    % generate channels & blockages exactly as in simulateBlockagePoint,
  % but with p_blk = 0.5
    %[H_Bs_RIS,H_RIS_user,H_BS_user] = generateChannels(M,N,K);
    %% BS-IRS channel matrix
    H_BS_RIS = H_BS_RIS_all{it,1};

    %% IRS-user channel matrix
    H_RIS_user = H_RIS_user_all{it,1};

    %% BS-user channel matrix
    H_BS_user = H_BS_user_all{it,1};
    I_d = rand(K,1)>p_d;
    I_r = rand(K,1)>p_r;
    for k=1:K
        if ~I_d(k), H_BS_user(:,k)  = 0; end
        if ~I_r(k), H_RIS_user(:,k) = 0; end
    end
    tmpWSR = NaN;
    % run your AO‐WMMSE + Manopt (or whatever scheme you want
  % to pick τ for—often you choose the *best* scheme here)
    %Theta0 = diag(exp(1i*2*pi*rand(N,1)));  % e.g. random SC
    %W0     = update_W_WMMSE_closed_form(H_Bs_RIS,H_RIS_user,H_BS_user,Theta0,...
    %weights,P_tx,noise);
    [tmpWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(H_BS_user,H_RIS_user,H_BS_RIS,tolerance,P_tx, noise, weights);
    wsr_pilot(it) = tmpWSR;
    %wsr_pilot(it) = computeWSR(H_BS_user,H_RIS_user,H_Bs_RIS,Theta0,W0,...
    %noise,weights);
end

%% 2) Choose tau as the median of that distribution
tau = median(wsr_pilot, 'omitnan');
fprintf('Chosen outage threshold \\tau = %.3f bps/Hz\n', tau);

% now plug this tau into your main simulate_blockage_effects.m
%params.tau = tau;
%
%tau = 5.679;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create channels here and load them inside the loop
if exist('channels.mat', 'file')
    load('channels.mat', 'H_BS_RIS_all', 'H_RIS_user_all', 'H_BS_user_all');
else
    H_BS_RIS_all   = cell(nMonte, 1);
    H_RIS_user_all = cell(nMonte, 1);
    H_BS_user_all  = cell(nMonte, 1);

    for iMonte = 1:nMonte
        H_BS_RIS_all{iMonte,1}   = BS_IRS(Nx, Ny(1), Mx, My, L, d_BR);
        H_RIS_user_all{iMonte,1} = FF_User(Nx, Ny(1), K, d_RU);
        H_BS_user_all{iMonte,1}  = BS_UE_v2(Mx, My, K, L, d_BU);
    end
    save('channels.mat', 'H_BS_RIS_all', 'H_RIS_user_all', 'H_BS_user_all');
end


% Pre-allocate
WSR_avg   = NaN(numel(eta_list),numel(schemes));
Outage    = NaN(size(WSR_avg));

% Main loop over blockage probability
if isempty(gcp('nocreate'))
    parpool('local');
end

for ip = 1:numel(eta_list)
  eta = eta_list(ip);
  fprintf('Simulating p=%.1f …\n', eta);
  for is = 1:numel(schemes)
    scheme = schemes{is};
    %[rWSR, rOut] = simulateBlockagePoint(p_blk,scheme,params,channels,nMonte);
    [rWSR, rOut] = simulateBlockageDensity(eta, scheme, M, N, K, P_tx, noise, weights, tolerance, Ng, tau, d_BU, d_RU, H_BS_RIS_all, H_RIS_user_all, H_BS_user_all, nMonte);
    WSR_avg(ip,is) = rWSR;
    Outage(ip,is)  = rOut;
  end
end

%% save the result
save('plot_Outage_FC__N_16', 'WSR_avg', 'Outage', 'eta_list', 'schemes');
%load('plot_Outage_FC__N_16', 'WSR_avg', 'Outage', 'eta_list', 'schemes');

% --- Plot Figure 7: Average WSR vs. Blockage Probability ---
figure; hold on; box on;
for is = 1:numel(schemes)
  plot(eta_list, WSR_avg(:,is),'LineWidth',1.5);
end
xlabel('Link blockage probability, p');
ylabel('Average WSR (bit/s/Hz)');
legend(schemes,'Location','southwest');
title('Fig.7: WSR vs. Blockage');
ylim([0, max(WSR_avg(:))*1.1]);

% --- Plot Figure 8: Outage vs. Blockage Probability ---
figure; hold on; box on;
for is = 1:numel(schemes)
  plot(eta_list, Outage(:,is),'LineWidth',1.5);
end
xlabel('Link blockage probability, p');
ylabel('Outage probability (R < 0.5 bps/Hz)');
legend(schemes,'Location','northwest');
title('Fig.8: Outage vs. Blockage');
ylim([0,1]);
