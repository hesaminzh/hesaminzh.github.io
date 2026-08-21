%% Outage_3D.m
% Generates Fig. 9 (3D outage)

clear; close all; clc;
rng(3,'twister');

% --- System & Simulation Parameters ---

%% BS parameters settings
Mx = 8;    % TX antenna nums
My = 4;    % TX antenna nums
K = 2;     % user nums
%L = 2;     % NLOS path nums
%SNR = 30; % SNR in dB
M = Mx*My;

%% surface parameters
Nx = 4;                      % x-axis surface element nums
Ny = [1,2,4,8,16];   % y-axis surface element nums
%Ny = [1,8];   % y-axis surface element nums
%Ny = 8;
N  = Nx*Ny;
Ng = 4; % BD-RIS group size

BS_loc = [0,0]; % location of BS
RIS_loc = [60,0]; % location of RIS in m
user_center = [70, 0];
user_range = [0, 10]; % range of user location from RIS in m



%% optimization parameters
%tolerance = 1e-5;
tolerance = 1e-3;

P_max_dBm = 30;
P_tx = 10^(P_max_dBm/10)*1e-3;
%P_max = 1; % W (30dBm)
%y = 10^(xDbm/10)*1e-3;
noise_dBm = -85;
noise = 10^(noise_dBm/10)*1e-3;
weights = ones(K,1); % user weights

nMonte = 500;      % # Monte Carlo trials
%nMonte = 1000;      % # Monte Carlo trials
%tau         = 0.5;      % outage threshold in bps/Hz



%% distances
if exist(['distances__','K_', num2str(K),'.mat'], 'file')
    load(['distances__','K_', num2str(K),'.mat'], 'user_loc', 'user_azimuths', 'd_BR', 'd_RU', 'd_BU');
else
    [user_loc, user_azimuths, d_BR, d_RU, d_BU] = generate_user_location(user_range, user_center, RIS_loc, BS_loc, K);
    save(['distances__','K_', num2str(K),'.mat'], 'user_loc', 'user_azimuths', 'd_BR', 'd_RU', 'd_BU');
end
%
%




%% channels

% Check if the channel matrices have already been computed and saved
if exist(['channels__','K_', num2str(K),'__', num2str(nMonte),'.mat'],'file')
    load(['channels__','K_', num2str(K),'__', num2str(nMonte),'.mat'], 'H_BS_RIS_all', 'H_RIS_user_all', 'H_BS_user_all');
else
    % Precompute channel matrices for all iterations and Ny values
    H_BS_RIS_all = cell(nMonte, length(Ny));
    H_RIS_user_all = cell(nMonte, length(Ny));
    H_BS_user_all = cell(nMonte, length(Ny));

    for iter = 1:nMonte
        for mi = 1:length(Ny)
            [H_BS_user_all{iter,mi}, H_BS_RIS_all{iter,mi}, H_RIS_user_all{iter,mi}] = generate_SV_channels(M, Nx, Ny(mi), K, d_BR, d_RU, d_BU, user_azimuths);
        end
    end
    % Save the precomputed channels to a .mat file
    save(['channels__','K_', num2str(K),'__', num2str(nMonte),'.mat'], 'H_BS_RIS_all', 'H_RIS_user_all', 'H_BS_user_all');
end









%p_list         = 0:0.1:1;  % blockage probabilities
%schemes = { ...
%  'FC' ...               % fully-connected RIS (AO)
%  'GC', ...              % group-connected RIS (AO, Ng=4)
%  'SC', ...              % single-connected RIS (optimized AO)
%  'NoRIS', ...           % direct path only
%};
schemes = {'FC', 'GC', 'SC'};


p_min = 0;    % 0% blocked at d_max
p_max = 0.99;    % 99% blocked at d_max
numPts = 10;

p_RU_max_list = linspace(p_min, p_max, numPts);  
% = [0.1, 0.3, 0.5, 0.7, 0.9]


%d_RU = [10.1161; 2.9741];    % your two distances
d_max = max(d_RU);

% For each target p at d_max, solve  p = 1 - exp(-η * d_max )
eta_list = -log(1 - p_RU_max_list) / d_max;  


% rows: the 5 η’s;  columns: the two users
p_RU = (1 - exp(-eta_list .* d_RU)).';  


%%%%

%nMonte = 500;
%nMonte = 1000;

numEta     = numel(eta_list);
numSchemes = numel(schemes);

% 1) Define your RIS sizes
%N_list    = [4, 8, 16, 32, 64];
N_list    = N;
numN      = numel(N_list);
%R_thresh = tau;


%%% find R_thresh
%nPilot = 200;
nPilot = 500;

% This will hold every single run’s WSR for each (scheme,N,r)
% Dimensions: (scheme, N, trial)
WSR_all = NaN(numN, numSchemes, nPilot);

% wrap your precomputed channels in constants
Hbsr_const = parallel.pool.Constant(H_BS_RIS_all);
Hrus_const = parallel.pool.Constant(H_RIS_user_all);
Hbuu_const = parallel.pool.Constant(H_BS_user_all);

if isempty(gcp('nocreate'))
    parpool('local');  % use up to 64 workers if available
end

%% 2) Parallelize over the Monte-Carlo index
parfor r = 1:nPilot
    % local buffer for this trial
    wsr_run = NaN(numN, numSchemes);

    % pull the r-th snapshot of all channel sets
    for iN = 1:numN
        H_BS_RIS_template = Hbsr_const.Value{r,iN};
        for is = 1:numSchemes
            sch = schemes{is};
            % copy and (if you want) sample blockages here,
            % or set eta=0 for the no-blockage pilot
            H_BS_RIS   = H_BS_RIS_template;
            H_RIS_user = Hrus_const.Value{r,iN};
            H_BS_user  = Hbuu_const.Value{r,iN};

            % run your AO or MRT/WMMSE
            switch sch
                case 'FC'
                    [wsr_run(iN,is), ~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(...
                        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights);
                case 'GC'
                    [wsr_run(iN,is), ~] = Alternating_Optimization_WMMSE_CF_Manopt_GC(...
                        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights, Ng);
                %case 'GC_2'
                %    [wsr_run(iN,is), ~] = Alternating_Optimization_WMMSE_CF_Manopt_GC(...
                %        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights, 2);
                case 'SC'
                    [wsr_run(iN,is), ~] = Alternating_Optimization_WMMSE_CF_Manopt_SC(...
                        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights);
                case 'NoRIS'
                    % MRT on (possibly blocked) direct channel only
                    %W = computeMRT(H_BS_user, P_tx);
                    %Theta = zeros(N);
                    %wsr_run(iN,is) = computeWSR(H_BS_user, H_RIS_user, H_BS_RIS, Theta, W, noise, weights);
                    W_initial = randn(M, K) + 1i * randn(M, K);
                    W = W_initial / norm(W_initial, 'fro') * sqrt(P_tx);
                    [~, wsr_run(iN,is)] = WMMSE_algorithm_closed_form(H_BS_user, W, tolerance, P_tx, K, noise, weights);

                case 'RandPhase'
                    % Random-phase RIS (SC)
                    Theta = diag(exp(1i*2*pi*rand(N,1)));
                    %W     = update_W_WMMSE_closed_form(H_BS_RIS, H_RIS_user, H_BS_user, Theta, ...
                    %        weights, P_tx, noise);
                    %wsr_run(iN,is) = computeWSR(H_BS_user, H_RIS_user, H_BS_RIS, Theta, W, noise, weights);
                    W_initial = randn(M, K) + 1i * randn(M, K);
                    W = W_initial / norm(W_initial, 'fro') * sqrt(P_tx);
                    H_eff = H_BS_user' + H_RIS_user' * Theta * H_BS_RIS;  % K x M effective channel
                    C = H_eff';
                    [~, wsr_run(iN,is)] = WMMSE_algorithm_closed_form(C, W, tolerance, P_tx, K, noise, weights);


                otherwise
                    error('Unknown scheme: %s', scheme);
            end

        end
    end

    % store this trial’s results
    WSR_all(:,:,r) = wsr_run;
end

%% 3) Compute the average WSR over Monte-Carlo
avgWSR = mean(WSR_all, 3, 'omitnan');   % size = [numN x numSchemes]

global_min = min(avgWSR(:));
global_max = max(avgWSR(:));
R_thresh = (global_min + global_max)/2;










% 2) Build 4-D index vectors
[etaIdx, schemeIdx, Nidx, monteIdx] = ndgrid(1:numEta, 1:numSchemes, 1:numN, 1:nMonte);

etaIdx    = etaIdx   (:);
schemeIdx = schemeIdx(:);
Nidx      = Nidx     (:);
monteIdx  = monteIdx (:);

T = numel(etaIdx);    % = numEta*numSchemes*numN*nMonte

% 3) Preallocate an outage flag
isOutAll = false(T,1);

% (Optional) constant channel storage if you’ve precomputed for each N 
% H_BS_RIS_all should then be a {numN × nMonte} cell, etc.
Hbsr_const = parallel.pool.Constant(H_BS_RIS_all);
Hrus_const = parallel.pool.Constant(H_RIS_user_all);
Hbuu_const = parallel.pool.Constant(H_BS_user_all);

if isempty(gcp('nocreate')), parpool; end

% 4) Run the big parfor
parfor t = 1:T
    ip    = etaIdx(t);
    is    = schemeIdx(t);
    iN    = Nidx(t);
    im    = monteIdx(t);

    eta    = eta_list(ip);
    scheme = schemes{is};
    N      = N_list(iN);

    % pull channels for this N and trial
    H_BS_RIS   = Hbsr_const.Value{im,iN};
    H_RIS_user = Hrus_const.Value{im,iN};
    H_BS_user  = Hbuu_const.Value{im,iN};

    % sample blockages
    p_BU = 1 - exp(-eta * d_BU);
    p_RU = 1 - exp(-eta * d_RU);
    I_d = rand(K,1) > p_BU;
    I_r = rand(K,1) > p_RU;
    H_BS_user(:,~I_d)  = 0;
    H_RIS_user(:,~I_r) = 0;

    % compute WSR for this scheme
    switch scheme
        case 'FC'
            [wsr, ~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(H_BS_user,H_RIS_user,H_BS_RIS,tolerance,P_tx,noise,weights);
        case 'GC'
            [wsr, ~] = Alternating_Optimization_WMMSE_CF_Manopt_GC(H_BS_user,H_RIS_user,H_BS_RIS,tolerance,P_tx,noise,weights,Ng);
        case 'SC'
            [wsr, ~] = Alternating_Optimization_WMMSE_CF_Manopt_SC(H_BS_user,H_RIS_user,H_BS_RIS,tolerance,P_tx,noise,weights);
        case 'NoRIS'
            % your existing MRT/WMMSE fallback
            W0 = randn(M,K)+1i*randn(M,K);
            W  = W0/norm(W0,'fro')*sqrt(P_tx);
            [~, wsr] = WMMSE_algorithm_closed_form(H_BS_user,W,tolerance,P_tx,K,noise,weights);
    end

    % 5) flag outage if below your threshold R_thresh
    isOutAll(t) = (wsr < R_thresh);
end

% reshape into [η × schemes × N × trials]
out4d = reshape(isOutAll, [numEta, numSchemes, numN, nMonte]);

% average along the 4th dim → [η × schemes × N]
outProb = squeeze( mean(out4d, 4, 'omitnan') );


save(['plot_Outage_3D','__K_',num2str(K),'.mat'], 'R_thresh', 'outProb', 'N_list' , 'eta_list', 'schemes');

% -------------------------------
% 1) Build the grid for plotting
% -------------------------------
[ETA, Ngrid] = meshgrid(eta_list, N_list);
%   ETA(i,j)   = eta_list(j)
%   Ngrid(i,j) = N_list(i)

% -------------------------------
% 2) Predefine colors & handles
% -------------------------------
colors = [ ...
    0.0000 0.4470 0.7410;   % FC    (blue)
    0.8500 0.3250 0.0980;   % GC    (orange)
    0.4660 0.6740 0.1880;   % SC    (green)
    0.4940 0.1840 0.5560    % NoRIS (purple)
];
h = gobjects(1,numSchemes);

% -------------------------------
% 3) Plot each scheme as a surface
% -------------------------------
figure; 
hold on; 
grid on;

for is = 1:numSchemes
    % extract the [numEta × numN] slice and transpose so that
    % rows→N_list,  cols→eta_list
    Z = squeeze( outProb(:, is, :) ).';  

    h(is) = surf( ...
        ETA, Ngrid, Z, ...
        'FaceColor',   colors(is,:), ...
        'EdgeColor',   'w', ...
        'LineWidth',   0.5, ...
        'FaceAlpha',   0.8, ...
        'DisplayName', schemes{is} ...
    );
end

% -------------------------------
% 4) Tidy up the axes & legend
% -------------------------------
colormap([])    % disable automatic colormap
%shading flat    % preserve your solid FaceColor
view(45,30)     % tilt for 3D perspective

xlabel('\eta (blockage density)','FontSize',12)
ylabel('RIS elements, N','FontSize',12)
zlabel('Outage probability','FontSize',12)
title('Outage vs. Blockage Density and RIS Size','FontSize',14)

legend(h, schemes, 'Location','northeast')
rotate3d on     % enable interactive rotation
