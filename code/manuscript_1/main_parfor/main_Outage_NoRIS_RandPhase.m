%% main_simulate_blockage_effects.m
% Generates Fig.7 (WSR vs. blockage) and Fig.8 (Outage vs. blockage)

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
%Ny = [1,2,4,8,16];   % y-axis surface element nums
%Ny = [1,8];   % y-axis surface element nums
%Ny = 8;
Ny = 16;
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

%d_BR = 60;
%d_RU = [5;15];
%d_BU = d_BR + d_RU;


%%% tau
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% find the correct value of tau so that outage vs blockage shows the S-curve in the figure.
%nPilot   = 200;          % enough to get a decent median
%
%%% Create channels here then load wherever needed
%H_BS_RIS_pilot   = cell(nPilot, 1);
%H_RIS_user_pilot = cell(nPilot, 1);
%H_BS_user_pilot  = cell(nPilot, 1);
%
%for iPilot = 1:nPilot
%    [H_BS_user_pilot{iPilot,1}, H_BS_RIS_pilot{iPilot,1}, H_RIS_user_pilot{iPilot,1}] = generate_SV_channels(M, Nx, Ny(1), K, d_BR, d_RU, d_BU, user_azimuths);
%end
%
%%
%%
%%% 1) Pilot simulation at p=0.5 to estimate the WSR distribution
%%p_pilot  = 0.5;
%%nPilot   = 200;          % enough to get a decent median
%eta_pilot = log(2)/mean(d_RU);
%p_d = 1 - exp(-eta_pilot*d_BU);
%p_r = 1 - exp(-eta_pilot*d_RU);
%
%if isempty(gcp('nocreate'))
%    parpool('local');
%end
%wsr_pilot = NaN(nPilot,1);
%
%parfor it=1:nPilot
%    % generate channels & blockages exactly as in simulateBlockagePoint,
%  % but with p_blk = 0.5
%    %[H_Bs_RIS,H_RIS_user,H_BS_user] = generateChannels(M,N,K);
%    %% BS-IRS channel matrix
%    H_BS_RIS = H_BS_RIS_pilot{it,1};
%
%    %% IRS-user channel matrix
%    H_RIS_user = H_RIS_user_pilot{it,1};
%
%    %% BS-user channel matrix
%    H_BS_user = H_BS_user_pilot{it,1};
%    % sample blockage masks
%    I_d = rand(K,1) > p_d;
%    I_r = rand(K,1) > p_r;
%    H_BS_user(:,~I_d)   = 0;
%    H_RIS_user(:,~I_r)  = 0;
%    tmpWSR = NaN;
%    % run your AO‐WMMSE + Manopt (or whatever scheme you want
%  % to pick τ for—often you choose the *best* scheme here)
%    %Theta0 = diag(exp(1i*2*pi*rand(N,1)));  % e.g. random SC
%    %W0     = update_W_WMMSE_closed_form(H_Bs_RIS,H_RIS_user,H_BS_user,Theta0,...
%    %weights,P_tx,noise);
%    [tmpWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(H_BS_user,H_RIS_user,H_BS_RIS,tolerance,P_tx, noise, weights);
%    wsr_pilot(it) = tmpWSR;
%    %wsr_pilot(it) = computeWSR(H_BS_user,H_RIS_user,H_Bs_RIS,Theta0,W0,...
%    %noise,weights);
%end
%
%%% 2) Choose tau as the median of that distribution
%tau = median(wsr_pilot, 'omitnan');
%fprintf('Chosen outage threshold \\tau = %.3f bps/Hz\n', tau);
%
%%tau = 2.306;
%%tau = tau + 0.5;
tau = 5.679;


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

%schemes = {'FC', 'GC', 'SC'};
schemes = {'NoRIS','RandPhase'};
%schemes = {'GC_2', 'NoRIS'};

%p_min = 0.1;   % 10 % blocked
%p_max = 0.9;   % 90 % blocked
%%d_RU = [10.1161,  2.9741];   % in metres
%d_RU_max = max(d_RU);  % the longest RIS–user path
%d_RU_min = min(d_RU);  % the shortest RIS–user path
%%
%eta_low  = -log(1 - p_min) / d_RU_max;  % ensures at η_low, p_RU≈0.1
%eta_high = -log(1 - p_max) / d_RU_min;  % ensures at η_high, p_RU≈0.9
%%numPoints = 9;                                   % however many you like
%numPoints = 5;                                   % however many you like
%eta_list  = linspace(eta_low, eta_high, numPoints);
%eta_list = eta_list.';
%
%%% 5) Check the actual RIS–user blockage sweep
%p_RU = 1 - exp(-eta_list .* max(d_RU));
%p_BU = 1 - exp(-eta_list .* max(d_BU));


%%%% test the eta

%p_min = 0.1;    % 10% blocked at d_max
p_min = 0;    % 0% blocked at d_max
%p_max = 0.9;    % 90% blocked at d_max
p_max = 0.99;    % 99% blocked at d_max
numPts = 10;

p_RU_max_list = linspace(p_min, p_max, numPts);  
% = [0.1, 0.3, 0.5, 0.7, 0.9]


%d_RU = [10.1161; 2.9741];    % your two distances
d_max = max(d_RU);

% For each target p at d_max, solve  p = 1 - exp(-η * d_max )
% eta_list = -log(1 - p_RU_max_list) / d_max;  
eta_list = 0:0.025:0.45;
%last_point = -log(1 - 0.99)/d_max;
% eta_list is now a 1×5 vector
%eta_list = [0 eta_list];
%eta_list = [eta_list last_point];

% rows: the 5 η’s;  columns: the two users
p_RU = 1 - exp(-eta_list.' * d_RU.');  

% e.g. p_RU(i,1) = 1-exp(-eta_i * d_RU(1)), 
%      p_RU(i,2) = 1-exp(-eta_i * d_RU(2))


%%%%

%nMonte = 500;
%nMonte = 1000;

% --- Initialization ---
numEta     = numel(eta_list);         % e.g., 7 blockage levels
numSchemes = numel(schemes);         % e.g., {'FC','GC','SC'}

% Preallocate full 3D result tensor
WSR_all = NaN(numEta, numSchemes, nMonte);   % [η × scheme × trial]

% Convert channels to constants for efficient access
Hbsr_const = parallel.pool.Constant(H_BS_RIS_all);
Hrus_const = parallel.pool.Constant(H_RIS_user_all);
Hbuu_const = parallel.pool.Constant(H_BS_user_all);

% Start parallel pool if not active
if isempty(gcp('nocreate'))
    parpool('local', 64);  % use up to 64 workers if available
end

% --- Main parallel loop ---
parfor im = 1:nMonte

    % Pull channel instances (1×numEta cell arrays)
    %Hbsr_all = Hbsr_const.Value{im, :};  % 1×numEta
    %Hrus_all = Hrus_const.Value{im, :};  % 1×numEta
    %Hbuu_all = Hbuu_const.Value{im, :};  % 1×numEta
    H_BS_RIS = Hbsr_const.Value{im, :};  % 1×numEta
    H_RIS_user = Hrus_const.Value{im, :};  % 1×numEta
    H_BS_user = Hbuu_const.Value{im, :};  % 1×numEta

    % Local result buffer for this trial
    local_wsr = NaN(numEta, numSchemes);

    for ip = 1:numEta
        eta = eta_list(ip);

        % Fetch channel matrices for this eta and trial
        %H_BS_RIS   = Hbsr_all{ip};
        %H_RIS_user = Hrus_all{ip};
        %H_BS_user  = Hbuu_all{ip};

        % Apply blockage
        p_BU = 1 - exp(-eta * d_BU);
        p_RU = 1 - exp(-eta * d_RU);
        I_d = rand(K,1) > p_BU;
        I_r = rand(K,1) > p_RU;
        H_BS_user(:,~I_d)  = 0;
        H_RIS_user(:,~I_r) = 0;

        for is = 1:numSchemes
            scheme = schemes{is};
            tempWSR = NaN;

            switch scheme
                %case 'FC'
                %    [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(...
                %        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights);
                %case 'GC'
                %    [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_GC(...
                %        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights, Ng);
                %%case 'GC_2'
                %%    [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_GC(...
                %%        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights, 2);
                %case 'SC'
                %    [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_SC(...
                %        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights);
                case 'NoRIS'
                    % MRT on (possibly blocked) direct channel only
                    %W = computeMRT(H_BS_user, P_tx);
                    %Theta = zeros(N);
                    %tempWSR = computeWSR(H_BS_user, H_RIS_user, H_BS_RIS, Theta, W, noise, weights);
                    W_initial = randn(M, K) + 1i * randn(M, K);
                    W = W_initial / norm(W_initial, 'fro') * sqrt(P_tx);
                    [~, tempWSR] = WMMSE_algorithm_closed_form(H_BS_user, W, tolerance, P_tx, K, noise, weights);

                case 'RandPhase'
                    % Random-phase RIS (SC)
                    Theta = diag(exp(1i*2*pi*rand(N,1)));
                    %W     = update_W_WMMSE_closed_form(H_BS_RIS, H_RIS_user, H_BS_user, Theta, ...
                    %        weights, P_tx, noise);
                    %tempWSR = computeWSR(H_BS_user, H_RIS_user, H_BS_RIS, Theta, W, noise, weights);
                    W_initial = randn(M, K) + 1i * randn(M, K);
                    W = W_initial / norm(W_initial, 'fro') * sqrt(P_tx);
                    H_eff = H_BS_user' + H_RIS_user' * Theta * H_BS_RIS;  % K x M effective channel
                    C = H_eff';
                    [~, tempWSR] = WMMSE_algorithm_closed_form(C, W, tolerance, P_tx, K, noise, weights);


                otherwise
                    error('Unknown scheme: %s', scheme);
            end

            % Store result for current eta/scheme/trial
            local_wsr(ip, is) = tempWSR;
        end
    end

    % Assign slice into global tensor (sliced by 'im')
    WSR_all(:, :, im) = local_wsr;
end


WSR_avg = mean(WSR_all, 3, 'omitnan');         % η × scheme
Outage  = sum(WSR_all < tau, 3, 'omitnan') ./ ...
          sum(~isnan(WSR_all), 3, 'omitnan');

%save(['plot_Outage_FC_GC_SC__N_', num2str(N),'__K_',num2str(K),'.mat'],'WSR_all' ,'WSR_avg', 'Outage', 'eta_list', 'schemes', 'numSchemes');
save(['plot_Outage_NoRIS_RandPhase__N_', num2str(N),'__K_',num2str(K),'.mat'],'WSR_all' ,'WSR_avg', 'Outage', 'eta_list', 'schemes', 'numSchemes');




figure; hold on
for is = 1:numel(schemes)
  plot( eta_list, WSR_avg(:,is), '-o', 'DisplayName', schemes{is},'LineWidth',1.5 );
end
xlabel('$\eta$ (blockage density)', 'Interpreter', 'LaTeX');
ylabel('Average WSR');
legend('Location','best');
grid on


figure; hold on;
for is = 1:numel(schemes)
  plot( eta_list, Outage(:,is), '-o', 'DisplayName', schemes{is}, 'LineWidth',1.5 );
end
xlabel('$\eta$ (blockage density)', 'Interpreter', 'LaTeX');
ylabel('Outage Probability');
legend('Location','best');
grid on
set(gca, 'YScale', 'log');
