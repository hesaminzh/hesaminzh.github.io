%% main_simulate_blockage_effects.m
% Generates Fig.7 (WSR vs. blockage) and Fig.8 (Outage vs. blockage)

clear; close all; clc;
rng(3,'twister');

% --- System & Simulation Parameters ---

%% BS parameters settings
Mx = 8;    % TX antenna nums
My = 4;    % TX antenna nums
K = 2;     % user nums
%K = 4;     % user nums
%K = 6;     % user nums
%L = 2;     % NLOS path nums
%SNR = 30; % SNR in dB
M = Mx*My;

%% surface parameters
Nx = 4;                      % x-axis surface element nums
%Ny = [1,2,4,8,16];   % y-axis surface element nums
%Ny = [1,8];   % y-axis surface element nums
Ny = 8;
N  = Nx*Ny;
Ng = 4; % BD-RIS group size

BS_loc = [0,0]; % location of BS
RIS_loc = [60,0]; % location of RIS in m
user_center = [70, 0];
user_range = [0, 10]; % range of user location from RIS in m



%% optimization parameters
%tolerance = 1e-5;
tolerance = 1e-3;

%P_max_dBm = 30;
%P_tx = 10^(P_max_dBm/10)*1e-3;
%P_max = 1; % W (30dBm)
%y = 10^(xDbm/10)*1e-3;
P_dBm = [20,25,30,35,40,45];
P_watts = 10 .^ ((P_dBm) / 10) * 1e-3;
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

%tau = 2.263;
%tau = 2.306;
%tau = 3.263;


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
%schemes = {'GC_2', 'NoRIS'};


%nMonte = 500;
%nMonte = 1000;

% --- Initialization ---
numPower     = numel(P_watts);         % e.g., different power values 
numSchemes = numel(schemes);         % e.g., {'FC','GC','SC'}

% Preallocate full 3D result tensor
WSR_all = NaN(numPower, numSchemes, nMonte);
%WSR_all = NaN(numEta, numSchemes, nMonte);   % [η × scheme × trial]

% Convert channels to constants for efficient access
Hbsr_const = parallel.pool.Constant(H_BS_RIS_all);
Hrus_const = parallel.pool.Constant(H_RIS_user_all);
Hbuu_const = parallel.pool.Constant(H_BS_user_all);

% Start parallel pool if not active
if isempty(gcp('nocreate'))
    parpool('local');  % use up to 64 workers if available
end

% --- Main parallel loop ---
parfor im = 1:nMonte

    % Pull channel instances (1×numEta cell arrays)
    %Hbsr_all = Hbsr_const.Value{im, :};
    %Hrus_all = Hrus_const.Value{im, :};
    %Hbuu_all = Hbuu_const.Value{im, :};
    H_BS_RIS = Hbsr_const.Value{im, :}; 
    H_RIS_user = Hrus_const.Value{im, :};
    H_BS_user = Hbuu_const.Value{im, :};

    % Local result buffer for this trial
    local_wsr = NaN(numPower, numSchemes);

    for ip = 1:numPower
        P_tx = P_watts(ip);

        for is = 1:numSchemes
            scheme = schemes{is};
            tempWSR = NaN;

            switch scheme
                case 'FC'
                    [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(...
                        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights);
                case 'GC'
                    [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_GC(...
                        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights, Ng);
                %case 'GC_2'
                %    [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_GC(...
                %        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights, 2);
                case 'SC'
                    [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_SC(...
                        H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights);
                %case 'NoRIS'
                %    % MRT on (possibly blocked) direct channel only
                %    %W = computeMRT(H_BS_user, P_tx);
                %    %Theta = zeros(N);
                %    %tempWSR = computeWSR(H_BS_user, H_RIS_user, H_BS_RIS, Theta, W, noise, weights);
                %    W_initial = randn(M, K) + 1i * randn(M, K);
                %    W = W_initial / norm(W_initial, 'fro') * sqrt(P_tx);
                %    [~, tempWSR] = WMMSE_algorithm_closed_form(H_BS_user, W, tolerance, P_tx, K, noise, weights);

                %case 'RandPhase'
                %    % Random-phase RIS (SC)
                %    Theta = diag(exp(1i*2*pi*rand(N,1)));
                %    %W     = update_W_WMMSE_closed_form(H_BS_RIS, H_RIS_user, H_BS_user, Theta, ...
                %    %        weights, P_tx, noise);
                %    %tempWSR = computeWSR(H_BS_user, H_RIS_user, H_BS_RIS, Theta, W, noise, weights);
                %    W_initial = randn(M, K) + 1i * randn(M, K);
                %    W = W_initial / norm(W_initial, 'fro') * sqrt(P_tx);
                %    H_eff = H_BS_user' + H_RIS_user' * Theta * H_BS_RIS;  % K x M effective channel
                %    C = H_eff';
                %    [~, tempWSR] = WMMSE_algorithm_closed_form(C, W, tolerance, P_tx, K, noise, weights);


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


WSR_avg = mean(WSR_all, 3, 'omitnan');         % power × scheme

save(['power__', num2str(N),'__K_',num2str(K),'.mat'],'WSR_all' ,'WSR_avg', 'P_watts', 'schemes', 'numSchemes');




figure('DefaultAxesFontSize',12); hold on;
LineW = 1.5;
MarkS = 8;
for is = 1:numel(schemes)
  plot( P_dBm, WSR_avg(:,is), '-o','LineWidth',LineW,'MarkerSize',MarkS,'DisplayName', schemes{is} );
end
xlabel('Transmit power (dBm)', 'Interpreter', 'LaTeX');
ylabel('WSR (bps/Hz)');
legend('Location','best');
grid on

