 UPA_reflecting_surface
clear;close all;clc;
% Set the seed for reproducibility
rng(3, 'twister');

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
P_max = 10^(P_max_dBm/10)*1e-3;
%P_max = 1; % W (30dBm)
%y = 10^(xDbm/10)*1e-3;
noise_dBm = -85;
noise = 10^(noise_dBm/10)*1e-3;
%sigma = 3.16e-12; % dBm-->W (-85dBm)
P_tx = P_max;
%noise = sigma;
weights = ones(K,1); % user weights




%% distances
if exist(['distances__','K_', num2str(K),'.mat'], 'file')
    load(['distances__','K_', num2str(K),'.mat'], 'user_loc', 'user_azimuths', 'd_BR', 'd_RU', 'd_BU');
else
    [user_loc, user_azimuths, d_BR, d_RU, d_BU] = generate_user_location(user_range, user_center, RIS_loc, BS_loc, K);
    save(['distances__','K_', num2str(K),'.mat'], 'user_loc', 'user_azimuths', 'd_BR', 'd_RU', 'd_BU');
end

%num_iterations = 500;
num_iterations = 1000;


% Check if the channel matrices have already been computed and saved
if exist(['channels__','K_', num2str(K),'__', num2str(num_iterations),'.mat'],'file')
    load(['channels__','K_', num2str(K),'__', num2str(num_iterations),'.mat'], 'H_BS_RIS_all', 'H_RIS_user_all', 'H_BS_user_all');
else
    % Precompute channel matrices for all iterations and Ny values
    H_BS_RIS_all = cell(num_iterations, length(Ny));
    H_RIS_user_all = cell(num_iterations, length(Ny));
    H_BS_user_all = cell(num_iterations, length(Ny));

    for iter = 1:num_iterations
        for mi = 1:length(Ny)
            [H_BS_user_all{iter,mi}, H_BS_RIS_all{iter,mi}, H_RIS_user_all{iter,mi}] = generate_SV_channels(M, Nx, Ny(mi), K, d_BR, d_RU, d_BU, user_azimuths);
        end
    end
    % Save the precomputed channels to a .mat file
    save(['channels__','K_', num2str(K),'__', num2str(num_iterations),'.mat'], 'H_BS_RIS_all', 'H_RIS_user_all', 'H_BS_user_all');
end


% outside parfor
Hbsr_const  = parallel.pool.Constant(H_BS_RIS_all);
Hrus_const  = parallel.pool.Constant(H_RIS_user_all);
Hbuu_const  = parallel.pool.Constant(H_BS_user_all);



%global_wsr = NaN(num_iterations, length(Ny), 3); % 3 corresponds to the number of scenarios (FC, SC, GC_4)
%global_time = NaN(num_iterations, length(Ny), 3); % 3 corresponds to the number of scenarios (FC, SC, GC_4)

%maxCores = feature('numCores');
%parpool('local', maxCores);
%parpool; % initialize parallel pool with all available cores

schemes = {'NoRIS', 'RandPhase'};
numSchemes = numel(schemes);
numElements = length(Ny);
numIter = num_iterations;

[iterIdx, miIdx, schIdx] = ndgrid(1:numIter, 1:numElements, 1:numSchemes);
iterIdx = iterIdx(:);
miIdx   = miIdx(:);
schIdx  = schIdx(:);




if isempty(gcp('nocreate'))
    parpool('local')
end

res      = NaN(numIter, numElements, numSchemes);
res_time = NaN(size(res));

parfor it = 1:numIter
    localRes      = NaN(numElements, numSchemes);
    localResTime  = NaN(numElements, numSchemes);
    for mi = 1:numElements
        % pull constants once per (it,mi)
        H_BS_RIS   = Hbsr_const.Value{it, mi};
        H_RIS_user = Hrus_const.Value{it, mi};
        H_BS_user  = Hbuu_const.Value{it, mi};
        for sch = 1:numSchemes
            tStart = tic;
            switch schemes{sch}
                case 'NoRIS'
                    W_initial = randn(M, K) + 1i * randn(M, K);
                    W = W_initial / norm(W_initial, 'fro') * sqrt(P_tx);
                    noise_NoRIS=ones(K,1) .* noise;
                    [~, wsr] = WMMSE_algorithm_closed_form(H_BS_user, W, tolerance, P_tx, K, noise_NoRIS, weights);
                case 'RandPhase'
                    N_RandPhase = N(mi);
                    Theta = diag(exp(1i*2*pi*rand(N_RandPhase,1)));
                    W_initial = randn(M, K) + 1i * randn(M, K);
                    W = W_initial / norm(W_initial, 'fro') * sqrt(P_tx);
                    noise_RandPhase=ones(K,1) .* noise;
                    H_eff = H_BS_user' + H_RIS_user' * Theta * H_BS_RIS;  % K x M effective channel
                    C = H_eff';
                    [~, wsr] = WMMSE_algorithm_closed_form(C, W, tolerance, P_tx, K, noise_RandPhase, weights);
            end
            localRes(mi, sch)     = wsr;
            if ~isnan(wsr)
                    localResTime(mi, sch) = toc(tStart);
            end
        end
    end
    res(it,:,:)      = localRes;
    res_time(it,:,:) = localResTime;
end

%save('fig1_other_benchmarks.mat', 'res', 'res_time', 'schemes','numSchemes', 'Ny');
save('fig1_NoRIS_RandPhase.mat', 'res', 'res_time', 'schemes','numSchemes', 'Ny');

%schemes = {'FC', 'GC', 'SC'};
%numSchemes = numel(schemes);

% Compute mean over the first (iteration) dimension
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numMi × numSchemes]
mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numMi × numSchemes]

figure; hold on;
for sch = 1:numSchemes
    plot(Ny, mean_wsr(:,sch), '-o', ...
         'DisplayName', schemes{sch}, 'LineWidth',1.5);
end
xlabel('Number of RIS elements, $N$', 'Interpreter', 'LaTeX');
ylabel('Weighted sum-rate (bps/Hz)');
legend('Location','best');
grid on;

figure; hold on;
for sch = 1:numSchemes
    plot(Ny, mean_time(:,sch), '-o', ...
         'DisplayName', schemes{sch}, 'LineWidth',1.5);
end
xlabel('Number of RIS elements, $N$', 'Interpreter', 'LaTeX');
ylabel('Average CPT time (sec)');
legend('Location','best');
grid on;

