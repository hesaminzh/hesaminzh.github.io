% UPA_reflecting_surface
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
SNR = P_max;
%noise = sigma;
weights = ones(K,1); % user weights




%% distances
if exist(['distances__','K_', num2str(K),'.mat'], 'file')
    load(['distances__','K_', num2str(K),'.mat'], 'user_loc', 'user_azimuths', 'd_BR', 'd_RU', 'd_BU');
else
    [user_loc, user_azimuths, d_BR, d_RU, d_BU] = generate_user_location(user_range, user_center, RIS_loc, BS_loc, K);
    save(['distances__','K_', num2str(K),'.mat'], 'user_loc', 'user_azimuths', 'd_BR', 'd_RU', 'd_BU');
end

num_iterations = 500;


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





%global_wsr = NaN(num_iterations, length(Ny), 3); % 3 corresponds to the number of scenarios (FC, SC, GC_4)
%global_time = NaN(num_iterations, length(Ny), 3); % 3 corresponds to the number of scenarios (FC, SC, GC_4)

%maxCores = feature('numCores');
%parpool('local', maxCores);
%parpool; % initialize parallel pool with all available cores

schemes = {'FC', 'GC', 'SC'};
numSchemes = numel(schemes);
numElements = length(Ny);
numIter = num_iterations;
totalJobs = numIter * numElements * numSchemes;

res = NaN(totalJobs,1);
res_time = NaN(totalJobs,1);



if isempty(gcp('nocreate'))
    parpool('local')
end

parfor job = 1:totalJobs
    if mod(job,100) == 0
        fprintf('Iteration: %d\n', job);
    end
    % 1) map linear index to 3D-subscript
    [it, mi, sch] = ind2sub([numIter, numElements, numSchemes], job);

    % 2) pull out the channel matrices for this iteration & Mi
    H_BS_RIS  = H_BS_RIS_all{ it, mi };
    H_RIS_user = H_RIS_user_all{ it, mi };
    H_BS_user  = H_BS_user_all{ it, mi };

    % 3) run exactly one of the three A-O routines
    tStart = tic;
    switch schemes{sch}
        case 'FC'
            [wsr,~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(H_BS_user, H_RIS_user, H_BS_RIS, tolerance, SNR, noise, weights);
        case 'SC'
            [wsr,~] = Alternating_Optimization_WMMSE_CF_Manopt_SC(H_BS_user, H_RIS_user, H_BS_RIS, tolerance, SNR, noise, weights);
        case 'GC'
            arch_param = Ng;
            [wsr,~] = Alternating_Optimization_WMMSE_CF_Manopt_GC(H_BS_user, H_RIS_user, H_BS_RIS, tolerance, SNR, noise, weights, arch_param);
    end
    res(job) = wsr;

    if ~isnan(wsr)
        res_time(job) = toc(tStart);
    end
end

global_wsr = reshape(res, [numIter, numMi, numSchemes]);
global_time = reshape(res_time, [numIter, numMi, numSchemes]);

% Save results
%save('WSR_vs_N__FC_SC.mat')
save('fig1_proposed_closed_form_K_2.mat', 'global_wsr', 'global_time', 'schemes', 'Ny');
load('fig1_proposed_closed_form_K_2.mat', 'global_wsr', 'global_time', 'schemes', 'Ny');

schemes = {'FC', 'GC', 'SC'};
numSchemes = numel(schemes);

% Compute mean over the first (iteration) dimension
mean_wsr = squeeze( mean(global_wsr, 1, 'omitnan') );  % size = [numMi × numSchemes]
mean_time = squeeze( mean(global_time, 1, 'omitnan') );  % size = [numMi × numSchemes]

figure; hold on;
for sch = 1:numSchemes
    plot(Ny, mean_wsr(:,sch), '-o', ...
         'DisplayName', schemes{sch}, 'LineWidth',1.5);
end
xlabel('Number of IRS elements, N_y');
ylabel('Average sum-rate (bps/Hz)');
legend('Location','best');
grid on;

figure; hold on;
for sch = 1:numSchemes
    plot(Ny, mean_time(:,sch), '-o', ...
         'DisplayName', schemes{sch}, 'LineWidth',1.5);
end
xlabel('Number of IRS elements, N_y');
ylabel('Average sum-rate (bps/Hz)');
legend('Location','best');
grid on;

