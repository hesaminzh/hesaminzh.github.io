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
%[user_loc, user_azimuths, d_BR, d_RU, d_BU] = generate_user_location(user_range, user_center, RIS_loc, BS_loc, K);





%% function to create distances

%d_BR = 50;
%cord = [(2*rand(K,1)-1)+4, 2*rand(K,1)-1]*10;
%for ii=1:K
%    c_IU(ii,:) = cord(ii,:)-[40, 30];
%end

%d_RU = sqrt(sum(c_IU.^2, 2));
%d_BU = sqrt(sum(cord.^2, 2));

%% distances
if exist('distances.mat', 'file')
    load('distances.mat', 'user_loc', 'user_azimuths', 'd_BR', 'd_RU', 'd_BU');
else
    [user_loc, user_azimuths, d_BR, d_RU, d_BU] = generate_user_location(user_range, user_center, RIS_loc, BS_loc, K);
    save('distances.mat', 'user_loc', 'user_azimuths', 'd_BR', 'd_RU', 'd_BU');
end
%dist = struct();
%dist.d_BR = d_BR;
%dist.d_RU = d_RU;
%dist.d_BU = d_BU;
%
%

num_iterations = 500;
%[H_BS_user, H_BS_RIS, H_RIS_user] = generate_SV_channels(M, Nx, Ny(mi), K, d_BR, d_RU, d_BU, user_azimuths);
%[H_BS_user_all{iter,mi}, H_BS_RIS_all{iter,mi}, H_RIS_user_all{iter,mi}] = generate_SV_channels(M, Nx, Ny(mi), K, d_BR, d_RU, d_BU, user_azimuths);

% Check if the channel matrices have already been computed and saved
if exist('channels.mat', 'file')
    load('channels.mat', 'H_BS_RIS_all', 'H_RIS_user_all', 'H_BS_user_all');
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
    save('channels.mat', 'H_BS_RIS_all', 'H_RIS_user_all', 'H_BS_user_all');
end





%% optimization parameters
%tolerance = 1e-5;
tolerance = 1e-3;
P_max = 1; % W (30dBm)
sigma = 3.16e-12; % dBm-->W (-85dBm)
SNR = P_max;
noise = sigma;
weights = ones(K,1); % user weights


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




parfor job = 1:totalJobs
    if mod(job,100) == 0
        fprintf('Iteration: %d\n', job);
    end
    % 1) map linear index to 3D-subscript
    [it, mi, sch] = ind2sub([numIter, numElements, numSchemes], job);

    % 2) pull out the channel matrices for this iteration & Mi
    H_BS_RIS  = H_BS_RIS_all{ it, mi };
    H_RIS_usr = H_RIS_user_all{ it, mi };
    H_BS_usr  = H_BS_user_all{ it, mi };

    % 3) run exactly one of the three A-O routines
    tStart = tic;
    switch schemes{sch}
        case 'FC'
            [wsr,~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(H_BS_usr, H_RIS_usr, H_BS_RIS, tolerance, SNR, noise, weights);
        case 'SC'
            [wsr,~] = Alternating_Optimization_WMMSE_CF_Manopt_SC(H_BS_usr, H_RIS_usr, H_BS_RIS, tolerance, SNR, noise, weights);
        case 'GC'
            arch_param = Ng;
            [wsr,~] = Alternating_Optimization_WMMSE_CF_Manopt_GC(H_BS_usr, H_RIS_usr, H_BS_RIS, tolerance, SNR, noise, weights, arch_param);
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

