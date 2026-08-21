% UPA_reflecting_surface
clear;close all;clc;
% Set the seed for reproducibility
rng(3, 'twister');

%% BS parameters settings
% Mx = 8;    % TX antenna nums
Mx = 4;    % TX antenna nums
My = 4;    % TX antenna nums
K = 2;     % user nums
L = 2;     % NLOS path nums
%SNR = 30; % SNR in dB
M = Mx*My;

%% surface parameters
Nx = 4;                      % x-axis surface element nums
Ny = [1,2,4,8,16];   % y-axis surface element nums
%Ny = [1,8];   % y-axis surface element nums
N  = Nx*Ny;

BS_loc = [0,0]; % location of BS
RIS_loc = [60,0]; % location of RIS in m
user_center = [70, 0];
user_range = [0, 10]; % range of user location from RIS in m
[user_loc, user_azimuths, d_BR, d_RU, d_BU] = generate_user_location(user_range, user_center, RIS_loc, BS_loc, K);





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
    load('distances.mat', 'd_BR', 'd_RU', 'd_BU');
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

% num_iterations = 1000;
num_iterations = 500;
%[H_BS_user, H_BS_RIS, H_RIS_user] = generate_SV_channels(M, Nx, Ny(mi), K, d_BR, d_RU, d_BU, user_azimuths);
[H_BS_user_all{iter,mi}, H_BS_RIS_all{iter,mi}, H_RIS_user_all{iter,mi}] = generate_SV_channels(M, Nx, Ny(mi), K, d_BR, d_RU, d_BU, user_azimuths);

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


global_wsr = NaN(num_iterations, length(Ny), 3); % 3 corresponds to the number of scenarios (FC, SC, GC_4)
%global_time = NaN(num_iterations, length(Ny), 3); % 3 corresponds to the number of scenarios (FC, SC, GC_4)

%maxCores = feature('numCores');
%parpool('local', maxCores);
%parpool; % initialize parallel pool with all available cores
parfor iter = 1:num_iterations
    if mod(iter,100) == 0
        fprintf('Iteration: %d\n', iter);
    end

    local_wsr = NaN(length(Ny), 3);
    %local_time = NaN(length(Ny), 3);
    for mi = 1:length(Ny)

        % Index precomputed channel matrices
        %% BS-IRS channel matrix
        H_BS_RIS = H_BS_RIS_all{iter, mi};

        %% IRS-user channel matrix
        H_RIS_user = H_RIS_user_all{iter, mi};

        %% BS-user channel matrix
        H_BS_user = H_BS_user_all{iter, mi};

        %% AO algorithm

        tic;
        [local_wsr(mi, 1), ~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(H_BS_user,H_RIS_user,H_BS_RIS,tolerance,SNR, noise, weights);
        %local_time(mi, 1) = toc;

        %tic;
        %[local_wsr(mi, 2), ~] = Alternating_Optimization_WMMSE_CF_Manopt_SC(H_BS_user,H_RIS_user,H_BS_RIS,tolerance,SNR, noise, weights);
        %local_time(mi, 2) = toc;

        %tic;
        %[local_wsr(mi, 3), ~] = Alternating_Optimization_WMMSE_CF_Manopt_GC_4(H_BS_user,H_RIS_user,H_BS_RIS,tolerance,SNR, noise, weights);
        %local_time(mi, 3) = toc;




        %% Check if wsr is NaN
        %if isnan(wsr)
        %    warning('NaN result encountered for N = %d in iteration %d\n. Skipping this iteration', Nx*Ny(mi), iter);
        %    continue; % Skip this iteration and move to the next one
        %end

        
        %fprintf('WSR: %d ---- N: %d ---- Iteration: %d\n', wsr, Nx * Ny(mi), iter);
        %fprintf('WSR_FC: %d ---- N: %d ---- Iteration: %d\n', local_wsr(mi, 1), Nx * Ny(mi), iter);
        %fprintf('WSR_SC: %d ---- N: %d ---- Iteration: %d\n', local_wsr(mi, 2), Nx * Ny(mi), iter);
        %fprintf('WSR_GC_4: %d ---- N: %d ---- Iteration: %d\n', local_wsr(mi, 3), Nx * Ny(mi), iter);
        %disp('-----------------------------------------------');

    end

    % Store the local results in the global variable
    global_wsr(iter,:,:) = local_wsr;
    %global_time(iter,:,:) = local_time;

end

% Save results
%save('WSR_vs_N__FC_SC.mat')
save('WSR_vs_N__FC_4_32.mat')

%delete(gcp('nocreate')); % Shut down the parallel pool

% Compute the mean sum-rate for each Ny value, ignoring NaNs
mean_wsr_FC = mean(global_wsr(:,:,1), 1, 'omitnan');
%mean_wsr_SC = mean(global_wsr(:,:,2), 1, 'omitnan');
%mean_wsr_GC_4 = mean(global_wsr(:,:,3), 1, 'omitnan');

%mean_time_FC = mean(global_time(:,:,1), 1, 'omitnan');
%mean_time_SC = mean(global_time(:,:,2), 1, 'omitnan');
%mean_time_GC_4 = mean(global_time(:,:,3), 1, 'omitnan');



figure(1);
%plot(N, mean_wsr_FC, 'b-.o', 'LineWidth', 1, 'MarkerSize', 8);
hold on;
plot(N, mean_wsr_FC, 'b--s', 'LineWidth', 1, 'MarkerSize', 8);
%plot(N, mean_wsr_SC, 'r--s', 'LineWidth', 1, 'MarkerSize', 8);
%plot(N, mean_wsr_GC_4, 'g-^', 'LineWidth', 1, 'MarkerSize', 8);


set(gca,'xtick',N);
xlim([N(1), N(end)]);
xlabel('N'),ylabel('WSR'); %xlim([1, iter]);
grid on;
%handle = legend('SC_Manopt_FP', 'GC_4_Manopt_FP');
%handle = legend('FC', 'SC');
%handle=legend('Proposed, $\mathcal{F}_{\rm c}$','Proposed, $\mathcal{F}_{\rm d}$, $b=1$',...
    %'Proposed, $\mathcal{F}_{\rm d}$, $b=2$', 'ZF+Random PBF');
set(handle, 'Interpreter', 'LaTex');
set(handle, 'FontSize', 13, 'Location', 'best');
set(gca, 'FontSize', 12);



%figure(2);
%%plot(N, mean_wsr_FC, 'b-.o', 'LineWidth', 1, 'MarkerSize', 8);
%hold on;
%plot(N, mean_time_FC, 'b--s', 'LineWidth', 1, 'MarkerSize', 8);
%plot(N, mean_time_SC, 'r--s', 'LineWidth', 1, 'MarkerSize', 8);
%%plot(N, mean_time_GC_4, 'g-^', 'LineWidth', 1, 'MarkerSize', 8);
