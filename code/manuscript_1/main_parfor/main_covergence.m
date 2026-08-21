% UPA_reflecting_surface
clear;close all;clc;
% Set the seed for reproducibility
rng(3, 'twister');

%% BS parameters settings
Mx = 8;    % TX antenna nums
My = 4;    % TX antenna nums
K = 2;     % user nums
%K = 4;     % user nums
%K = 6;     % user nums
%L = 2;     % NLOS path nums
%SNR = 30; % SNR in dB
M = Mx*My;

BS_loc = [0,0]; % location of BS
RIS_loc = [60,0]; % location of RIS in m
user_center = [70, 0];
user_range = [0, 10]; % range of user location from RIS in m


%% surface parameters
Nx = 4;                      % x-axis surface element nums
%Ny = [1,2,4,8,16];   % y-axis surface element nums
Ny = 8;
N  = Nx*Ny;
Ng = 4;

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







%% optimization parameters
%tolerance = 1e-5;
tolerance = 1e-3;
P_max_dBm = 30;
P_tx = 10^(P_max_dBm/10)*1e-3;
%sigma = 3.16e-12; % dBm-->W (-85dBm)
%P_tx = P_max;
noise_dBm = -85;
noise = 10^(noise_dBm/10)*1e-3;
weights = ones(K,1); % user weights

%num_iterations = 500;
%fixedIter = 50;  % fixed number of AO iterations for plotting
fixedIter = 100;  % fixed number of AO iterations for plotting
global_convergence = NaN(num_iterations, fixedIter, 3); % 3 scenarios (FC, SC, GC_4)

% Convert channels to constants for efficient access
Hbsr_const = parallel.pool.Constant(H_BS_RIS_all);
Hrus_const = parallel.pool.Constant(H_RIS_user_all);
Hbuu_const = parallel.pool.Constant(H_BS_user_all);

% Start parallel pool if not active
if isempty(gcp('nocreate'))
    parpool('local');  % use up to 64 workers if available
end


parfor iMonte = 1:num_iterations
    %if mod(iMonte, 100) == 0
    %    fprintf('Monte Carlo Run: %d\n', iMonte);
    %end
    
    % Fixed RIS configuration: use mi = 4 (32 elements)
    H_BS_RIS = Hbsr_const.Value{iMonte, 1};
    H_RIS_user = Hrus_const.Value{iMonte, 1};
    H_BS_user = Hbuu_const.Value{iMonte, 1};
    
    % Pre-allocate local convergence traces for each scenario
    temp_conv = NaN(fixedIter, 4);
    % local_conv_GC4 = NaN(1, fixedIter); % Uncomment if using GC_4
    
    %% Run AO algorithms and post-process their convergence trace
    % Fully Cooperative (FC) scenario
    [~, wsr_all_FC] = Alternating_Optimization_WMMSE_CF_Manopt_FC(H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights);
    len_FC = length(wsr_all_FC);
    if len_FC < fixedIter
        wsr_all_FC = [wsr_all_FC, repmat(wsr_all_FC(end), 1, fixedIter - len_FC)];
    else
        wsr_all_FC = wsr_all_FC(1:fixedIter);
    end
    temp_conv(:, 1) = wsr_all_FC;
    
    
    [~, wsr_all_GC] = Alternating_Optimization_WMMSE_CF_Manopt_GC(H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights,Ng);
    len_GC = length(wsr_all_GC);
    if len_GC < fixedIter
        wsr_all_GC = [wsr_all_GC, repmat(wsr_all_GC(end), 1, fixedIter - len_GC)];
    else
        wsr_all_GC = wsr_all_GC(1:fixedIter);
    end
    temp_conv(:, 2) = wsr_all_GC;
    
    [~, wsr_all_GC_2] = Alternating_Optimization_WMMSE_CF_Manopt_GC(H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights, 2);
    len_GC = length(wsr_all_GC_2);
    if len_GC < fixedIter
        wsr_all_GC_2 = [wsr_all_GC_2, repmat(wsr_all_GC_2(end), 1, fixedIter - len_GC)];
    else
        wsr_all_GC_2 = wsr_all_GC_2(1:fixedIter);
    end
    temp_conv(:, 3) = wsr_all_GC_2;

    % Single connected case (SC) scenario
    [~, wsr_all_SC] = Alternating_Optimization_WMMSE_CF_Manopt_SC(H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights);
    len_SC = length(wsr_all_SC);
    if len_SC < fixedIter
        wsr_all_SC = [wsr_all_SC, repmat(wsr_all_SC(end), 1, fixedIter - len_SC)];
    else
        wsr_all_SC = wsr_all_SC(1:fixedIter);
    end
    temp_conv(:, 4) = wsr_all_SC;
    
    
    % Store the convergence data into the global array
    global_convergence(iMonte, :, :) = temp_conv;
end


%save('convergence__FC_SC.mat', 'fixedIter', 'global_convergence');
save(['convergence__','K_',num2str(K),'__',num2str(num_iterations),'.mat'], 'fixedIter', 'global_convergence');

%% Post-Processing: Average Convergence Behavior
mean_convergence = squeeze(mean(global_convergence, 1, 'omitnan'));
%std_convergence  = squeeze(std(global_convergence, 1));

%% Plotting the Convergence Behavior
iterations = 1:fixedIter;
figure;
hold on;
plot(iterations, mean_convergence(:, 1), 'b', 'LineWidth', 1.5);
plot(iterations, mean_convergence(:, 2), 'g', 'LineWidth', 1.5);
plot(iterations, mean_convergence(:, 3), 'y', 'LineWidth', 1.5);
plot(iterations, mean_convergence(:, 4), 'r', 'LineWidth', 1.5);
%errorbar(iterations, mean_convergence(:, 1), std_convergence(:, 1), 'LineWidth', 1.5);
%errorbar(iterations, mean_convergence(:, 2), std_convergence(:, 2), 'LineWidth', 1.5);
% Uncomment for GC_4:
% errorbar(iterations, mean_convergence(:, 3), std_convergence(:, 3), 'LineWidth', 1.5);
xlabel('AO Iteration');
ylabel('Weighted Sum-Rate (WSR)');
title('Convergence Behavior of AO Algorithms');
legend('FC', 'GC (group size: 4)','GC (group size: 2)', 'SC'); % Add 'GC_4' if needed
grid on;
hold off;
