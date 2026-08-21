addpath('/Users/hesami/Documents/Research/ris_isac/Hesamilab/matlab_gpt/bdris_mmWave/bdris_new/data/figure_1/nMonte_1000/');
load('fig1_proposed_closed_form_K_2.mat');
%% Proposed
% Compute mean over the first (iteration) dimension
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numMi × numSchemes]
%mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numMi × numSchemes]

N = [4,8,16,32,64];

%% colors
color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
colors = {color1, color2, color3};

figure('DefaultAxesFontSize',12); hold on;
LineW = 1.5;
MarkS = 8;
plot(N, mean_wsr(:,1), '-o', 'Color', colors{1},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);
plot(N, mean_wsr(:,2), '-o', 'Color', colors{2},'DisplayName', schemes{2}, 'LineWidth', LineW, 'MarkerSize', MarkS);
plot(N, mean_wsr(:,3), '-o', 'Color', colors{3},'DisplayName', schemes{3}, 'LineWidth', LineW, 'MarkerSize', MarkS);

clearvars mean_wsr res res_time schemes Ny numSchemes;
load('cvx_FC.mat');
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numElements × numSchemes]
%mean_time = squeeze( mean(global_time, 1, 'omitnan') );  % size = [numElements × numSchemes]
plot(N, mean_wsr, ':^', 'Color', colors{1},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);

clearvars mean_wsr res res_time schemes Ny numSchemes;
load('cvx_GC.mat');
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numElements × numSchemes]
plot(N, mean_wsr, ':^', 'Color', colors{2},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);

clearvars mean_wsr res res_time schemes Ny numSchemes;
load('cvx_SC.mat');
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numElements × numSchemes]
plot(N, mean_wsr, ':^', 'Color', colors{3},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);


%% Two-stage
clearvars mean_wsr res res_time schemes Ny numSchemes;
load('fig1_TwoStage.mat');

% Compute mean over the first (iteration) dimension
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numMi × numSchemes]
%mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numMi × numSchemes]

plot(N, mean_wsr(:,1), '-.s', 'Color', colors{1},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);
plot(N, mean_wsr(:,2), '-.s', 'Color', colors{2},'DisplayName', schemes{2}, 'LineWidth', LineW, 'MarkerSize', MarkS);
plot(N, mean_wsr(:,3), '-.s', 'Color', colors{3},'DisplayName', schemes{3}, 'LineWidth', LineW, 'MarkerSize', MarkS);

%% RZF
clearvars mean_wsr res res_time schemes Ny numSchemes;
load('fig1_other_benchmarks.mat');
% Compute mean over the first (iteration) dimension
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numMi × numSchemes]
%mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numMi × numSchemes]
plot(N, mean_wsr(:,4), '-v', 'Color', colors{1}, 'DisplayName', schemes{4}, 'LineWidth', LineW, 'MarkerSize', MarkS);
plot(N, mean_wsr(:,5), '-v', 'Color', colors{2}, 'DisplayName', schemes{5}, 'LineWidth', LineW, 'MarkerSize', MarkS);
plot(N, mean_wsr(:,6), '-v', 'Color', colors{3}, 'DisplayName', schemes{6}, 'LineWidth', LineW, 'MarkerSize', MarkS);


%% NoRIS
%plot(N, mean_wsr(:,7), '-k', 'DisplayName', schemes{7}, 'LineWidth', LineW, 'MarkerSize', MarkS);

%% RandPhase
%plot(N, mean_wsr(:,8), '-m', 'DisplayName', schemes{8}, 'LineWidth', LineW, 'MarkerSize', MarkS);

xlabel('Number of RIS elements', 'Interpreter', 'LaTeX');
ylabel('Weighted sum-rate (bps/Hz)');
legend('Location','best');
set(gca,'xtick',N);
xlim([N(1), N(end)]);
grid on;
