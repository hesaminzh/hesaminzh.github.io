addpath('/Users/hesami/Documents/Research/ris_isac/Hesamilab/matlab_gpt/bdris_mmWave/bdris_new/data/figure_1/nMonte_1000/');
load('fig1_proposed_closed_form_K_2.mat');

%% colors
color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
colors = {color1, color2, color3};


%% Proposed
% Compute mean over the first (iteration) dimension
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numMi × numSchemes]
%mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numMi × numSchemes]

%N_orig = [0,4,8,16,32,64];
N_orig = [4,8,16,32,64];
%firstval = 1.906788539005997;
y_orig1 = mean_wsr(:, 1);        % 5×1 vector for scheme sch
%y_orig1 = [firstval;y_orig1];
y_orig2 = mean_wsr(:, 2);        % 5×1 vector for scheme sch
%y_orig2 = [firstval;y_orig2];
y_orig3 = mean_wsr(:, 3);        % 5×1 vector for scheme sch
%y_orig3 = [firstval;y_orig3];

% new x‐axis with “in‐between” points
N_new  = [4, 8, 16, 24, 32, 40, 48, 56, 64];

% interpolate (choose 'pchip' or 'spline'), here we chose 'spline'
y_new1 = interp1(N_orig, y_orig1, N_new, 'spline');
y_new2 = interp1(N_orig, y_orig2, N_new, 'spline');
y_new3 = interp1(N_orig, y_orig3, N_new, 'spline');

figure('DefaultAxesFontSize',12); %hold on;
ax1 = axes;
hold(ax1, 'on');
LineW = 1.5;
MarkS = 8;
p1 = plot(ax1, N_new, y_new1, '-o', 'Color', colors{1},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);
p2 = plot(ax1, N_new, y_new2, '-o', 'Color', colors{2},'DisplayName', schemes{2}, 'LineWidth', LineW, 'MarkerSize', MarkS);
p3 = plot(ax1, N_new, y_new3, '-o', 'Color', colors{3},'DisplayName', schemes{3}, 'LineWidth', LineW, 'MarkerSize', MarkS);




%% CVX
clearvars mean_wsr res res_time schemes Ny numSchemes y_new1 y_new2 y_new3;
load('cvx_FC.mat');
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numElements × numSchemes]
y_orig1 = mean_wsr;
y_new1 = interp1(N_orig, y_orig1, N_new, 'spline');
%mean_time = squeeze( mean(global_time, 1, 'omitnan') );  % size = [numElements × numSchemes]
p4 = plot(ax1, N_new, y_new1, ':^', 'Color', colors{1},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);

clearvars mean_wsr res res_time schemes Ny numSchemes y_new1 y_new2 y_new3;
load('cvx_GC.mat');
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numElements × numSchemes]
y_orig2 = mean_wsr;
y_new2 = interp1(N_orig, y_orig2, N_new, 'spline');
p5 = plot(ax1, N_new, y_new2, ':^', 'Color', colors{2},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);

clearvars mean_wsr res res_time schemes Ny numSchemes y_new1 y_new2 y_new3;
load('cvx_SC.mat');
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numElements × numSchemes]
y_orig3 = mean_wsr;
y_new3 = interp1(N_orig, y_orig3, N_new, 'spline');
p6 = plot(ax1, N_new, y_new3, ':^', 'Color', colors{3},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);


%% Two-stage
clearvars mean_wsr res res_time schemes Ny numSchemes y_new1 y_new2 y_new3;
load('fig1_TwoStage.mat');

% Compute mean over the first (iteration) dimension
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numMi × numSchemes]
%mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numMi × numSchemes]

y_orig1 = mean_wsr(:,1);
y_orig2 = mean_wsr(:,2);
y_orig3 = mean_wsr(:,3);

y_new1 = interp1(N_orig, y_orig1, N_new, 'spline');
y_new2 = interp1(N_orig, y_orig2, N_new, 'spline');
y_new3 = interp1(N_orig, y_orig3, N_new, 'spline');

p7 = plot(ax1, N_new, y_new1, '-.s', 'Color', colors{1},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);
p8 = plot(ax1, N_new, y_new2, '-.s', 'Color', colors{2},'DisplayName', schemes{2}, 'LineWidth', LineW, 'MarkerSize', MarkS);
p9 = plot(ax1, N_new, y_new3, '-.s', 'Color', colors{3},'DisplayName', schemes{3}, 'LineWidth', LineW, 'MarkerSize', MarkS);

%% RZF
clearvars mean_wsr res res_time schemes Ny numSchemes y_new1 y_new2 y_new3;
load('fig1_other_benchmarks.mat');
% Compute mean over the first (iteration) dimension
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numMi × numSchemes]
%mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numMi × numSchemes]

y_orig1 = mean_wsr(:,4);
y_orig2 = mean_wsr(:,5);
y_orig3 = mean_wsr(:,6);

y_new1 = interp1(N_orig, y_orig1, N_new, 'spline');
y_new2 = interp1(N_orig, y_orig2, N_new, 'spline');
y_new3 = interp1(N_orig, y_orig3, N_new, 'spline');

p10 = plot(ax1, N_new, y_new1, '-v', 'Color', colors{1}, 'DisplayName', schemes{4}, 'LineWidth', LineW, 'MarkerSize', MarkS);
p11 = plot(ax1, N_new, y_new2, '-v', 'Color', colors{2}, 'DisplayName', schemes{5}, 'LineWidth', LineW, 'MarkerSize', MarkS);
p12 = plot(ax1, N_new, y_new3, '-v', 'Color', colors{3}, 'DisplayName', schemes{6}, 'LineWidth', LineW, 'MarkerSize', MarkS);


%% NoRIS
%plot(N_new, mean_wsr(:,7), '-k', 'DisplayName', schemes{7}, 'LineWidth', LineW, 'MarkerSize', MarkS);

%% RandPhase
%plot(N_new, mean_wsr(:,8), '-m', 'DisplayName', schemes{8}, 'LineWidth', LineW, 'MarkerSize', MarkS);




%% second axis for color legneds
ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'off');
hold(ax2, 'on');
% dummy plots
c1 = plot(ax2, NaN, NaN, '-', 'Color', colors{1}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'FC');
c2 = plot(ax2, NaN, NaN, '-', 'Color', colors{2}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'GC');
c3 = plot(ax2, NaN, NaN, '-', 'Color', colors{3}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'SC');

lg2 = legend(ax2, [c1, c2, c3], {'FC', 'GC', 'SC'}, 'Location', 'best');
set(lg2, 'Color', 'white');

%% third axis linestyle legends
ax3 = axes('Position', get(ax1, 'Position'), 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'off');
hold(ax3, 'on');
%h1 = plot(ax2, NaN, NaN, 'k-', 'LineWidth', 1.5, 'DisplayName', '-');
s1 = plot(ax3, NaN, NaN, 'k-o', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'Proposed' );
s2 = plot(ax3, NaN, NaN, 'k:^', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'WMMSE-CVX MO' );
s3 = plot(ax3, NaN, NaN, 'k-.s', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'two-stage');
s4 = plot(ax3, NaN, NaN, '-v', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'RZF');

lg3 = legend(ax3, [s1, s2, s3, s4], {'Proposed', 'WMMSE-CVX MO', 'two-stage', 'RZF'}, 'Location', 'best');
set(lg3, 'Color', 'white');

% Link the axes to synchronize zooming and panning
linkaxes([ax1, ax2, ax3]);

% Ensure the second axes does not interfere with the first
set(ax2, 'HitTest', 'off', 'HandleVisibility', 'off');
set(ax3, 'HitTest', 'off', 'HandleVisibility', 'off');
xlabel('Number of RIS elements $N$', 'Interpreter', 'LaTeX');
ylabel('Weighted sum-rate (bps/Hz)');
legend('Location','best');
set(gca,'xtick',N_new);
xlim([N_new(1), N_new(end)]);
grid on;
