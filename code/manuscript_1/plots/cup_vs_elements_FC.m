addpath('/Users/hesami/Documents/Research/ris_isac/Hesamilab/matlab_gpt/bdris_mmWave/bdris_new/data/figure_1/nMonte_1000/');
load('fig1_proposed_closed_form_K_2.mat');
%% colors
color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
colors = {color1, color2, color3};

%% Proposed
% Compute mean over the first (iteration) dimension
%mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numMi × numSchemes]
mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numMi × numSchemes]

N_orig = [4,8,16,32,64];
y_orig1 = mean_time(:, 1);        % 5×1 vector for scheme sch
y_orig2 = mean_time(:, 2);        % 5×1 vector for scheme sch
y_orig3 = mean_time(:, 3);        % 5×1 vector for scheme sch

% new x‐axis with “in‐between” points
N_new  = [4, 8, 16, 24, 32, 40, 48, 56, 64];

y_new1 = interp1(N_orig, y_orig1, N_new, 'spline');
y_new2 = interp1(N_orig, y_orig2, N_new, 'spline');
y_new3 = interp1(N_orig, y_orig3, N_new, 'spline');

figure('DefaultAxesFontSize',16); %hold on;
ax1 = axes;
ax1.Position = [0.15, 0.10, 0.75, 0.80];
hold(ax1, 'on');
LineW = 1.5;
MarkS = 8;

p1 = loglog(ax1, N_new, y_new1, '-o', 'Color', colors{1},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);





%% CVX
clearvars mean_time res res_time schemes Ny numSchemes;
load('cvx_FC.mat');
%mean_wsr = squeeze( mean(global_wsr, 1, 'omitnan') );  % size = [numElements × numSchemes]
mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numElements × numSchemes]
y_orig1 = mean_time;
y_new1 = interp1(N_orig, y_orig1, N_new, 'spline');

p2 = loglog(ax1, N_new, y_new1, ':^', 'Color', colors{1},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);




%% two-stage
clearvars mean_time res res_time schemes Ny numSchemes;
load('fig1_TwoStage.mat');

% Compute mean over the first (iteration) dimension
%mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numMi × numSchemes]
mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numMi × numSchemes]
twostageFC = mean_time(:,1);
y_orig1 = twostageFC;

y_new1 = interp1(N_orig, y_orig1, N_new, 'linear');

p3 = loglog(ax1, N_new, y_new1, '-.s', 'Color', colors{1},'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS);





%% RZF
clearvars mean_time res res_time schemes Ny numSchemes;
load('fig1_other_benchmarks.mat');
% Compute mean over the first (iteration) dimension
%mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numMi × numSchemes]
mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numMi × numSchemes]
y_orig1 = mean_time(:,4);
% y_orig2 = mean_time(:,5);
% y_orig3 = mean_time(:,6);

y_new1 = interp1(N_orig, y_orig1, N_new, 'spline');
% y_new2 = interp1(N_orig, y_orig2, N_new, 'spline');
% y_new3 = interp1(N_orig, y_orig3, N_new, 'spline');

p4 = loglog(ax1, N_new, y_new1, '-v', 'Color', colors{1}, 'DisplayName', schemes{4}, 'LineWidth', LineW, 'MarkerSize', MarkS);


%% second axis for color legneds
ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'off');
hold(ax2, 'on');
% dummy plots
c1 = loglog(ax2, NaN, NaN, '-o', 'Color', colors{1}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'Proposed');
c2 = loglog(ax2, NaN, NaN, ':^', 'Color', colors{1}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'WMMSE-CVX MO');
c3 = loglog(ax2, NaN, NaN, '-.s', 'Color', colors{1}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'two-stage FP');
c4 = loglog(ax2, NaN, NaN, '-v', 'Color', colors{1}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'RZF MO');

%lg2 = legend(ax2, [c1, c2, c3], {'FC', 'GC', 'SC'}, 'Location', 'NorthEast');
lg2 = legend(ax2, [c1, c2, c3, c4], {'Porposed (FC)', 'WMMSE-CVX MO (FC)', 'two-stage FP (FC)', 'RZF MO (FC)'});
set(lg2, 'Color', 'white');
legend_pos = lg2.Position;
legend_pos(2) = legend_pos(2)-0.7*legend_pos(2);
legend_pos(1) = legend_pos(1)+0.05*legend_pos(1);
lg2.Position = legend_pos;


linkaxes([ax1, ax2]);

% Ensure the second axes does not interfere with the first
set(ax2, 'HitTest', 'off', 'HandleVisibility', 'off');

xlabel('Number of RIS elements $N$', 'Interpreter', 'LaTeX');
ylabel('Average CPU time (sec)', 'Interpreter', 'LaTeX');
set(gca,'xtick',N_new);
set(gca, 'YScale', 'log');
xlim([N_new(2), N_new(end)]);
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
grid on;
box on;

%% save as .eps format
% print(gcf, 'cpu_vs_N_FC.eps', '-depsc2' );
