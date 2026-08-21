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
y_orig1(end) = 8.264284078551960;
%8.264284078551960
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

% figure('DefaultAxesFontSize',12); %hold on;
figure('DefaultAxesFontSize',16); %hold on;
ax1 = axes;
ax1.Position = [0.15, 0.10, 0.75, 0.80];
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

%lg2 = legend(ax2, [c1, c2, c3], {'FC', 'GC', 'SC'}, 'Location', 'NorthEast');
lg2 = legend(ax2, [c1, c2, c3], {'FC', 'GC', 'SC'});
set(lg2, 'Color', 'white');

%% third axis linestyle legends
ax3 = axes('Position', get(ax1, 'Position'), 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'off');
hold(ax3, 'on');
%h1 = plot(ax2, NaN, NaN, 'k-', 'LineWidth', 1.5, 'DisplayName', '-');
s1 = plot(ax3, NaN, NaN, 'k-o', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'Proposed' );
s2 = plot(ax3, NaN, NaN, 'k:^', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'WMMSE-CVX MO' );
s3 = plot(ax3, NaN, NaN, 'k-.s', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'two-stage FP');
s4 = plot(ax3, NaN, NaN, 'k-v', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'RZF MO');

lg3 = legend(ax3, [s1, s2, s3, s4], {'Proposed', 'WMMSE-CVX MO', 'two-stage FP', 'RZF MO'}, 'Location', 'NorthWest');
set(lg3, 'Color', 'white');

pos3 = lg3.Position;
pos2 = lg2.Position;
%pos2 = pos3;
pos2(1) = pos3(1) + pos3(4) + 0.115;
%pos2(3) = pos(3) * 0.6;
set(lg2, 'Position', pos2);
% set the position of the legends
%pos2 = lg2.Position;
%pos3 = pos2;
%pos3(2) = pos2(2) - pos2(4) - 0.01; % shift down
%set(lg3, 'Position', pos3 )

% Link the axes to synchronize zooming and panning
linkaxes([ax1, ax2, ax3]);

% Ensure the second axes does not interfere with the first
set(ax2, 'HitTest', 'off', 'HandleVisibility', 'off');
set(ax3, 'HitTest', 'off', 'HandleVisibility', 'off');
%xlabel('\textrm{Number of RIS elements} $N$','FontSize', 12, 'Interpreter', 'LaTeX');
%ylabel('\textrm{Weighted sum-rate (bps/Hz)}', 'FontSize', 12, 'Interpreter', 'LaTeX');
xlabel('Number of RIS elements $N$', 'Interpreter', 'LaTeX');
ylabel('Weighted sum-rate (bps/Hz)', 'Interpreter', 'LaTeX');
set(gca,'xtick',N_new);
xlim([N_new(2), N_new(end)]);
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
grid on;
box on;


%% save as .eps format
%print(gcf, 'wsr_vs_elements.eps', '-depsc2' );
