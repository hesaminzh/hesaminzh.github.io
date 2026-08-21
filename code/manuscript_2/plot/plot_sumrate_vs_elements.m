color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
color4 = [251,144,98] / 255; %(251,144,98)
color5 = [106,13,131] / 255; % (106,13,131)
colors = {color1, color2, color3, color4, color5};

LineW = 1.5;
MarkS = 8;
figure('DefaultAxesFontSize',16); %hold on;
ax1 = axes;
% ax1.Position = [0.15, 0.10, 0.75, 0.80];
hold(ax1, 'on');

%% Proposed
load('ppadmm_sumrate_w_optimal_spacing__K_2__nMonte_500.mat');
clearvars arch_list
arch_list = {'Proposed (stem-conn.)','Proposed (tree-conn.)','Proposed (single-conn.)'};

p1 = plot(ax1, M_RIS, avg_sumrate(:,1), '-o', 'Color', colors{1}, 'DisplayName', arch_list{1}, 'LineWidth', LineW, 'MarkerSize', MarkS );
p2 = plot(ax1, M_RIS, avg_sumrate(:,2), '-o', 'Color', colors{2}, 'DisplayName', arch_list{2}, 'LineWidth', LineW, 'MarkerSize', MarkS );
p3 = plot(ax1, M_RIS, avg_sumrate(:,3), '-o', 'Color', colors{3}, 'DisplayName', arch_list{3}, 'LineWidth', LineW, 'MarkerSize', MarkS );

clearvars avg_sumrate arch_list



%% pp-ADMM
load('ppadmm_sumrate__K_2__nMonte_500.mat');
clearvars arch_list
arch_list = {'pp-ADMM with fixed $d_s$ (stem-conn.) [9]', 'pp-ADMM with fixed $d_s$ (tree-conn.) [9]', 'pp-ADMM with fixed $d_s$ (single-conn.) [9]' };
p4 = plot(ax1, M_RIS, avg_sumrate(:,1), '-.*', 'Color', colors{1}, 'DisplayName', arch_list{1}, 'LineWidth', LineW, 'MarkerSize', MarkS );
% Tridiagonal
p5 = plot(ax1, M_RIS, avg_sumrate(:,2), '-.*', 'Color', colors{2}, 'DisplayName', arch_list{2}, 'LineWidth', LineW, 'MarkerSize', MarkS );
% Single-conn.
p6 = plot(ax1, M_RIS, avg_sumrate(:,3), '-.*', 'Color', colors{3}, 'DisplayName', arch_list{3}, 'LineWidth', LineW, 'MarkerSize', MarkS );


lgd = legend(ax1, 'show', 'Interpreter', 'LaTeX');
lgd.Location = 'southeast';
xlabel('Number of RIS elements per panel, $M$', 'Interpreter', 'LaTeX');
% xlabel('Number of RIS elements per panel');
ylabel('Weighted sum-rate (bps/Hz)', 'Interpreter', 'LaTeX');
% ylabel('Weighted sum-rate (bps/Hz)');
xlim([0, M_RIS(end)]);
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
% set(gca, 'xtick', M_RIS);
grid on;
box on;

% print(gcf, 'wsr_vs_elements_K_2.eps', '-depsc2');


%% K = 4
LineW = 1.5;
MarkS = 8;
figure('DefaultAxesFontSize',16); %hold on;
ax1 = axes;
% ax1.Position = [0.15, 0.10, 0.75, 0.80];
hold(ax1, 'on');

load('ppadmm_sumrate_w_optimal_spacing__K_4__nMonte_500.mat');
arch_list = {'Proposed (stem-conn.)','Proposed (tree-conn.)','Proposed (single-conn.)'};

% Stem-conn.
p1 = plot(ax1, M_RIS, avg_sumrate(:,1), '-o', 'Color', colors{1}, 'DisplayName', arch_list{1}, 'LineWidth', LineW, 'MarkerSize', MarkS );
% Tridiagonal
p2 = plot(ax1, M_RIS, avg_sumrate(:,2), '-o', 'Color', colors{2}, 'DisplayName', arch_list{2}, 'LineWidth', LineW, 'MarkerSize', MarkS );
% Single-conn.
p3 = plot(ax1, M_RIS, avg_sumrate(:,3), '-o', 'Color', colors{3}, 'DisplayName', arch_list{3}, 'LineWidth', LineW, 'MarkerSize', MarkS );

%% pp-ADMM
% Stem-conn.
load('ppadmm_sumrate__K_4__nMonte_500.mat');
arch_list = {'pp-ADMM with fixed $d_s$ (stem-conn.) [9]', 'pp-ADMM with fixed $d_s$ (tree-conn.) [9]', 'pp-ADMM with fixed $d_s$ (single-conn.) [9]' };

p4 = plot(ax1, M_RIS, avg_sumrate(:,1), '-.*', 'Color', colors{1}, 'DisplayName', arch_list{1}, 'LineWidth', LineW, 'MarkerSize', MarkS );
% Tridiagonal
p5 = plot(ax1, M_RIS, avg_sumrate(:,2), '-.*', 'Color', colors{2}, 'DisplayName', arch_list{2}, 'LineWidth', LineW, 'MarkerSize', MarkS );
% Single-conn.
p6 = plot(ax1, M_RIS, avg_sumrate(:,3), '-.*', 'Color', colors{3}, 'DisplayName', arch_list{3}, 'LineWidth', LineW, 'MarkerSize', MarkS );


lgd = legend(ax1, 'show', 'Interpreter', 'LaTeX');
lgd.Location = 'southeast';
xlabel('Number of RIS elements per panel, $M$', 'Interpreter', 'LaTeX');
% xlabel('Number of RIS elements per panel');
ylabel('Weighted sum-rate (bps/Hz)', 'Interpreter', 'LaTeX');
% ylabel('Weighted sum-rate (bps/Hz)');
xlim([0, M_RIS(end)]);
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
% set(gca, 'xtick', M_RIS);
grid on;
box on;

% print(gcf, 'wsr_vs_elements_K_4.eps', '-depsc2');
