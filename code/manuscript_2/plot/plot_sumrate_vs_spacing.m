color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
color4 = [251,144,98] / 255; %(251,144,98)
color5 = [106,13,131] / 255; % (106,13,131)
colors = {color1, color2, color3, color4, color5};

%% M = 16
load('Spacings__K_4__M_16__nMonte_500.mat');
clearvars arch_list
arch_list = {'Stem-conn.', 'Tree-conn.', 'Single-conn.'};

figure('DefaultAxesFontSize',20); %hold on;
ax1 = axes;
hold(ax1, 'on');
LineW = 1.5;
MarkS = 8;
p1 = plot(ax1, panel_spacing_list, avg_sumrate_spacing(:,1), '-o', 'Color', colors{1}, 'DisplayName', arch_list{1}, 'LineWidth', LineW, 'MarkerSize', MarkS );
p2 = plot(ax1, panel_spacing_list, avg_sumrate_spacing(:,2), '-s', 'Color', colors{2}, 'DisplayName', arch_list{2}, 'LineWidth', LineW, 'MarkerSize', MarkS );
p3 = plot(ax1, panel_spacing_list, avg_sumrate_spacing(:,3), '-*', 'Color', colors{3}, 'DisplayName', arch_list{3}, 'LineWidth', LineW, 'MarkerSize', MarkS );

lgd = legend(ax1, 'show');
lgd.Location = 'northeast';

xlabel('Inter-panel spacing $d_s$ (m)', 'Interpreter', 'LaTeX');
ylabel('Weighted sum-rate (bps/Hz)', 'Interpreter', 'LaTeX');
xlim([panel_spacing_list(1), panel_spacing_list(end)]);
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
% set(gca, 'xtick', panel_spacing_list);
grid on;
box on;

% save as .eps format
% print(gcf, 'wsr_vs_spacing_M_16.eps', '-depsc2' );

%% M = 32
clearvars avg_sumrate_spacing
load('Spacings__K_4__M_32__nMonte_500.mat');
arch_list = {'Stem-conn.', 'Tree-conn.', 'Single-conn.'};

figure('DefaultAxesFontSize',20); %hold on;
ax1 = axes;
hold(ax1, 'on');
LineW = 1.5;
MarkS = 8;
p1 = plot(ax1, panel_spacing_list, avg_sumrate_spacing(:,1), '-o', 'Color', colors{1}, 'DisplayName', arch_list{1}, 'LineWidth', LineW, 'MarkerSize', MarkS );
p2 = plot(ax1, panel_spacing_list, avg_sumrate_spacing(:,2), '-s', 'Color', colors{2}, 'DisplayName', arch_list{2}, 'LineWidth', LineW, 'MarkerSize', MarkS );
p3 = plot(ax1, panel_spacing_list, avg_sumrate_spacing(:,3), '-*', 'Color', colors{3}, 'DisplayName', arch_list{3}, 'LineWidth', LineW, 'MarkerSize', MarkS );

lgd = legend(ax1, 'show');
lgd.Location = 'southeast';

xlabel('Inter-panel spacing $d_s$ (m)', 'Interpreter', 'LaTeX');
ylabel('Weighted sum-rate (bps/Hz)', 'Interpreter', 'LaTeX');
xlim([panel_spacing_list(1), panel_spacing_list(end)]);
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
% set(gca, 'xtick', panel_spacing_list);
grid on;
box on;
%print(gcf, 'wsr_vs_spacing_M_32.eps', '-depsc2' );

%% M = 64
clearvars avg_sumrate_spacing
load('Spacings__K_4__M_64__nMonte_500.mat');
arch_list = {'Stem-conn.', 'Tree-conn.', 'Single-conn.'};

figure('DefaultAxesFontSize',20); %hold on;
ax1 = axes;
hold(ax1, 'on');
LineW = 1.5;
MarkS = 8;
p1 = plot(ax1, panel_spacing_list, avg_sumrate_spacing(:,1), '-o', 'Color', colors{1}, 'DisplayName', arch_list{1}, 'LineWidth', LineW, 'MarkerSize', MarkS );
p2 = plot(ax1, panel_spacing_list, avg_sumrate_spacing(:,2), '-s', 'Color', colors{2}, 'DisplayName', arch_list{2}, 'LineWidth', LineW, 'MarkerSize', MarkS );
p3 = plot(ax1, panel_spacing_list, avg_sumrate_spacing(:,3), '-*', 'Color', colors{3}, 'DisplayName', arch_list{3}, 'LineWidth', LineW, 'MarkerSize', MarkS );

lgd = legend(ax1, 'show');
lgd.Location = 'southeast';

xlabel('Inter-panel spacing $d_s$ (m)', 'Interpreter', 'LaTeX');
ylabel('Weighted sum-rate (bps/Hz)', 'Interpreter', 'LaTeX');
xlim([panel_spacing_list(1), panel_spacing_list(end)]);
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
% set(gca, 'xtick', panel_spacing_list);
grid on;
box on;
%print(gcf, 'wsr_vs_spacing_M_64.eps', '-depsc2' );
