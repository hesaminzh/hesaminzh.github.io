addpath('/Users/hesami/Documents/Research/ris_isac/Hesamilab/matlab_gpt/bdris_mmWave/bdris_new/data/outage/new/');
%% WSR vs Blockage Density
load('plot_Outage_FC_GC_SC__N_64__K_2.mat');

color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
color4 = [251,144,98] / 255; %(251,144,98) % RandPhase
color5 = [106,13,131] / 255; % (106,13,131) % NoRIS
colors = {color1, color2, color3, color4, color5};


figure('DefaultAxesFontSize',16); %hold on;
ax1 = axes;
hold(ax1, 'on');
LineW = 1.5;
MarkS = 8;
p1 = plot(ax1, eta_list, WSR_avg(:,1), '-o', 'Color', colors{1}, 'DisplayName', schemes{1},'LineWidth', LineW, 'MarkerSize', MarkS );
p2 = plot(ax1, eta_list, WSR_avg(:,2), '-s', 'Color', colors{2}, 'DisplayName', schemes{2},'LineWidth', LineW, 'MarkerSize', MarkS );
p3 = plot(ax1, eta_list, WSR_avg(:,3), '-*', 'Color', colors{3}, 'DisplayName', schemes{3},'LineWidth', LineW, 'MarkerSize', MarkS );

% NoRIS and RandPhase
clearvars WSR_all WSR_avg Outage
load('plot_Outage_NoRIS_RandPhase__N_64__K_2.mat')
p4 = plot(ax1, eta_list, WSR_avg(:,1), '-x', 'Color', colors{5}, 'DisplayName', schemes{1},'LineWidth', LineW, 'MarkerSize', MarkS );
p5 = plot(ax1, eta_list, WSR_avg(:,2), '-^', 'Color', colors{4}, 'DisplayName', schemes{2},'LineWidth', LineW, 'MarkerSize', MarkS );




%% second axis for color legneds
ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'off');
hold(ax2, 'on');
% dummy plots
c1 = plot(ax2, NaN, NaN, '-o', 'Color', colors{1}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'FC');
c2 = plot(ax2, NaN, NaN, '-s', 'Color', colors{2}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'GC');
c3 = plot(ax2, NaN, NaN, '-*', 'Color', colors{3}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'SC');
c4 = plot(ax2, NaN, NaN, '-^', 'Color', colors{4}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'Random Phase');
c5 = plot(ax2, NaN, NaN, '-x', 'Color', colors{5}, 'LineWidth', LineW, 'MarkerSize', MarkS,'DisplayName', 'No RIS');

%lg2 = legend(ax2, [c1, c2, c3], {'FC', 'GC', 'SC'}, 'Location', 'NorthEast');
lg2 = legend(ax2, [c1, c2, c3, c4, c5], {'FC', 'GC', 'SC', 'Random Phase', 'No RIS'});
set(lg2, 'Color', 'white');

%%% third axis linestyle legends
%ax3 = axes('Position', get(ax1, 'Position'), 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'off');
%hold(ax3, 'on');
%%h1 = plot(ax2, NaN, NaN, 'k-', 'LineWidth', 1.5, 'DisplayName', '-');
%s1 = plot(ax3, NaN, NaN, 'k-o', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'Proposed' );
%s2 = plot(ax3, NaN, NaN, 'k:^', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'WMMSE-CVX MO' );
%s3 = plot(ax3, NaN, NaN, 'k-.s', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'two-stage FP');
%s4 = plot(ax3, NaN, NaN, 'k-v', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'RZF MO');
%
%lg3 = legend(ax3, [s1, s2, s3, s4], {'Proposed', 'WMMSE-CVX MO', 'two-stage FP', 'RZF MO'}, 'Location', 'NorthWest');
%set(lg3, 'Color', 'white');
%
%pos3 = lg3.Position;
%pos2 = lg2.Position;
%%pos2 = pos3;
%pos2(1) = pos3(1) + pos3(4) + 0.05;
%%pos2(3) = pos(3) * 0.6;
%set(lg2, 'Position', pos2);
%% set the position of the legends
%%pos2 = lg2.Position;
%%pos3 = pos2;
%%pos3(2) = pos2(2) - pos2(4) - 0.01; % shift down
%%set(lg3, 'Position', pos3 )
%
%% Link the axes to synchronize zooming and panning
%linkaxes([ax1, ax2, ax3]);
linkaxes([ax1, ax2]);

% Ensure the second axes does not interfere with the first
set(ax2, 'HitTest', 'off', 'HandleVisibility', 'off');
%set(ax3, 'HitTest', 'off', 'HandleVisibility', 'off');
%xlabel('\textrm{Number of RIS elements} $N$','FontSize', 12, 'Interpreter', 'LaTeX');
%ylabel('\textrm{Weighted sum-rate (bps/Hz)}', 'FontSize', 12, 'Interpreter', 'LaTeX');
xlabel('\textrm{Blockage Density} $\eta$', 'Interpreter', 'LaTeX');
ylabel('\textrm{Weighted sum-rate (bps/Hz)}', 'Interpreter', 'LaTeX');
%set(gca,'xtick',N_new);
xlim([eta_list(1), eta_list(end)]);
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
grid on;
box on;


%% save as .eps format
% print(gcf, 'wsr_vs_blockage.eps', '-depsc2' );




