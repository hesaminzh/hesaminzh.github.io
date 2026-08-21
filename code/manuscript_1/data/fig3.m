addpath('/Users/hesami/Documents/Research/ris_isac/Hesamilab/matlab_gpt/bdris_mmWave/bdris_new/data/power/')

%% colors
color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
colors = {color1, color2, color3};


%% K = 6
load('power__32__K_4_6.mat');
figure('DefaultAxesFontSize',12); %hold on;
ax1 = axes;
hold(ax1, 'on');
LineW = 1.5;
MarkS = 8;
plot(ax1, P_dBm, WSR_avg(2,:,1), '-o','Color',colors{1},'LineWidth',LineW,'MarkerSize',MarkS,'DisplayName', schemes{1} );
plot(ax1, P_dBm, WSR_avg(2,:,2), '-o','Color',colors{2},'LineWidth',LineW,'MarkerSize',MarkS,'DisplayName', schemes{2} );
plot(ax1, P_dBm, WSR_avg(2,:,3), '-o','Color',colors{3},'LineWidth',LineW,'MarkerSize',MarkS,'DisplayName', schemes{3} );




%% K = 4

plot(ax1, P_dBm, WSR_avg(1,:,1), '-x','Color',colors{1},'LineWidth',LineW,'MarkerSize',MarkS,'DisplayName', schemes{1} );
plot(ax1, P_dBm, WSR_avg(1,:,2), '-x','Color',colors{2},'LineWidth',LineW,'MarkerSize',MarkS,'DisplayName', schemes{1} );
plot(ax1, P_dBm, WSR_avg(1,:,3), '-x','Color',colors{3},'LineWidth',LineW,'MarkerSize',MarkS,'DisplayName', schemes{1} );


%% K = 2
clearvars WSR_all WSR_avg
P_dBm = [20,25,30,35,40,45];
load('power__32__K_2.mat');
plot(ax1, P_dBm, WSR_avg(:,1), '-s','Color',colors{1},'LineWidth',LineW,'MarkerSize',MarkS,'DisplayName', schemes{1} );
plot(ax1, P_dBm, WSR_avg(:,2), '-s','Color',colors{2},'LineWidth',LineW,'MarkerSize',MarkS,'DisplayName', schemes{2} );
plot(ax1, P_dBm, WSR_avg(:,3), '-s','Color',colors{3},'LineWidth',LineW,'MarkerSize',MarkS,'DisplayName', schemes{3} );




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
s1 = plot(ax3, NaN, NaN, 'k-o', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'K = 6' );
s2 = plot(ax3, NaN, NaN, 'k-x', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'K = 4' );
s3 = plot(ax3, NaN, NaN, 'k-s', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'K = 2');

lg3 = legend(ax3, [s1, s2, s3], {'$K = 6$', '$K = 4$', '$K = 2$'}, 'Location', 'NorthWest', 'Interpreter', 'LaTeX');
set(lg3, 'Color', 'white');

pos3 = lg3.Position;
pos2 = lg2.Position;
%pos2 = pos3;
pos2(1) = pos3(1) + pos3(4) + 0.005;
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
xlabel('\textrm{Transmit power (dBm)}', 'Interpreter', 'LaTeX');
ylabel('\textrm{Weighted sum-rate (bps/Hz)}', 'Interpreter', 'LaTeX');
set(gca,'xtick',P_dBm);
xlim([P_dBm(1), P_dBm(end)]);
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
grid on;

%% save as .eps format
print(gcf, 'figure3.eps', '-depsc2' );
