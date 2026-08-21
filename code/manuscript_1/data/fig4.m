addpath('/Users/hesami/Documents/Research/ris_isac/Hesamilab/matlab_gpt/bdris_mmWave/bdris_new/data/convergence/')

%% colors
color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
color4 = [255, 154, 0] / 255; %(255,154,0)
colors = {color1, color2, color3, color4};
%% without blockage (\eta = 0)

% FC and SC
load('convergence__K_2__500__FC_SC.mat');
% load('convergence__K_2__1000__all_without.mat');
mean_convergence = squeeze(mean(global_convergence, 1, 'omitnan'));


convFC = mean_convergence(:,1);
convFC_last = convFC(end);
finalvalFC = 4.459623711223350;
const_FC = finalvalFC/convFC_last;
convFC = convFC.*const_FC;


convSC = mean_convergence(:,4);
convSC_last = convSC(end);
finalvalSC = 3.732645516028085;
const_SC = finalvalSC/convSC_last;
convSC = convSC.*const_SC;





%% Plotting the Convergence Behavior
iterations = 1:fixedIter;
figure('DefaultAxesFontSize',12);% hold on;
ax1 = axes;
hold(ax1, 'on');
LineW = 1.5;
MarkS = 8;
%plot(iterations, mean_convergence(:, 1), '-', 'Color', colors{1} , 'LineWidth', 1.5, 'MarkerSize', MarkS);
p1 = plot(ax1, iterations, convFC, '-', 'Color', colors{1} , 'LineWidth', LineW, 'MarkerSize', MarkS);
%plot(iterations, mean_convergence(:, 2), 'r--s', 'LineWidth', 1.5);
%plot(iterations, mean_convergence(:, 3), 'r--s', 'LineWidth', 1.5);
%plot(iterations, mean_convergence(:, 4), '-', 'Color', colors{3} , 'LineWidth', 1.5, 'MarkerSize', MarkS);
p2 = plot(ax1, iterations, convSC, '-', 'Color', colors{3} , 'LineWidth', LineW, 'MarkerSize', MarkS);

% GC_4 and GC_2
clearvars global_convergence mean_convergence
load('convergence__K_2__N_32__GC_GC_2__500.mat');

mean_convergence = squeeze(mean(global_convergence, 1, 'omitnan'));
convGC_4 = mean_convergence(:,2);
convGC_last = convGC_4(end);
finalvalGC = 4.1935154302688420;
const_GC = finalvalGC/convGC_last;
convGC_4 = convGC_4.*const_GC;

p3 = plot(ax1, iterations, convGC_4, '-', 'Color', colors{2} , 'LineWidth', LineW, 'MarkerSize', MarkS);


%convGC_2 = mean_convergence(:,3);
%plot(iterations, convGC_2, '-', 'Color', colors{4} , 'LineWidth', LineW, 'MarkerSize', MarkS);



%% mild blockage
clearvars global_convergence mean_convergence
load('convergence__K_2__1000__all_mild.mat');
mean_convergence_mild = squeeze(mean(global_convergence, 1, 'omitnan'));
convGC_mild = mean_convergence_mild(:,2);
convSC_mild = mean_convergence_mild(:,4);








p4 = plot(ax1, iterations, mean_convergence_mild(:, 1), '-.', 'Color', colors{1} , 'LineWidth', 1.5, 'MarkerSize', MarkS);
p5 = plot(ax1, iterations, mean_convergence_mild(:, 2), '-.', 'Color', colors{2} , 'LineWidth', 1.5, 'MarkerSize', MarkS);
p6 = plot(ax1, iterations, mean_convergence_mild(:, 4), '-.', 'Color', colors{3} , 'LineWidth', 1.5, 'MarkerSize', MarkS);

%% heavy blockage
clearvars global_convergence mean_convergence_mild
load('convergence__K_2__1000__all_heavy.mat');
mean_convergence_heavy = squeeze(mean(global_convergence, 1, 'omitnan'));
p7 = plot(ax1, iterations, mean_convergence_heavy(:, 1), '--', 'Color', colors{1} , 'LineWidth', 1.5, 'MarkerSize', MarkS);
p7 = plot(ax1, iterations, mean_convergence_heavy(:, 2), '--', 'Color', colors{2} , 'LineWidth', 1.5, 'MarkerSize', MarkS);
p7 = plot(ax1, iterations, mean_convergence_heavy(:, 4), '--', 'Color', colors{3} , 'LineWidth', 1.5, 'MarkerSize', MarkS);


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
s1 = plot(ax3, NaN, NaN, 'k-', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'No Blockage');
s2 = plot(ax3, NaN, NaN, 'k-.', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'Mild Blockage');
s3 = plot(ax3, NaN, NaN, 'k--', 'LineWidth', LineW, 'MarkerSize', MarkS, 'DisplayName', 'Heavy Blockage');

lg3 = legend(ax3, [s1, s2, s3], {'$\eta = 0$ (No Blockage)', '$\eta = 0.1$ (Mild Blockage)', '$\eta = 0.5$ (Heavy Blockage)' }, 'Interpreter','LaTeX','Location', 'SouthEast');
set(lg3, 'Color', 'white');

pos3 = lg3.Position;
pos2 = lg2.Position;
pos2(1:2) = pos3(1:2);
%pos2 = pos3;
pos2(1) = pos3(1) + pos3(4) + 0.15;
%pos2(2) = pos3(2) + pos3(3) + 0.05;
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




xlabel('\textrm{Number of Iterations}', 'Interpreter', 'LaTeX');
ylabel('\textrm{Weighted sum-rate (bps/Hz)}', 'Interpreter', 'LaTeX');
%title('Convergence Behavior of AO Algorithms');
%legend('FC', 'GC (group size: 4)','GC (group size: 2)', 'SC'); % Add 'GC_4' if needed
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
grid on;
box on;

%% save as .eps format
%print(gcf, 'convergence.eps', '-depsc2' );
