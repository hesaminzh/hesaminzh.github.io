%% Outage
N = 32;
K = 2;
addpath('/Users/hesami/Documents/Research/ris_isac/Hesamilab/matlab_gpt/bdris_mmWave/bdris_new/data/outage/new/')
load(['plot_Outage_FC_GC_SC__N_', num2str(N),'__K_',num2str(K),'.mat'],'WSR_all' ,'WSR_avg', 'Outage', 'eta_list', 'schemes', 'numSchemes');
%load(['plot_Outage_NoRIS__N_', num2str(N),'__K_',num2str(K),'.mat'],'WSR_all' ,'WSR_avg', 'Outage', 'eta_list', 'schemes', 'numSchemes');





figure('DefaultAxesFontSize',12); hold on;
LineW = 1.5;
MarkS = 8;
for is = 1:numel(schemes)
  plot( eta_list, WSR_avg(:,is), '-o', 'DisplayName', schemes{is}, 'LineWidth', LineW, 'MarkerSize', MarkS );
end
xlabel('$\eta$ (blockage density)', 'Interpreter', 'LaTeX');
ylabel('WSR (bps/Hz)');
legend('Location','best');
grid on
xlim([min(eta_list), max(eta_list)]);

% Step 1: Select xtick positions — keep first, then skip every other
reduced_eta = eta_list([1, 3:2:end]);

% Step 2: Create custom labels
labels = [{'0'}, arrayfun(@(x) sprintf('%.2f', x), reduced_eta(2:end), 'UniformOutput', false)];

% Step 3: Apply to plot
xticks(reduced_eta);
xticklabels(labels);

clear all; N = 32; K = 2; LineW = 1.5; MarkS = 8;
load(['plot_Outage_NoRIS__N_', num2str(N),'__K_',num2str(K),'.mat'],'WSR_all' ,'WSR_avg', 'Outage', 'eta_list', 'schemes', 'numSchemes');
for is = 1:numel(schemes)
  plot( eta_list, WSR_avg(:,is), '-o', 'DisplayName', schemes{is}, 'LineWidth', LineW, 'MarkerSize', MarkS );
end
print(gcf, 'figure7.eps', '-depsc2' );


%%% figure 8

figure('DefaultAxesFontSize',12); hold on;
LineW = 1.5;
MarkS = 8;
for is = 1:numel(schemes)
  plot( eta_list, Outage(:,is), '-o', 'DisplayName', schemes{is}, 'LineWidth', LineW, 'MarkerSize', MarkS );
end
xlabel('$\eta$ (blockage density)', 'Interpreter', 'LaTeX');
ylabel('Outage Probability');
legend('Location','best');
grid on
set(gca, 'YScale', 'log');
reduced_eta = eta_list([1, 3:2:end]);

% Step 2: Create custom labels
labels = [{'0'}, arrayfun(@(x) sprintf('%.2f', x), reduced_eta(2:end), 'UniformOutput', false)];

% Step 3: Apply to plot
xticks(reduced_eta);
xticklabels(labels);

clear all; N = 32; K = 2; LineW = 1.5; MarkS = 8;
load(['plot_Outage_NoRIS__N_', num2str(N),'__K_',num2str(K),'.mat'],'WSR_all' ,'WSR_avg', 'Outage', 'eta_list', 'schemes', 'numSchemes');
for is = 1:numel(schemes)
  plot( eta_list, Outage(:,is), '-o', 'DisplayName', schemes{is}, 'LineWidth', LineW, 'MarkerSize', MarkS );
end

print(gcf, 'figure8.eps', '-depsc2' );
