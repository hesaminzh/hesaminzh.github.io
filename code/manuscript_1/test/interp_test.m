addpath('/Users/hesami/Documents/Research/ris_isac/Hesamilab/matlab_gpt/bdris_mmWave/bdris_new/data/figure_1/nMonte_1000/');
load('fig1_proposed_closed_form_K_2.mat');
%% Proposed
% Compute mean over the first (iteration) dimension
mean_wsr = squeeze( mean(res, 1, 'omitnan') );  % size = [numMi × numSchemes]
mean_time = squeeze( mean(res_time, 1, 'omitnan') );  % size = [numMi × numSchemes]

Nx = 4;                      % x-axis surface element nums
Ny = [1,2,4,8,16];   % y-axis surface element nums
%Ny = [1,8];   % y-axis surface element nums
N  = Nx*Ny;

color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
colors = {color1, color2, color3};

figure('DefaultAxesFontSize',12); hold on;
LineW = 1.5;
MarkS = 8;
plot(N, mean_wsr(:,1), '-o', 'Color', color1,'DisplayName', schemes{1}, 'LineWidth', LineW, 'MarkerSize', MarkS); hold on;
plot(N, mean_wsr(:,2), '-o', 'Color', color2,'DisplayName', schemes{2}, 'LineWidth', LineW, 'MarkerSize', MarkS);
plot(N, mean_wsr(:,3), '-o', 'Color', color3,'DisplayName', schemes{3}, 'LineWidth', LineW, 'MarkerSize', MarkS);
%
% original x‐ and y‐values
N_orig = [4, 8, 16, 32, 64];
y_orig1 = mean_wsr(:, 1);        % 5×1 vector for scheme sch
y_orig2 = mean_wsr(:, 2);        % 5×1 vector for scheme sch
y_orig3 = mean_wsr(:, 3);        % 5×1 vector for scheme sch

% new x‐axis with “in‐between” points
N_new  = [4, 8, 16, 24, 32, 40, 48, 56, 64];

% interpolate (choose 'pchip' or 'spline')
y_new1 = interp1(N_orig, y_orig1, N_new, 'spline');
y_new2 = interp1(N_orig, y_orig2, N_new, 'spline');
y_new3 = interp1(N_orig, y_orig3, N_new, 'spline');

% now plot the smooth curve
figure;
plot( N_new, y_new1, '-o','DisplayName', 'New', 'LineWidth', LineW, 'MarkerSize', MarkS ); hold on;
plot( N_new, y_new2, '-o','DisplayName', 'New', 'LineWidth', LineW, 'MarkerSize', MarkS );
plot( N_new, y_new3, '-o','DisplayName', 'New', 'LineWidth', LineW, 'MarkerSize', MarkS );
xlabel('N');
ylabel('Mean WSR');
grid on;
legend('Location','best');
title(['Interpolated performance for ', schemes{sch}]);
%for sch = 1:numSchemes
%    plot(N, mean_wsr(:,sch), '-o', 'Color', colors{sch},'DisplayName', schemes{sch}, 'LineWidth', LineW, 'MarkerSize', MarkS);
%end
xlabel('Number of RIS elements, $N$', 'Interpreter', 'LaTeX');
ylabel('Weighted sum-rate (bps/Hz)');
legend('Location','best');
grid on;
