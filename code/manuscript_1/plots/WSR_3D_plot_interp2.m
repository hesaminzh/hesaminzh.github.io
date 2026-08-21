clear all; close all; %setup_figures;

load('plot_WSR_3D_FC_GC_SC__K_2.mat');
% schemes = {'FC','GC','SC'};
schemes = {'FC','GC','SC','Random Phase','No RIS'};
numSchemes = numel(schemes);
% load('plot_WSR_3D__NoRIS_RandPhase__K_2.mat');
% -------------------------------
% 1) Build the grid for plotting
% -------------------------------
[ETA, Ngrid] = meshgrid(eta_list, N_list);
%   ETA(i,j)   = eta_list(j)
%   Ngrid(i,j) = N_list(i)

% Finer eta and N vectors
% eta_fine = linspace(min(eta_list), max(eta_list), 100);  % 100 pts
% eta_fine = eta_list;
eta_fine = 0:0.025:0.45;
% N_fine   = [4, 8, 16, 24, 32, 40, 48, 56, 64];
N_fine   = [8, 16, 24, 32, 40, 48, 56, 64];

% Create new fine meshgrid
[ETA_fine, N_fine_grid] = meshgrid(eta_fine, N_fine);
% -------------------------------
% 2) Predefine colors & handles
% -------------------------------
color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
color4 = [251,144,98] / 255; % Orange
color5 = [106,13,131] / 255; % Purple

colors = {color1, color2, color3, color4, color5};
h = gobjects(1,numSchemes);
Z_vals = cell(1,numSchemes);

% -------------------------------
% 3) Plot each scheme as a surface
% -------------------------------
figure('DefaultAxesFontSize',18); %hold on;
hold on; 
grid on;

for is = 1:numSchemes
    % extract the [numEta × numN] slice and transpose so that
    % rows→N_list,  cols→eta_list
    if is <= 3
        Z = squeeze( wsr3d(:, is, :) ).';  
        Z_smooth = interp2(ETA, Ngrid, Z, ETA_fine, N_fine_grid, 'makima');
        Z_vals{is} = Z;
    else
        clearvars wsr4d wsr3d schemes numSchemes
        load('plot_WSR_3D__NoRIS_RandPhase__K_2.mat');
        clearvars schemes
        schemes = {'FC','GC','SC','Random Phase','No RIS'};
        if strcmp(schemes{is},'No RIS')
            Z = squeeze( wsr3d(:, 1, :) ).';  
            Z_smooth = interp2(ETA, Ngrid, Z, ETA_fine, N_fine_grid, 'makima');
            Z_vals{is} = Z;
        elseif strcmp(schemes{is},'Random Phase')
            Z = squeeze( wsr3d(:, 2, :) ).';  
            Z_smooth = interp2(ETA, Ngrid, Z, ETA_fine, N_fine_grid, 'makima');
            Z_vals{is} = Z;
        end


    end
    % Interpolate using spline (smoothest) or 'pchip' (shape-preserving)
    % Z_smooth = interp2(ETA, Ngrid, Z, ETA_fine, N_fine_grid, 'makima');
    % Z = squeeze( wsr3d(:, is, :) ).';  

    hold on;
    h(is) = mesh( ...
        ETA_fine, N_fine_grid, Z_smooth, ...
        'FaceColor',   colors{is}, ...
        'EdgeColor',   'w', ...
        'LineWidth',   0.5, ...
        'FaceAlpha',   1, ...
        'DisplayName', schemes{is} ...
    );
end

% -------------------------------
% 4) Tidy up the axes & legend
% -------------------------------
colormap([])    % disable automatic colormap
%shading flat    % preserve your solid FaceColor
% view(45,30)     % tilt for 3D perspective
view(164,20);
% axis tight

legend(h, schemes, 'Location','nw');
rotate3d on;     % enable interactive rotation
view(163,20);
set(gcf, 'Color', [1,1,1]);
axis tight;

txtZ = text(0.61, 60, 7.5, ...
     '\textrm{Weighted sum-rate (bps/Hz)}', ...
     'Interpreter', 'LaTeX', ...
     'FontSize', 18, ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom');

txtZ.Rotation = 90;

txtX = text(0.20, 70, -1.4, '\textrm{Blockage density $\eta$}', ...
    'Interpreter', 'LaTeX', 'FontSize', 18, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

txtX.Rotation = 2.6;


txtY = text(0.56, 55, 0.3, '\textrm{Number of RIS elements $N$}', ...
    'Interpreter', 'LaTeX', 'FontSize', 18, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

txtY.Rotation = -38;

set(gca, 'YTick', 0:16:64);

%%% Datatips
%dt = datatip(h(1)); % adds datatip to the 1st surface
%dt.Position = [0.1, 64, Z_value];

% print(gcf, 'wsr_3D.eps', '-depsc2' );

