clear all; close all; setup_figures;

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

% -------------------------------
% 3) Plot each scheme as a surface
% -------------------------------
figure('DefaultAxesFontSize',14); %hold on;
hold on; 
grid on;

for is = 1:numSchemes
    % extract the [numEta × numN] slice and transpose so that
    % rows→N_list,  cols→eta_list
    if is <= 3
        Z = squeeze( wsr3d(:, is, :) ).';  
    else
        clearvars wsr4d wsr3d Z schemes numSchemes
        load('plot_WSR_3D__NoRIS_RandPhase__K_2.mat');
        clearvars schemes
        schemes = {'FC','GC','SC','Random Phase','No RIS'};
        if strcmp(schemes{is},'No RIS')
            Z = squeeze( wsr3d(:, 1, :) ).';  
        elseif strcmp(schemes{is},'Random Phase')
            Z = squeeze( wsr3d(:, 2, :) ).';  
        end


    end
    % Z = squeeze( wsr3d(:, is, :) ).';  

    hold on;
    h(is) = surf( ...
        ETA, Ngrid, Z, ...
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
view(160,20)
% axis tight

legend(h, schemes, 'Location','nw');
rotate3d on;     % enable interactive rotation
view(160,20);
set(gcf, 'Color', [1,1,1]);
axis tight;

txtZ = text(0.64, 60, 7.5, ...
     '\textrm{Weighted sum-rate (bps/Hz)}', ...
     'Interpreter', 'LaTeX', ...
     'FontSize', 14, ...
     'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom');

txtZ.Rotation = 90;

txtX = text(0.20, 70, -1.2, '\textrm{Blockage density $\eta$}', ...
    'Interpreter', 'LaTeX', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

txtX.Rotation = 3.5;


txtY = text(0.56, 48, 0.2, '\textrm{Number of RIS elements $N$}', ...
    'Interpreter', 'LaTeX', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

txtY.Rotation = -35;




