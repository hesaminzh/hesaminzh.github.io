%--- Define parameters --------------------------------------------
schemes = {'FC','GC','SC'};        % 3 schemes
Kvals    = [2, 4, 6];              % values of K
styles   = {':','--','-'};         % line styles for K=2,4,6
%% colors
color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
colors = {color1, color2, color3};

%--- Common P_dBm vector (in dBm) -------------------------------
%P_dBm = [20,25,30,35,40,45];       % length = 6

%--- Create figure ------------------------------------------------
%figure('DefaultAxesFontSize',12);
%ax1 = axes; hold(ax1,'on');
%--- Parameters ---------------------------------------------------
schemes     = {'FC','GC','SC'};
Kvals       = [2,4,6];
nK          = numel(Kvals);
clusterWidth= 0.65;    % total normalized width for all 3 bars
gap         = 0.02;   % small space between bars
P_dBm       = [20,25,30,35,40,45];

%--- Compute bar width and offsets -------------------------------
barWidth = (clusterWidth - (nK-1)*gap)/nK;
%offsets  = ((1:nK) - (nK+1)/2) * (barWidth + gap);
offsets = [-clusterWidth-clusterWidth/1, 0, clusterWidth+clusterWidth/1]; %% it works

% e.g. offsets ≈ [-0.2867, 0, +0.2867] for clusterWidth=0.8, gap=0.02
%drawOrder = [1,3,2];
drawOrder = [1,2,3];


LineW = 1.5;
%LineW = 1.65;
MarkS = 8;
f = figure('Color', [1 1 1], 'DefaultAxesFontSize', 16);
ax1 = axes; hold(ax1,'on');
set(f, 'Units', 'pixels', 'Position', [100, 100, 500, 500]);  % 500×500 square figure

%--- Loop over each K and each scheme ----------------------------
for idx = 1:numel(Kvals)
    ik = drawOrder(idx);
    K = Kvals(ik);
    
    %--- Load the right .mat file --------------------------------
    if K==6
        datafile = 'power__64__K_6_.mat';
        load(datafile,'WSR_avg');     % WSR_avg is 6×3
        W = WSR_avg;                  % M×3
    else
        datafile = 'power__64__K_2_4.mat';
        load(datafile,'WSR_avg');     % WSR_avg is 2×6×3
        % extract the slice for this K
        idx_k = find(Kvals==K,1,'first');   % idx_k=1 for K=2, idx_k=2 for K=4
        W = squeeze( WSR_avg(idx_k,:,:) );   % M×3
    end
    
    %--- Compute horizontal offsets to separate the K‐groups --------
    %xOffset = (ik - (numel(Kvals)+1)/2) * (groupWidth/numel(Kvals));
    
    %--- Plot one bar‐group per scheme -----------------------------
    for is = 1:numel(schemes)
        %x = P_dBm + xOffset;
        x = P_dBm + offsets(ik);
        y = W(:,is);  % M×1
        
        bar(ax1, x, y, 'grouped',...
            'BarWidth', barWidth, ...
            'FaceColor', colors{is}, ...
            'EdgeColor','k', ...
            'LineStyle', styles{ik}, ...
            'LineWidth',LineW, ...
            'DisplayName', sprintf('K=%d, %s',K,schemes{is}) );
    end
end

%--- Axes labels and grid -----------------------------------------
xlabel('Transmit power (dBm)');
ylabel('Weighted sum-rate (bps/Hz)');
set(ax1,'XTick',P_dBm);
%grid on;



%% second axis for color legneds
ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'off');
hold(ax2, 'on');
% dummy plots
c1 = bar(ax2, NaN, NaN, 'LineStyle', 'none', 'FaceColor', colors{1}, 'LineWidth', LineW, 'DisplayName', 'FC');
c2 = bar(ax2, NaN, NaN, 'LineStyle', 'none', 'FaceColor', colors{2}, 'LineWidth', LineW, 'DisplayName', 'GC');
c3 = bar(ax2, NaN, NaN, 'LineStyle', 'none', 'FaceColor', colors{3}, 'LineWidth', LineW, 'DisplayName', 'SC');

%lg2 = legend(ax2, [c1, c2, c3], {'FC', 'GC', 'SC'}, 'Location', 'NorthEast');
lg2 = legend(ax2, [c1, c2, c3], {'FC', 'GC', 'SC'});
set(lg2, 'Color', 'white');

%% third axis linestyle legends
ax3 = axes('Position', get(ax1, 'Position'), 'Color', 'none', 'XTick', [], 'YTick', [], 'Box', 'off');
hold(ax3, 'on');
%h1 = bar(ax2, NaN, NaN, 'k-', 'LineWidth', 1.5, 'DisplayName', '-');
s1 = bar(ax3, NaN, NaN, 'LineStyle', '-', 'FaceColor', [1,1,1],'LineWidth', LineW, 'DisplayName', 'K = 6' );
s2 = bar(ax3, NaN, NaN, 'LineStyle', '--', 'FaceColor', [1,1,1],'LineWidth', LineW, 'DisplayName', 'K = 4' );
s3 = bar(ax3, NaN, NaN, 'LineStyle', ':', 'FaceColor', [1,1,1], 'LineWidth', LineW, 'DisplayName', 'K = 2');

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
%xlabel('\textrm{Transmit power (dBm)}', 'Interpreter', 'LaTeX');
%ylabel('\textrm{Weighted sum-rate (bps/Hz)}', 'Interpreter', 'LaTeX');
xlabel('Transmit power (dBm)');
ylabel('Weighted sum-rate (bps/Hz)');
set(gca,'xtick',P_dBm);
%xlim([P_dBm(1), P_dBm(end)]);
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
%grid on;
box on;

xlim_step = mean(diff(P_dBm));
xlim([P_dBm(1)-xlim_step/2, P_dBm(end)+xlim_step/2]);


%set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');
ax = gca;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridColor = [0.7 0.7 0.7]; % Light gray
ax.MinorGridAlpha = 0.5;           % Semi-transparent


% save as .eps format
%print(gcf, 'barplot_64.eps', '-depsc2' );
