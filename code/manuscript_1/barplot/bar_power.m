%--- Define parameters --------------------------------------------
schemes = {'FC','GC','SC'};        % 3 schemes
Kvals    = [2, 4, 6];              % values of K
styles   = {'-','--',':'};         % line styles for K=2,4,6
%% colors
color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
colors = {color1, color2, color3};

%--- Common P_dBm vector (in dBm) -------------------------------
P_dBm = [20,25,30,35,40,45];       % length = 6

%--- Create figure ------------------------------------------------
figure('DefaultAxesFontSize',12);
ax1 = axes; hold(ax1,'on');
%groupWidth = 0.8;  % total width of each P_dBm cluster
%groupWidth = 0.95;  % total width of each P_dBm cluster
groupWidth = 1.0;  % total width of each P_dBm cluster
offsetScale = 1.3;
barOffsets = [-offsetScale, 0, +offsetScale];
barWidth = groupWidth/(2*offsetScale);
drawOrder = [1,3,2];

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
        x = P_dBm + barOffsets(ik)*barWidth;
        y = W(:,is);  % M×1
        
        bar(ax1, x, y, ...
            'BarWidth', groupWidth/numel(Kvals), ...
            'FaceColor', colors{is}, ...
            'EdgeColor','k', ...
            'LineStyle', styles{ik}, ...
            'LineWidth',1.2, ...
            'DisplayName', sprintf('K=%d, %s',K,schemes{is}) );
    end
end

%--- Axes labels and grid -----------------------------------------
xlabel('P (dBm)');
ylabel('Average WSR');
set(ax1,'XTick',P_dBm);
xlim_step = mean(diff(P_dBm));
grid on;

%--- First legend: schemes (colors) ------------------------------
% Create dummy patches for the color legend
for is = 1:numel(schemes)
    p(is) = patch(NaN, NaN, colors{is}, 'DisplayName', schemes{is});
end
lgd1 = legend(ax1, p, schemes, 'Location','NorthEast');
lgd1.Title.String = 'Scheme';
lgd1.Color = 'white';

%--- Second legend: K (line styles) ------------------------------
% Overlay a transparent axes for the style legend
ax2 = axes('Position', get(ax1,'Position'), ...
           'Color','none', 'XTick',[], 'YTick',[], 'Box','off');
hold(ax2,'on');

% Dummy lines to represent each K
for ik = 1:numel(Kvals)
    h(ik) = plot(ax2, NaN, NaN, 'k', ...
                 'LineStyle', styles{ik}, ...
                 'LineWidth',1.5, ...
                 'DisplayName', sprintf('K = %d',Kvals(ik)));
end
lgd2 = legend(ax2, h, arrayfun(@(k)sprintf('K = %d',k),Kvals,'uni',0), ...
              'Location','SouthEast');
lgd2.Title.String = 'K value';
lgd2.Color = 'white';

% Prevent the overlay axes from capturing clicks/visibility
set(ax2,'HitTest','off','HandleVisibility','off');
linkaxes([ax1 ax2]);

xlim([P_dBm(1)-xlim_step/2, P_dBm(end)+xlim_step/2]);
