addpath('/Users/hesami/Documents/Research/ris_isac/Hesamilab/matlab_gpt/bdris_mmWave/bdris_new/data/power/')

%% colors
color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)
colors = {color1, color2, color3};


% Suppose you have:
schemes = {'FC', 'GC', 'SC'};       % nSchemes = 3
%K_list       = [2, 4, 6];             % nK       = 3
%P_dBm   = [0, 5, 10, 15];        % example power vector

% 1a) Fix colors per scheme
%colors = lines(numel(schemes));  

% 1b) Define distinct line styles (for bar edges) per K
styles = {'-', '--', ':'};       
%assert(numel(styles)==numel(K_list));




%% N = 64

% K = 6
load('power__64__K_6_.mat');

figure('DefaultAxesFontSize',12); 
ax1 = axes; hold(ax1,'on');

groupWidth = 0.8;  % total cluster width
nSchemes   = numel(schemes);
nK         = numel(K_list);

for ik = 1:nK
    % horizontal offset of this K-group within each power
    xOffset = (ik-(nK+1)/2) * (groupWidth/nK);
    
    for is = 1:nSchemes
        % x-positions for this (K,scheme)
        x = P_dBm + xOffset;
        y = squeeze( WSR_avg(:,is) );   % M×1
        
        % draw bars
        h = bar(ax1, x, y, groupWidth/nK,'stacked', 'FaceColor', colors{is}, ...
                'EdgeColor','k', 'LineStyle', styles{ik}, ...
                'LineWidth',1.2, 'DisplayName', ...
                sprintf('K=%d, %s', K_list(ik), schemes{is}));
        % we will not use this legend entry directly
    end
end

xlabel('P (dBm)'); ylabel('Average WSR');
grid on;





