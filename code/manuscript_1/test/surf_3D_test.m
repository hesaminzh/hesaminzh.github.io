% --- toy grid & data ---
eta_list = linspace(0.1,0.9,5);
N_list   = [4, 8, 16, 32, 64];
[ETA, Ngrid] = meshgrid(eta_list, N_list);

% synthetic “outage” surfaces
Z_FC = (1-ETA).*log2(1 + Ngrid/8)*4 + 2;
Z_GC = Z_FC - 1;
Z_SC = Z_GC - 1;

schemes = {'FC','GC','SC'};

% fixed RGB for each scheme
colors = [ ...
  0.0000 0.4470 0.7410;   % blue for FC
  0.8500 0.3250 0.0980;   % orange for GC
  0.4660 0.6740 0.1880    % green for SC
];

% --- plotting ---
figure; 
hold on; 
grid on;

% preallocate handles
h = gobjects(1,3);

for i = 1:3
  switch schemes{i}
    case 'FC', Z = Z_FC;
    case 'GC', Z = Z_GC;
    case 'SC', Z = Z_SC;
  end
  h(i) = surf( ...
    ETA, Ngrid, Z, ...
    'FaceColor',   colors(i,:), ...   % solid color
    'EdgeColor',   'w',       ...     % white mesh lines
    'LineWidth',   0.5,       ...     % thin mesh
    'FaceAlpha',   0.8        ...     % slight transparency
  );
end

% set the view and labels
view(45,30)
xlabel('\eta (blockage density)')
ylabel('RIS elements, N')
zlabel('Outage probability')
title('Outage vs. \eta and RIS size for FC / GC / SC')

% explicit legend using the handles
legend(h, schemes, 'Location','northeast')

% ensure we’re not accidentally using a colormap
colormap([])    
%shading flat    % flat shading preserves the FaceColor
