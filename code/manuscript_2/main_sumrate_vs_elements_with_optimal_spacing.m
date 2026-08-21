% close all
clear all
clc

% addpath(genpath('./function/'))

%% Simulation parameters

seed = 91;
rng(seed,'twister');

% BS antenna dimensions (ULA of N elements)
N_size = [8, 2];
% N_size = [4, 2];
Nx = N_size(1);
Ny = N_size(2);
N  = prod(N_size);

% Number of users (increase K to stress spatial separation)
% K  = 2;
K  = 4;
% K  = 6;

% Number of panels (RIS sub‑surfaces)
% I_panels = 2;
I_panels = 4;
% I_panels = 6;

% M_RIS_y  = [2,4,8,12,16]; % RIS elements
Mx = 4;
M_RIS_y  = [2,4,8,16]; % RIS elements
M_RIS = Mx*M_RIS_y;

Gs = 4;              % group size for GC

sigma2_dbm= -85;
sigma2    = 10^((sigma2_dbm - 30)/10);
PT_dbm    = 30;
PT        = 10^((PT_dbm - 30)/10);

% BS and initial RIS reference location (meters)
BS_loc  = [0, 0];
RIS_ref = [50, 0];

% User cluster center and range. A larger radius improves angular
% diversity at the RIS. We choose a disc of radius up to 10 m centered at
% (60, 0).  Users are placed uniformly at random within this disc.
% user_center = [55, 5];
user_center = [60, 0];
user_range  = [0, 10];

% Generate random user locations;
if exist('generate_user_location','file') == 2
    [user_loc, user_azimuths, d_BR, d_RU, ~] = generate_user_location(user_range, user_center, RIS_ref, BS_loc, K);
else
    % Sample K users uniformly from a disc defined by user_center and user_range
    radii = user_range(1) + (user_range(2) - user_range(1)) * sqrt(rand(K,1));
    angles = 2*pi*rand(K,1);
    user_loc = user_center + [radii .* cos(angles), radii .* sin(angles)];
end


% pp-ADMM params
rho  = 0.2;
tau  = 0.1;
xi   = 0.1;
Z0   = 50;

% Convergence controls
eps_rate  = 1e-3;
max_iters = 800;
tol_W     = 1e-3;

% Normalization controls
use_gain_norm       = true;   % <== TURN ON normalization during optimization
keep_sinr_invariant  = true;  % if false, also scale noise so noise power stays identical
target_entry_power  = 1;      % target mean |entry|^2 after scaling

% ---------------- Parameter packing for objective function ----------------
% Many simulation parameters are reused inside the objective function that
% computes the sum rate for a given spacing. To avoid repeatedly
% capturing these variables with nested functions, we pack them into a structure.
% The objective function `objective_sum_rate_for_spacing` accepts this structure together with the spacing.
params = struct();
params.I_panels = I_panels;
params.Gs       = Gs;
params.RIS_ref  = RIS_ref;
params.user_loc = user_loc;
params.BS_loc   = BS_loc;
params.N        = N;
params.Mx       = Mx;
% params.My       = My;
params.K        = K;
params.sigma2   = sigma2;
params.PT       = PT;
% params.arch     = arch;
% params.M_list   = M_list;
% params.Mtot     = Mtot;
params.eps_rate = eps_rate;
params.max_iters= max_iters;
params.tol_W    = tol_W;
params.rho      = rho;
params.tau      = tau;
params.xi       = xi;
params.Z0       = Z0;
params.use_gain_norm       = use_gain_norm;
params.keep_sinr_invariant = keep_sinr_invariant;
params.target_entry_power  = target_entry_power;

% fixed_ds = 2; % fixed panel spacing for ppADMM
if K == 2
    load('neldermead_ds__stem_tridiagonal_single__K_2__nMonte_500.mat');
elseif K == 4
    load('neldermead_ds__stem_tridiagonal_single__K_4__nMonte_500.mat');
end

nMonte = 500;
numElements = numel(M_RIS_y);
arch_list = {'stem', 'tridiagonal', 'single'};   % or 'fully', 'tridiagonal', 'single'
% arch_list = {'single'};
numArch = numel(arch_list);
sum_rate = NaN(nMonte,numElements,numArch);

% Start parallel pool if not active
if isempty(gcp('nocreate'))
    parpool('local', 2);  % use up to 4 workers if available
%     parpool('local', 40);  % use up to 40 workers if available
end
%delete(gcp('nocreate')); % Shut down the parallel pool


parfor iMonte = 1:nMonte
    if mod(iMonte,10) == 0
        disp(['Iteration:', num2str(iMonte)]);
    end

    sum_rate_local = NaN(numElements,numArch);


    for iElement = 1:numElements
        rng(seed + 100000*iMonte + iElement,'twister');

        My  = M_RIS_y(iElement);
        M = Mx*My;
        M_list = M*ones(1,I_panels);
        Mtot = sum(M_list);
        for iArch = 1:numArch
            arch = arch_list{iArch};
            optimal_ds = avg_ds(iElement,iArch);
            val = objective_sum_rate_for_spacing(optimal_ds, My, M_list, Mtot, arch, params);
            ppadmm_rate = -val; % val was intended for minimization so we negate that.
            sum_rate_local(iElement,iArch) = ppadmm_rate;

        end % end of architecture loop

    end % end of element loop
    sum_rate(iMonte,:,:) = sum_rate_local;

end % end of monte loop


% Average and plotting
avg_sumrate = squeeze(mean(sum_rate,1,'omitnan')); % numElements x numArch

    


% Choose a new style for the multi-panel curve
styles = {'r-s','g-x','b-*','m-o'};

figure; hold on; grid on; box on;
for iArch = 1:numArch
    plot(M_RIS, avg_sumrate(:,iArch), styles{iArch}, 'LineWidth', 1.8);
end

legend(arch_list, 'Location', 'best', 'Interpreter', 'LaTeX');

save(['ppadmm_sumrate_w_optimal_spacing__K_',num2str(K),'__nMonte_',num2str(nMonte),'.mat'], 'M_RIS', 'avg_sumrate', 'styles', 'arch_list', 'numArch');

% ========================= End of script =========================
