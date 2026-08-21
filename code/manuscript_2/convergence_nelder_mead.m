close all
clear all
clc

%CONVERGENCE_NELDER_MEAD  Run a two‑step algorithm that jointly
% optimizes the RIS susceptance matrix, BS precoder, and inter‑panel
% spacing in a multi‑panel mmWave system. The first step uses the
% partially-proximal ADMM (ppADMM) optimizer to update the susceptance
% matrix B and beamformer W given a fixed panel spacing. The second
% step updates the panel spacing using a derivative‑free Nelder–Mead
% (NM) algorithm. This script iterates between these steps until
% convergence and plots the convergence of the weighted sum rate.

% The NM algorithm implemented in the helper function
% NELDER_MEAD_PANEL_SPACING follows the standard procedure of
% reflection, expansion, contraction and shrink operations with
% coefficients alpha=1, gamma=2, rho=1/2 and delta=1/2. It searches for
% the spacing that maximizes the sum rate by minimizing the
% negative of the sum rate. The sum rate objective is obtained
% by evaluating the ppADMM optimizer and computing the true
% (unscaled) sum rate for the multi‑panel RIS configuration.

% ---------------- Simulation parameters ----------------
% seed = 106;
seed = 91;
rng(seed,'twister');

% BS antenna dimensions (ULA of N elements)
N_size = [8, 2];
Nx = N_size(1);
Ny = N_size(2);
N  = prod(N_size);

% Number of users
% K  = 2;
K  = 4;

% Number of panels (RIS sub‑surfaces)
% I_panels = 2;
I_panels = 4;

% Number of elements per panel
Mx = 4;
My = 8;
M  = Mx*My;

% Group size for GC (not used directly here)
Gs = 4;

% Noise and power parameters
sigma2_dbm= -85;
sigma2    = 10^((sigma2_dbm - 30)/10);
PT_dbm    = 30;
PT        = 10^((PT_dbm - 30)/10);

% BS and RIS reference locations (meters)
BS_loc  = [0, 0];
RIS_ref = [50, 0];

% User cluster centre and range
user_center = [60, 0];
user_range  = [0, 10];

% Generate random user locations if helper is unavailable
if exist('generate_user_location','file') == 2
    [user_loc, user_azimuths, d_BR, d_RU, ~] = generate_user_location(user_range, user_center, RIS_ref, BS_loc, K);
else
    radii = user_range(1) + (user_range(2) - user_range(1)) * sqrt(rand(K,1));
    angles = 2*pi*rand(K,1);
    user_loc = user_center + [radii .* cos(angles), radii .* sin(angles)];
end

% pp-ADMM parameters
rho  = 0.2;
tau  = 0.1;
xi   = 0.1;
Z0   = 50;

% Convergence controls for pp-ADMM
eps_rate  = 1e-3;
max_iters = 800;
tol_W     = 1e-3;

% Normalization controls
use_gain_norm       = true;
keep_sinr_invariant = true;
target_entry_power  = 1;

% Number of elements per RIS panel (vector) and total elements
M_list = M*ones(1,I_panels);
Mtot   = sum(M_list);

% Architecture for each panel
% arch = 'single';
arch = 'tridiagonal';

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
params.My       = My;
params.K        = K;
params.sigma2   = sigma2;
params.PT       = PT;
params.arch     = arch;
params.M_list   = M_list;
params.Mtot     = Mtot;
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

% Initial guess
% ds0 = [0.5, 1, 2, 4, 6, 8, 10]; % list of initial panel spacing
% ds0 = [0.5, 1, 2, 4, 6, 8]; % list of initial panel spacing
% ds0 = [1, 2, 4, 6, 8]; % list of initial panel spacing
Ns = 6; % in paper we have Ns = 5, here Ns is actually Ns+1
ds0_min = 0.5;
ds0_max = 8;
ds0 = ds0_min + (ds0_max - ds0_min).*rand(1,Ns);
ds0 = sort(ds0);

% Convergence tolerance and maximum iterations for NM
tol_ds    = 1e-3;
% max_iter_nm = 1e2;
max_iter_nm = 8e2;

obj_handle = @(d_s) objective_sum_rate_for_spacing(d_s, My, M_list, Mtot, arch, params);

[opt_ds, opt_val, mean_vals, std_vals] = nelder_mead_panel_spacing(ds0, obj_handle, tol_ds, max_iter_nm);
best_ds = opt_ds;
best_rate = -opt_val;
fprintf('Optimum distance is %g\n', best_ds);
fprintf('Best rate is %g\n', best_rate);

valid_indices = find(~isnan(mean_vals));
last_iter = valid_indices(end);
mean_to_plot = mean_vals(1:last_iter);
std_to_plot = std_vals(1:last_iter);

figure('DefaultAxesFontSize',16); %hold on;
hold on;

% Calculate upper and lower bounds for the shaded region
upper_bound = mean_to_plot + std_to_plot;
lower_bound = mean_to_plot - std_to_plot;

plot(mean_to_plot, 'b-', 'LineWidth', 2, 'DisplayName', 'Mean');

% Shade the standard deviation area
% We create vertices for the shaded area by going forward along the upper bound
% and backwards along the lower bound to close the shape.
x_values = 1:length(mean_to_plot);
% fill([x_values,fliplr(x_values)], [upper_bound, fliplr(lower_bound)], [0.8, 0.8, 0.8], ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Std')
fill([x_values,fliplr(x_values)], [upper_bound, fliplr(lower_bound)], [150, 206, 180]/255, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Std')
xlabel('Number of iterations');
ylabel('$D_q$', 'Interpreter', 'LaTeX');
legend;
set(gcf, 'Color', [1,1,1]);
set(gca, 'LineWidth',1);
xlim([0, last_iter]);
% ylim([0, max(upper_bound)+mean(upper_bound)/2]);
ylim([0, 3]);
% set(gca, 'xtick', M_RIS);
grid on;
box on;

print(gcf, 'convergence_nm__tridiagonal__K_4', '-dsvg' );
% print(gcf, 'convergence_nm__tridiagonal__K_2', '-dsvg' );


