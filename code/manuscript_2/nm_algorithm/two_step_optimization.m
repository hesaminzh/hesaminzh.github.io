% close all
clear all
clc

%TWO_STEP_OPTIMIZATION  Run a two‑step algorithm that jointly
% optimizes the RIS susceptance matrix, BS precoder, and inter‑panel
% spacing in a multi‑panel mmWave system. The first step uses the
% partially-proximal ADMM (ppADMM) optimizer to update the susceptance
% matrix B and beamformer W given a fixed panel spacing. The second
% step updates the panel spacing using a derivative‑free Nelder–Mead
% (NM) algorithm. This script iterates between these steps until
% convergence and plots the convergence of t>he weighted sum rate.

% The NM algorithm implemented in the helper function
% NELDER_MEAD_PANEL_SPACING follows the standard procedure of
% reflection, expansion, contraction and shrink operations with
% coefficients alpha=1, gamma=2, rho=1/2 and delta=1/2. It searches for
% the spacing that maximizes the sum rate by minimizing the
% negative of the sum rate. The sum rate objective is obtained
% by evaluating the ppADMM optimizer and computing the true
% (unscaled) sum rate for the multi‑panel RIS configuration.

% ---------------- Simulation parameters ----------------
seed = 98;
rng(seed,'twister');

% BS antenna dimensions (ULA of N elements)
N_size = [8, 2];
Nx = N_size(1);
Ny = N_size(2);
N  = prod(N_size);

% Number of users
K  = 2;
% K  = 4;

% Number of panels (RIS sub‑surfaces)
I_panels = 2;
% I_panels = 4;

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
arch = 'single';

% ---------------- Parameter packing for objective function ----------------
% Many simulation parameters are reused inside the objective function that
% computes the sum rate for a given spacing. To avoid repeatedly
% capturing these variables with nested functions, we pack them into a structure.
% The objective function `objective_sum_rate_for_spacing` accepts this structure together with the spacing.
params = struct();
params.I_panels = I_panels;
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

% ---------------- Two‑step optimization loop ----------------
% Initial guess and step size for the panel spacing
% ds0  = 2;       % initial spacing (e.g., 2 metres)
% ds0  = 4;       % initial spacing (e.g., 4 metres)
% ds0  = 6;       % initial spacing (e.g., 6 metres)
ds0 = [0.5, 1, 2, 4, 6, 8, 10]; % list of initial panel spacing
step = 1;       % initial simplex step size

% Convergence tolerance and maximum iterations for NM
tol_ds    = 1e-3;
max_iter_nm = 2e2;

% Outer loop convergence settings
% max_outer   = 5;
max_outer   = 1e2;
outer_tol   = 1e-3;

% ds_history = zeros(max_outer,1);
ds_history = zeros(max_outer,length(ds0));
% rate_history = zeros(max_outer,1);
rate_history = zeros(max_outer,length(ds0));
best_ds_list = zeros(max_outer,1);

prev_ds = ds0;
outer_idx = 0;
while true
    outer_idx = outer_idx+1;
    fprintf('Outer iteration %d: optimizing panel spacing...\n', outer_idx);
    % Run NM optimizer for the spacing. The objective takes the
    % spacing and returns the negative sum rate using the packed
    % parameters.
    obj_handle = @(d_s) objective_sum_rate_for_spacing(d_s, params);
    [opt_ds, opt_val] = nelder_mead_panel_spacing(prev_ds, step, obj_handle, tol_ds, max_iter_nm);
    if isnan(opt_val)
        continue;
    end
    % Record best spacing and corresponding sum rate (note opt_val is negative sum rate)
    ds_history(outer_idx,:)   = opt_ds;
    best_ds_list(outer_idx) = opt_ds(1); % for mean and std convergence
    rate_history(outer_idx) = -opt_val;
    fprintf('  Optimal spacing = %.4f m, sum rate = %.6f bps/Hz\n', opt_ds, -opt_val);

    % Check convergence based on change in spacing     [[  Isn't it redundant?!  ]]
    if outer_idx > 1
        first_val = ds_history(outer_idx);
        first_val = first_val(1);
        second_val = ds_history(outer_idx-1);
        second_val = second_val(1);
        if abs( first_val - second_val ) < outer_tol
            fprintf('  Converged: change in spacing below tolerance (%.3e)\n', outer_tol);
            break;
        end
    end
    if outer_idx >= max_outer
        ds_history = NaN(max_outer,1);
        fprintf('  Reached max_outer iter (%g)\n', max_outer);
        break; % return in function !!!!!!!!!!!!!!!!                        return;!!!
    end

    % if ( outer_idx > 1 && abs(ds_history(outer_idx) - ds_history(outer_idx-1)) < outer_tol ) || outer_idx >= max_outer
    %     if outer_idx >= max_outer
    %         ds_history = NaN(max_outer,1);
    %         break; % return in function !!!!!!!!!!!!!!!!
    %     else
    %         fprintf('  Converged: change in spacing below tolerance (%.3e)\n', outer_tol);
    %         break;
    %     end
    % end

    % Update initial guess for next outer iteration
    prev_ds = opt_ds;
    % Optionally reduce the step size for the next NM call
    step = max(step/2, tol_ds*10);
end

% best_ds = ds_history(find(ds_history,1,'last'));
% Find the linear index of the last non-zero element in the entire matrix
last_nonzero_linear_idx = find(ds_history(:), 1, 'last');

% Check if any non-zero element was found
if ~isempty(last_nonzero_linear_idx)
    % Convert the linear index to row and column indices
    [row_idx, ~] = ind2sub(size(ds_history), last_nonzero_linear_idx);

    % best_ds is the first element of that specific row
    best_ds = ds_history(row_idx, 1);
else
    % If the entire matrix is zeros, best_ds remains undefined
    best_ds = NaN;
end
fprintf('Running ppADMM with final spacing %.4f m...\n', best_ds);

% Build panel positions for final spacing
panel_positions = zeros(I_panels,2);
center_offset  = -(I_panels-1)/2;
for ii = 1:I_panels
    panel_positions(ii,:) = [RIS_ref(1), RIS_ref(2) + (center_offset + (ii-1)) * best_ds];
end
% Generate channels and normalise
G_cell = cell(I_panels,1);
H_cell = cell(I_panels,1);
for ii = 1:I_panels
    panel_loc = panel_positions(ii,:);
    [G_cell{ii}, H_cell{ii}] = generate_SV_channel(N, Mx, My, user_loc, BS_loc, panel_loc);
end
G_agg = vertcat(G_cell{:});
H_agg = vertcat(H_cell{:});

powG = mean(abs(G_agg(:)).^2);
powH = mean(abs(H_agg(:)).^2);
sG = sqrt(target_entry_power / powG);
sH = sqrt(target_entry_power / powH);
if use_gain_norm
    Gw = sG * G_agg;
    Hw = sH * H_agg;
    if keep_sinr_invariant
        q = (sG * sH)^2;
        sigma2w = q * sigma2;
    else
        sigma2w = sigma2;
    end
else
    Gw      = G_agg;
    Hw      = H_agg;
    sigma2w = sigma2;
end

[B_final, W_final, rate_all, iters] = ppADMM_optimizer(Gw, Hw, sigma2w, PT, arch, I_panels, Gs, M_list, Mtot, N, K, eps_rate, max_iters, tol_W, rho, tau, xi, Z0);

Theta_final = theta_from_B(B_final, Z0);
final_true_rate = compute_sum_rate_theta(Theta_final, H_agg, G_agg, W_final, sigma2);
fprintf('Final TRUE unscaled sum rate (multi-BD-RIS, I=%d): %.6f bps/Hz\n', I_panels, final_true_rate);

% --------------- Plotting convergence ---------------
figure; hold on; grid on;
% Plot the ppADMM convergence for the final spacing
T = iters;
plot(0:T, rate_all(1:T+1), 'k--', 'LineWidth', 1.5, 'DisplayName', 'ppADMM convergence with optimized spacing');
% Plot the outer loop convergence of the sum rate versus iteration
% outer_iters = find(rate_history);
outer_iters = find(~isnan(rate_history));
plot(outer_iters-1, rate_history(outer_iters), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Two‑step outer convergence');
xlabel('Iteration index');
ylabel('Sum rate (bps/Hz)');
legend('Location','best');
title(sprintf('Convergence of Two‑Step Optimization (final spacing = %.2f m)', best_ds));
