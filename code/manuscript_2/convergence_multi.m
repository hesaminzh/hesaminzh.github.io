close all
clear all
clc

% addpath(genpath('./function/'))

%% Simulation parameters

seed = 93;
rng(seed,'twister');

% BS antenna dimensions (ULA of N elements)
N_size = [8, 2];
% N_size = [4, 2];
Nx = N_size(1);
Ny = N_size(2);
N  = prod(N_size);

% Number of users (increase K to stress spatial separation)
K  = 2;
% K  = 4;

% Number of panels (RIS sub‑surfaces)
I_panels = 2;
% I_panels = 4;
% I_panels = 6;

Mx = 4;
%M_RIS_y  = [2,4,8,16]; % RIS elements
My  = 8; % RIS elements
M = Mx*My;

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
% if generate_user_location is unavailable, draw users manually.
if exist('generate_user_location','file') == 2
    [user_loc, user_azimuths, d_BR, d_RU, ~] = generate_user_location(user_range, user_center, RIS_ref, BS_loc, K);
else
    % Sample K users uniformly from a disc defined by user_center and user_range
    radii = user_range(1) + (user_range(2) - user_range(1)) * sqrt(rand(K,1));
    angles = 2*pi*rand(K,1);
    user_loc = user_center + [radii .* cos(angles), radii .* sin(angles)];
end

%% Place RIS panels along the wall to create panel diversity

% Panels lie at x = RIS_ref(1).  Spread them along the y-axis with
% spacing of 10 metres. The central panel retains the reference y
% coordinate of RIS_ref; other panels are offset relative to this.
% panel_spacing = 10;
panel_spacing = 2;
% panel_spacing = 1;
panel_positions = zeros(I_panels,2);
center_offset  = -(I_panels-1)/2;
for i = 1:I_panels
    panel_positions(i,:) = [RIS_ref(1), RIS_ref(2) + (center_offset + (i-1)) * panel_spacing];
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




%% -------- Normalization controls --------
use_gain_norm       = true;   % <== TURN ON normalization during optimization
keep_sinr_invariant  = true;  % if false, also scale noise so noise power stays identical
target_entry_power  = 1;      % target mean |entry|^2 after scaling








%% ===========================
%  Multi-BD-RIS extension (I panels)
%  ===========================

% Number of BD-RIS panels
% M_list   = [16, 16];            % elements per panel (sum can differ from single-panel M)
M_list = M*ones(1,I_panels);
Mtot = sum(M_list);


% ---- Generate channels per panel ----
G_cell = cell(I_panels,1);
H_cell = cell(I_panels,1);
for i = 1:I_panels
    panel_loc = panel_positions(i,:);
    % Use the extended SV generator with increased multipath richness (L=6)
    [G_cell{i}, H_cell{i}] = generate_SV_channel(N, Mx, My, user_loc, BS_loc, panel_loc);
end

% ---- Aggregate channels ----
G_agg = vertcat(G_cell{:});   % Mtot x N
H_agg = vertcat(H_cell{:});   % Mtot x K

%% -------- Normalization for multi-panel case --------

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
    % PTw = PT;
else
    Gw      = G_agg;
    Hw      = H_agg;
    sigma2w = sigma2;
    % PTw     = PT;
end

%% -------- Build block-diagonal architecture for multi-BD-RIS --------
% Here, for simplicity, we choose the same architecture for EACH panel as
% the single-panel architecture 'arch' used above (e.g., 'fully').
% You can easily generalize to per-panel architectures if needed.

arch = 'single';   % or 'group', 'tridiagonal', 'single', 'stem', 'band', ...
% arch = 'group';   % or 'group', 'tridiagonal', 'single', 'stem', 'band', ...

% First, build local S_sets for each panel (using their local sizes)
S_local_cell  = cell(I_panels,1);
for i = 1:I_panels
    Mi  = M_list(i);
    switch arch
        case 'fully'
            S_local_cell{i} = arch_sparsity_sets_all(Mi, 'fully', []);
        case 'group'
            S_local_cell{i} = arch_sparsity_sets_all(Mi, 'group', struct('groups',[],'group_size', Gs));
        case 'tridiagonal'
            S_local_cell{i} = arch_sparsity_sets_all(Mi, 'tridiagonal', []);
        case 'single'
            S_local_cell{i} = arch_sparsity_sets_all(Mi, 'single', []);
        case 'stem'
            q = compute_q_band_stem(Mi, N, ones(K,1));
            S_local_cell{i} = arch_sparsity_sets_all(Mi, 'stem', struct('q',q));
        case 'band'
            q = compute_q_band_stem(Mi, N, ones(K,1));
            S_local_cell{i} = arch_sparsity_sets_all(Mi, 'band', struct('q',q));
        otherwise
            error('Unsupported architecture for multi-BD-RIS: %s', arch);
    end
end

% Now, build GLOBAL S_sets_big over indices 1..Mtot by offsetting
% each panel's S_local indices into the global index space.
S_sets_big = cell(Mtot,1);
row_offset = 0;
for i = 1:I_panels
    Mi   = M_list(i);
    Sloc = S_local_cell{i};
    for r = 1:Mi
        global_row = row_offset + r;
        S_local_row = Sloc{r};
        S_sets_big{global_row} = row_offset + S_local_row;  % shift indices by row_offset
    end
    row_offset = row_offset + Mi;
end

[offsets_big, num_vars_big] = pack_offsets(S_sets_big);

%% -------- Initialization for multi-panel pp-ADMM --------
PTw = PT;
W      = (randn(N,K) + 1i*randn(N,K))/sqrt(2);
W      = W/norm(W,'fro')*sqrt(PTw);
U      = Hw;                       % start with U = H (normalized)
Lambda = zeros(Mtot, K);           % dual variable
B      = zeros(Mtot);              % susceptance
x      = zeros(num_vars_big, 1);         % packed upper-tri vars
I  = eye(Mtot);

% Tracking convergence
rate_all   = zeros(max_iters+1, 1);
Theta      = theta_from_B(B, Z0);
rate_prev  = compute_sum_rate_theta(Theta, Hw, Gw, W, sigma2w);
rate_all(1)= rate_prev;

fprintf('========== Multi-BD-RIS (I = %d, arch = %s) ==========\n', I_panels, arch);
fprintf('it=%3d  R=%.6f bps/Hz (init)\n', 0, rate_prev);

it = 0;
while true
    it = it + 1;

    % (1) y, gamma update (FP block)
    [y, gamma] = updateYGamma_ppadmm(U, Gw, W, sigma2w);

    % (2) W update
    W = update_W_ppadmm(W, U, y, gamma, Gw, tau, PTw, tol_W);

    % (3) B update (on the BIG block-diagonal B)
    [B, x] = update_B_ppadmm(U, Hw, Lambda, S_sets_big, offsets_big, x, rho, xi, Z0);

    % (4) U update
    U = update_U_ppadmm(B, Gw, Hw, W, y, gamma, Lambda, rho, Z0);

    % (5) Dual ascent
    residual = (I - 1i*Z0*B)*U - (I + 1i*Z0*B)*Hw;
    Lambda   = Lambda + rho * residual;

    % ---- Rate and stopping check ----
    Theta  = theta_from_B(B, Z0);
    rate_curr = compute_sum_rate_theta(Theta, Hw, Gw, W, sigma2w);
    rate_all(it+1) = rate_curr;

    % if rate_prev == 0
    %     dRate = abs(rate_curr - rate_prev);
    % else
    %     dRate = abs(rate_curr - rate_prev) / max(1, abs(rate_prev));
    % end

    dRate = abs(rate_curr - rate_prev);

    if (dRate < eps_rate) || (it >= max_iters)
        if it >= max_iters
            fprintf('Multi-BD-RIS: Reached max iters (%d). Final R=%.6f\n', it, rate_curr);
        else
            fprintf('Multi-BD-RIS: Converged at it=%d. Final R=%.6f\n', it, rate_curr);
        end
        fprintf('\n');
        break;
    end

    rate_prev = rate_curr;
end

% Inspect sparsity pattern
% figure; spy(B ~= 0);

% ---- Final TRUE unscaled rate for multi-BD-RIS ----
W_true = W;
Theta_true = theta_from_B(B, Z0);
final_true_rate = compute_sum_rate_theta(Theta_true, H_agg, G_agg, W_true, sigma2);

fprintf('Final TRUE unscaled (multi-BD-RIS, I=%d): %.6f bps/Hz\n', I_panels, final_true_rate);


% --- Add multi-BD-RIS convergence curve ---
T = it;   % number of iterations in the multi-BD-RIS run

% Choose a new style for the multi-panel curve
style  = 'k--'; % black dashed
label  = sprintf('multi-%s, I = %d', arch, I_panels);

figure; hold on; grid on;
plot(0:T, rate_all(1:T+1), ...
     style, 'LineWidth', 1.5, ...
     'DisplayName', label);

legend('Location','best', 'Interpreter', 'LaTeX');

% ========================= End of script =========================



