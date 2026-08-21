function val = objective_sum_rate_for_spacing(d_s, My, M_list, Mtot, arch, params)
%OBJECTIVE_SUM_RATE_FOR_SPACING  Objective function for spacing optimization.
%
%   VAL = OBJECTIVE_SUM_RATE_FOR_SPACING(D_S, PARAMS) returns the
%   negative of the weighted sum rate (bps/Hz) obtained after
%   performing Step 1 of the two‑step algorithm with the panel
%   spacing D_S. This function is intended to be used with
%   derivative‑free optimizers (e.g. Nelder–Mead) which minimize
%   objective values. To maximize the sum rate we return the
%   negative of the sum rate.
%
%   PARAMS is a structure containing all required simulation
%   parameters and function handles. The fields of PARAMS are:
%     .I_panels, .RIS_ref, .user_loc, .BS_loc, .N, .Mx, .My,
%     .K, .sigma2, .PT, .arch, .M_list, .Mtot, .eps_rate,
%     .max_iters, .tol_W, .rho, .tau, .xi, .Z0,
%     .use_gain_norm, .keep_sinr_invariant, .target_entry_power,
%     .generate_SV_channel_handle, .theta_from_B_handle,
%     .compute_sum_rate_theta_handle, .ppADMM_optimizer_handle.
%
%   The returned value VAL = -final_sum_rate.

% Enforce a minimum positive spacing to avoid degenerate panel placement
min_spacing = 0.1;
if d_s < min_spacing
    val = inf;
    return;
end

% Unpack parameters
I_panels = params.I_panels;
Gs       = params.Gs;
RIS_ref  = params.RIS_ref;
user_loc = params.user_loc;
BS_loc   = params.BS_loc;
N        = params.N;
Mx       = params.Mx;
% My       = params.My;
K        = params.K;
sigma2   = params.sigma2;
PT       = params.PT;
% arch     = params.arch;
% M_list   = params.M_list;
% Mtot     = params.Mtot;
eps_rate = params.eps_rate;
max_iters= params.max_iters;
tol_W    = params.tol_W;
rho      = params.rho;
tau      = params.tau;
xi       = params.xi;
Z0       = params.Z0;
use_gain_norm       = params.use_gain_norm;
keep_sinr_invariant = params.keep_sinr_invariant;
target_entry_power  = params.target_entry_power;

% Compute panel positions for spacing d_s
panel_positions = zeros(I_panels,2);
center_offset  = -(I_panels-1)/2;
for ii = 1:I_panels
    panel_positions(ii,:) = [RIS_ref(1), RIS_ref(2) + (center_offset + (ii-1)) * d_s];
end

% Generate channel matrices per panel and aggregate them
G_cell = cell(I_panels,1);
H_cell = cell(I_panels,1);
for ii = 1:I_panels
    panel_loc = panel_positions(ii,:);
    [G_cell{ii}, H_cell{ii}] = generate_SV_channel(N, Mx, My, user_loc, BS_loc, panel_loc);
end
G_agg = vertcat(G_cell{:});
H_agg = vertcat(H_cell{:});

% Normalize channels
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

% Run Step 1: ppADMM optimizer to update B and W
[B_opt, W_opt, ~, ~] = ppADMM_optimizer(Gw, Hw, sigma2w, PT, arch, I_panels, Gs, M_list, Mtot, N, K, eps_rate, max_iters, tol_W, rho, tau, xi, Z0);

% Compute the true (unscaled) sum rate for the candidate spacing
Theta_opt = theta_from_B(B_opt, Z0);
final_rate = compute_sum_rate_theta(Theta_opt, H_agg, G_agg, W_opt, sigma2);

% Negative sign because Nelder–Mead performs minimization
val = -final_rate;
end
