%close all
clear all
clc
% addpath(genpath('./function/'))

%% Simulation parameters

seed = 98;
rng(seed,'twister');

% BS antenna dimensions (ULA of N elements)
N_size = [8, 2];
Nx = N_size(1);
Ny = N_size(2);
N  = prod(N_size);

% Number of users (increase K to stress spatial separation)
% K  = 2;
K  = 4;

% Number of panels (RIS sub‑surfaces)
% I_panels = 2;
I_panels = 4;
% I_panels = 6;

Mx = 4;
%M_RIS_y  = [2,4,8,16]; % RIS elements
My  = 8; % RIS elements
M = Mx*My;

M_list = M*ones(1,I_panels);
Mtot = sum(M_list);

Gs = 4;              % group size for GC
Z0 = 50;

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
% =========================
% AO + spacing update (robust while-loop version)
% =========================

% Outer-loop settings
t = 0;
converged = false;
maxOuterIter = 10;

% Tolerances
tol_R = 1e-4;            % relative WSR improvement tolerance
tol_g = 1e-6;            % gradient tolerance (optional)

% pp-ADMM params
rho  = 0.2;
tau  = 0.1;
xi   = 0.1;
Z0   = 50;

% Convergence controls
eps_rate  = 1e-3;
max_iters = 800;
tol_W     = 1e-3;




% Armijo backtracking settings for maximization
d_min = 0.5; d_max = 8; % feasible spacing interval
alpha0 = 0.5; % step size
beta = 0.5;              % shrink factor in (0,1)
c1 = 1e-4;               % Armijo sufficient increase constant
delta = 0.2; % finite difference step
tol_d = 1e-2;
maxOuterIter = 10;
maxBacktrack = 20;       % safeguard
% d = 2; % initial d
d = 2; % initial d

% Track previous accepted objective (outer objective at current d)
R_prev_outer = -Inf;

% Precompute panel index offsets (avoids repeated for-loop logic)
panel_idx = (0:I_panels-1).';
center_offset = -(I_panels-1)/2;

% panel coordinates
build_panel_positions = @(dval) [RIS_ref(1) * ones(I_panels,1), RIS_ref(2) + (center_offset + panel_idx) * dval];
% for i = 1:I_panels
%     panel_positions(i,:) = [RIS_ref(1), RIS_ref(2) + (center_offset + (i-1)) * d];
% end

arch = 'single';   % or 'group', 'tridiagonal', 'single', 'stem', 'band', ...

while ~converged && (t < maxOuterIter)

    t = t + 1;
    d_old = d;

    % ----- Step 1: optimize (B,W) for fixed spacing d -----
    panel_positions = build_panel_positions(d_old);

    % Generate channels G(d), H(d)
    G_cell = cell(I_panels,1);
    H_cell = cell(I_panels,1);

    % IMPORTANT:
    % Use a fixed user/channel realization within this AO iteration when
    % evaluating d+delta, d-delta, and Armijo trial points (common random numbers).
    [user_loc, ~, ~, ~, ~] = generate_user_location(user_range, user_center, RIS_ref, BS_loc, K);


    for i = 1:I_panels
        panel_loc = panel_positions(i,:);
        [G_cell{i}, H_cell{i}] = generate_SV_channel(N, Mx, My, user_loc, BS_loc, panel_loc);
    end

    G_agg = vertcat(G_cell{:});
    H_agg = vertcat(H_cell{:});

    % Normalization (keep consistent across objective evaluations)
    sG = sqrt(1 / mean(abs(G_agg(:)).^2));
    sH = sqrt(1 / mean(abs(H_agg(:)).^2));
    Gw = sG * G_agg;
    Hw = sH * H_agg;
    sigma2w = (sG * sH)^2 * sigma2;

    % pp-ADMM block update
    [Bopt, Wopt, ~, ~] = ppADMM_optimizer(Gw, Hw, sigma2w, PT, arch, I_panels, Gs, M_list, Mtot, N, K, eps_rate, max_iters, tol_W, rho, tau, xi, Z0);


    % Current objective at d_old
    Theta = theta_from_B(Bopt, Z0);
    R_curr = compute_sum_rate_theta(Theta, Hw, Gw, Wopt, sigma2w);

    % ----- Step 2: finite-difference gradient for spacing d -----
    d_plus = min(d_max, d_old + delta);
    d_minus = max(d_min, d_old - delta);

    denom_fd = d_plus - d_minus;
    if denom_fd <= eps
        warning('Finite-difference interval collapsed at boundary. Stopping AO.'); %!!!
        break;
    end

    % Evaluate WSR at d_plus and d_minus using SAME user_loc / random realization
    panel_positions_plus = build_panel_positions(d_plus);
    panel_positions_minus = build_panel_positions(d_minus);

    % NOTE:
    % Modify evaluate_WSR_fixedBW(...) so it reuses the same user_loc (and same
    % channel-generation randomness, if applicable) to reduce Monte Carlo noise.
    R_plus = evaluate_WSR_fixedBW(panel_positions_plus, Bopt, Wopt, N, Mx, My, I_panels, user_loc, BS_loc, sigma2, Z0, sG, sH);

    R_minus = evaluate_WSR_fixedBW(panel_positions_minus, Bopt, Wopt, N, Mx, My, I_panels, user_loc, BS_loc, sigma2, Z0, sG, sH);

    g = (R_plus - R_minus) / denom_fd;

    if ~isfinite(g)
        warning('Non-finite spacing gradient encountered. Stopping AO.');
        break;
    end

    % Optional early stop if gradient is tiny
    if abs(g) < tol_g
        fprintf('Iter %d: tiny gradient |g|=%.3e, stopping.\n', t, abs(g));
        break;
    end

    % ----- Step 3: Armijo backtracking (while loop) for spacing update -----
    alpha_try = alpha0;
    bt = 0;
    accepted = false;

    % Ascent direction for scalar maximization
    p = g;

    % If p is extremely small, no meaningful move
    if abs(p) < eps
        d_new = d_old;
        R_trial = R_curr;
        accepted = true;
    else
        while (bt < maxBacktrack)

            % Projected trial update
            d_trial = d_old + alpha_try * p;
            d_trial = max(d_min, min(d_max, d_trial));

            % If projection kills the step, stop line search
            if abs(d_trial - d_old) < tol_d
                d_new = d_old;
                R_trial = R_curr;
                accepted = true;
                break;
            end

            % Evaluate objective at trial spacing using fixed (B,W)
            panel_positions_trial = build_panel_positions(d_trial);

            R_trial = evaluate_WSR_fixedBW(panel_positions_trial, Bopt, Wopt, N, Mx, My, I_panels, user_loc, BS_loc, sigma2, Z0, sG, sH);

            % Armijo sufficient increase for maximization:
            % R(d_trial) >= R(d_old) + c1 * alpha * g * p
            % Here p = g, so RHS = R_curr + c1 * alpha * g^2
            if isfinite(R_trial) && (R_trial >= R_curr + c1 * alpha_try * g * p)
                d_new = d_trial;
                accepted = true;
                break;
            end

            % Backtrack
            alpha_try = beta * alpha_try;
            bt = bt + 1;
        end
    end

    % Fallback if line search fails
    if ~accepted
        d_new = d_old;
        R_trial = R_curr;
        alpha_try = 0;
    end

    % Update spacing
    d = d_new;

    % Outer stopping criteria
    % rel_improve = abs(R_curr - R_prev_outer);
    rel_improve = abs(R_curr - R_prev_outer) / max(1, abs(R_prev_outer));
    step_small = abs(d - d_old) < tol_d;
    obj_stall = (t > 1) && (rel_improve < tol_R);

    fprintf(['Iter %d: d=%.4f m, WSR=%.4f bps/Hz, g=%.4e, alpha=%.3e, BT=%d, relImp=%.3e\n'], t, d_old, R_curr, g, alpha_try, bt, rel_improve);

    % Optional: warm-start next line search with accepted alpha
    if accepted && (alpha_try > 0)
        alpha0 = min(2 * alpha_try, 1.0);   % mild growth cap
    end

    % Convergence decision
    converged = step_small || obj_stall;
    R_prev_outer = R_curr;
end
