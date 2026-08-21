% spacing_optimization_main.m

% --- System parameters (see report for detailed description) ---
seed = 91; rng(seed,'twister');
N_size = [8,2]; N = prod(N_size);
K = 4; I_panels = 4; Mx = 4; My = 4; % example: M=16
sigma2_dbm= -85; sigma2 = 10^((sigma2_dbm - 30)/10);
PT_dbm    = 30; PT    = 10^((PT_dbm - 30)/10);
BS_loc  = [0, 0]; RIS_ref = [50, 0];
user_center = [60,0]; user_range = [0,10];

% Outer optimisation parameters
d_min = 0.5; d_max = 8;            % feasible spacing interval
alpha = 0.5; delta = 0.2;          % step size and finite difference step
tol_d = 1e-2; maxOuterIter = 10;
d = 2;                             % initial spacing guess

for t = 1:maxOuterIter
    % ----- Step 1: optimise (B,W) for fixed spacing -----
    panel_positions = zeros(I_panels,2);
    center_offset = -(I_panels-1)/2;
    for i = 1:I_panels
        panel_positions(i,:) = [RIS_ref(1), RIS_ref(2) + (center_offset + (i-1)) * d];
    end
    % Generate channels G(d) and H(d)
    G_cell = cell(I_panels,1); H_cell = cell(I_panels,1);
    user_loc = generate_random_users(K,user_range,user_center);  % helper
    for i = 1:I_panels
        [G_cell{i}, H_cell{i}] = generate_SV_channel(N, Mx, My, user_loc, BS_loc, panel_positions(i,:));
    end
    G_agg = vertcat(G_cell{:}); H_agg = vertcat(H_cell{:});
    % Normalisation
    sG = sqrt(1/mean(abs(G_agg(:)).^2));
    sH = sqrt(1/mean(abs(H_agg(:)).^2));
    Gw = sG*G_agg; Hw = sH*H_agg; sigma2w = (sG*sH)^2*sigma2;
    % Run pp-ADMM to get B and W (assume single architecture for simplicity)
    [Bopt,Wopt] = ppADMM_optimizer(Gw,Hw,sigma2w,PT); % user-defined
    % Compute WSR for current spacing
    R_curr = compute_sum_rate_theta(theta_from_B(Bopt,50),Hw,Gw,Wopt,sigma2w);

    % ----- Step 2: update spacing via finite difference -----
    % Evaluate WSR at d+delta
    d_plus = min(d_max, d+delta);
    panel_positions_plus = panel_positions;
    for i = 1:I_panels
        panel_positions_plus(i,2) = RIS_ref(2) + (center_offset + (i-1)) * d_plus;
    end
    R_plus = evaluate_WSR_fixedBW(panel_positions_plus,Bopt,Wopt,N,Mx,My,K,user_range,user_center,BS_loc,sigma2);
    % Evaluate WSR at d-delta
    d_minus = max(d_min, d-delta);
    panel_positions_minus = panel_positions;
    for i = 1:I_panels
        panel_positions_minus(i,2) = RIS_ref(2) + (center_offset + (i-1)) * d_minus;
    end
    R_minus = evaluate_WSR_fixedBW(panel_positions_minus,Bopt,Wopt,N,Mx,My,K,user_range,user_center,BS_loc,sigma2);
    % Finite difference gradient
    g = (R_plus - R_minus)/(d_plus - d_minus);
    % Update spacing
    d_new = d + alpha * g;
    d_new = max(d_min, min(d_max, d_new));
    fprintf('Iter %d: d=%.2f m, WSR=%.3f bps/Hz, grad=%.3f\n', t, d, R_curr, g);
    if abs(d_new - d) < tol_d
        d = d_new; break;
    end
    d = d_new;
end

fprintf('Optimal spacing estimate: d*=%.2f m\n', d);

