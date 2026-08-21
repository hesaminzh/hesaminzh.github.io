function [B, W, rate_all, it] = ppADMM_optimizer(Gw, Hw, sigma2w, PT, arch, I_panels, Gs, M_list, Mtot, N, K, eps_rate, max_iters, tol_W, rho, tau, xi, Z0)
% PPADMM_OPTIMIZER implements pp-ADMM algorithm that optimizes susceptance matrix B and beamforming matrix W in sum-rate maximization problem
% Detailed input arguments go here
%
% Hesam Nezhadmohammad (last revised: 20/2/2026)


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

% fprintf('========== Multi-BD-RIS (I = %d, arch = %s) ==========\n', I_panels, arch);
% fprintf('it=%3d  R=%.6f bps/Hz (init)\n', 0, rate_prev);

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

    if rate_prev == 0
        dRate = abs(rate_curr - rate_prev);
    else
        dRate = abs(rate_curr - rate_prev) / max(1, abs(rate_prev));
    end

    % dRate = abs(rate_curr - rate_prev);

    if (dRate < eps_rate) || (it >= max_iters)
        if it >= max_iters
            fprintf('Multi-BD-RIS: Reached max iters (%d). Final R=%.6f\n', it, rate_curr);
            B = NaN;
            W = NaN;
        % else
        %     fprintf('Multi-BD-RIS: Converged at it=%d. Final R=%.6f\n', it, rate_curr);
        end
        % fprintf('\n');
        break;
    end

    rate_prev = rate_curr;
end

% Inspect sparsity pattern
% figure; spy(B ~= 0);

end
