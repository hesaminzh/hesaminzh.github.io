function W_new = update_W_ppadmm(W_prev, U, y, gamma, G, tau, PT, tol)
%UPDATE_W_PPADMM  Proximal pp-ADMM W-update with total power constraint.
%   W_prev (N x K), U (M x K), y,gamma (K x 1), G (M x N), tau>0, PT>0.
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)

    if nargin < 8 || isempty(tol), tol = 1e-6; end
    maxBisIter = 800;
    K  = size(U,2);
    N  = size(G,2);
    th = tau/2;

    % Build A = U' * G   (K x N), rows are u_k^H G
    A = U' * G;

    % Weights (paper): Q = A' diag((1+gamma)|y|^2) A,  c_k = (1+gamma_k) y_k (G^H u_k)
    % FPw = (1 + gamma(:)) .* (abs(y(:)).^2);   % Kx1
    FPw = (1 + gamma(:)) .* (pow_abs(y(:),2));   % Kx1
    t          = (1 + gamma(:)) .* y(:);          % Kx1
    % t          = (1 + gamma(:)) .* conj(y(:));          % Kx1

    % Numerics guards (soft)
    FPw(~isfinite(FPw)) = 0;
    t(~isfinite(t))     = 0;

    % Weighted Q and RHS (vectorized; avoids diag)
    Q   = A' * (FPw .* A);          % (N x N)
    % Q = zeros(N, N);
    % for k = 1:K
    %     uk = U(:,k);                  % M x 1
    %     v  = G' * uk;                 % N x 1  (equals G^H u_k)
    %     Q  = Q + FPw(k) * (v * v');   % v v^H = G^H u_k u_k^H G
    % end
    Q   = (Q + Q')/2;                      % Hermitian hygiene
    RHS = A' .* t.' + th * W_prev;         % (N x K) : columns c_k + (tau/2)w_k^t

    % First try eta=0: (Q + tau/2 I) W = RHS
    LHS0 = Q + th * eye(N);
    W0   = LHS0 \ RHS;
    p0   = norm(W0,'fro')^2;
    if p0 <= PT
        W_new = W0;
        return;
    end

    % EVD path for 1-D bisection on eta
    % [Uq, d] = eig(Q + th*eye(N), 'vector');      % d: Mx1 eigenvalues of Q+tau/2 I
    % d  = max(real(d), 0);                        % clip tiny negatives
    % % d  = d + 1e-12*max(1, norm(d,inf));          % small floor
    [Uq, D] = eig(Q + th*eye(N));
    lambda = real(diag(D));                     % eigenvalues as vector
    [lambda_sorted, idx] = sort(lambda, 'descend');
    Uq = Uq(:, idx);                      % reorder eigenvectors to match
    d = lambda_sorted;
    F_tilde = Uq' * RHS;                         % (N x K)
    % s = sum(abs(F_tilde).^2, 2);                 % s_n = ||row n||^2
    s = sum(pow_abs(F_tilde,2), 2);                 % s_n = ||row n||^2

    % Power p(eta) = sum_n s_n / (d_n + eta)^2
% --- after building Q, RHS, and EVD giving Uq, d (vector), F_tilde, s ---

    p_of = @(eta) sum( s ./ (d + eta).^2 );    % power vs eta (strictly decreasing)

    % Bracket (unchanged)
    eta_low = 0;
    p_low   = p_of(eta_low);
    if p_low <= PT
        W_new = (Q + (tau/2)*eye(N)) \ RHS;    % inactive: unconstrained proximal
        return;
    end
    eta_high = 1; maxGrow = 60;
    while p_of(eta_high) > PT && maxGrow > 0
        eta_high = 2*eta_high; maxGrow = maxGrow-1;
    end
    if p_of(eta_high) > PT
        eta_high = eta_high * 2^10;            % last resort
    end

    % === power-based bisection (key change) ===
    p_tol_rel = 1e-8;                     % or tie to 'tol'
    p_tol_abs = p_tol_rel * max(1, PT);
    for it = 1:maxBisIter
        eta = 0.5*(eta_low + eta_high);
        p   = p_of(eta);
        if abs(p - PT) <= p_tol_abs
            break;
        end
        if p > PT
            eta_low = eta;
        else
            eta_high = eta;
        end
    end

    % Build W at the eta that meets the power equality (not the last feasible)
    W_new = Uq * ((1./(d + eta)) .* F_tilde);

    % Keep only a tiny safety shrink if roundoff pushed us over
    p_final = norm(W_new,'fro')^2;
    if p_final > PT*(1+1e-10)
        W_new = W_new * sqrt(PT / p_final);
    end
end

