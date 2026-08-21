function [P, out] = FP_update_P_v2(H, Pt, alpha_k, beta_k, tolerance, res_k, noise, weights)

[Nt, K] = size(H);
noise = ones(K, 1) .* noise;

% Regularization term for numerical stability
%regularization = 1e-8;

% Initial minimum value of mu
mu_min = 0;

% Maximum iterations for inner loops to prevent infinite loops
%max_inner_iterations = 800;
max_inner_iterations = 1500;

% Turn off warnings for nearly singular matrices
warning('off', 'MATLAB:nearlySingularMatrix')

Q = (abs(beta_k').*H)*(abs(beta_k').*H)' + mu_min * eye(size(H, 1));
% 1) build the “true” KKT matrix
A0        = Q + mu_min*eye(Nt);

% 2) decide how well-conditioned you need A to be
rcond_min = 1e-12;          % you want cond(A) <= 1e12

% 3) start with the tiniest relative loading
qmax      = norm(Q,2);      
regularization       = eps(qmax);      % ≃ 2e-19 · qmax

% 4) if Q==0 then eps(qmax)==0, so force a tiny absolute floor
abs_floor = 1e-16;          % 1e-14 is still 100× smaller than 1e-12
regularization       = max(regularization, abs_floor);

% 5) now grow only if needed
maxIter   = 50;             % you’ll rarely need more than a few
iter      = 0;
while true
    A    = A0 + regularization*eye(Nt);
    rc   = rcond(A);
    if rc >= rcond_min
        break;
    end
    if iter >= maxIter
        warning("Regularization search hit %d iters; rcond=%.3e", iter, rc);
        W = NaN;
        out = NaN;
        return;
    end
    regularization  = regularization * 10;
    iter = iter + 1;
end


% Calculate initial P with regularization term
P = ((abs(beta_k').*H)*(abs(beta_k').*H)' + mu_min * eye(size(H, 1)) + regularization * eye(size(H, 1))) \ ...
    (transpose(sqrt(weights.*(1 + alpha_k)) .* beta_k) .* H);

out = NaN;

% Check if power constraint is satisfied
if sum(sum(abs(P).^2)) - Pt > 0
    % Binary search for mu
    muleft = 0;
    muright = 1;

    % Increase muright until the power constraint is satisfied
    inner_count = 0;
    while inner_count < max_inner_iterations
        P = ((abs(beta_k').*H)*(abs(beta_k').*H)' + eye(Nt) * muright + regularization * eye(Nt)) \ ...
            (transpose(sqrt(weights.*(1 + alpha_k)) .* beta_k) .* H);

        if sum(sum(abs(P).^2)) - Pt <= 0
            break;
        end

        muright = muright * 10;
        inner_count = inner_count + 1;
    end

    if inner_count >= max_inner_iterations
        warning('Maximum inner iterations reached in first inner loop.');
        return;
    end

    % Perform binary search to find the optimal mu
    inner_count = 0;
    while inner_count < max_inner_iterations
        mu = (muleft + muright) / 2;

        P = ((abs(beta_k').*H)*(abs(beta_k').*H)' + eye(Nt) * mu + regularization * eye(Nt)) \ ...
            (transpose(sqrt(weights.*(1 + alpha_k)) .* beta_k) .* H);

        distance = sum(sum(abs(P).^2)) - Pt;

        if abs(distance) <= tolerance
            break;
        end

        if distance > 0
            muleft = mu;
        else
            muright = mu;
        end
        inner_count = inner_count + 1;
    end

    if inner_count >= max_inner_iterations
        warning('Maximum inner iterations reached in second inner loop.');
        return;
    end
end

% Calculate T_kcvx and E_k
T_kcvx = sum(abs(H' * P).^2, 2) + noise;
E_k = -abs(beta_k).^2 .* T_kcvx + 2 * sqrt(weights.*(1 + alpha_k)) .* real(conj(beta_k) .* diag(H' * P)) + res_k;
out = sum(E_k);

% Turn warnings back on
warning('on', 'MATLAB:nearlySingularMatrix')

end

