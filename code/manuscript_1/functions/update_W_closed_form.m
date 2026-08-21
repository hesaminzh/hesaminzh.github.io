function [W, out] = update_W_closed_form(Pt, H, noise, w, g, weights)

%converged_WMMSE = true;
[Nt, K] = size(H);
noise = ones(K, 1) .* noise;

% Regularization term for numerical stability
%regularization = 1e-8;
%regularization = 0;

%tolerance = epsilon;
tolerance = 1e-6;

% Initial minimum value of mu
mu_min = 0; % initial guess

% Maximum iterations for inner loops to prevent infinite loops
%max_inner_iterations = 800;
%max_inner_iterations = 1000;
max_inner_iterations = 1500;

% Turn off warnings for nearly singular matrices
warning('off', 'MATLAB:nearlySingularMatrix')


% closed-from solution of W
Q = ( abs( (sqrt(w).*g)' ) .* H ) * ( abs( (sqrt(w).*g)' ) .* H )';

%%q_norm = norm(Q,2);
%%regularization = max( 1e-12, 1e-6*q_norm);
%%regularization = max(eps(norm(Q)), 1e-8);
%
%%%%%
%% --- 1) estimate the largest eigenvalue of Q by its 2-norm ---
%qnorm = norm(Q,2);   % ≈ max eigenvalue of Q
%
%% --- 2) choose a target max condition number, e.g. 1e10 or 1e12 ---
%cond_max = 1e10;     
%
%% --- 3) compute the minimum regularization that guarantees cond(Q+regularization I) <= cond_max ---
%regularization = qnorm/cond_max;
%
%% --- 4) also enforce at least machine–eps relative to qnorm ---
%regularization = max(regularization, eps*qnorm);
%
%% now regularization is typically in the 1e–13…1e–15 range for your Q
%%%%

%A0         = Q + mu_min*eye(Nt);
%rcond_min  = 1e-12;
%cond_max   = 1/rcond_min;
%qmax       = norm(Q,2);
%
%reg_cond   = qmax/(cond_max - 1);
%reg_floor  = max(1e-12, 1e-7*qmax);
%reg_eps    = eps(qmax);
%
%regularization = max([reg_cond, reg_floor, reg_eps]);
%
%
%%A0       = Q + mu_min*eye(Nt);
%%rcond_min= 1e-12;        % must be < 1
%%%rcond_min = 1e-11;
%%%rcond_min = 1e-10;
%%regularization      = eps*norm(Q,2);
%%regularization      = max(regularization, 1e-8);
%%maxIter  = 50;           % safety cap
%%iter     = 0;
%
%while true
%  A    = A0 + regularization*eye(Nt);
%  rc   = rcond(A);
%  if rc >= rcond_min
%    break;
%  end
%  if iter >= maxIter
%    warning("Regularization search hit %d iters; rcond=%.3e", iter, rc);
%    W = NaN;
%    out = NaN;
%    return;
%  end
%  regularization  = regularization * 10;
%  iter = iter + 1;
%end

%P = A \ B;


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

% 6) safe, one-shot inversion
%P = (A0 + reg*eye(Nt)) \ B;



W = ( Q + mu_min * eye(Nt) + regularization * eye(Nt)  ) \ ...
    ( transpose( w .* conj(g)) .* H );

out = NaN;


% Check if power constraint is satisfied
if sum(sum(abs(W).^2)) - Pt > 0
    % Binary search for mu
    muleft = 0;
    muright = 1;

    % Increase muright until the power constraint is satisfied
    inner_count = 0;
    while inner_count < max_inner_iterations % Adjust the upper bound (muright)
        Q = ( abs( (sqrt(w).*g)' ) .* H ) * ( abs( (sqrt(w).*g)' ) .* H )';
        W = ( Q + muright * eye(Nt) + regularization * eye(Nt)  ) \ ...
            ( transpose( w .* conj(g)) .* H );

        if sum(sum(abs(W).^2)) - Pt <= 0
            break;
        end

        muright = muright * 10;
        inner_count = inner_count + 1;
    end

    if inner_count >= max_inner_iterations
        warning('Maximum inner iterations reached in first inner loop.');
        %converged_WMMSE = false;
        return;
    end

    % Perform binary search to find the optimal mu
    inner_count = 0;
    while inner_count < max_inner_iterations
        mu = (muleft + muright) / 2;

        Q = ( abs( (sqrt(w).*g)' ) .* H ) * ( abs( (sqrt(w).*g)' ) .* H )';
        W = ( Q + mu * eye(Nt) + regularization * eye(Nt)  ) \ ...
            ( transpose( w .* conj(g)) .* H );

        distance = sum(sum(abs(W).^2)) - Pt;

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
        %converged_WMMSE = false;
        return;
    end
end



out = WSR(H, W, K, noise, weights);

% Turn warnings back on
warning('on', 'MATLAB:nearlySingularMatrix')


end
