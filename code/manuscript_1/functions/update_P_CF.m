function [P, out] = update_P_CF(Pt, H, noise, w, g, weights, epsilon)

%converged_WMMSE = true;
[Nt, K] = size(H);
noise = ones(K, 1) .* noise;

% Regularization term for numerical stability
regularization = 1e-8;

tolerance = epsilon;

% Initial minimum value of mu
mu_min = 0; % initial guess

% Maximum iterations for inner loops to prevent infinite loops
max_inner_iterations = 800;

% Turn off warnings for nearly singular matrices
warning('off', 'MATLAB:nearlySingularMatrix')


% closed-from solution of P
Q = ( abs( (sqrt(w).*g)' ) .* H ) * ( abs( (sqrt(w).*g)' ) .* H )';
P = ( Q + mu_min * eye(Nt) + regularization * eye(Nt)  ) \ ...
    ( transpose( w .* conj(g)) .* H );

out = NaN;


% Check if power constraint is satisfied
if sum(sum(abs(P).^2)) - Pt > 0
    % Binary search for mu
    muleft = 0;
    muright = 1;

    % Increase muright until the power constraint is satisfied
    inner_count = 0;
    while inner_count < max_inner_iterations % Adjust the upper bound (muright)
        Q = ( abs( (sqrt(w).*g)' ) .* H ) * ( abs( (sqrt(w).*g)' ) .* H )';
        P = ( Q + muright * eye(Nt) + regularization * eye(Nt)  ) \ ...
            ( transpose( w .* conj(g)) .* H );

        if sum(sum(abs(P).^2)) - Pt <= 0
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
        P = ( Q + mu * eye(Nt) + regularization * eye(Nt)  ) \ ...
            ( transpose( w .* conj(g)) .* H );

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
        %converged_WMMSE = false;
        return;
    end
end



out = WSR(H, P, K, noise, weights);

% Turn warnings back on
warning('on', 'MATLAB:nearlySingularMatrix')


end
