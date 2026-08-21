function [P, Rate, converged] = FP_algorithm_v2(H, transmit_SNR, tolerance, noise)
    % INPUTS:
    % H: equivalent channel --- NtxK
    % transmit_SNR: transmit_SNR at the BS
    % tolerance: accuracy of convergence 
    % OUTPUTS:
    % P: beamformer matrix
    % Rate: sum-rate of the users
    % converged: flag indicating whether the algorithm converged

    [Nt, K] = size(H);
    %noise = ones(K, 1) .* dBm2pow(noise);
    noise = ones(K, 1) .* noise;

    %Pt = dBm2pow(transmit_SNR);
    Pt = transmit_SNR;
    P_p = (0.4 * Pt) / K;

    P = H ./ vecnorm(H) * sqrt(P_p);

    flag = 1;
    obj_past = 0;
    count = 0;
    max_iterations = 800;
    max_inner_iterations = 800; % Limit for inner loops
    converged = 0; % initialized the converged flag
    Rate = NaN; % initialize the Rate variable 

    while flag
        T_k = sum(square_abs(H' * P), 2) + noise;
        alpha_k = T_k ./ (T_k - square_abs(diag(H' * P))) - 1;
        beta_k = sqrt(1 + alpha_k) .* diag(H' * P) ./ T_k;

        mu_min = 0; % test whether mu = 0
        warning('off', 'MATLAB:nearlySingularMatrix');
        P = ((abs(beta_k') .* H) * (abs(beta_k') .* H)' + mu_min * eye(size(H, 1))) \ (transpose(sqrt(1 + alpha_k) .* beta_k) .* H);

        if sum(sum(pow_abs(P, 2))) - Pt <= 0
            continue;
        end

        % Bisection search to find optimal dual variable mu
        muleft = 0;
        muright = 1;
        inner_count = 0;

        while inner_count < max_inner_iterations
            P = ((abs(beta_k') .* H) * (abs(beta_k') .* H)' + eye(Nt) * muright) \ (transpose(sqrt(1 + alpha_k) .* beta_k) .* H);

            if sum(sum(pow_abs(P, 2))) - Pt <= 0
                break;
            end

            muright = muright * 10;
            inner_count = inner_count + 1;
        end

        if inner_count >= max_inner_iterations
            warning('Maximum inner iterations reached in first inner loop.');
            return;
        end

        %disp(['assert_1:', num2str(count)]);
        inner_count = 0;

        while inner_count < max_inner_iterations
            mu = (muleft + muright) / 2;

            P = ((abs(beta_k') .* H) * (abs(beta_k') .* H)' + eye(Nt) * mu) \ (transpose(sqrt(1 + alpha_k) .* beta_k) .* H);

            distance = sum(sum(pow_abs(P, 2))) - Pt;

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

        %disp(['assert_2:', num2str(count)]);

        T_ktest = sum(square_abs(H' * P), 2) + noise;
        rate_p = log2(T_ktest ./ (T_ktest - square_abs(diag(H' * P))));
        obj = sum(rate_p);

        if abs(obj - obj_past) <= tolerance
            flag = 0;
            converged = 1; % mark as converged
        else
            obj_past = obj;
            count = count + 1;
        end

        if count >= max_iterations
            warning('Maximum iterations reached. Convergence not achieved.');
            flag = 0;
        end
    end

    Rate = obj;
end

