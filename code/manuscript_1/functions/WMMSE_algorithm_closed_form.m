function [P, wsr] = WMMSE_algorithm_closed_form(C, P, epsilon, Pt, K, noise, weights)
% Input
% P: Initial beamforming matrix (M x K)
% C: Equivalent channel matrix (M x K)
% Pt: Total transmit power
% M: Number of transmit antennas
% K: Number of users
% noise: Noise variance for each user (K x 1)
% weights: Weights for each user (K x 1)

% Input
%P = exp(1i .* rand(M, M) .* 2 .* pi);
%P = sqrt(Pt/M) * (P ./ vecnorm(P));
%P = P;
%converged_WMMSE = false;

% Step 1: calculate WSR form P and Theta
noise=ones(K,1) .* noise;
wsr = WSR(C, P, K, noise, weights);
%wsr_all = wsr;
diff = 10;
step = 0;
while diff > epsilon
    wsr_pre = wsr;
    step = step + 1;

    % Step 3: update gk
    g_MMSE = MMSE_receiver(C, P, K, noise);

    % Step 4: update wk (updated weights)
    e_MMSE = MSE(C, P, K, noise, g_MMSE);
    w = weights ./ e_MMSE;
    
    % Step 5: update P with updated wk
    %[P, out] = update_P_CF(Pt, C, noise, w, g_MMSE, weights, epsilon);
    [P, out] = update_W_closed_form(Pt, C, noise, w, g_MMSE, weights);
    % check if CVX solved with no problem
    if isnan(out)
        %fprintf('Closed form solution didn''t work at step: %d\n', step);
        wsr = NaN;
        return;
    end

    wsr = WSR(C, P, K, noise, weights);
    diff = abs(wsr - wsr_pre);
    %wsr_all = [wsr_all wsr];
    %fprintf('WMMSE step-- %d, WSR-- %d\n\n', step, wsr);
    if step >= 8000%5000
        wsr = NaN;
        return;
    end
end

end


