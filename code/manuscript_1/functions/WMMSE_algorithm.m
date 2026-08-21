function [P, wsr, converged_WMMSE] = WMMSE_algorithm(C, P, epsilon, Pt, M, K, noise, weights)
% Input
%P = exp(1i .* rand(M, M) .* 2 .* pi);
%P = sqrt(Pt/M) * (P ./ vecnorm(P));

% Step 1: calculate WSR form P and Theta
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
    [P, ~, converged_WMMSE] = update_P(M, Pt, C, P, K, noise, w, g_MMSE);
    % check if CVX solved with no problem
    if converged_WMMSE == false
        fprintf('CVX did not solve properly at step %d\n', step);
        break;
    end

    wsr = WSR(C, P, K, noise, weights);
    diff = abs(wsr - wsr_pre);
    %wsr_all = [wsr_all wsr];
    %fprintf('WMMSE step-- %d, WSR-- %d\n\n', step, wsr);
end

end
