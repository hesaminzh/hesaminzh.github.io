function R = compute_sum_rate_theta(Theta, Hru, G, W, sigma2)
%COMPUTE_SUM_RATE_THETA  Sum-rate using the Eq.(3) (Theta-form) SINR.
%
% Inputs:
%   Theta : M x M  (BD-RIS scattering matrix)
%   Hru   : M x K  (RIS-to-user channels, columns h_k)
%   G     : M x N  (BS-to-RIS channel)
%   W     : N x K  (BS beamforming matrix, columns w_k)
%   sigma2: scalar noise power
%
% Output:
%   R     : scalar sum-rate in bps/Hz
%
% Hesam Nezhadmohammad (19/10/2025)

    % K x K with entries S(k,j) = h_k^H * Theta * G * w_j
    S    = (Hru') * (Theta * (G * W));
    sig  = pow_abs(diag(S), 2);            % |h_k^H Θ G w_k|^2
    tot  = sum(pow_abs(S, 2), 2);          % ∑_j |h_k^H Θ G w_j|^2
    intf = tot - sig;                      % ∑_{j≠k} |·|^2
    R    = sum(log2(1 + (sig ./ (intf + sigma2))));
end

