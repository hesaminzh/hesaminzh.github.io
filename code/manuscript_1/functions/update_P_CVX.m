function [P, f, convergence_WMMSE] = update_P_CVX(M, Pt, c, P_init, K, noise, w, g_MMSE)
% Update P using cvx toolbox
%  [P, f] = update_P(M, Pt, c, P, K, noise, w, g_MMSE)
%INPUTS:
%   M: number of transmit antennas at BS
%   Pt: transmit power at BS
%   c: equivalent channel ---- M x K
%   P: beamfomer matrix ---- M x K
%   K: number of users
%   noise: AWGN at users (watt)
%   w: weight of MSE
%   g_MMSE: MMSE receiver
%OUTPUTS:
%   P: updated P
%   f: optimal value of objective function
% Date: 21-Jul-24
% Author: Hessam Nezhadmohammad

cvx_begin quiet
    variable P(M,K) complex

    % constraint
    quad_over_lin(vec(P),1) <= Pt;

    % objective function
    f = objective_func(c, P, K, noise, g_MMSE, w);

    minimize(f);
cvx_end
% Check if P contains NaN values
if any(isnan(P(:)))
    convergence_WMMSE = false;
    P = P_init; % Optionally you can reset P to the inital value if NaN occurs
    f = NaN; % set f to NaN or any other appropriate value to indicate failure
else
    convergence_WMMSE = true;
    f = objective_func(c, P, K, noise, g_MMSE, w);
end

end

function [f] = objective_func(c, P, K, noise, g_MMSE, w)
    term_1 = w' * MSE(c, P, K, noise, g_MMSE);
    f = term_1;
end
