function [e] = MSE(c, P, K, noise, g)
% The Mean Square Error at the MMSE receiver
%  [e] = MSE(c, P, K, noise, g)
%INPUTS:
%   c: equivalent channel ---- M x K
%   P: beamformer matrix at BS ---- M x K
%   K: number of users
%   noise: AWGN noise at users (watt)
%   g: MMSE receiver
%OUTPUTS:
%   e: mean square error
% Date: 21-Jul-24
% Author: Hessam Nezhadmohammad

% MSE at the MMSE receiver
e = [];
for k = 1:K
    pk = P(:,k);
    ck = c(:,k);
    gk = g(k);

    ek = pow_abs(gk, 2) * ( sum(pow_abs(ck' * P, 2)) + noise(k)) - 2*real(gk * ck' * pk) + 1;
    e = [e;ek]; % ---- K x 1
end
end
