function [gamma] = SINR(c, P, K, noise)
%  [gamma] = SINR(c, P, K, noise)
% INPUTS:
%   c: equivalent channel ---- M x K
%   P: beamformer at BS ---- M x K
%   K: number of users
%   noise: AWGN at user (watt) ---- K x 1
% OUTPUTS:
%   gamma: SINR
% Date: 21-Jul-24
% Author: Hessam Nezhadmohammad

gamma = zeros(K, 1);
for k = 1:K
    pk = P(:,k);
    P_inter = P; P_inter(:,k) = [];
    ck = c(:,k);

    % SINR
    gamma(k) = pow_abs(ck' * pk, 2) / ( noise(k) + sum(pow_abs(ck' * P_inter, 2)) );
end
end
