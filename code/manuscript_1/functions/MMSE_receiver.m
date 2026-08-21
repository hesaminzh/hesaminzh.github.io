function [g] = MMSE_receiver(c, P, K, noise)
% The MMSE receiver
%  [g] = MMSE_receiver(c, P, K, noise)
%INPUTS:
% c: equivalent channel ---- M x K
% P: beafomer matrix at BS ---- M x K
% K: number of users
% noise: AWGN at users (watt)
%OUTPUTS:
% g: MMSE receiver   
% Date: 21-Jul-24
% Author: Hessam Nezhadmohammad
g = zeros(K, 1);
for k = 1:K
    pk = P(:,k);
    ck = c(:,k);
    g(k) = (pk' * ck) / ( sum(pow_abs(ck'*P, 2)) + noise(k) );
end
end
