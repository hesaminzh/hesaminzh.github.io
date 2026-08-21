function [wsr] = WSR(c, P, K, noise, weights)
% The Weighted Sum Rate of the system
%  [wsr] = WSR(c, P, K, noise, weights)
% INPUTS:
%   c: equivalent channel ---- M x K
%   P: beamformer matrix at BS ---- M x K
%   K: number of users
%   noise: AWGN at users (watt)
%   weights: user weights
% OUTPUTS:
%   wsr: weighted sum rate
% M: number of BS antennas
% K: number of users
%Date: 21-Jul-24
%Author: Hessam Nezhadmohammad
gamma = SINR(c, P, K, noise);
R = log2(1+gamma);

wsr = weights' * R;
end
