function [y,gamma] = updateYGamma_ppadmm(U,G,W,sigma2)
%UPDATEYGAMMA_PPADMM  FP/QT auxiliary updates for sum-rate.
%   U (M x K): columns u_k = Θ^H h_k
%   G (M x N): RIS->BS
%   W (N x K): BS precoder
%   sigma2:    noise variance (scalar)
%
%   y_k   = S(k,k) / (sum_j |S(k,j)|^2 + sigma2)
%   gamma = |S(k,k)|^2 / (sum_j |S(k,j)|^2 + sigma2 - |S(k,k)|^2)
%
%   where S = (U^H G) W.
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)

    S      = (U' * G) * W;             % K x K
    rowpow = sum(abs(S).^2, 2);         % K x 1
    diagS  = diag(S);                   % K x 1
    den_y  = rowpow + sigma2;           % K x 1
    num2   = abs(diagS).^2;             % K x 1

    y      = diagS ./ den_y;            % complex K x 1
    den_inr= max(den_y - num2, eps);    % avoid div-by-zero
    gamma  = num2 ./ den_inr;           % real K x 1


    % Or

    % K = size(W,2); y=zeros(1,K); gamma=zeros(1,K);
    % GW  = G*W; 
    % UGW = U' * GW;             % K x K
    % for k = 1:K
    %     num      = UGW(k,k);                   % u_k^H G w_k
    %     denom_y  = sum(abs(UGW(k,:)).^2) + sigma2;
    %     y(k)     = num / denom_y;              % complex
    %     denom_g  = denom_y - abs(num)^2;       % sum_{j≠k}|⋅|^2 + σ^2
    %     gamma(k) = abs(num)^2 / max(denom_g, eps);
    % end
end

