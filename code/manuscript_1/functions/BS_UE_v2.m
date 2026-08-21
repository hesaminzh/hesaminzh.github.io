function [B] = BS_UE_v2(Mx, My, K, L, d)
% BS-UE channel matrix
M = Mx * My;

lambda_gain = 9.60; % \varrho_t = 9.82 dBi
a = 61.4; % \varrho_a or PL0 in O2I
b = 2.92; % \varrho_b for NLoS
sigma = 8.7; % 8.7 dB for NLoS
kappa = a + 10 * b * log10(d) + normrnd(0, sigma, [K, 1]); % PL(r) for each user, size: K x 1

% Channel attenuation for each user and each path
var = 0.01 * ones(K, L) .* (sqrt(1/2 * 10.^(-0.1 * kappa))) .* (randn(K, L) + 1j * randn(K, L)); % size: K x L
% L --> NLoS paths

% ULA array steering vector
%theta = rand(K, L) - 0.5; % Random angles for each user and each path
%Arr_BS = exp(-2 * 1i * pi * theta .* ((0:M-1) - (M-1)/2)'); % size: M x L x K

% Path gains
beta_BL = lambda_gain * var; % size: K x L

% Preallocate B
B = zeros(M, K);
for k = 1:K
    % Random distribution of angles for user k
    angU = rand(1, L) - 0.5;

    % Array steering vector for user k
    arrU = exp(-2 * 1i * pi * angU' * ((0:M-1) - (M-1)/2))'; % size: M x L

    % Construct channel matrix for each user
    B(:, k) = arrU * beta_BL(k, :)';
end
end
