function [B] = BS_UE(Nx, Ny, K, L, d)
    % BS-UE channel matrix
    M = Nx * Ny;
    
    lambda_gain = 9.60; % \varrho_t = 9.82 dBi
    a = 61.4; % \varrho_a or PL0 in O2I
    b = 2.92; % \varrho_b for NLoS
    sigma = 8.7; % 8.7 dB for NLoS
    kappa = a + 10 * b * log10(d) + normrnd(0, sigma, 1); % PL(r)
    % channel attenuation
    var = 0.01 * ones(1, L) .* (sqrt(1/2 * 10 .^ (-0.1 * kappa))) * (randn(1) + 1j * randn(1)); % size: 1 x L
    % L --> NLoS paths

    % ULA array steering vector
    %theta = rand(1, L) - 0.5;
    %Arr_BS = exp(-2 * 1i * pi * theta' * ((0:M-1) - (M-1)/2))'; % size: M x L

    % Path gains
    beta_BL = lambda_gain * var;

    % Preallocate B
    B = zeros(M, K);
    for k = 1:K
        % Random distribution of angles
        angU = rand(1, L) - 0.5;

        % Array steering vector
        arrU = exp(-2 * 1i * pi * angU' * ((0:M-1) - (M-1)/2))'; % size: M x L

        % Construct channel matrix for each user
        B(:, k) = arrU * beta_BL';
    end
end
