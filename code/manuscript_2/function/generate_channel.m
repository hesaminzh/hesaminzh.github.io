function H = generate_channel(Rx,Tx,amp,kappa)
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)
    LoS = generate_LoS(Rx,Tx);
    NLoS = 1/sqrt(2) .* (randn(Rx,Tx)+1i*randn(Rx,Tx));
    H = amp.*(sqrt(kappa/(kappa+1))*LoS + sqrt(1/(kappa+1))*NLoS);
    % H = (sqrt(kappa/(kappa+1)))*LoS + amp.*(sqrt(1/(kappa+1))*NLoS);
end

function LoS = generate_LoS(Rx, Tx)
%GENERATE_LOS  LoS BS->RIS channel (Rx x Tx) using ULA steering vectors.

    % Random AoD and AoA in [0, 2*pi)
    theta_t = 2*pi*rand;   % AoD from BS to RIS
    theta_r = 2*pi*rand;   % AoA at RIS from BS

    % Indices start at 0
    n_bs  = (0:Tx-1).';     % Tx x 1
    n_ris = (0:Rx-1).';     % Rx x 1

    % ULA steering vectors (half-wavelength spacing)
    a_t = exp(1i * pi * sin(theta_t) * n_bs);   % Tx x 1
    a_r = exp(1i * pi * sin(theta_r) * n_ris);  % Rx x 1

    % % Normalize
    % a_t = a_t / sqrt(Tx);
    % a_r = a_r / sqrt(Rx);

    % Complex path gain (set as needed, here unit gain)
    alpha = 1;

    % LoS channel: Rx x Tx (RIS rows, BS columns)
    LoS = alpha * (a_r * a_t'); 
end
