function Theta = theta_from_B(B, Z0)
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)
    M  = size(B,1);
    Lm = eye(M) - 1i*Z0*B;   % (I - j Z0 B)
    Lp = eye(M) + 1i*Z0*B;   % (I + j Z0 B)
    % Correct Cayley transform: Theta = (I + j Z0 B)^(-1) * (I - j Z0 B)
    Theta = Lp \ Lm;          % (I + j Z0 B)^{-1} * (I - j Z0 B)
end


