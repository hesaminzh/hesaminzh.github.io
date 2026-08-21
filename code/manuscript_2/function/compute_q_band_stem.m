function q = compute_q_band_stem(NI, NT, Nk_vec)
% Nk_vec: (K x 1) number of Rx antennas per user (Nk) or streams if you use dk instead
% in MU-MISO we have Nk_vec = ones(K,1)
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)

D = min(sum(Nk_vec), NT); % multiuser MIMO DoF
L = min(D, floor(NI/2));  % L as in the paper
q = 2*L - 1;              % band/stem width ensuring fully-conn. performance
end

