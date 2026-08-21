function [sum_rate, converged] = Two_stage_beamforming_FP_v2(G,H,E,tolerance,SNR, noise)
%TWO_STAGE_BEAMFORMING_FP 


Theta=symuni(H*G'*E');
F=G+E'*Theta'*H;

[~, sum_rate, converged]=FP_algorithm_v2(F,SNR,tolerance, noise);

if ~converged
    sum_rate = NaN; % or some other indicator for non-convergence
end

