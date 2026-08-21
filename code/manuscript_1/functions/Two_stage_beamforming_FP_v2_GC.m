function [sum_rate] = Two_stage_beamforming_FP_v2_GC(G,H,E,tolerance,Pt, noise, weights, arch_param)
%TWO_STAGE_BEAMFORMING_FP 


Theta=group_symuni(H*G'*E', arch_param);
F=G+E'*Theta'*H;

[Nt, K] = size(G);

P_initial = randn(Nt, K) + 1i * randn(Nt, K);
P = P_initial / norm(P_initial, 'fro') * sqrt(Pt);

T_k=sum(square_abs(F'*P),2)+noise;
alpha_k=T_k./(T_k-square_abs(diag(F'*P)))-1; % first auxiliary variable ---> SINR
beta_k=sqrt(weights.*(1+alpha_k)).*diag(F'*P)./T_k; % second auxiliary variable
res_k=weights.*log2(1+alpha_k)-weights.*alpha_k; % constant terms ---> const(alpha_k)
%[~, sum_rate, converged]=FP_algorithm_v2(F,Pt,tolerance, noise);
[~, sum_rate] = FP_update_P_v2(F, Pt, alpha_k, beta_k, tolerance, res_k, noise, weights);

end
