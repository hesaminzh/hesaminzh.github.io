function wsr = RZF_update_P(G, SNR, noise, weights)
%CLASSICAL_RZF 

K=size(G,2);
%noise=ones(K,1) .* dBm2pow(noise);
noise = ones(K,1) .* noise;
F=G;

%Pt=dBm2pow(SNR);
Pt = SNR;
D=F/(F'*F+1/Pt*eye(K));
D=D./vecnorm(D);
P=diag(sqrt(Pt/K)*ones(K,1));
W=D*P;


T_k=sum(square_abs(F'*W),2)+noise;
rate_p=log2(T_k./(T_k-square_abs(diag(F'*W))));
wsr_p = weights .* rate_p;
%sum_rate=sum(rate_p);
wsr = sum(wsr_p);
end

