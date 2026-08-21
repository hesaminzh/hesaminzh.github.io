function wsr = Alternative_optimization_Manifold_PIMRC(G,H,E,tolerance,SNR, noise, weights)
% ALTERNATIVE_OPTIMIZATION_MANIFOLD here is the summary of the function
% Detailed description is displayed here
H_d=G;
H_r=H;
G=E;
%Pt=db2pow(SNR);
Pt = SNR;
[Nr,K]=size(H_r);
%Nr=size(H_r,1);
Nt = size(G, 2);
noise=ones(K,1) .* noise;

% Group connected IRS model
num_groups = 1; % FC
Theta = zeros(Nr, Nr);
Ng = Nr / num_groups; 
for i = 1:num_groups
    Theta(Ng*(i-1)+1:Ng*i, Ng*(i-1)+1:Ng*i) = exp(1i.*rand(Ng,Ng).*2.*pi);   
end


%initial Theta as indentity matrix
%Theta=eye(Nr);
P_initial = randn(Nt, K) + 1i * randn(Nt, K);
P = P_initial / norm(P_initial, 'fro') * sqrt(Pt);

H=H_d+G'*Theta*H_r;
%P_p=(0.4* Pt)/K;
%P=H_d./vecnorm(H_d)*sqrt(P_p);
 
flag=1;
t_past=0;
count=0;
wsr = NaN;

max_retries = 5; % Maximum number of retries for Manopt_update_Theta
retry_count = 0;

while(flag)

    T_k=sum(square_abs(H'*P),2)+noise;
    alpha_k=T_k./(T_k-square_abs(diag(H'*P)))-1; % first auxiliary variable ---> SINR
    beta_k=sqrt(weights.*(1+alpha_k)).*diag(H'*P)./T_k; % second auxiliary variable

    res_k=weights.*log2(1+alpha_k)-weights.*alpha_k; % constant terms ---> const(alpha_k)

    [P,out]=FP_update_P_v2(H,Pt,alpha_k,beta_k,tolerance,res_k, noise, weights);
    if isnan(out)
        return;
    end

    [Theta,t]=Manopt_update_Theta(H_d,H_r,G,P,alpha_k,beta_k,Theta,noise,weights);
    % Check the validity of the optimized Theta
    if isnan(t) || isinf(t) || retry_count >= max_retries
        warning('Manopt_update_Theta failed to produce a valid Theta. Retrying...')
    for i = 1:num_groups
        Theta(Ng*(i-1)+1:Ng*i, Ng*(i-1)+1:Ng*i) = exp(1i.*rand(Ng,Ng).*2.*pi); 
    end
        retry_count = retry_count + 1;
        continue;
    end


    retry_count = 0; % Reset retry count after a successful update
    H=H_d+G'*Theta'*H_r;

    if abs(t-t_past)<=tolerance
        flag=0;
    else
        t_past=t;
        count=count+1;
        fprintf('count: %d\n', count);
    end
    if count>=2000
        return;
    end
end

%project Theta onto the symmetric unitary matrix
Theta=symuni(Theta); % FC
F=H_d+G'*Theta'*H_r;

%calculate the final sum_rate
%[~,sum_rate]=FP_algorithm_v2(F,SNR,tolerance, noise);
wsr = WSR(F, P, K, noise, weights);

end

