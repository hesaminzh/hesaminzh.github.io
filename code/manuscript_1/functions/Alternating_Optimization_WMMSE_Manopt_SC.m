function wsr = Alternating_Optimization_WMMSE_Manopt_SC(G,H,E,tolerance,SNR, noise, weights)
% ALTERNATING_OPTIMIZATION_MANOPT
% INPUTS:
%   G: BS-Users channel matrix          ---- Nt x K
%   H: RIS-Users channel matrix         ---- Nr x K
%   E: BS-RIS channel matrix            ---- Nr x Nt
%   SNR: maximum transmit SNR (watt)
%   noise: AWGN noise at the users (watt) 
%   weights: User weights               ---- K x 1
% OUTPUTS:
%   wsr: weighted sum-rate at the users ---- K x 1

H_d=G;
H_r=H;
G=E;
%Pt=db2pow(SNR);
Pt = SNR;
[Nr,K]=size(H_r);
%Nr=size(H_r,1);
Nt = size(G, 2);
noise=ones(K,1) .* noise;

%% Initialization

P_initial = randn(Nt, K) + 1i * randn(Nt, K);
P = P_initial / norm(P_initial, 'fro') * sqrt(Pt);
% Group connected IRS model
%num_groups = 1; % FC
%Theta = zeros(Nr, Nr);
%Ng = Nr / num_groups; 
%for i = 1:num_groups
%    Theta(Ng*(i-1)+1:Ng*i, Ng*(i-1)+1:Ng*i) = exp(1i.*rand(Ng,Ng).*2.*pi);   
%end

theta=exp(1i.*rand(Nr,1).*2.*pi);
Theta = diag(theta);



%initial Theta as indentity matrix
%Theta=eye(Nr);

C = H_d+G'*Theta*H_r; % equivalent channel
%P_p=(0.4* Pt)/K;
%P=H_d./vecnorm(H_d)*sqrt(P_p);

wsr = WSR(C, P, K, noise, weights);
%wsr_all = wsr; % for convergence behavior

%flag=1;
%t_past=0;
epsilon = tolerance;
diff = 10;
step = 0;

max_retries = 5; % Maximum number of retries for Manopt_update_Theta
retry_count = 0;

while diff > epsilon
    wsr_pre = wsr;
    step = step + 1;

    [P, ~, converged_WMMSE] = WMMSE_algorithm(C, P, epsilon, Pt, Nt, K, noise, weights);
    if converged_WMMSE == false
        wsr = NaN;
        return;
    end
    T_k=sum(square_abs(C'*P),2)+noise;
    alpha_k=T_k./(T_k-square_abs(diag(C'*P)))-1; % first auxiliary variable ---> SINR
    beta_k=sqrt(weights.*(1+alpha_k)).*diag(C'*P)./T_k; % second auxiliary variable

    %res_k=weights.*log2(1+alpha_k)-weights.*alpha_k; % constant terms ---> const(alpha_k)

    %[P,out]=FP_update_P_v2(C,Pt,alpha_k,beta_k,epsilon,res_k, noise, weights);
    %if isnan(out)
    %    wsr = NaN;
    %    return;
    %end

    [Theta,t]=Manopt_update_Theta(H_d,H_r,G,P,alpha_k,beta_k,Theta,noise,weights);
    % Check the validity of the optimized Theta
    if isnan(t) || isinf(t) || retry_count >= max_retries
        warning('Manopt_update_Theta failed to produce a valid Theta. Retrying...')
        theta=exp(1i.*rand(Nr,1).*2.*pi);
        Theta = diag(theta);
        retry_count = retry_count + 1;
        continue;
    end

    if isnan(t) || isinf(t)
        warning('Manopt_update_Theta didn''t work, continue to the next iteration...')
        wsr = NaN;
        return;
    end
    retry_count = 0; % Reset retry count after a successful update

    % Project Theta here
    % project Theta onto the symmetric unitary matrix
    %Theta=symuni(Theta); % FC
    Theta = group_symuni(Theta, 1); % SC
    C = H_d+G'*Theta'*H_r;

    %calculate the final sum_rate
    %[~,sum_rate]=FP_algorithm_v2(C,SNR,tolerance, noise);
    wsr = WSR(C, P, K, noise, weights);
    diff = abs(wsr - wsr_pre);
    %wsr_all = [wsr_all wsr]; % for convergence behavior
    %disp(['step--' num2str(step) ', WSR--' num2str(wsr)]);


    %C=H_d+G'*Theta'*H_r;
    %Theta=symuni(Theta); % FC
    %C=H_d+G'*Theta'*H_r;
    %wsr = WSR(C, P, K, noise, weights);
    %diff = abs(wsr - wsr_pre);
    %wsr_all = [wsr_all wsr];
    %disp(['step--' num2str(step) ', WSR--' num2str(wsr)]);


    %if abs(t-t_past)<=tolerance
    %    flag=0;
    %else
    %    t_past=t;
    %    step=step+1;
    %    %fprintf('AO step: %d\n', step);
    %end

    if step>=1500
        wsr = NaN;
        return;
    end
end

%project Theta onto the symmetric unitary matrix
%Theta=symuni(Theta); % FC
%F=H_d+G'*Theta'*H_r;

%calculate the final sum_rate
%[~,sum_rate]=FP_algorithm_v2(F,SNR,tolerance, noise);
%wsr = WSR(F, P, K, noise, weights);

end

