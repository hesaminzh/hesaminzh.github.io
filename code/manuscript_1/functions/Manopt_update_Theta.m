function [Theta,out]=Manopt_update_Theta(H_d,H_r,G,P,alpha_k,beta_k,Theta_init,noise,weights)
[Nr,K]=size(H_r);

[A_k,C_k]=deal(zeros(Nr,Nr));

e_k=G*P;
t_k=H_d'*P;

eta_k=weights.*log2((1+alpha_k))-weights.*alpha_k-square_abs(beta_k).*(sum(square_abs(t_k),2)+noise); % constants w.r.t Theta
lambda_k=sqrt(weights.*(1+alpha_k)).*conj(beta_k).*diag(t_k);

for k=1:K
    j_user=1:K;

    z_k=conj(t_k(k,j_user)).*e_k*repmat(H_r(:,k)',K,1); % repmat(H_r(:,k)',K,1) is KxN matrix of h_k' (ris-user) channels

    A_k=A_k+sqrt(weights(k).*(1+alpha_k(k)))*conj(beta_k(k))*e_k(:,k)*H_r(:,k)'-square_abs(beta_k(k))*z_k;
    
    C_k=C_k+square_abs(beta_k(k))*H_r(:,k)*H_r(:,k)';
   
end
B_k=e_k*e_k';

% -2*real(trace(Theta_init*A_k))+real(trace(Theta_init*B_k*Theta_init'*C_k))

manifold=stiefelcomplexfactory(Nr, Nr);
problem.M=manifold;

problem.cost=@(Theta) -2*real(trace(Theta*A_k))+real(trace(Theta*B_k*Theta'*C_k)); % it's the remaining part by removing the constants w.r.t Theta
problem.egrad=@(Theta) 2*C_k*Theta*B_k-2*A_k';

warning('off', 'manopt:getHessian:approx');
options = struct();
options.verbosity = 0;

[Theta, xcost] = conjugategradient(problem,Theta_init,options);

out=-xcost+sum(eta_k+2*real(lambda_k)); % Final WSR is (constants w.r.t Theta) + (remaining part by removing the constants w.r.t Theta )

end

