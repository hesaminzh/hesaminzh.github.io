function U_new = update_U_ppadmm(B, G, H, W, y, gamma, Lambda, rho, Z0)
%UPDATE_U_PPADMM  Update the auxiliary variable U in the pp‑ADMM algorithm.
%   This function implements the per‑user update for the columns of U as
%   described in the closed-form update of U in the paper.  For each user k, we solve
%
%   u_k^{t+1} = [ (1+gamma_k)|y_k|^2 G W W^H G^H + (rho/2)(I + Z0^2 B^2) ]^{-1}
%                 * [ (1+gamma_k) conj(y_k) G w_k
%                     - 0.5 * (I + i Z0 B) lambda_k^t
%                     + 0.5 * rho * (I + i Z0 B)^2 h_k ]
%
%   where h_k and lambda_k denote the k‑th columns of H and Lambda, respectively.
%
% In the original pp-ADMM algorithm (see Eq. (11e) and closed-form update of U in the
% paper), the update of U is not a simple multiplication of the
% scattering matrix Theta; instead, each column u_k solves an
% unconstrained convex quadratic problem that incorporates the
% current precoder W, the FP variables y and gamma, the proximal
% parameter rho, and the Lagrange multiplier Lambda. A simplified
% update U = Theta^H H fails to couple the BD-RIS susceptance B
% with the beamforming and FP variables, causing B to remain
% zero. Here we implement the full update derived in the
% paper
%
% To avoid redundant computations, we precompute the common terms
% outside the user loop.
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)


% Precompute G*W once: GW = G * W (M x K)
%

% Dimensions
[M, ~] = size(G);
K = size(H, 2);

% Precompute common matrices
% G*W (N×K) used to build G W W^H G^H
GW = G * W;
S = GW*GW'; % M x M Hermitian
I = eye(M);
% Real matrix M1 = I + Z0^2 B^2
M1 = I + (Z0^2) * (B*B);
% Complex matrices for RHS
ThetaFactor  = I + 1j * Z0 * B;      % I + i Z0 B
ThetaFactor2 = ThetaFactor * ThetaFactor; % (I + i Z0 B)^2

% Initialise U
U_new = zeros(M, K);

% Loop over users
for k = 1:K
    % Left‑hand side matrix
    alpha = (1 + gamma(k)) * abs(y(k))^2;
    Rk = alpha*S + (rho/2) * M1;

    % Right‑hand side vector
    term1 = (1 + gamma(k)) * conj(y(k)) * GW(:,k);
    term2 = -0.5 * (ThetaFactor * Lambda(:,k));
    term3 = 0.5 * rho * (ThetaFactor2 * H(:,k));
    fk = term1 + term2 + term3;

    % Solve for u_k
    U_new(:,k) = Rk \ fk;
end

end
