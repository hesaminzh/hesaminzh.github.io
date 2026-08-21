function [B, x_new] = update_B_ppadmm(U, H, Lambda, S_list, offsets, x_prev, rho, xi, Z0)
%UPDATE_B_PPADMM  pp-ADMM B-update (paper-exact, realified LS with proximal term).
%   Inputs:
%       U,H        : M x K complex
%       Lambda     : M x K complex (dual)
%       S_list     : {S_i} upper-triangular index sets per row i
%       offsets    : 0-based cumulative offsets for x packing (length M)
%       x_prev     : previous x (sum_i |S_i| x 1 real)
%       rho, xi    : ADMM penalty and proximal weight (xi > 0)
%       Z0         : real scalar (>0)
%   Outputs:
%       B          : M x M real symmetric, obeying architecture
%       x_new      : updated parameter vector
%       stats      : diagnostic residuals, method flag
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)


    M = size(U,1);
    K = size(U,2);
    Ptot = offsets(end) + numel(S_list{end});

    % Conversion to real domain: Mmat = Re/Im{i Z0 (U+H)}, Gamma = Re/Im{U - H + Lambda/rho}
    MU  = [real(1i*Z0*(U + H)),  imag(1i*Z0*(U + H))];                      % M x 2K
    GammaU  = [real((U - H) + (Lambda/rho)), imag((U - H) + (Lambda/rho))];                   % M x 2K

    % b = vec(Gamma^T)
    b = reshape(GammaU.', [], 1);                             % (2MK) x 1

    % Architecture-encoded regression matrix
    A = build_bigA_from_M_and_S(MU, S_list, offsets);     % (2MK) x Ptot

    % Normal equations with proximal center
    At = A.';                                                % Ptot x (2MK)
    AtA = At*A;                                              % Ptot x Ptot
    Atb = At*b;                                              % Ptot x 1

    use_dual = (Ptot >= 2*M*K);

    if use_dual
        % Dual / Woodbury form
        t    = Atb + xi * x_prev;                            % Ptot x 1
        G    = (A*A.' + xi * eye(2*M*K)) \ (A * t);          % (2MK) x 1
        x_new = (t - At * G) / xi;                           % Ptot x 1
        method = 'dual';


    else
        % Primal form
        x_new = (AtA + xi * eye(Ptot)) \ (Atb + xi * x_prev);
        % x_new = (rho*AtA + xi * eye(Ptot)) \ (rho*Atb + xi * x_prev); % if xi is not passed as xi/rho in the main script
        method = 'primal';
    end

    % Map x -> B (upper triangle per S_list, then mirror)
    B = zeros(M);
    for i = 1:M
        Si = S_list{i};
        for q = 1:numel(Si) % q: local index inside row i's block
            j   = Si(q); % actual column index
            idx = offsets(i) + q; % position in x_new
            B(i,j) = x_new(idx);
            B(j,i) = x_new(idx); % symmetry
        end
    end

    % % Diagnostics
    % stationarity = norm((AtA + xi*eye(Ptot))*x_new - (Atb + xi*x_prev)) / max(1,norm(Atb + xi*x_prev));
    % primal_res   = norm(reshape((B*MU - GammaU).', [], 1)) / max(1,norm(b));
    % stats = struct('method', method, 'stationarity', stationarity, 'primal_residual', primal_res);
end

