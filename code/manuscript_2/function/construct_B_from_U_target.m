function B_arch = construct_B_from_U_target(U_target, H_used, S_sets, Z0)
%CONSTRUCT_B_FROM_U_TARGET  Theorem 2 construction of B for a given U.
%   Given:
%       U_target  : M x K complex matrix (effective RIS-side channel we want)
%       H_used    : M x K complex RIS->users channel (must be consistent with U_target)
%       S_sets    : cell array {S_i} from arch_sparsity_sets_all (architecture G)
%       Z0        : reference impedance
%
%   Returns:
%       B_arch    : M x M real symmetric susceptance matrix, with support
%                   constrained by S_sets, such that
%                   (I - j Z0 B_arch) U_target ≈ (I + j Z0 B_arch) H_used
%                   (up to numerical precision).
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)

    [M, K] = size(U_target);

    % === 1) Build MU and GammaU as in the TIT paper (single-antenna users) ===
    % MU = [Re(j Z0 (H + U)), Im(j Z0 (H + U))]
    % GammaU = [Re(U - H), Im(U - H)]
    MU      = [ real(1i * Z0 * (H_used + U_target)), ...
                imag(1i * Z0 * (H_used + U_target)) ];   % M x (2K)
    GammaU  = [ real(U_target - H_used), ...
                imag(U_target - H_used) ];               % M x (2K)

    % Vectorize right-hand side: b = vec(GammaU.')
    % Size: (2K*M) x 1
    b = reshape(GammaU.', [], 1);

    % === 2) Prepare upper-tri index sets and count unknowns ===
    S_upper = cell(M,1);
    num_vars = 0;
    for i = 1:M
        Si = S_sets{i};
        Si = Si(Si >= i);       % keep only upper-tri (j >= i)
        S_upper{i} = Si;
        num_vars = num_vars + numel(Si);
    end

    % === 3) Build the big matrix A such that vec(B MU) = A x ===
    % We want BMU = GammaU  => vec(BMU) = vec(GammaU) = b.
    % Dimensions:
    %   MU: M x (2K)
    %   B:  M x M (unknown, but structured)
    %   => BMU: M x (2K)
    %   => vec(BMU): (2K*M) x 1
    %
    % x stacks all independent upper-tri elements:
    %   x = [ B(1,S_upper{1}), B(2,S_upper{2}), ..., B(M,S_upper{M}) ]^T

    rowsA = 2*K*M;
    A = zeros(rowsA, num_vars);

    col = 0;
    for i = 1:M
        Si = S_upper{i};
        for t = 1:numel(Si)
            j = Si(t);          % j >= i
            col = col + 1;

            % This unknown corresponds to b_ij = B(i,j) with j >= i.
            % Because B is symmetric, B(j,i) = B(i,j) = x(col).
            % Contribution to row i of BMU is b_ij * MU(j,:).
            row_i = (i-1)*2*K + 1 : (i*2*K);
            A(row_i, col) = MU(j, :).';

            % If off-diagonal, also contributes to row j as b_ij * MU(i,:).
            if j ~= i
                row_j = (j-1)*2*K + 1 : (j*2*K);
                A(row_j, col) = MU(i, :).';
            end
        end
    end

    % === 4) Solve least-squares system A x = b ===
    % Theorem 2 says there exists (at least one) exact solution for generic U.
    % Numerically we just do LS; residual should be ~ machine precision.
    x = A \ b;
    % x = pinv(A)*b;

    % === 5) Map x back to a symmetric B_arch with the desired sparsity ===
    B_arch = zeros(M);
    col = 0;
    for i = 1:M
        Si = S_upper{i};
        for t = 1:numel(Si)
            j = Si(t);
            col = col + 1;
            val = x(col);
            B_arch(i,j) = val;
            B_arch(j,i) = val;  % enforce symmetry
        end
    end
end

