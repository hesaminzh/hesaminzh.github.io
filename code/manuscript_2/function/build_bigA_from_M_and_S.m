function A = build_bigA_from_M_and_S(M_real, S_sets, offsets)
% BUILD_BIGA_FROM_M_AND_S  Construct A ∈ R^{(2KM) x (sum_i |S_i|)}.
%   Paper identity:
%       || B M - Gamma ||_F^2   ≡   || A x - b ||_2^2,
%     with b = vec(Gamma^T), M = real-expanded matrix (M x 2K), and
%     x stacks the *upper-tri* (symmetric) entries {B_{i,ell}} with ell \in S_i.
%
%   A’s structure (per variable x_{i,ell} --> B_{i,ell} = B_{ell,i}):
%     - contributes a_ell to block-row i,
%     - contributes a_i to block-row ell,
%     - zeros elsewhere,
%     where a_r^T is the r-th row of M (i.e., M_real(r,:)), and each block-row
%     has height 2K (since M_real is M x 2K).
%
%   Inputs:
%     M_real    : M x (2K) real-expanded matrix used in the LS (paper’s M)
%     S_sets    : cell(M,1), S_sets{i} lists allowed upper-tri columns in row i
%     offsets   : 1 x M zero-based start positions for each row-block in x
%
%   Output:
%     A         : (M*2K) x (sum_i |S_i|) real matrix encoding B→BM mapping
%
%   Remarks:
%     - This function does not build b. In code: b = reshape(Gamma.', [], 1).
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)


    M    = size(M_real,1);            % rows of B and M_real
    K2   = size(M_real,2);            % equals 2K (real/imag stacking)
    use_sparse = (M*K2 > 2e4); % heuristic
    Ptot = offsets(end) + numel(S_sets{end});  % total columns in A

    % Optionally allocate sparse with a rough nnz estimate:
    %   Each off-diagonal variable adds TWO (2K)-blocks; diagonal only ONE.
    if use_sparse
        % upper bound: every column has two blocks
        nnz_est = 2 * K2 * Ptot;
        A = spalloc(M*K2, Ptot, nnz_est);
    else
        A = zeros(M*K2, Ptot);
    end

    % Precompute a_r = M_real(r,:).' (column vectors of length 2K)
    a = cell(1,M);
    for r = 1:M
        a{r} = M_real(r,:).';
    end

    % Fill A according to the Zheyu's paper (pp-ADMM and TIT25)
    for i = 1:M
        r1 = (i-1)*K2 + 1;
        r2 = i*K2;
        for j = 1:M
            Sj = S_sets{j};               % columns ell in row j
            base = offsets(j);            % zero-based start for row j’s block
            for q = 1:numel(Sj)
                col  = base + q;          % index column in A
                ell  = Sj(q);             % ell = S_j(q)

                if j == i
                    % block-row i gets a_ell
                    A(r1:r2, col) = A(r1:r2, col) + a{ell};
                elseif (ell == i)
                    % block-row i == ell gets a_j
                    A(r1:r2, col) = A(r1:r2, col) + a{j};
                end
                % else: zero contribution
            end
        end
    end
end


