function S_sets = arch_sparsity_sets(M, architecture, group_size)
% ARCH_SPARSITY_SETS  Build {S_i} index sets encoding the BD-RIS architecture.
%   Paper context (pp-ADMM, B-update linearization):
%     We vectorize the symmetric scattering matrix B by stacking its
%     upper-triangular independent variables row-by-row. For each row i,
%     S_i \in {i,i+1,...,M} lists the allowed nonzeros.
%
%   Inputs:
%     M            : number of RIS ports (size of B is M x M)
%     architecture : 'fully' or 'group'
%     group_size   : G_s (only used for 'group')
%
%   Output:
%     S_sets : cell(M,1), where S_sets{i} is a row-vector of column indices
%              (\ge i) permitted in row i of the *upper triangle* of B.
%
%   Notes:
%     - 'fully'   -> S_i = i:M
%     - 'group'   -> S_i limited to the current group’s end (in-group upper-tri)
%
%   See paper: pp-ADMM B-update -> LS form ||A x - b||_2^2 with x = stacked B_ij
%   over j \in S_i (upper-tri only).
% Hesam Nezhadmohammad (last revised: 19/10/2025)

    S_sets = cell(M,1);
    Adj = zeros(M);  % symmetric 0-1 adjacency, zero diagonal

    switch lower(architecture)
        case 'fully'

            % for i = 1:M
            %     S_sets{i} = i:M;  % all upper-tri positions incl. diagonal
            % end

            % Complete graph (fully-connected).
            Adj = ones(M) - eye(M);

        case 'group'

            % if nargin < 3 || isempty(group_size)
            %     error('arch_sparsity_sets: group_size must be provided for ''group''.');
            % end
            % Gs = group_size;
            % for i = 1:M
            %     g         = floor((i-1)/Gs);
            %     start_idx = g*Gs + 1;
            %     end_idx   = min((g+1)*Gs, M);
            %     S_sets{i} = i:end_idx;  % in-group upper-tri
            % end

            for i = 1:M
                g = floor((i-1)/group_size);
                idx = (g*group_size+1):min((g+1)*group_size, M);
                Adj(idx,idx) = 1; Adj(idx,idx) = Adj(idx,idx) - eye(numel(idx));
            end



        case 'single'   % single-connected (diagonal only)

            % for i = 1:M
            %     S_sets{i} = i;   % only diagonal
            % end


        case 'tridiagonal'
            % % % Tree-connected (path graph / tridiagonal): |i-j|>1 are zero off-diagonals.
            % % % Upper triangle keeps {i} and {i+1} if exists.
            % for i = 1:M
            %     if i < M
            %         S_sets{i} = [i, i+1];
            %     else
            %         S_sets{i} = M;
            %     end
            % end

            for i = 1:M-1
                Adj(i,i+1)=1; Adj(i+1,i)=1;
            end

        otherwise
            error('arch_sparsity_sets: unknown architecture "%s".', architecture);
    end
    Adj = (Adj>0); Adj = triu(Adj,1) + triu(Adj,1).';  % symmetric, 0/1, no diag
    for i = 1:M
        upper_neighbors = find(Adj(i, i+1:M)) + i;   % neighbors with index > i
        S_sets{i} = [i, upper_neighbors];          % include diagonal by construction
    end
end

