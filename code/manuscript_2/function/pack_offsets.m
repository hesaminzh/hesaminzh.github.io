function [offsets, num_vars] = pack_offsets(S_sets)
% PACK_OFFSETS  Compute column offsets for stacking upper-tri variables into x.
%   Paper mapping:
%     x = [ B_{1,S_1}, B_{2,S_2}, ..., B_{M,S_M} ]^T  \in R^{sum_i |S_i|}
%     This function returns zero-based column starts for each row-block.
%
%   Input:
%     S_sets : cell(M,1) with S_sets{i} = indices in row i (upper-tri)
%
%   Outputs:
%     offsets : 1 x M, zero-based start position of row i's block in x
%     num_vars: total number of stacked variables = sum_i |S_i|
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)

    M = numel(S_sets);
    offsets  = zeros(1,M);
    current   = 0;
    for i = 1:M
        offsets(i)  = current;
        current     = current + numel(S_sets{i});
    end
    num_vars = current;
end

