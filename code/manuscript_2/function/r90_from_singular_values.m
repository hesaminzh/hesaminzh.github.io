function r90 = r90_from_singular_values(X)
%R90_FROM_SINGULAR_VALUES computes the effective rank of a matrix -> number of strong spatial modes
% NOTE: If you truly want to serve K = 4 or K = 6 users simultaneously with good separation, you typically need r90(H_eff) closer to 4 or 6 (and N >= K)
%
% Hesam Nezhadmohammad
    s = svd(X);
    E = cumsum(s.^2) / sum(s.^2);
    r90 = find(E >= 0.90, 1, 'first');
end

