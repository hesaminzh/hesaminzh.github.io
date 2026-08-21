function R = compute_sum_rate(U, W, G, sigma2)
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)
    % U (M x K), W (N x K), G (M x N); RIS-only path
    H_eff = U' * G;           % K x N, rows: u_k^H G
    S     = H_eff * W;        % K x K, entries: u_k^H G w_j
    sig   = abs(diag(S)).^2;  % |u_k^H G w_k|^2
    tot   = sum(abs(S).^2,2); % sum_j |u_k^H G w_j|^2
    intf  = tot - sig;        % sum_{j!=k} ...
    R     = sum(log2(1 + sig ./ (intf + sigma2)));
end

