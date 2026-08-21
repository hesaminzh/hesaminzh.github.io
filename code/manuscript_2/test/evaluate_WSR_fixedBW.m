function R = evaluate_WSR_fixedBW(panel_positions, B, W, N, Mx, My, I_panels, user_loc, BS_loc, sigma2, Z0, sG, sH, rng_state);

%EVALUATE_WSR_FIXEDBW
%    Evaluate WSR for a given panel placement using fixed (Theta, W),
%    by regenerating channels and calling compute_sum_rate_theta(...).
%
% Inputs:
%   panel_positions : I_panels x 2
%   B               : M x M
%   W               : N x K
%   N, Mx, My       : system dimensions
%   I_panels        : number of panels
%   user_loc        : fixed user locations (common random numbers)
%   BS_loc          : BS coordinates
%   sigma2          : raw noise power
%   sG, sH          : normalization factors from current AO iteration
%
% Output:
%   R               : scalar WSR (bps/Hz)
%
% Hesam Nezhadmohammad (last revised: 20/2/2026)

    rng(rng_state); % force common random numbers
    Theta = theta_from_B(B, Z0);
    G_cell = cell(I_panels, 1);
    H_cell = cell(I_panels, 1);

    for i = 1:I_panels
        panel_loc = panel_positions(i,:);
        [G_cell{i}, H_cell{i}] = generate_SV_channel(N, Mx, My, user_loc, BS_loc, panel_loc);
    end

    G_agg = vertcat(G_cell{:});
    H_agg = vertcat(H_cell{:});

    % Use the SAME normalization factors as the current point (important)
    Gw = sG * G_agg;
    Hw = sH * H_agg;
    sigma2w = (sG * sH)^2 * sigma2;

    % Evaluate using your existing function
    R = compute_sum_rate_theta(Theta, Hw, Gw, W, sigma2w);
end
