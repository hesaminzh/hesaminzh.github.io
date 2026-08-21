function [avgWSR, outageProb] = simulateBlockageDensity(eta, scheme, M, N, K, P_tx, noise, weights, tolerance, Ng, tau, d_BU, d_RU, H_BS_RIS_all, H_RIS_user_all, H_BS_user_all, nMonte)
% eta: blockage density
wsr_all = NaN(nMonte,1);

%% use parfor instead of for!
%
%channels_const = parallel.pool.Constant(channels); % send the entire channels to each worker only once at the starts
parfor iMonte = 1:nMonte
    % --- 1) Generate random channels (Saleh-Valenzuela or your own) ---
    % Here we use your BS-RIS, RIS-user, BS-user generators (load channels.mat or call BS_IRS, FF_User, BS_UE_v2)
    %[H_BS_RIS, H_RIS_user, H_BS_user] = generateChannels(M, N, K);
    %% BS-IRS channel matrix
    H_BS_RIS   = H_BS_RIS_all{iMonte,1};

    %% IRS-user channel matrix
    H_RIS_user = H_RIS_user_all{iMonte,1};

    %% BS-user channel matrix
    H_BS_user  = H_BS_user_all{iMonte,1};

    %%%--- distance-dependent blockage probabilities for each user k
    p_BU = 1 - exp(-eta * d_BU);
    p_RU = 1 - exp(-eta * d_RU);
    % --- 2) Sample blockages for each user k ---
    I_d = rand(K,1) > p_BU;   % BS–user link available?
    I_r = rand(K,1) > p_RU;   % RIS–user link available?
    for k = 1:K
        if ~I_d(k), H_BS_user(:,k) = 0; end
        if ~I_r(k), H_RIS_user(:,k) = 0; end
    end

    %%% Equivalent way to remove blocked channels, no loop!
    %% 1 means “blocked,” 0 means “unblocked”
    %isBlocked = rand(1,K) < p_blk;
    %
    %% zero blocked channels in one line, no loop
    %H_BS_user(:, isBlocked) = 0;
    %H_RIS_user(:, isBlocked) = 0;

    % --- 4) Compute WSR for the chosen scheme ---
    tempWSR = NaN;
    switch scheme
        case 'NoRIS'
            % MRT on (possibly blocked) direct channel only
            %W = computeMRT(H_BS_user, P_tx);
            %Theta = zeros(N);
            %tempWSR = computeWSR(H_BS_user, H_RIS_user, H_BS_RIS, Theta, W, noise, weights);
            W_initial = randn(M, K) + 1i * randn(M, K);
            W = W_initial / norm(W_initial, 'fro') * sqrt(P_tx);
            [~, tempWSR] = WMMSE_algorithm_closed_form(H_BS_user, W, tolerance, P_tx, K, noise, weights);

        case 'RandPhase'
            % Random-phase RIS (SC)
            Theta = diag(exp(1i*2*pi*rand(N,1)));
            %W     = update_W_WMMSE_closed_form(H_BS_RIS, H_RIS_user, H_BS_user, Theta, ...
            %        weights, P_tx, noise);
            %tempWSR = computeWSR(H_BS_user, H_RIS_user, H_BS_RIS, Theta, W, noise, weights);
            W_initial = randn(M, K) + 1i * randn(M, K);
            W = W_initial / norm(W_initial, 'fro') * sqrt(P_tx);
            H_eff = H_BS_user' + H_RIS_user' * Theta * H_BS_RIS;  % K x M effective channel
            C = H_eff';
            [~, tempWSR] = WMMSE_algorithm_closed_form(C, W, tolerance, P_tx, K, noise, weights);

        case {'SC','GC','FC'}
            % AO WMMSE + Manopt RIS
            if strcmpi(scheme,'GC')
                arch_param = Ng;
                [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_GC(H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights, arch_param);
            elseif strcmpi(scheme,'FC')
                arch_param = [];
                [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights);
            elseif strcmpi(scheme,'SC')
                [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_SC(H_BS_user, H_RIS_user, H_BS_RIS, tolerance, P_tx, noise, weights);
            end

        otherwise
            error('Unknown scheme %s', scheme);
    end
    wsr_all(iMonte) = tempWSR;
end

% --- 5) Compute statistics ---
avgWSR     = mean(wsr_all, 'omitnan');
outageProb = sum(wsr_all < tau, 'omitnan') / (sum(~isnan(wsr_all)));
end
