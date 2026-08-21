function [avgWSR, outageProb] = simulateBlockagePoint(p_blk, scheme, params, channels, nMonte)
% Runs nMonte trials at blockage probability p_blk for one scheme.
wsr_all = NaN(nMonte,1);

%% use parfor instead of for!
%
channels_const = parallel.pool.Constant(channels); % send the entire channels to each worker only once at the starts
parfor iMonte = 1:nMonte
    % --- 1) Generate random channels (Saleh-Valenzuela or your own) ---
    % Here we use your BS-RIS, RIS-user, BS-user generators (load channels.mat or call BS_IRS, FF_User, BS_UE_v2)
    %[H_BS_RIS, H_RIS_user, H_BS_user] = generateChannels(params.M, params.N, params.K);
    %% BS-IRS channel matrix
    H_BS_RIS   = channels_const.Value.H_BS_RIS_all{iMonte,1};

    %% IRS-user channel matrix
    H_RIS_user = channels_const.Value.H_RIS_user_all{iMonte,1};

    %% BS-user channel matrix
    H_BS_user  = channels_const.Value.H_BS_user_all{iMonte,1};

    % --- 2) Sample blockages for each user k ---
    I_d = rand(params.K,1) > p_blk;   % BS–user link available?
    I_r = rand(params.K,1) > p_blk;   % RIS–user link available?
    for k = 1:params.K
        if ~I_d(k), H_BS_user(:,k) = 0; end
        if ~I_r(k), H_RIS_user(:,k) = 0; end
    end

    %%% Equivalent way to remove blocked channels, no loop!
    %% 1 means “blocked,” 0 means “unblocked”
    %isBlocked = rand(1,params.K) < p_blk;
    %
    %% zero blocked channels in one line, no loop
    %H_BS_user(:, isBlocked) = 0;
    %H_RIS_user(:, isBlocked) = 0;
    % --- 3) Pack into channels struct ---
    %channels.H_BS_RIS  = H_BS_RIS;
    %channels.H_RIS_user= H_RIS_user;
    %channels.H_BS_user = H_BS_user;
    ch = struct();
    ch.H_BS_RIS = H_BS_RIS;
    ch.H_RIS_user = H_RIS_user;
    ch.H_BS_user = H_BS_user;

    % --- 4) Compute WSR for the chosen scheme ---
    tempWSR = NaN;
    switch scheme
        case 'NoRIS'
            % MRT on (possibly blocked) direct channel only
            %W = computeMRT(ch.H_BS_user, params.P_tx);
            %Theta = zeros(params.N);
            %tempWSR = computeWSR(ch.H_BS_user, ch.H_RIS_user, ch.H_BS_RIS, Theta, W, params.noise, params.weights);
            W_initial = randn(params.M, params.K) + 1i * randn(params.M, params.K);
            W = W_initial / norm(W_initial, 'fro') * sqrt(params.P_tx);
            [~, tempWSR] = WMMSE_algorithm_closed_form(ch.H_BS_user, W, params.tolerance, params.P_tx, params.K, params.noise, params.weights);

        case 'RandPhase'
            % Random-phase RIS (SC)
            Theta = diag(exp(1i*2*pi*rand(params.N,1)));
            %W     = update_W_WMMSE_closed_form(ch.H_BS_RIS, ch.H_RIS_user, ch.H_BS_user, Theta, ...
            %        params.weights, params.P_tx, params.noise);
            %tempWSR = computeWSR(ch.H_BS_user, ch.H_RIS_user, ch.H_BS_RIS, Theta, W, params.noise, params.weights);
            W_initial = randn(params.M, params.K) + 1i * randn(params.M, params.K);
            W = W_initial / norm(W_initial, 'fro') * sqrt(params.P_tx);
            H_eff = ch.H_BS_user' + ch.H_RIS_user' * Theta * ch.H_BS_RIS;  % K x M effective channel
            C = H_eff';
            [~, tempWSR] = WMMSE_algorithm_closed_form(C, W, params.tolerance, params.P_tx, params.K, params.noise, params.weights);

        case {'SC','GC','FC'}
            % AO WMMSE + Manopt RIS
            if strcmpi(scheme,'GC')
                arch_param = params.Ng;
                [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_GC(ch.H_BS_user, ch.H_RIS_user, ch.H_BS_RIS, params.tolerance, params.P_tx, params.noise, params.weights, arch_param);
            elseif strcmpi(scheme,'FC')
                arch_param = [];
                [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_FC(ch.H_BS_user, ch.H_RIS_user, ch.H_BS_RIS, params.tolerance, params.P_tx, params.noise, params.weights);
            elseif strcmpi(scheme,'SC')
                [tempWSR, ~] = Alternating_Optimization_WMMSE_CF_Manopt_SC(ch.H_BS_user, ch.H_RIS_user, ch.H_BS_RIS, params.tolerance, params.P_tx, params.noise, params.weights);
            end

        otherwise
            error('Unknown scheme %s', scheme);
    end
    wsr_all(iMonte) = tempWSR;
end

% --- 5) Compute statistics ---
avgWSR     = mean(wsr_all, 'omitnan');
outageProb = sum(wsr_all < params.tau, 'omitnan') / (sum(~isnan(wsr_all)));
end
