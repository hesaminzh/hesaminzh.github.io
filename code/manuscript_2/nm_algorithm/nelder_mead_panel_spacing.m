function [best_ds,best_rate,mean_history,std_history] = nelder_mead_panel_spacing(ds_init, obj_fun, tol, max_iter)
%NELDER_MEAD_PANEL_SPACING  Optimize a scalar spacing parameter using the
%Nelder–Mead simplex search.
%
%   [BEST_DS,BEST_RATE,HISTORY] = NELDER_MEAD_PANEL_SPACING(DS_INIT,OBJ_FUN,TOL,MAX_ITER)
%   performs a derivative‑free optimization of a one–dimensional parameter
%   using the Nelder–Mead algorithm. DS_INIT specifies the starting
%   guess for the parameter (e.g., an initial inter‑panel spacing), STEP
%   is the initial simplex step size, OBJ_FUN is a function handle
%   returning a scalar cost (to be minimized) for a given spacing,
%   TOL sets the termination threshold on the simplex size, and
%   MAX_ITER is the maximum number of iterations allowed. The
%   function returns the best spacing found (BEST_DS), the
%   corresponding objective value (BEST_RATE), the history of mean
%   values at each iteration (MEAN_HISTORY), and the history of
%   standard deviation (std) values (STD_HISTORY).
%
%
%   The Nelder–Mead method maintains a simplex of n+1 points in n
%   dimensions (here n=1). At each iteration the points are ordered
%   according to their objective values, and the worst point is
%   replaced by a new point obtained by reflection, expansion,
%   contraction or shrink steps. Standard coefficients are used for
%   these operations, as described in the Nelder–Mead literature: the
%   reflection coefficient alpha=1, expansion coefficient gamma=2, contraction
%   coefficient rho=1/2 and shrink coefficient delta=1/2.
%
% Source 1:
% X. Long, Y. Li, M. Pauli, B. Nuss, and T. Zwick, “Joint BS deployment
% and power optimization for minimum EMF exposure with RL in real-
% world based urban scenario,” 2026, arXiv:2504.05017.
%
% Source 2:
% Multi-RIS Deployment Optimization for mmWave
% ISAC Systems in Real-World Environments
% I didn't have access to internet for the exact url :')
% Date: 10/4/2026
% Hesam Nezhadmohammad


% Standard coefficients
alpha = 1;    % reflection
gamma = 2;    % expansion
rho   = 0.5;  % contraction
sigma = 0.5;  % shrink

% Initialize simplex (two points in 1D)
% ds = [ds_init, ds_init + step];
ds = ds_init; % row vector
% if size(ds_init,1) == 1
%     ds = ds_init; % we work with row vector so later pdist goes as pdist(ds.')
% elseif size(ds_init,2) == 1
%     ds = ds_init.';
% end

fvals = zeros(size(ds));
for i = 1:length(ds)
    fvals(i) = obj_fun(ds(i));
end

% History storage: each row is [ds(1) fvals(1) ds(2) fvals(2) ... ds(end) fvals(end)]
% history = zeros(max_iter,length(ds)*2);

mean_history = NaN(1, max_iter);
std_history = NaN(1, max_iter);
iter = 0;
while iter < max_iter % q
    iter = iter + 1;
    % Order simplex by objective value (ascending)
    [fvals, idx] = sort(fvals);
    ds = ds(idx);

    % Save history
    % history(iter,:) = reshape([ds(:) fvals(:)].', 1, [] );

    % Compute the centroid of the best point(s) excluding worst
    % centroid = ds(1);
    x_c = mean(ds(1:end-1), 2); % centroid exluding the worst

    % Reflection step
    % xr = x_c + alpha*(x_c - ds(2));
    xr = x_c + alpha*(x_c - ds(end));
    fr = obj_fun(xr);

    % Determine region of reflected point
    if fr < fvals(1)
        % Expansion
        xe = x_c + gamma*(xr - x_c);
        fe = obj_fun(xe);
        if fe < fr
            % ds(2)   = xe;
            ds(end)   = xe;
            % fvals(2)= fe;
            fvals(end)= fe;
        else
            % ds(2)   = xr;
            ds(end)   = xr;
            % fvals(2)= fr;
            fvals(end)= fr;
        end
    % elseif fr < fvals(2)
    elseif fvals(1) <= fr && fr < fvals(end-1)
        % Accept reflection
        % ds(2)   = xr;
        ds(end)   = xr;
        % fvals(2)= fr;
        fvals(end)= fr;
    else
        % Contraction
        % if fr < fvals(2)
        if fvals(end-1) <= fr && fr < fvals(end)                           % ----- outside contraction
            % Outside contraction
            x_oc = x_c + rho*(xr - x_c);
            fc = obj_fun(x_oc);
            if fc <= fr
                ds(end)   = x_oc;
                fvals(end)= fc;
            else
                % Shrink
                % ---- shrink whole simplex towards best vertex ----
                for j = 2:length(ds)
                    ds(j) = ds(1) + sigma * (ds(j) - ds(1));
                    fvals(j)= obj_fun(ds(j));
                end
            end
        else
        % elseif fvals(end) <= fr
            % Inside contraction
            x_ic = x_c - rho*(xr - x_c);
            fc = obj_fun(x_ic);
            if fc < fvals(end)
                ds(end)   = x_ic;
                fvals(end)= fc;
            else
                % Shrink
                % ---- shrink whole simplex towards best vertex ----
                for j = 2:length(ds)
                    ds(j) = ds(1) + sigma * (ds(j) - ds(1));
                    fvals(j)= obj_fun(ds(j));
                end
            end
        end
    end

    if iter >= max_iter
        best_ds = NaN;
        best_rate = NaN;
        % history = NaN(max_iter,length(ds)*2);
        fprintf('Reached max_iter in NM (%g)\n', max_iter );
        return;
    end


    pairwise_dist = pdist( ds.' ); % ds is a row vector so we transpose it
    mean_history(iter) = mean(pairwise_dist);
    std_history(iter) = std(pairwise_dist);
    simplexSize = max(pairwise_dist);      % largest edge length
    % fRange = max(fvals) - min(fvals);
    if simplexSize < tol% || fRange < tol
        break;
    end

end

% Order final simplex
[fvals, idx] = sort(fvals);
ds = ds(idx);

best_ds   = ds(1);
best_rate = fvals(1);
% Trim unused history rows
% history = history(1:iter,:);
end
