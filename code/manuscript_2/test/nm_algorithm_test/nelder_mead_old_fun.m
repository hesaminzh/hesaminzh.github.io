function [x_opt, f_opt, iter, history] = nelder_mead_old_fun(fun, x0, options)
%NELDER_MEAD_TEST_FUN  Simple implementation of the Nelder‑Mead simplex algorithm.
%
% INPUT
%   fun      – function handle,  f = fun(x)  (x is a column vector)
%   x0       – initial guess (column vector, dimension n)
%   options  – structure with fields (all optional)
%       .alpha   – reflection coefficient (default 1)
%       .gamma   – expansion coefficient (default 2)
%       .beta    – contraction coefficient (default 0.5)
%       .delta   – shrinkage coefficient (default 0.5)
%       .maxIter – maximum number of iterations (default 200)
%       .tolSize – simplex size tolerance (default 1e-6)
%       .tolFun  – function‑value tolerance (default 1e-8)
%       .initialSimplex – matrix n×(n+1) with user‑supplied initial simplex.
%                        If omitted, a regular simplex around x0 is built.
%
% OUTPUT
%   x_opt   – best vertex found (column vector)
%   f_opt   – function value at x_opt
%   iter    – number of iterations performed
%   history – struct with fields X (cell array of simplexes) and f (cell array
%             of function values) useful for plotting or debugging.
%
% The algorithm follows the classic steps:
%   1) order vertices,
%   2) compute centroid of the n best points (excluding the worst vertex),
%   3) reflect,
%   4) possibly expand,
%   5) possibly contract (outside or inside),
%   6) possibly shrink.
%
%--------------------------------------------------------------------------

% % ---------- 1.  Parse options & set defaults ----------
if nargin < 3
    options = struct();
end
alpha = getOption(options, 'alpha', 1);
gamma = getOption(options, 'gamma', 2);
beta  = getOption(options, 'beta' , 0.5);
delta = getOption(options, 'delta', 0.5);
maxIter = getOption(options, 'maxIter', 200);
tolSize = getOption(options, 'tolSize', 1e-6);
tolFun  = getOption(options, 'tolFun' , 1e-8);

% alpha = options.alpha;
% gamma = options.gamma;
% beta = options.beta;
% delta = options.delta;
% maxIter = options.maxIter;
% tolSize = options.tolSize;
% tolFun = options.tolFun;

% ---------- 2.  Initialise simplex ----------
n = length(x0);                     % problem dimension
if isfield(options, 'initialSimplex')
    X = options.initialSimplex;     % user supplied (n × (n+1))
    if size(X,1)~=n || size(X,2)~=n+1
        error('initialSimplex must be of size n × (n+1).');
    end
else
    % Build a regular simplex around x0 (Nelder & Mead, 1965)
    %   x0 is the first vertex, the other n vertices are displaced
    %   along each coordinate direction by a factor of 0.05 * max(|x0|,1).
    scale = 0.05 * max(abs(x0),1);
    X = zeros(n, n+1);
    X(:,1) = x0;
    for i = 1:n
        ei = zeros(n,1);
        ei(i) = 1;
        X(:,i+1) = x0 + scale(i) * ei;
    end
end

% ---------- 3.  Allocate history storage ----------
history.X = cell(maxIter+1,1);
history.f = cell(maxIter+1,1);
history.X{1} = X;

% ---------- 4.  Main loop ----------
iter = 0;
while iter < maxIter
    iter = iter + 1;

    % ---- 4.1  Evaluate function at all vertices -----------------
    fVals = zeros(1, n+1);
    for j = 1:n+1
        fVals(j) = fun( X(:,j) );
    end

    % ---- 4.2  Sort vertices (ascending) -------------------------
    [fVals, idx] = sort(fVals);
    X = X(:, idx);               % reorder columns to match sorted fVals

    % ---- 4.3  Check termination criteria ------------------------
    simplexSize = max(pdist( X.' ));      % largest edge length
    fRange = max(fVals) - min(fVals);
    if simplexSize < tolSize || fRange < tolFun
        break;
    end

    % ---- 4.4  Centroid of the best n points (exclude worst) -----
    x_c = mean( X(:,1:n), 2 );    % column vector

    % ---- 4.5  Reflection ----------------------------------------
    x_r = x_c + alpha * (x_c - X(:,end));
    f_r = fun(x_r);

    if f_r < fVals(1)                                 % (a) reflected is best -> expansion
        x_e = x_c + gamma * (x_r - x_c);
        f_e = fun(x_e);
        if f_e < f_r
            X(:,end) = x_e;       % replace worst by expanded point
        else
            X(:,end) = x_r;       % replace worst by reflected point
        end

    elseif f_r < fVals(end-1)                        % (b) reflected is good enough
        X(:,end) = x_r;           % simply accept the reflected point

    else                                             % (c) reflected is poor -> contraction
        if f_r < fVals(end)                           % ----- outside contraction
            x_oc = x_c + beta * (x_r - x_c);
            f_oc = fun(x_oc);
            if f_oc <= f_r
                X(:,end) = x_oc;   % accept outside‑contracted point
            else
                % ---- shrink whole simplex towards best vertex ----
                for j = 2:n+1
                    X(:,j) = X(:,1) + delta * (X(:,j) - X(:,1));
                end
            end

        else                                            % ----- inside contraction
            x_ic = x_c - beta * (x_c - X(:,end));
            f_ic = fun(x_ic);
            if f_ic < fVals(end)
                X(:,end) = x_ic;   % accept inside‑contracted point
            else
                % ---- shrink whole simplex towards best vertex ----
                for j = 2:n+1
                    X(:,j) = X(:,1) + delta * (X(:,j) - X(:,1));
                end
            end
        end
    end

    % ---- store current simplex for post‑processing ----------------
    history.X{iter+1} = X;
    history.f{iter+1} = fVals;
end

% ---------- 5.  Return results ----------
x_opt = X(:,1);          % best vertex (lowest f)
f_opt = fVals(1);
history.X = history.X(1:iter+1);
history.f = history.f(1:iter+1);
end

% -------------------------------------------------------------------------
function val = getOption(opt, name, default)
% Helper that returns opt.(name) if it exists, otherwise default.
if isfield(opt, name)
    val = opt.(name);
else
    val = default;
end
end

