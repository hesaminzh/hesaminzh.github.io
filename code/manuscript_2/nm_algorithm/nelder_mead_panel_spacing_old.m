function [best_ds, best_rate, history] = nelder_mead_panel_spacing_old(ds_init, step, obj_fun, tol, max_iter)
%NELDER_MEAD_PANEL_SPACING Optimise spacing variables using the Nelder–Mead algorithm.
%
%   [BEST_DS, BEST_RATE, HISTORY] = NELDER_MEAD_PANEL_SPACING(DS_INIT, STEP, OBJ_FUN, TOL, MAX_ITER)
%   performs an unconstrained minimisation of the objective function OBJ_FUN
%   starting from an initial guess DS_INIT using the Nelder–Mead simplex
%   algorithm.  DS_INIT should be a column vector of length N containing
%   the initial values of the optimisation variables.  STEP can be either
%   a scalar or a vector of length N specifying the size of the initial
%   simplex along each coordinate direction.  OBJ_FUN must be a function
%   handle accepting a column vector of length N and returning a scalar
%   objective value.  TOL is the convergence tolerance used when
%   determining whether the simplex has collapsed sufficiently or the
%   function values no longer change significantly.  MAX_ITER sets a
%   hard upper bound on the number of iterations; a robust `while` loop is
%   used to enforce this without prematurely breaking out of the loop.
%
%   This implementation generalises the classic Nelder–Mead method to
%   multi‑dimensional search spaces.  For an N‑dimensional problem, the
%   simplex contains N+1 vertices, each of which is a length‑N column
%   vector.  The initial simplex is formed by adding STEP to each
%   coordinate of DS_INIT in turn.  Throughout the iterations the
%   algorithm orders the vertices according to the objective values,
%   computes the centroid of the best N vertices, and then attempts
%   reflection, expansion, contraction or shrink operations in order to
%   progress towards a minimum.  A robust inside contraction formula,
%   XC = CENTROID − RHO*(XR − CENTROID), is used as recommended in
%   recent literature.
%
%   HISTORY records the state of the simplex at each iteration.  Each
%   row contains the concatenated simplex coordinates followed by their
%   corresponding objective values.  For an N‑dimensional problem the
%   history array therefore has (N+1)² columns.
%
%   Example:
%       % Minimise the function f(x) = (x(1)-1)^2 + (x(2)+3)^2
%       obj = @(x) (x(1)-1)^2 + (x(2)+3)^2;
%       init = [0; 0];
%       step = 1;
%       [x_opt, f_opt, hist] = nelder_mead_panel_spacing(init, step, obj, 1e-6, 200);
%       fprintf('Optimal x: %g %g\n', x_opt(1), x_opt(2));
%       fprintf('Optimal f: %g\n', f_opt);

% Check and reshape the initial guess
ds_init = ds_init(:); % ensure column vector
n = numel(ds_init);    % dimensionality

% Handle the step input
if isscalar(step)
    % uniform step size along each dimension
    step_vec = repmat(step, n, 1);
else
    % ensure step is a column vector of length n
    step_vec = step(:);
    if numel(step_vec) ~= n
        error('STEP must be either a scalar or a vector of the same length as DS_INIT.');
    end
end

% Set the standard Nelder–Mead coefficients
alpha = 1.0;  % reflection
gamma = 2.0;  % expansion
rho   = 0.5;  % contraction
sigma = 0.5;  % shrink

% Construct the initial simplex: N+1 vertices in an N‑dimensional space
% The first vertex is the initial guess.  Subsequent vertices add
% STEP_VEC(j) to the j‑th coordinate.
ds = repmat(ds_init, 1, n+1);
for j = 1:n
    ds(:, j+1) = ds(:, j+1) + step_vec(j);
end

% Evaluate the objective function at each vertex
fvals = zeros(1, n+1);
for j = 1:(n+1)
    fvals(j) = obj_fun(ds(:, j));
end

% Preallocate the history matrix.  Each row stores (N+1) sets of
% coordinates followed by their function values: (N+1)*(N+1) columns.
history = zeros(max_iter, (n+1)^2);

% Iteration counter
iter = 0;

% Main optimisation loop
while iter < max_iter
    iter = iter + 1;

    % Order the simplex according to the objective values (ascending)
    [fvals, order] = sort(fvals);
    ds = ds(:, order);

    % Save current simplex and objective values into history
    % Flattened as [x1', x2', ..., x_{n+1}', f1, f2, ..., f_{n+1}]
    row = zeros(1, (n+1)^2);
    idx = 1;
    for j = 1:(n+1)
        row(idx:(idx+n-1)) = ds(:, j).';
        idx = idx + n;
    end
    row(idx:end) = fvals;
    history(iter, :) = row;

    % Compute the centroid of the best N vertices (all except the worst)
    centroid = mean(ds(:, 1:n), 2);

    % Reflection step
    xr = centroid + alpha * (centroid - ds(:, n+1));
    fr = obj_fun(xr);

    % Case 1: the reflected point is better than the best point
    if fr < fvals(1)
        % Expansion step
        xe = centroid + gamma * (xr - centroid);
        fe = obj_fun(xe);
        if fe < fr
            ds(:, n+1) = xe;
            fvals(n+1) = fe;
        else
            ds(:, n+1) = xr;
            fvals(n+1) = fr;
        end

    % Case 2: the reflected point is better than the second worst
    elseif fr < fvals(n)
        % Accept the reflected point
        ds(:, n+1) = xr;
        fvals(n+1) = fr;

    else
        % The reflected point is not an improvement over most of the simplex
        % Decide between outside contraction, inside contraction or shrink
        if fr < fvals(n+1)
            % Outside contraction
            xc = centroid + rho * (xr - centroid);
            fc = obj_fun(xc);
            if fc <= fr
                ds(:, n+1) = xc;
                fvals(n+1) = fc;
            else
                % Shrink the simplex towards the best vertex
                for j = 2:(n+1)
                    ds(:, j) = ds(:, 1) + sigma * (ds(:, j) - ds(:, 1));
                    fvals(j) = obj_fun(ds(:, j));
                end
            end
        else
            % Inside contraction (use robust formulation)
            xc = centroid - rho * (xr - centroid);
            fc = obj_fun(xc);
            if fc < fvals(n+1)
                ds(:, n+1) = xc;
                fvals(n+1) = fc;
            else
                % Shrink the simplex towards the best vertex
                for j = 2:(n+1)
                    ds(:, j) = ds(:, 1) + sigma * (ds(:, j) - ds(:, 1));
                    fvals(j) = obj_fun(ds(:, j));
                end
            end
        end
    end

    % Convergence check based on simplex size and function range
    simplex_size = 0;
    for j = 2:(n+1)
        simplex_size = max(simplex_size, norm(ds(:, j) - ds(:, 1)));
    end
    f_range = max(fvals) - min(fvals);
    if simplex_size < tol || f_range < tol
        % Trim unused history rows before exiting
        history = history(1:iter, :);
        break;
    end
end

% If maximum iterations reached, trim unused history rows
if iter >= max_iter
    history = history(1:iter, :);
end

% Order the final simplex and return the best point
[fvals, order] = sort(fvals);
ds = ds(:, order);
best_ds = ds(:, 1);
best_rate = fvals(1);
end
