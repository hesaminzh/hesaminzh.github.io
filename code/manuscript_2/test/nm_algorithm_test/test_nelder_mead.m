% DEMONELDERMEAD  Demonstrates the NM implementation on a simple 2‑D problem.
% The objective is
%      f(x) = (x1-2)^2 + (x2+3)^2,
% whose global minimum is at (2, -3) with f_min = 0.

clear;
%close all;
clc;

% ---------- 2.1  Define the objective function ----------
% fun = @(x) (x(1)-2)^2 + (x(2)+3)^2;   % anonymous function, expects column vector
fun = @(x) (x(1)-0.2)^2 + (x(2)+22.5)^2;   % anonymous function, expects column vector
% fun = @(x) (x(1)-17)^2 + (x(2)+2)^2;   % anonymous function, expects column vector

% ---------- 2.2  Choose an initial guess ----------
x0 = [0; 0];          % far away from the optimum on purpose

% ---------- 2.3  (optional) tweak NM parameters ----------
opts.alpha  = 1;      % reflection
opts.gamma  = 2;      % expansion
opts.beta   = 0.5;    % contraction
opts.delta  = 0.5;    % shrink
opts.maxIter= 200;
opts.tolSize= 1e-8;
opts.tolFun = 1e-12;

% ---------- 2.4  Call the Nelder‑Mead routine ----------
% [x_opt, f_opt, it, hist] = nelder_mead_old_fun(fun, x0, opts);
[x_opt, f_opt, it, hist] = my_NM_algorithm(fun, x0, opts);

% ---------- 2.5  Report the solution ----------
fprintf('Nelder‑Mead finished after %d iterations.\n', it);
fprintf('Estimated optimum    : (%.6f , %.6f)\n', x_opt(1), x_opt(2));
fprintf('Objective at optimum : %.6e\n', f_opt);
fprintf('True optimum         : (2 , -3)  with f = 0\n');

% ---------- 2.6  Visualise the simplex evolution (optional) ----------
% Plot the contour of the quadratic and superimpose the simplex at each iteration.
[X1, X2] = meshgrid(linspace(-2,4,200), linspace(-7,1,200));
Z = (X1-2).^2 + (X2+3).^2;

figure; hold on; box on;
contour(X1, X2, Z, 30, 'LineColor', [0.6 0.6 0.6]);   % iso‑lines of the objective
plot(2, -3, 'kp', 'MarkerSize',12, 'MarkerFaceColor','k'); % true optimum

colormap(parula);
for k = 1:length(hist.X)
    simplex = hist.X{k};
    % close the triangle for nicer drawing
    xs = [simplex(1,:) simplex(1,1)];
    ys = [simplex(2,:) simplex(2,1)];
    plot(xs, ys, '-o', 'Color', [0 0.6 0], 'LineWidth',1.2);
    pause(0.05);   % slow‑down to see the motion
end
title('Nelder‑Mead simplex moving toward the minimum');
xlabel('x_1'); ylabel('x_2');
legend('objective contours','true optimum','simplex trajectory','Location','best');
hold off;

