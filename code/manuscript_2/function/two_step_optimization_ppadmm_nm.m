function [best_ds, best_rate] = two_step_optimization_ppadmm_nm(My, M_list, Mtot, arch, Ns, ds0_min, ds0_max, tol_ds, max_iter_nm, params)
%TWO_STEP_OPTIMIZATION_PPADMM_NM  two step optimization with ppadmm and nelder-mead algorithm
% Detailed explanation goes here
%
% Hesam Nezhadmohammad


ds0 = ds0_min + (ds0_max - ds0_min).*rand(1,Ns);
ds0 = sort(ds0);


obj_handle = @(d_s) objective_sum_rate_for_spacing(d_s, My, M_list, Mtot, arch, params);
% [opt_ds, opt_val, mean_vals, std_vals] = nelder_mead_panel_spacing(ds0, obj_handle, tol_ds, max_iter_nm);
[opt_ds, opt_val, ~, ~] = nelder_mead_panel_spacing(ds0, obj_handle, tol_ds, max_iter_nm);
best_ds = opt_ds;
best_rate = -opt_val;

% fprintf('Optimum distance is %g\n', best_ds);
% fprintf('Best rate is %g\n', best_rate);
%
% valid_indices = find(~isnan(mean_vals));
% last_iter = valid_indices(end);
% mean_to_plot = mean_vals(1:last_iter);
% std_to_plot = std_vals(1:last_iter);
%
% figure('DefaultAxesFontSize',16); %hold on;
% hold on;
%
% % Calculate upper and lower bounds for the shaded region
% upper_bound = mean_to_plot + std_to_plot;
% lower_bound = mean_to_plot - std_to_plot;
%
% plot(mean_to_plot, 'b-', 'LineWidth', 2, 'DisplayName', 'Mean');
%
% % Shade the standard deviation area
% % We create vertices for the shaded area by going forward along the upper bound
% % and backwards along the lower bound to close the shape.
% x_values = 1:length(mean_to_plot);
% % fill([x_values,fliplr(x_values)], [upper_bound, fliplr(lower_bound)], [0.8, 0.8, 0.8], ...
% %     'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Std')
% fill([x_values,fliplr(x_values)], [upper_bound, fliplr(lower_bound)], [150, 206, 180]/255, ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Std')
% xlabel('Number of iterations');
% ylabel('$D_q$', 'Interpreter', 'LaTeX');
% legend;
% set(gcf, 'Color', [1,1,1]);
% set(gca, 'LineWidth',1);
% xlim([0, last_iter]);
% % ylim([0, max(upper_bound)+mean(upper_bound)/2]);
% ylim([0, 3]);
% % set(gca, 'xtick', M_RIS);
% grid on;
% box on;

% print(gcf, 'convergence_nm', '-depsc2' )

end
