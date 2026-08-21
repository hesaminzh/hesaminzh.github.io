% In your main script or command window:

% Define fixed parameters
target_area = 50;
my_pi = 3.14159;

% Create the anonymous function handle
% This anonymous function 'objective_func' takes only one argument 'r' (radius)
% and calls 'calculate_circle_area' with 'r' and the fixed 'my_pi'.
objective_func = @(r) calculate_circle_area(r, my_pi);

% Define the function we want to minimize: the difference between calculated
% area and the target area. We want to find 'r' where this difference is zero.
% fminbnd finds the minimum, so we'll minimize the *absolute difference*.
error_func = @(r) abs(objective_func(r) - target_area);

% Use fminbnd to find the radius that minimizes the error.
% We'll search for the radius between 0 and 10.
initial_guess_radius = 5; % A starting point for the search
radius_result = fminbnd(error_func, 0, 10);

% Display the results
fprintf('The radius that gives an area closest to %f is: %f\n', target_area, radius_result);

% Verify the area with the found radius
calculated_area = calculate_circle_area(radius_result, my_pi);
fprintf('The calculated area with this radius is: %f\n', calculated_area);

% Let's see what happens if we tried to use calculate_circle_area directly
% This would fail if the optimizer expects a single-input function:
% fminbnd(@calculate_circle_area, 0, 10); % This would give an error!

