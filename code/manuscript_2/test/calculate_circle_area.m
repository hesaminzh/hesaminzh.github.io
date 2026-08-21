% calculate_circle_area.m
function area = calculate_circle_area(radius, pi_value)
    % Calculates the area of a circle.
    % radius: the radius of the circle
    % pi_value: the value of pi to use (e.g., built-in pi or an approximation)
    area = pi_value * radius.^2;
end

