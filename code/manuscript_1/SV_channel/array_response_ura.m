function a = array_response_ura(My, Mz, az, el)
%ARRAY_RESPONSE_URA Compute array response for a URA of size My x Mz at given azimuth and elevation.
% Elements arranged with My along horizontal (y-axis) and Mz along vertical (z-axis).
% az: azimuth angle (rad) from array broadside (0 = broadside direction in horizontal plane)
% el: elevation angle (rad) from array broadside (0 = same horizontal plane, positive = upward)
    % Element index grid
    my_idx = (0:My-1).';   % horizontal index (along y-axis)
    mz_idx = (0:Mz-1).';   % vertical index (along z-axis)
    % Wavenumber phase increments for half-wavelength spacing:
    % Horizontal phase increment = -pi * sin(az) * cos(el)
    % Vertical phase increment   = -pi * sin(el)
    ky = -pi * sin(az) * cos(el);
    kz = -pi * sin(el);
    % Construct array response as kronecker of horizontal and vertical steering vectors
    a_y = exp(1j * ky * my_idx);  % My x 1
    a_z = exp(1j * kz * mz_idx);  % Mz x 1
    A = a_y * a_z.';              % My x Mz matrix of element phases
    a = reshape(A, My*Mz, 1);
    a = a / norm(a);
end
