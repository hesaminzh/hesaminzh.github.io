function a = array_response_ula(N, phi)
%ARRAY_RESPONSE_ULA Compute array response for a ULA of N elements at angle phi.
% Assumes half-wavelength antenna spacing and broadside reference.
% phi: angle of departure/arrival (in radians) measured from array broadside (0 = broadside).
    n = (0:N-1).';   % element indices
    % For half-wavelength spacing, phase shift = pi * sin(phi) * element_index
    a = exp(-1j * pi * sin(phi) * n);
    % Normalize the vector
    a = a / norm(a);
end
