function a = array_response_ula(N,phi)
%ARRAY_RESPONSE_ULA Compute the array response of a uniform linear array (ULA).
%
%   a = ARRAY_RESPONSE_ULA(N,phi) returns the steering vector for
%   a ULA with N elements and half‑wavelength spacing. The angle phi is the
%   direction of departure/arrival measured from broadside (in radians).
%
%   The k‑th element of the steering vector is exp(-j*pi*(k-1)*sin(phi)).
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)


    n = (0:N-1).';
    a = exp(-1j*pi*sin(phi)*n);
end
