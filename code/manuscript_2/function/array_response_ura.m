function a = array_response_ura(My,Mz,az,el)
%ARRAY_RESPONSE_URA Compute the array response for a uniform rectangular array (URA).
%
%   a = ARRAY_RESPONSE_URA(My,Mz,az,el) returns the steering vector
%   for a URA with My elements along the horizontal (y) axis and Mz elements
%   along the vertical (z) axis. The azimuth angle az and elevation angle el
%   are measured in radians. Half‑wavelength spacing is assumed between
%   neighbouring elements.
%
%   The horizontal and vertical steering vectors are given by
%     a_y(m) = exp(j * (-pi * sin(az) * cos(el)) * m)
%     a_z(n) = exp(j * (-pi * sin(el)) * n)
%   and the full URA response is the Kronecker product of a_y and a_z
% reference: Secrecy Capacity Maximization for a Hybrid Relay-RIS Scheme in mmWave MIMO Networks (https://arxiv.org/pdf/2205.13904)
%
% Hesam Nezhadmohammad (last revised: 20/10/2025)

    my_idx = (0:My-1).';
    mz_idx = (0:Mz-1).';
    ky = -pi * sin(az) * cos(el);
    kz = -pi * sin(el);
    a_y = exp(1j * ky * my_idx);
    a_z = exp(1j * kz * mz_idx);
    A = a_y * a_z.';
    a = reshape(A, My*Mz, 1);
end
