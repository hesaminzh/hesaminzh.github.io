function a = array_response_UPA(a1,a2,N,p0)
%ARRAY_RESPONSE_UPA Compute the array response of a uniform planar array (UPA).
%
%   a = array_response_UPA(a1,a2,N,p0) returns an N×1 vector containing the
%   array steering vector for a uniform planar array (UPA) with N elements.
%   The array is assumed to have p0 elements along the vertical dimension and
%   N/p0 elements along the horizontal dimension. The azimuth angle is a1 and
%   the elevation angle is a2, both expressed in radians. The element spacing
%   is assumed to be half a wavelength, which results in a phase shift of
%   π for each unit distance in the direction cosines.
%
%   Inputs:
%     a1  - azimuth angle in radians
%     a2  - elevation angle in radians
%     N   - total number of antenna elements (product of vertical and horizontal counts)
%     p0  - number of elements along the vertical (row) dimension
%
%   Output:
%     a   - steering vector of size N×1
%
%   This implementation follows the expression used in many mmWave channel
%   models where the steering vector entry corresponding to the (m,n)-th
%   element is exp(j π (m sin(a1) sin(a2) + n cos(a2))).
%
%
% Hesam Nezhadmohammad (last revised: 19/10/2025)

% Preallocate the output vector
a = zeros(N,1);

% Compute the number of elements in each dimension
num_vert = p0;
num_horiz = N/p0;

% Fill in the response vector
idx = 1;
for m = 0:(num_vert-1)
    for n = 0:(num_horiz-1)
        a(idx) = exp(1i * pi * ( m * sin(a1) * sin(a2) + n * cos(a2) ));
        idx = idx + 1;
    end
end

end


