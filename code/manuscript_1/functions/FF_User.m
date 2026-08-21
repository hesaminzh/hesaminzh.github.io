function [H2] = FF_User(Nx,Ny,K,d)

% RIS-user channel matrix

N = Nx*Ny;
lambda_gain = 9.60; % 9.82 dBi
a = 61.4;
b = 2;
sigma = 3.8;  % 5.8 dB
kappa = a+10*b*log10(d)+normrnd(0, sigma, 1);

% far-field array
angU_e = rand(1, K)-0.5;                         % elevation sin_e
angU_a = rand(1, K)-0.5;    % azimuth cos_a*cos_e

% array steering vector
arrU_a = (1/sqrt(Nx))*exp(-2*1i*pi*angU_a'*((0:Nx-1)-(Nx-1)/2))';
arrU_e = (1/sqrt(Ny))*exp(-2*1i*pi*angU_e'*((0:Ny-1)-(Ny-1)/2))';

% loss gain
beta_user = sqrt(N)*lambda_gain*(sqrt(1/2*10.^(-0.1*kappa))).*(randn(K, 1)+1j*randn(K, 1));

for k=1:K
    H2(:, k) = beta_user(k) * kron(arrU_a(:, k), arrU_e(:, k));
end


end

