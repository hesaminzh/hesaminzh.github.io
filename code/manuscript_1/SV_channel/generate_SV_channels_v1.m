function [H_BS_user, H_BS_RIS, H_RIS_user] = generate_SV_channels_v1(M, Nx, Ny, K, d_BR_horiz, d_BU_horiz)
% generate_SV_channels  SV mmWave channels for BS--RIS, RIS--users, BS--users
% Inputs:
%   M             – # BS antennas
%   Nx, Ny        – RIS array dims, so N = Nx*Ny
%   K             – # users
%   d_BR_horiz    – horizontal BS→RIS distance (m)
%   d_BU_horiz    – 1×K vector of horizontal BS→user distances (m)
% Outputs:
%   H_BS_user (M×K)       – BS→user channel
%   H_BS_RIS (N×M)       – BS→RIS channel
%   H_RIS_user (N×K)      – RIS→user channels

% heights (m)
h_BS = 6;    h_RIS = 3;    h_U = 1.5;
% vertical separations
dv_BR = abs(h_BS - h_RIS);
dv_RU = abs(h_RIS - h_U);
dv_BU = abs(h_BS - h_U);

N = Nx * Ny;
% SV parameters
a_LOS = 61.4;   b_LOS = 2.0;   sigma_LOS = 5.8;
a_NLOS = 72.0;  b_NLOS = 2.92; sigma_NLOS = 8.7;
K_dB   = 13.2;  Krician = 10^(K_dB/10);
%lambda_gain = 9.60;    % linear from dBi

%% 1) BS→RIS (LOS Rician)
d_BR = sqrt(d_BR_horiz^2 + dv_BR^2);
xi_LOS = sigma_LOS * randn();
PL_BR_dB = a_LOS + 10*b_LOS*log10(d_BR) + xi_LOS;
PL_BR = 10^(-PL_BR_dB/10);

% array responses
% BS ULA (N×1) broadside
phi_BS = 0;
a_t = array_response_ula(M, phi_BS);  % Tx
% RIS URA (N×1) broadside
theta_el = -atan(dv_BR / max(d_BR_horiz,1e-6));
theta_az = 0;
a_r = array_response_ura(Nx, Ny, theta_az, theta_el);  % Rx

G_LOS = a_r * a_t';                  % N×M
G_NLOS = (randn(N, M) + 1j*randn(N, M))/sqrt(2);
H_BS_RIS = sqrt(PL_BR) * ( ...
      sqrt(Krician/(Krician+1)) * sqrt(N*M) * G_LOS + ...
      sqrt(1/(Krician+1))      * sqrt(N*M) * (1/sqrt(N*M)) * G_NLOS );

%% 2) RIS→User (each user k)
H_RIS_user = zeros(N, K);
for k = 1:K
    d_RU = sqrt((abs(d_BU_horiz(k) - d_BR_horiz))^2 + dv_RU^2);
    xi_LOS = sigma_LOS * randn();
    PL_RU = 10^(-(a_LOS + 10*b_LOS*log10(d_RU) + xi_LOS)/10);
    
    % LOS
    if d_BU_horiz(k) <= d_BR_horiz
        az = 0; 
    else
        az = pi;
    end
    el = -atan(dv_RU / max(abs(d_BU_horiz(k)-d_BR_horiz),1e-6));
    a_t = array_response_ura(Nx, Ny, az, el);  % URA Tx
    
    h_NLOS = (randn(N,1) + 1j*randn(N,1))/sqrt(2);
    H_RIS_user(:,k) = sqrt(PL_RU) * ( ...
                  sqrt(Krician/(Krician+1)) * sqrt(N) * a_t + ...
                  sqrt(1/(Krician+1))      * sqrt(N) * (1/sqrt(N)) * h_NLOS );
end

%% 3) BS→User (NLOS SV with Ld paths)
Ld = 3;
H_BS_user = zeros(M, K);
for k = 1:K
    d_BU = sqrt(d_BU_horiz(k)^2 + dv_BU^2);
    xi_NLOS = sigma_NLOS * randn();
    PL_BU = 10^(-(a_NLOS + 10*b_NLOS*log10(d_BU) + xi_NLOS)/10);
    
    for l = 1:Ld
        alpha = sqrt(PL_BU) * (randn+1j*randn)/sqrt(2*Ld);
        phi   = 2*pi*rand();
        a_t   = array_response_ula(M, phi);
        H_BS_user(:,k) = H_BS_user(:,k) + sqrt(M) * alpha * a_t;
    end
end
end
