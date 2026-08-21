function [H_BS_RIS, H_RIS_user] = generate_SV_channel(N, Mx, My, user_loc, BS_loc, panel_loc)
%GENERATE_SV_CHANNEL Generate Saleh–Valenzuela (SV) mmWave channels for a single RIS panel.
%
%   [H_BS_RIS, H_RIS_user] = generate_SV_channel_extended(N, Mx, My, user_loc,
%   BS_loc, panel_loc) returns the BS->RIS channel matrix H_BS_RIS (size
%   M×N) and the RIS->user channel matrix H_RIS_user (size M×K), where K
%   is the number of users.  The function supports arbitrary panel
%   positions and orientations.  The input arguments are:
%
%     N        : number of antennas at the BS (scalar)
%     Mx, My   : number of elements along the RIS vertical and horizontal
%                dimensions (scalars), so that M = Mx*My
%     user_loc : K×2 matrix of user coordinates in the horizontal plane
%                (x,y) relative to a common origin (meters)
%     BS_loc   : 1×2 vector specifying the horizontal coordinates of the BS
%     panel_loc: 1×2 vector specifying the horizontal coordinates of this
%                RIS panel
%
%
%     'a_LOS','b_LOS','sigma_LOS'
%     'a_NLOS','b_NLOS','sigma_NLOS'
%                       : path-loss model parameters (defaults taken from
%                         the IEEE TVT'20 RIS paper)
%     'K_dB'         : Rician K-factor in dB (default 0 ⇒ Rayleigh)
%
%   The function implements the SV geometric model for mmWave channels.
%   NLOS contributions are normalised by sqrt(Ld) to keep the average power 
%   independent of the number of scatterers.
%   The path-loss exponents follow the log-distance model with log-normal shadowing.  
%   See the reference: Secrecy Capacity Maximization for a Hybrid Relay-RIS Scheme in mmWave MIMO Networks 
%   https://arxiv.org/pdf/2205.13904; https://ieeexplore.ieee.org/document/9860450
% 
%
%   Example usage:
%
%     [G,H] = generate_SV_channel_extended(8,4,8,user_loc,[0,0],[50,0]);
%
% Hesam Nezhadmohammad
% last revised: 01/25/2026

% Heights = [h_BS,h_RIS,h_U], specifying the heights (meters) of the BS, RIS and users (default [6,3,1.5])
Heights = [6,3,1.5];
h_BS  = Heights(1);
h_RIS = Heights(2);
h_U   = Heights(3);

% Number of users
K = size(user_loc,1);

% Number of RIS elements
M = Mx*My;

% Path-loss parameters
a_LOS      = 61.4;
b_LOS      = 2.0;
sigma_LOS  = 5.8;
a_NLOS     = 72.0;
b_NLOS     = 2.92;
sigma_NLOS = 8.7;

% Multipath parameters
L  = 6; % total number of propagation paths per link (default 6) 
Ld = L - 1;       % number of NLOS scatterers

% Rician K-factor (power ratio between LOS and NLOS, set to 0 as default)
K_dB     = 0;
% K_dB     = 13.2; % TVT20 paper used this!
Krician  = 10^(K_dB/10);

% Antenna and element gains (set to 1 for simplicity)
lambda_gain = 1;

% Compute vertical height differences
dv_BR = abs(h_BS - h_RIS);
dv_RU = abs(h_RIS - h_U);

% Compute horizontal distances and angles for BS->RIS link
dx_BR = panel_loc(1) - BS_loc(1);
dy_BR = panel_loc(2) - BS_loc(2);
d_BR_horiz = sqrt(dx_BR^2 + dy_BR^2);

% Azimuth angle from RIS to BS (this sets the RIS broadside orientation)
theta_az_BR = atan2(BS_loc(2) - panel_loc(2), BS_loc(1) - panel_loc(1));

% Elevation angle of arrival at RIS from BS (negative since RIS is below BS)
% Use a small epsilon to avoid division by zero
el_BR = -atan(dv_BR / max(d_BR_horiz, 1e-6));

% BS steering angle (angle of departure from BS to RIS)
phi_BS = atan2(panel_loc(2) - BS_loc(2), panel_loc(1) - BS_loc(1));

% Generate array responses (unit-norm) for LOS path
a_t = array_response_ula(N, phi_BS);               % BS transmit array response
a_r = array_response_ura(Mx, My, theta_az_BR, el_BR);  % RIS receive array response
% a_r = array_response_UPA(theta_az_BR,el_BR,Mx*My,4); % UPA version, commento above and uncomment this to use it.
G_LOS = a_r * a_t';  % M×N rank-one LOS component

% Generate NLOS components for BS->RIS link
G_NLOS_sum = zeros(M,N);
for ell = 1:Ld
    % Random azimuth and elevation spreads (uniform azimuth, small elevation spread)
    az = 2*pi*rand();
    el = (randn()*5)*pi/180;  % ±5 degrees std dev
    ar = array_response_ura(Mx, My, az, el);
    % ar = array_response_UPA(az,el,Mx*My,4);
    at = array_response_ula(N, 2*pi*rand());
    alpha = (randn() + 1j*randn())/sqrt(2);
    G_NLOS_sum = G_NLOS_sum + alpha * (ar * at');
end
% Normalize NLOS sum by sqrt(Ld) to maintain unit average power per scatterer
G_NLOS = G_NLOS_sum / sqrt(Ld);

% Path-loss and shadowing for BS->RIS link
% For the BS->RIS link, LOS and NLOS components experience different
% large-scale fading according to the mmWave path-loss model, reference: Secrecy Capacity Maximization for a Hybrid Relay-RIS Scheme in mmWave MIMO Networks (https://arxiv.org/pdf/2205.13904; https://ieeexplore.ieee.org/document/9860450)
% Draw independent shadowing terms for LOS and NLOS components
xi_LOS_BR  = sigma_LOS  * randn();
xi_NLOS_BR = sigma_NLOS * randn();
% Compute large-scale path-loss for LOS and NLOS links using log-distance
PL_LOS_BR  = 10^(-(a_LOS  + 10*b_LOS  * log10(d_BR_horiz) + xi_LOS_BR )/10);
PL_NLOS_BR = 10^(-(a_NLOS + 10*b_NLOS * log10(d_BR_horiz) + xi_NLOS_BR)/10);

% Compose BS->RIS channel: scale LOS and NLOS terms separately
H_BS_RIS = lambda_gain * ( sqrt(PL_LOS_BR)  * sqrt(Krician/(Krician+1)) * sqrt(M*N) * G_LOS + ...
                           sqrt(PL_NLOS_BR) * sqrt(1/(Krician+1))   * sqrt(M*N) * G_NLOS );

% Generate RIS->User links
H_RIS_user = zeros(M, K);
for k = 1:K
    % Horizontal distance from panel to user k
    dx_RU = user_loc(k,1) - panel_loc(1);
    dy_RU = user_loc(k,2) - panel_loc(2);
    d_RU_horiz = sqrt(dx_RU^2 + dy_RU^2);
    % Azimuth angle from RIS to user (for the LOS steering)
    az_ru = atan2(user_loc(k,2) - panel_loc(2), user_loc(k,1) - panel_loc(1));
    % Elevation angle (negative because RIS is above user)
    el_RU = -atan(dv_RU / max(d_RU_horiz, 1e-6));
    % LOS steering
    h_LOS = array_response_ura(Mx, My, az_ru, el_RU);
    % h_LOS = array_response_UPA(az_ru,el_RU,My*My,4);
    % Generate NLOS contributions
    h_NLOS_sum = zeros(M,1);
    for ell = 1:Ld
        az = 2*pi*rand();
        el = (randn()*5)*pi/180;
        at = array_response_ura(Mx, My, az, el);
        % at = array_response_UPA(az,el,Mx*My,4);
        alpha = (randn() + 1j*randn())/sqrt(2);
        h_NLOS_sum = h_NLOS_sum + alpha * at;
    end
    % Normalize NLOS cluster sum
    h_NLOS = h_NLOS_sum / sqrt(Ld);
    % Path-loss and shadowing for RIS->User link
    % Each user experiences distinct LOS and NLOS large-scale fading.
    % Draw independent shadowing for LOS and NLOS components.
    xi_LOS_RU  = sigma_LOS  * randn();
    xi_NLOS_RU = sigma_NLOS * randn();
    % Compute path-loss for LOS and NLOS parts
    PL_LOS_RU  = 10^(-(a_LOS  + 10*b_LOS  * log10(d_RU_horiz) + xi_LOS_RU )/10);
    PL_NLOS_RU = 10^(-(a_NLOS + 10*b_NLOS * log10(d_RU_horiz) + xi_NLOS_RU)/10);
    % Compose RIS->User channel: scale LOS and NLOS components separately
    H_RIS_user(:,k) = lambda_gain * ( sqrt(PL_LOS_RU)  * sqrt(Krician/(Krician+1)) * sqrt(M) * h_LOS + ...
                                      sqrt(PL_NLOS_RU) * sqrt(1/(Krician+1))   * sqrt(M) * h_NLOS );
end

end
