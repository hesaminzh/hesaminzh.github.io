Mx = 8;    % TX antenna nums
My = 4;    % TX antenna nums
K = 2;     % user nums
L = 2;     % NLOS path nums
%SNR = 30; % SNR in dB
M = Mx*My;

%% surface parameters
Nx = 4;                      % x-axis surface element nums
%Ny = [1,2,4,8,16];   % y-axis surface element nums
%Ny = [1,8];   % y-axis surface element nums
Ny = 1;
N  = Nx*Ny;

para.BS_loc = [0,0]; % location of BS
para.RIS_loc = [60,0]; % location of RIS in m
para.user_center = [70, 0];
para.user_range = [0, 10]; % range of user location from RIS in m
[user_loc, user_azimuths, d_BR, d_RU, d_BU] = generate_user_location(para,K);

%% previous channels (worked!)

H_BS_RIS_v1 = BS_IRS(Nx, Ny, Mx, My, L, d_BR);
H_RIS_user_v1 = FF_User(Nx, Ny, K, d_RU);
H_BS_user_v1 = BS_UE_v2(Mx, My, K, L, d_BU);


[H_BS_user_v2, H_BS_RIS_v2, H_RIS_user_v2] = generate_SV_channels(M, Nx, Ny, K, d_BR, d_BU);
[H_BS_user_v3, H_BS_RIS_v3, H_RIS_user_v3] = generate_SV_channels_v2(M, Nx, Ny, K, d_BR, d_BU);
[H_BS_user_v4, H_BS_RIS_v4, H_RIS_user_v4] = generate_SV_channels_v3(M, Nx, Ny, K, d_BR, d_BU);
[H_BS_user_v5, H_BS_RIS_v5, H_RIS_user_v5] = generate_SV_channels_v4(M, Nx, Ny, K, d_BR, d_RU, d_BU, user_azimuths);


