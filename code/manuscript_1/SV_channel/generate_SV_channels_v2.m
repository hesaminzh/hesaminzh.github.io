function [H_BS_user, H_BS_RIS, H_RIS_user] = generate_SV_channels_v2(M, Nx, Ny, K, d_BR_horiz, d_BU_horiz)
  % heights
  h_BS = 6;  h_RIS = 3;  h_U = 1.5;
  dv_BR = abs(h_BS - h_RIS);
  dv_RU = abs(h_RIS - h_U);
  dv_BU = abs(h_BS - h_U);

  N = Nx * Ny;
  % --- SV parameters (same as before) ---
  a_LOS = 61.4;   b_LOS = 2.0;   sigma_LOS = 5.8;
  a_NLOS = 72.0;  b_NLOS = 2.92; sigma_NLOS = 8.7;
  K_dB   = 13.2;  Krician = 10^(K_dB/10);
  lambda_gain_dBi = 9.60;
  lambda_gain = 10^(lambda_gain_dBi/10);    % <-- bring back your antenna gain

  %% 1) BS→RIS
  d_BR = sqrt(d_BR_horiz^2 + dv_BR^2);
  xi_LOS = sigma_LOS * randn();
  PL_BI = 10^(-(a_LOS + 10*b_LOS*log10(d_BR) + xi_LOS)/10);

  % array responses (assume your array_response_* functions already normalize by sqrt(#elements))
  a_t = array_response_ula(M, 0);  % BS broadside
  el  = -atan(dv_BR / max(d_BR_horiz,1e-6));
  a_r = array_response_ura(Nx, Ny, 0, el);  % RIS broadside toward BS

  G_LOS  = a_r * a_t';                            % N×M
  G_NLOS = (randn(N,M) + 1j*randn(N,M))/sqrt(2);  % un-normalized Rayleigh

  % **note: no 1/sqrt(N*M) on the NLOS term**
  H_BS_RIS = lambda_gain * sqrt(PL_BI) * ( ...
        sqrt(Krician/(Krician+1)) * sqrt(N*M) * G_LOS + ...
        sqrt(1/(Krician+1))      * sqrt(N*M)       * G_NLOS );

  %% 2) RIS→User
  H_RIS_user = zeros(N, K);
  for k = 1:K
    d_RU = sqrt((d_BU_horiz(k)-d_BR_horiz)^2 + dv_RU^2);
    xi_LOS = sigma_LOS * randn();
    PL_IU = 10^(-(a_LOS + 10*b_LOS*log10(d_RU) + xi_LOS)/10);

    if d_BU_horiz(k) <= d_BR_horiz
      az = 0;
    else
      az = pi;
    end
    el = -atan(dv_RU / max(abs(d_BU_horiz(k)-d_BR_horiz),1e-6));
    a_t = array_response_ura(Nx, Ny, az, el);

    h_NLOS = (randn(N,1) + 1j*randn(N,1))/sqrt(2);
    % **again no 1/sqrt(N) on the NLOS term**
    H_RIS_user(:,k) = lambda_gain * sqrt(PL_IU) * ( ...
                 sqrt(Krician/(Krician+1)) * sqrt(N) * a_t + ...
                 sqrt(1/(Krician+1))      * sqrt(N) * h_NLOS );
  end

  %% 3) BS→User (unchanged SV NLOS cluster)
  Ld = 3;
  H_BS_user = zeros(M, K);
  for k = 1:K
    d_BU = sqrt(d_BU_horiz(k)^2 + dv_BU^2);
    xi_NLOS = sigma_NLOS * randn();
    PL_BU = 10^(-(a_NLOS + 10*b_NLOS*log10(d_BU) + xi_NLOS)/10);

    for ell = 1:Ld
      alpha = lambda_gain * sqrt(PL_BU) * (randn+1j*randn)/sqrt(2*Ld);
      phi   = 2*pi*rand();
      a_t   = array_response_ula(M, phi);
      H_BS_user(:,k) = H_BS_user(:,k) + sqrt(M) * alpha * a_t;
    end
  end
end
