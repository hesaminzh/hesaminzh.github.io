close all; clear all; clc;
% Value for blockage density (eta) such that blockage probability would vary between 0 to 1
%eta_list = [0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.09];
%eta_list = 0:0.1:1;
d_BU = [41.2286; 44.1635];
d_RU = [34.1976; 30.0230];

%eta_max = -log(0.1)/max(d_BU);
%eta_list = linspace(0, eta_max, 9);

p_grid   = 0:0.1:0.9;                % desired blockage probs
d_BU_max = max(d_BU);
eta_list = -log(1 - p_grid)/d_BU_max;
p_BU = 1 - exp(-eta_list .* d_BU);
p_RU = 1 - exp(-eta_list .* d_RU);

exp_BU = eta_list.*d_BU; % exponent in blockage probability of BS-User link without minus sign
exp_RU = eta_list.*d_RU; % exponent in blockage probability of RIS-User link without minus sign

p_BU = exp(-exp_BU); % LoS probability in BS-User link
q_BU = 1-p_BU; % Blockage probability in BS-User link

p_RU = exp(-exp_RU); % LoS probability in RIS-User link
q_RU = 1-p_RU; % Blockage probability in RIS-User link

disp(q_BU);
disp(q_RU);







%%%%%%%% both paths having a reasonable blockage probabilities
%% 1) Desired window of RIS–user blockage probabilities
p_min = 0.1;   % 10 % blocked
p_max = 0.9;   % 90 % blocked

%% 2) Your RIS–user distances
%d_RU = [10.1161,  2.9741];   % in metres
d_RU = [19, 18];
d_BU = [69, 68];

%% 3) Invert p = 1−exp(−η⋅d) to get η_low and η_high
d_RU_max = max(d_RU);  % the longest RIS–user path
d_RU_min = min(d_RU);  % the shortest RIS–user path

eta_low  = -log(1 - p_min) / d_RU_max;  % ensures at η_low, p_RU≈0.1
eta_high = -log(1 - p_max) / d_RU_min;  % ensures at η_high, p_RU≈0.9

%% 4) Build your η grid
numPoints = 9;                                   % however many you like
eta_list  = linspace(eta_low, eta_high, numPoints);

%% 5) Check the actual RIS–user blockage sweep
p_RU = 1 - exp(-eta_list' * d_RU');
p_BU = 1 - exp(-eta_list' * d_BU');
disp([ eta_list'  p_RU ]);
