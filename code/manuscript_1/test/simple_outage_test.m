%nMonte   = 500;           % number of trials
nMonte   = 1000;           % number of trials
mu       = 7.5;           % mean rate (bps/Hz)
sigma    = 1.0;           % standard deviation
wsr_all  = mu + sigma*randn(nMonte,1);   % Gaussian scatter around 7.5

tau_vec  = linspace(0,15,150);           % thresholds from 0 to 15
outageP  = arrayfun(@(t) mean(wsr_all < t), tau_vec);

figure();
plot(tau_vec, outageP, 'LineWidth',1.5);
xlabel('$\tau\;(bps/Hz)$', 'Interpreter', 'LaTeX');
ylabel('$P_{\rm out}(\tau)$', 'Interpreter', 'LaTex');
grid on;
