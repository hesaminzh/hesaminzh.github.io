% existing data
%eta_list = [η1, η2, η3, η4, η5];   % non-uniform
%out_FC   = [...];  % your 1×5 vector
% same for out_GC, out_SC, out_No

% choose a uniform η grid (still 5 points here)
eta_uniform = linspace(eta_list(1), eta_list(end), 15);

% pchip interpolation for a smooth match
wsrF_uni = interp1(eta_list, WSR_avg(:,1),   eta_uniform, 'pchip');
wsrG_uni = interp1(eta_list, WSR_avg(:,2),   eta_uniform, 'pchip');
wsrS_uni = interp1(eta_list, WSR_avg(:,3),   eta_uniform, 'pchip');
wsrN_uni = interp1(eta_list, out_No,   eta_uniform, 'pchip');

% plot
figure; hold on; grid on;
plot(eta_uniform, wsrF_uni, '-o','LineWidth',1.5,'DisplayName','FC');
plot(eta_uniform, wsrG_uni, '-s','LineWidth',1.5,'DisplayName','GC');
plot(eta_uniform, wsrS_uni, '-d','LineWidth',1.5,'DisplayName','SC');
%plot(eta_uniform, wsrN_uni, '-^','LineWidth',1.5,'DisplayName','NoRIS');
xlabel('\eta (blockage density)')
ylabel('Outage Probability')
title('Outage vs. \eta (interpolated)')
legend('Location','northwest')




outF_uni = interp1(eta_list, Outage(:,1),   eta_uniform, 'pchip');
outG_uni = interp1(eta_list, Outage(:,2),   eta_uniform, 'pchip');
outS_uni = interp1(eta_list, Outage(:,3),   eta_uniform, 'pchip');
%outN_uni = interp1(eta_list, out_No,   eta_uniform, 'pchip');

% plot
figure; hold on; grid on;
plot(eta_uniform, outF_uni, '-o','LineWidth',1.5,'DisplayName','FC');
plot(eta_uniform, outG_uni, '-s','LineWidth',1.5,'DisplayName','GC');
plot(eta_uniform, outS_uni, '-d','LineWidth',1.5,'DisplayName','SC');
plot(eta_uniform, outN_uni, '-^','LineWidth',1.5,'DisplayName','NoRIS');
xlabel('\eta (blockage density)')
ylabel('Outage Probability')
title('Outage vs. \eta (interpolated)')
legend('Location','northwest')

