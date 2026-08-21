color1 = [0, 91, 150] / 255;   % (0,91,150)
color2 = [0, 159, 136] / 255;  % (0,159,136)
color3 = [208, 44, 129] / 255; % (208,44,129)

x = 0:0.1:10;
y1 = sin(x);
y2 = cos(x);
y3 = sin(x).*cos(x);

% color1 = [0, 91, 150] / 255;
% color2 = [0, 159, 136] / 255;
% color3 = [208, 44, 129] / 255;

figure;
plot(x, y1, 'Color', color1, 'LineWidth', 2); hold on;
plot(x, y2, 'Color', color2, 'LineWidth', 2);
plot(x, y3, 'Color', color3, 'LineWidth', 2);
legend('sin(x)', 'cos(x)', 'sin(x)*cos(x)');
xlabel('x'); ylabel('y');
grid on;
title('Plot with Custom RGB Colors');
