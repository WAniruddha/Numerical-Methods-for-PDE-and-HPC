% Clear previous data
clear; clc; close all;

% Load data
analytical_data = load('Analytical_Solution.dat');
numerical_100 = load('IVP2_UPWIND_T1000_X100_CFL1.dat');
numerical_400 = load('IVP2_UPWIND_T1000_X400_CFL1.dat');
numerical_2000 = load('IVP2_UPWIND_T1000_X2000_CFL1.dat');

% Extract x and f values
x_analytical = analytical_data(:,1);
f_analytical = analytical_data(:,2);

x_num_100 = numerical_100(:,1);
f_num_100 = numerical_100(:,2);

x_num_400 = numerical_400(:,1);
f_num_400 = numerical_400(:,2);

x_num_2000 = numerical_2000(:,1);
f_num_1600 = numerical_2000(:,2);

% Interpolate analytical solution at numerical points
f_analytical_interp_100 = interp1(x_analytical, f_analytical, x_num_100, 'linear', 'extrap');
f_analytical_interp_400 = interp1(x_analytical, f_analytical, x_num_400, 'linear', 'extrap');
f_analytical_interp_1600 = interp1(x_analytical, f_analytical, x_num_2000, 'linear', 'extrap');

% Reduce the number of plotted points
skip_100 = ceil(length(x_num_100) / 100);
skip_400 = ceil(length(x_num_400) / 400);
skip_1600 = ceil(length(x_num_2000) / 1600);

% Define colors
color_100 = '#D95319';  % Orange
color_400 = '#EDB120';  % Yellow
color_1600 = '#0072BD'; % Blue

% Define line style
lineStyle = ':';  % Dotted line

% Define line width
lineWidth = 1.5;

% Create figure and adjust width
figure;
set(gcf, 'Position', [100, 100, 1300, 600]); % Adjust figure width

% -----------------------------------
% Subplot 1: X100 Grid
% -----------------------------------
ax1 = subplot('Position', [0.05, 0.15, 0.28, 0.7]); % Custom width spacing
plot(x_analytical, f_analytical, 'k-', 'LineWidth', lineWidth); hold on;
plot(x_num_100(1:skip_100:end), f_num_100(1:skip_100:end), lineStyle, 'Color', color_100, 'LineWidth', lineWidth);
xlabel('x', 'FontSize', 14);
ylabel('f', 'FontSize', 14);
title('X100', 'FontSize', 14);
legend('Analytical', 'Numerical (X100)', 'Location', 'northeast', 'FontSize', 12);
grid on;
text(min(x_analytical), max(f_analytical), '(a)', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

% Add Zoomed Plot for Subplot 1
zp1 = BaseZoom();
zp1.run;

% -----------------------------------
% Subplot 2: X400 Grid
% -----------------------------------
ax2 = subplot('Position', [0.38, 0.15, 0.28, 0.7]); % Custom width spacing
plot(x_analytical, f_analytical, 'k-', 'LineWidth', lineWidth); hold on;
plot(x_num_400(1:skip_400:end), f_num_400(1:skip_400:end), lineStyle, 'Color', color_400, 'LineWidth', lineWidth);
xlabel('x', 'FontSize', 14);
ylabel('f', 'FontSize', 14);
title('X400', 'FontSize', 14);
legend('Analytical', 'Numerical (X400)', 'Location', 'northeast', 'FontSize', 12);
grid on;
text(min(x_analytical), max(f_analytical), '(b)', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

% Add Zoomed Plot for Subplot 2
zp2 = BaseZoom();
zp2.run;

% -----------------------------------
% Subplot 3: X2000 Grid
% -----------------------------------
ax3 = subplot('Position', [0.71, 0.15, 0.28, 0.7]); % Custom width spacing
plot(x_analytical, f_analytical, 'k-', 'LineWidth', lineWidth); hold on;
plot(x_num_2000(1:skip_1600:end), f_num_1600(1:skip_1600:end), lineStyle, 'Color', color_1600, 'LineWidth', lineWidth);
xlabel('x', 'FontSize', 14);
ylabel('f', 'FontSize', 14);
title('X2000', 'FontSize', 14);
legend('Analytical', 'Numerical (X2000)', 'Location', 'northeast', 'FontSize', 12);
grid on;
text(min(x_analytical), max(f_analytical), '(c)', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

% Add Zoomed Plot for Subplot 3
zp3 = BaseZoom();
zp3.run;

% Title for the entire figure
sgtitle('IVP2 UPWIND T=1000 CFL=1', 'FontSize', 16);
