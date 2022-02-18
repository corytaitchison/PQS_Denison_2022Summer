% f01_220124_compare.m
%
% ------------------
% Created: 2022-01-22 2045
% Author: Cory
% Title: Comparison
% Description:
%     Comparison of the fit of u0 and u0 + u1 with the 
%     numeric data, using the new functions (Logbook_03) 
% ------------------
%  

dat = readtable('PureQuartic.csv');

s = 6^(1/4); 
Gamma = -24; %%%
phi = Gamma / (1280 * s^4);
% A = 0.4384; B = 1.0170; 
A = 0.4454085; B = 1.0281799;
cs = [-2, -9, -6, -9, 0.25, 0; ...
    9, 6, 9, -6, 0, -0.25; ...
    -6, -9, 6, -9, 0.25, 0; ...
    9, -6, 9, -2, 0, -0.25] * ...
    [A^3 * phi; A^2 * B * phi; A * B^2 * phi; B^3 * phi; ...
    A; B]; 

% ds = [0.0903392, -0.128515, 0.0271901, 0.00377054, -0.0630546, 0.132237];
% ds = [0.11333, -0.104456, 0.0700616, 0.0315229, -0.043174, 0.135931];

u0 = @(x) A * cos(s .* x) ./ cosh(s .* x) + B * sin(s .* x) ./ sinh(s .* x); 
u1 = @(xs) arrayfun(@(x) sum(arrayfun(@(k) cs(k+1) * (cos(s.*x) ./ cosh(s .* x)) .^ (3-k) .* (sin(s.*x) ./ sinh(s .* x)) .^ (k), 0:3)), xs);
% u2 = @(xs) arrayfun(@(x) sum(arrayfun(@(k) ds(k+1) * (cos(s.*x) ./ cosh(s .* x)) .^ (5-k) .* (sin(s.*x) ./ sinh(s .* x)) .^ (k), 0:5)), xs);

% plot(dat.Var1, dat.Var2);
% hold on

figure('color', 'white');
plot(dat.Var1, dat.Var2, 'k', 'LineWidth', 2);
hold on 
plot(dat.Var1, u0(dat.Var1), 'LineWidth', 2);
plot(dat.Var1, u0(dat.Var1) + u1(dat.Var1), 'LineWidth', 2);
% plot(dat.Var1, u0(dat.Var1) + u1(dat.Var1) + u2(dat.Var1), 'LineWidth', 2);
% plot(dat.Var1, xxx, 'LineWidth', 2);
hold off 
xlim([-2, 2])
legend('Data', 'u0', 'u0 + u1');
% legend('u0', 'u0 + u1 (new)', 'Data', 'u0 + u1 (old)');
xlabel('Time (ps)')
ylabel('Amplitude')
set(gca, 'FontSize', 16);

figure('color','white');
semilogy(dat.Var1, dat.Var2, 'k', 'LineWidth', 2);
hold on 
semilogy(dat.Var1, u0(dat.Var1), 'LineWidth', 2);
semilogy(dat.Var1, u0(dat.Var1) + u1(dat.Var1), 'LineWidth', 2);
% semilogy(dat.Var1, u0(dat.Var1) + u1(dat.Var1) + u2(dat.Var1), 'LineWidth', 2);
semilogy(dat.Var1, xxx, 'LineWidth', 2);
hold off 
% xlim([-1.27, -1.25])
xlim([-5.3, -5.25])
% legend('Data', 'u0', 'u0 + u1');
legend('u0', 'u0 + u1 (new)', 'Data', 'u0 + u1 (old)');
xlabel('Time (ps)')
ylabel('Amplitude')
set(gca, 'FontSize', 16);
