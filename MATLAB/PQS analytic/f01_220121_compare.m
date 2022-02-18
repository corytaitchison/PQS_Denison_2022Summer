dat = readtable('PureQuartic.csv');

s = 6^(1/4); 
Gamma = -24;
A = 0.4384;
B = 1.0170; 
phi = Gamma / (1280 * s ^ 4);
cs = [9 -6 9 -2; -6 -9 6 -9; 9 6 9 -6; -2 -9 -6 -9] * [A^3; A^2 * B; A * B^2; B^3] * phi; 

u0 = @(x) A * cos(s .* x) ./ cosh(s .* x) + B * sin(s .* x) ./ sinh(s .* x); 
u1 = @(xs) arrayfun(@(x) sum(arrayfun(@(k) cs(k+1) * (cos(s.*x) ./ cosh(s .* x)) .^ (3-k) .* (sin(s.*x) ./ sinh(s .* x)) .^ (k), 0:3)), xs);

% plot(dat.Var1, dat.Var2);
% hold on

figure('color', 'white');
plot(dat.Var1, u0(dat.Var1), 'LineWidth', 2);
hold on 
plot(dat.Var1, u0(dat.Var1) + u1(dat.Var1), 'LineWidth', 2);
plot(dat.Var1, dat.Var2, 'LineWidth', 2);
hold off 
xlim([-2, 2])
legend('u0', 'u0 + u1', 'Data');
xlabel('Time (ps)')
ylabel('Amplitude')
set(gca, 'FontSize', 16);

figure('color','white');
semilogy(dat.Var1, u0(dat.Var1), 'LineWidth', 2);
hold on 
semilogy(dat.Var1, u0(dat.Var1) + u1(dat.Var1), 'LineWidth', 2);
semilogy(dat.Var1, dat.Var2, 'LineWidth', 2);
hold off 
xlim([-1.27, -1.25])
legend('u0', 'u0 + u1', 'Data');
xlabel('Time (ps)')
ylabel('Amplitude')
set(gca, 'FontSize', 16);
