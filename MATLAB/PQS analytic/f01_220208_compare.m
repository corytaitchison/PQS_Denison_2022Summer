% f01_220208_compare.m
%
% ------------------
% Created: 2022-02-08 20:30
% Author: Cory
% Title: Comparison
% Description:
%     Interactive comparison of the numerical data and our ansatz
% ------------------
% 

close all

% ------------------
global dat s Gamma Xi theta phi cs ds es u0 un

dat = readtable('PureQuartic.csv');

s = 6^(1/4); 
Gamma = -24; 
% A = 0.4384; B = 1.0170; 
% A = 0.44446666; B = 1.0295334;
Xi = 1.474; 
XiLimits = [1.3, 1.7];
theta = 0.98948;
thetaLimits = [0.9, 1.1];
% ------------------

% ------------------
phi = Gamma / (1280 * s^4);
cs = @(A, B) [-2, -9, -6, -9, 0.25, 0; ...
    9, 6, 9, -6, 0, -0.25; ...
    -6, -9, 6, -9, 0.25, 0; ...
    9, -6, 9, -2, 0, -0.25] * ...
    [A^3 * phi; A^2 * B * phi; A * B^2 * phi; B^3 * phi; ...
    A; B]; 

ds = @(A, B) [A/8 + ((-95*A^5 + 86*A^4*B + 98*A^3*B^2 + 352*A^2*B^3 + 273*A*B^4 + 266*B^5)*Gamma^2)/(21299200*s^8) - (3 * (2*A^3 + 9*A^2*B + 6*A*B^2 + 9*B^3)*Gamma)/(5120*s^4); ...
    -(B/8) - ((86*A^5 + 671*A^4*B + 712*A^3*B^2 + 798*A^2*B^3 + 626*A*B^4 - 273*B^5)*Gamma^2)/(21299200*s^8) + (3 * (3*A^3 + 2*A^2*B + 3*A*B^2 - 2*B^3)*Gamma)/(5120*s^4); ...
    (Gamma * ((49*A^5 + 356*A^4*B + 50*A^3*B^2 + 532*A^2*B^3 - 399*A*B^4 + 176*B^5)*Gamma - 12480*B*(3*A^2 + 4*A*B + 3*B^2)*s^4))/(10649600*s^8); ...
    (Gamma * (-((176*A^5 + 399*A^4*B + 532*A^3*B^2 - 50*A^2*B^3 + 356*A*B^4 - 49*B^5)*Gamma) - 12480*A*(3*A^2 - 4*A*B + 3*B^2)*s^4))/(10649600*s^8); 
    -(A/8) + ((273*A^5 + 626*A^4*B - 798*A^3*B^2 + 712*A^2*B^3 - 671*A*B^4 + 86*B^5)*Gamma^2)/(21299200*s^8) + (3 * (2*A^3 + 3*A^2*B - 2*A*B^2 + 3*B^3)*Gamma)/(5120*s^4); ...
    B/8 - ((266*A^5 - 273*A^4*B + 352*A^3*B^2 - 98*A^2*B^3 + 86*A*B^4 + 95*B^5)*Gamma^2)/(21299200*s^8) + (3 * (-9*A^3 + 6*A^2*B - 9*A*B^2 + 2*B^3)*Gamma)/(5120*s^4)];

es = @(A, B) [(5*A)/ 64 + ((5032*A^7 + 436015*A^6*B - 346656*A^5*B^2 + 168605*A^4*B^3 - 1428120*A^3*B^4 - 1383459*A^2*B^5 - 1203728*A*B^6 - 1116049*B^7)*Gamma^3)/(30125588480000*s^12) + ((-95*A^5 + 86*A^4*B + 98*A^3*B^2 + 352*A^2*B^3 + 273*A*B^4 + 266*B^5)*Gamma^2)/(17039360*s^8) - (9*(2*A^3 + 9*A^2*B + 6*A*B^2 + 9*B^3)*Gamma)/(20480*s^4); ...
    -((5*B)/64) + ((-436015*A^7 + 728536*A^6*B + 2110275*A^5*B^2 + 3979200*A^4*B^3 + 7591715*A^3*B^4 + 2938008*A^2*B^5 + 5045425*A*B^6 - 1203728*B^7)*Gamma^3)/(30125588480000*s^12) - (3*(86*A^5 + 671*A^4*B + 712*A^3*B^2 + 798*A^2*B^3 + 626*A*B^4 - 273*B^5)*Gamma^2)/(85196800*s^8) + (3*(3*A^3 + 2*A^2*B + 3*A*B^2 - 2*B^3)*Gamma)/(20480*s^4); ... 
	-(A/64) - (3*(115552*A^7 + 703425*A^6*B + 1665496*A^5*B^2 + 3105875*A^4*B^3 + 1461680*A^3*B^4 + 2863603*A^2*B^5 - 979336*A*B^6 + 461153*B^7)*Gamma^3)/(30125588480000*s^12) + ((-377*A^5 + 1142*A^4*B + 590*A^3*B^2 + 2824*A^2*B^3 + 567*A*B^4 + 1682*B^5)*Gamma^2)/(85196800*s^8) + (3*(2*A^3 - 9*A^2*B - 18*A*B^2 - 9*B^3)*Gamma)/(20480*s^4); ... 
	B/64 + ((-33721*A^7 + 795840*A^6*B + 1863525*A^5*B^2 + 1095528*A^4*B^3 + 3415589*A^3*B^4 - 877008*A^2*B^5 + 1518343*A*B^6 - 285624*B^7)*Gamma^3)/(6025117696000*s^12) + ((94*A^5 - 1215*A^4*B - 1072*A^3*B^2 - 2494*A^2*B^3 - 1166*A*B^4 + 721*B^5)*Gamma^2)/(85196800*s^8) + (3*(-15*A^3 + 2*A^2*B - 15*A*B^2 + 6*B^3)*Gamma)/(20480*s^4); ... 
	 -(A/64) - ((285624*A^7 + 1518343*A^6*B + 877008*A^5*B^2 + 3415589*A^4*B^3 - 1095528*A^3*B^4 + 1863525*A^2*B^5 - 795840*A*B^6 - 33721*B^7)*Gamma^3)/(6025117696000*s^12) + ((-721*A^5 - 1166*A^4*B + 2494*A^3*B^2 - 1072*A^2*B^3 + 1215*A*B^4 + 94*B^5)*Gamma^2)/(85196800*s^8) + (3*(6*A^3 + 15*A^2*B + 2*A*B^2 + 15*B^3)*Gamma)/(20480*s^4); ... 
	B/64 + (3*(461153*A^7 + 979336*A^6*B + 2863603*A^5*B^2 - 1461680*A^4*B^3 + 3105875*A^3*B^4 - 1665496*A^2*B^5 + 703425*A*B^6 - 115552*B^7)*Gamma^3)/(30125588480000*s^12) + ((1682*A^5 - 567*A^4*B + 2824*A^3*B^2 - 590*A^2*B^3 + 1142*A*B^4 + 377*B^5)*Gamma^2)/(85196800*s^8) + (3*(9*A^3 - 18*A^2*B + 9*A*B^2 + 2*B^3)*Gamma)/(20480*s^4); ... 
	(5*A)/64 + ((-1203728*A^7 - 5045425*A^6*B + 2938008*A^5*B^2 - 7591715*A^4*B^3 + 3979200*A^3*B^4 - 2110275*A^2*B^5 + 728536*A*B^6 + 436015*B^7)*Gamma^3)/(30125588480000*s^12) - (3*(273*A^5 + 626*A^4*B - 798*A^3*B^2 + 712*A^2*B^3 - 671*A*B^4 + 86*B^5)*Gamma^2)/(85196800*s^8) - (3*(2*A^3 + 3*A^2*B - 2*A*B^2 + 3*B^3)*Gamma)/(20480*s^4); ... 
	-((5*B)/64) + ((1116049*A^7 - 1203728*A^6*B + 1383459*A^5*B^2 - 1428120*A^4*B^3 - 168605*A^3*B^4 - 346656*A^2*B^5 - 436015*A*B^6 + 5032*B^7)*Gamma^3)/(30125588480000*s^12) + ((266*A^5 - 273*A^4*B + 352*A^3*B^2 - 98*A^2*B^3 + 86*A*B^4 + 95*B^5)*Gamma^2)/(17039360*s^8) + (9*(9*A^3 - 6*A^2*B + 9*A*B^2 - 2*B^3)*Gamma)/(20480*s^4)];

u0 = @(x, A, B) A * cos(s .* x) ./ cosh(s .* x) + B * sin(s .* x) ./ sinh(s .* x); 
un = @(xs, n, coefs) arrayfun(@(x) sum(arrayfun(@(k) coefs(k+1) * (cos(s.*x) ./ cosh(s .* x)) .^ ((2*n+1)-k) .* (sin(s.*x) ./ sinh(s .* x)) .^ (k), 0:(2*n+1))), xs);
% ------------------

f = figure('color', 'white', 'Position', [10 10 1000 600]);
ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);

% y = u0(dat.Var1, Xi * cos(theta)^2, Xi * sin(theta)^2) + un(dat.Var1, 1, cs(Xi * cos(theta)^2, Xi * sin(theta)^2)) + ...
    % un(dat.Var1, 2, ds(Xi * cos(theta)^2, Xi * sin(theta)^2)) + un(dat.Var1, 3, es(Xi * cos(theta)^2, Xi * sin(theta)^2)); 
y = u0(dat.Var1, Xi * cos(theta)^2, Xi * sin(theta)^2);

subplot('Position', [0.1 0.39 0.35 0.54])
semilogy(dat.Var1, dat.Var2, 'LineWidth', 2);
hold on 
p1 = semilogy(dat.Var1, y, 'LineWidth', 2);
hold off
xlim([-14, -11.5])
xlabel('Time (ps)')
ylabel('Amplitude')
set(gca, 'FontSize', 14)

subplot('Position', [0.55 0.39 0.4 0.54])
p2 = semilogy(dat.Var1, dat.Var2 ./ y, 'LineWidth', 2);
yline(1, 'k--', 'LineWidth', 2);
xlim([-13.5, -12])
xlabel('Time (ps)')
ylabel('Relative Amplitude')
title(sprintf('A = %.7f, B = %.7f', Xi * cos(theta)^2, Xi * sin(theta)^2));
set(gca, 'FontSize', 14)

XiSlider = uicontrol('Parent', f, 'Style', 'slider', 'Position', [260, 115, 419, 23],...
    'value', Xi, 'min', XiLimits(1), 'max', XiLimits(2));
XiLabel = uicontrol('Parent', f, 'Style', 'text', 'Position', [450,150,100,23],...
    'String', sprintf('Xi: %.5f', Xi), 'FontSize', 16, 'BackgroundColor', 'white', ...
    'FontWeight', 'bold');

thetaSlider = uicontrol('Parent', f, 'Style', 'slider', 'Position', [260, 35, 419, 23],...
    'value', theta, 'min', thetaLimits(1), 'max', thetaLimits(2));
thetaLabel = uicontrol('Parent', f, 'Style', 'text', 'Position', [420,70,160,23],...
    'String', sprintf('theta: %.5f', theta), 'FontSize', 16, 'BackgroundColor', 'white', ...
    'FontWeight', 'bold');

XiSlider.Callback = @(es, ed) callback(es.Value, thetaSlider.Value, p1, p2, XiLabel, thetaLabel);
thetaSlider.Callback = @(es, ed) callback(XiSlider.Value, es.Value, p1, p2, XiLabel, thetaLabel);

function callback(Xi_, theta_, p1, p2, XiLabel, thetaLabel) 
    global dat s Gamma Xi theta phi cs ds es u0 un
    Xi = Xi_;
    theta = theta_;
    % y = u0(dat.Var1, Xi * cos(theta)^2, Xi * sin(theta)^2) + un(dat.Var1, 1, cs(Xi * cos(theta)^2, Xi * sin(theta)^2)) + ...
        % un(dat.Var1, 2, ds(Xi * cos(theta)^2, Xi * sin(theta)^2)) + un(dat.Var1, 3, es(Xi * cos(theta)^2, Xi * sin(theta)^2)); 
    y = u0(dat.Var1, Xi * cos(theta)^2, Xi * sin(theta)^2);
    p1.YData = y;
    title(sprintf('A = %.7f, B = %.7f', Xi * cos(theta)^2, Xi * sin(theta)^2));
    p2.YData = dat.Var2 ./ y;
    XiLabel.String = sprintf('Xi: %.5f', Xi);
    thetaLabel.String = sprintf('theta: %.5f',theta);
    drawnow
end
