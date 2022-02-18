% f02_220124_atheta.m
%
% ------------------
% Created: 2022-01-24 17:30
% Author: Cory
% Title: Finding A, theta
% Description:
%     Code for finding the variables Xi, theta that align with the numeric data
% ------------------
% 
close all

global s Gamma Xi theta dat
s = 6^(1/4); 
Gamma = 24;
dat = readtable('PureQuartic.csv');
Xi = 1.12051; 
XiLimits = [1.12, 1.121];
theta = 1.16200;
thetaLimits = [1.16, 1.165];

u0 = @(x) Xi * ( cos(theta) * cos(s .* x) ./ cosh(s .* x) + sin(theta) * sin(s .* x) ./ sinh(s .* x) ); 

f = figure('color', 'white', 'Position', [10 10 1000 600]);
ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);

subplot('Position', [0.1 0.39 0.35 0.54])
semilogy(dat.Var1, dat.Var2, 'LineWidth', 2);
hold on 
p1 = semilogy(dat.Var1, u0(dat.Var1), 'LineWidth', 2);
hold off
xlim([-14, -11.5])
xlabel('Time (ps)')
ylabel('Amplitude')
set(gca, 'FontSize', 14)

subplot('Position', [0.55 0.39 0.4 0.54])
p2 = semilogy(dat.Var1, dat.Var2 ./ u0(dat.Var1), 'LineWidth', 2);
yline(1, 'k--', 'LineWidth', 2);
xlim([-13.5, -12])
xlabel('Time (ps)')
ylabel('Relative Amplitude')
title(sprintf('A = %.7f, B = %.7f', Xi * cos(theta), Xi * sin(theta)));
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
    global s Gamma Xi theta dat
    Xi = Xi_;
    theta = theta_;
    u0 = @(x) Xi * ( cos(theta) * cos(s .* x) ./ cosh(s .* x) + sin(theta) * sin(s .* x) ./ sinh(s .* x) ); 
    p1.YData = u0(dat.Var1);
    title(sprintf('A = %.7f, B = %.7f', Xi * cos(theta), Xi * sin(theta)));
    p2.YData = dat.Var2 ./ u0(dat.Var1);
    XiLabel.String = sprintf('Xi: %.5f', Xi);
    thetaLabel.String = sprintf('theta: %.5f',theta);
    drawnow
    % disp([theta Xi])
end