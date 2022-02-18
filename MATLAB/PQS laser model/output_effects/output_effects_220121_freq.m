% output_effects_220121.m
%
% ------------------
% Created: 2022-01-21 11:37am
% Author: Cory
% Title: Output Effects (Frequency)
% Description:
%   Parametric sweep of frequency shift to determine the effect on the output of a PQS laser (SMF + WS).
%  Adapted from code by A. F. J. Runge
% ------------------

% ------------------
% beta coeff for WS, for quadratic dispersion 
beta2_WS = 0;               % beta2  (s^2/m) 
beta3_WS = -0.12e-39;       % beta3 with opposite sign for Waveshaper (s^3/m)
beta4_WS = 2.2e-54;         % beta4 with opposite sign for Waveshape (s^4/m)
% ------------------

% ------------------
freqMethod = "MULT";       % Frequency shift method
roundtrips = 30;           % Number of loops for the laser 
feedback = 0.3;            % percent output coupler
% ------------------

% ------------------
N = 20;
freqdeltas = linspace(-5e9, 5e9, N);  % Frequency shift per round trip (in Hz)
peaklambdas = zeros(1, N);   
% ------------------

for i = 1:N 

   close all;

   initialise;
   freqdelta = freqdeltas(i);

   laserLoop;

   % Get peak
   spect = 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2)); 
   [~, ind] = max(spect);
   peaklambdas(i) = lambdanm(ind);  % In nm

   disp({i, peaklambdas(i)});

end

% ------------------
lambda0 = mean(peaklambdas) * 1e-9;
f0 = c / lambda0;
titlestring = ['02__output_effects__' datestr(now, 'yyyymmddHHMM')]
% ------------------

% ------------------
figure('color', 'white');
lambdadeltas = c ./ (f0 + freqdeltas) - lambda0; 
subplot(211)
plot(lambdadeltas * 1e9, peaklambdas - lambda0 * 1e9, 'LineWidth', 2); 
xticks([-5:2:5] * 1e-2)
xlabel('Applied wavelength shift (nm)')
ylabel('Output shift (nm)')
title({sprintf('Output shift. Trips: %d, OC: %.2f', roundtrips, feedback), ... 
   titlestring})
set(gca, 'FontSize', 14)

subplot(212)
plot(freqdeltas * 1e-9, c ./ peaklambdas - c ./ lambda0 * 1e-9, 'LineWidth', 2); 
xticks([-5:1:5])
xlabel('Applied frequency shift (GHz)')
ylabel('Output shift (GHz)')
set(gca, 'FontSize', 14)
% ------------------

% ------------------
saveas(gcf, ['data/' titlestring '.fig'])
save(['data/' titlestring '.mat'], 'freqdeltas', 'peaklambdas', 'roundtrips', 'feedback', 'lambda0', 'f0')
% ------------------

return

% ------------------
% ------------------
% ------------------
titlestring = ['02__output_effects__' datestr(now, 'yyyymmddHHMM')];
close all;
open('data/02__output_effects__202201241210.fig')
subplot(211)

hold on 
load('data/02__output_effects__202201241116.mat')
plot(lambdadeltas * 1e9, peaklambdas - lambda0 * 1e9, 'LineWidth', 2); 
load('data/02__output_effects__202201241139.mat')
plot(lambdadeltas * 1e9, peaklambdas - lambda0 * 1e9, 'LineWidth', 2); 
% y = 30x
plot(lambdadeltas * 1e9, lambdadeltas * 1e9 * roundtrips, 'k--', 'LineWidth', 2);
hold off 

title({'Output shift. Trips: 30', titlestring})
legend('OC: 0.3', 'OC: 0.5', 'OC: 0.7', 'y = 30x', 'Location', 'SouthEast')
subplot(212)

hold on 
load('data/02__output_effects__202201241116.mat')
plot(freqdeltas * 1e-9, c ./ peaklambdas - c ./ lambda0 * 1e-9, 'LineWidth', 2); 
load('data/02__output_effects__202201241139.mat')
plot(freqdeltas * 1e-9, c ./ peaklambdas - c ./ lambda0 * 1e-9, 'LineWidth', 2); 
% y = 30x
plot(freqdeltas * 1e-9, freqdeltas * 1e-9 * roundtrips, 'k--', 'LineWidth', 2); 
hold off 

set(gca, 'FontSize', 14)

saveas(gcf, ['data/' titlestring '.fig'])