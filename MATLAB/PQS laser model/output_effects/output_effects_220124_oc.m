% output_effects_220121.m
%
% ------------------
% Created: 2022-01-21 11:37am
% Author: Cory
% Title: Output Effects (Output Coupling)
% Description:
%   Parametric sweep of output coupling to determine the effect on the output of a PQS laser (SMF + WS).
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
roundtrips = 250;           % Number of loops for the laser 
FREQDELTA = 5e9;          % Frequency shift per round trip (Hz)
% ------------------

% ------------------
N = 20;
feedbacks = linspace(0.1, 0.9, N);
peaklambdas = zeros(1, N); 
lambda0s = zeros(1, N);   
% ------------------

% ------------------
titlestring = ['02__output_effects__' datestr(now, 'yyyymmddHHMM') '__oc']
% ------------------

for i = 1:N 
   tic
   close all;

   initialise;
   feedback = feedbacks(i);

   freqdelta = FREQDELTA; 
   laserLoop;
   spect = 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2)); 
   [~, ind] = max(spect);
   peaklambdas(i) = lambdanm(ind);  % In nm

   freqdelta = 0;
   laserLoop;
   spect = 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2)); 
   [~, ind] = max(spect);
   lambda0s(i) = lambdanm(ind);  % In nm

   toc
   disp({i, peaklambdas(i), lambda0s(i)});

   save(['data/' titlestring '.mat'], 'FREQDELTA', 'peaklambdas', 'roundtrips', 'feedback', 'lambda0s')

end

% ------------------
figure('color', 'white');
subplot(211)
plot(feedbacks * 100, peaklambdas - lambda0s, '.', 'MarkerSize', 20); 
% xticks([-5:2:5] * 1e-2)
xlabel('Energy feedback per trip (%)')
ylabel('Output shift (nm)')
title({sprintf('Output shift. Trips: %d, Applied shift: %.2f GHz', roundtrips, FREQDELTA * 1e-9), ... 
   titlestring})
set(gca, 'FontSize', 14)

subplot(212)
plot(feedbacks * 100, c ./ peaklambdas - c ./ lambda0s, '.', 'MarkerSize', 20); 
% xticks([-5:1:5])
xlabel('Energy feedback per trip (%)')
ylabel('Output shift (GHz)')
set(gca, 'FontSize', 14)
% ------------------

% ------------------
saveas(gcf, ['data/' titlestring '.fig'])
save(['data/' titlestring '.mat'], 'FREQDELTA', 'peaklambdas', 'roundtrips', 'feedback', 'lambda0s')
% ------------------
