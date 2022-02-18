% output_effects_220121.m
%
% ------------------
% Created: 2022-01-26 11:00pm
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
roundtrips = 250;           % Number of loops for the laser 
% ------------------

% ------------------
N = 13;
freqdeltas = linspace(-5e9, 5e9, N);  % Frequency shift per round trip (in Hz)
M = 3;
feedbacks = [0.25, 0.5, 0.75];    % percent output coupler
peaklambdas = zeros(M, N);   
% ------------------

figure('color', 'white');

titlestring = ['02__output_effects__' datestr(now, 'yyyymmddHHMM')];

for m = 1:M 

   feedback = feedbacks(m);

   for n = 1:N 

      close all;

      tic

      initialise;
      freqdelta = freqdeltas(n);

      laserLoop;

      % Get peak
      spect = 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2)); 
      [~, ind] = max(spect);
      peaklambdas(m, n) = lambdanm(ind);  % In nm

      toc
      save(['data/' titlestring '.mat'], 'freqdeltas', 'peaklambdas', 'roundtrips', 'feedbacks', 'lambda0', 'f0')
      disp([m, n])

   end

   lambda0 = mean(peaklambdas(m, :)) * 1e-9;
   f0 = c ./ lambda0;
   lambdadeltas = c ./ (f0 + freqdeltas) - lambda0;

   subplot(211)
   hold on 
   plot(lambdadeltas * 1e9, peaklambdas(m, :) - lambda0 * 1e9, 'LineWidth', 2); 
   hold off 

   subplot(212) 
   hold on 
   plot(freqdeltas * 1e-9, c ./ peaklambdas(m, :) - c ./ lambda0 * 1e-9, 'LineWidth', 2); 
   hold off 

   disp(m)
   drawnow
end

subplot(211) 
hold on 
plot(lambdadeltas * 1e9, lambdadeltas * 1e9 * roundtrips, 'k--', 'LineWidth', 2);
hold off 
xticks([-5:2:5] * 1e-2)
xlabel('Applied wavelength shift (nm)')
ylabel('Output shift (nm)')
title({sprintf('Output shift. Trips: %d', roundtrips), ... 
   titlestring})
legend('OC: 0.25', 'OC: 0.50', 'OC: 0.75', 'y = 30x', 'Location', 'SouthEast')
set(gca, 'FontSize', 14)

subplot(212) 
hold on 
plot(freqdeltas * 1e-9, freqdeltas * 1e-9 * roundtrips, 'k--', 'LineWidth', 2); 
hold off 
xticks([-5:1:5])
xlabel('Applied frequency shift (GHz)')
ylabel('Output shift (GHz)')
set(gca, 'FontSize', 14)

% ------------------
saveas(gcf, ['data/' titlestring '.fig'])
save(['data/' titlestring '.mat'], 'freqdeltas', 'peaklambdas', 'roundtrips', 'feedbacks', 'lambda0', 'f0')
% ------------------

return 

% ------------------
% ------------------
dat = open 'data/02__output_effects__202201262341.mat';
titlestring = '02__output_effects__202201262341';

figure('color', 'white');

for (m = 1:length(dat.feedbacks))

   lambda0 = mean(dat.peaklambdas(m, :)) * 1e-9;
   f0 = c ./ lambda0;
   lambdadeltas = c ./ (f0 + dat.freqdeltas) - lambda0;

   subplot(211)
   hold on 
   plot(lambdadeltas * 1e9, dat.peaklambdas(m, :) - lambda0 * 1e9, '.-', 'LineWidth', 2, 'Markersize', 14); 
   hold off 

   subplot(212) 
   hold on 
   plot(dat.freqdeltas * 1e-9, c ./ dat.peaklambdas(m, :) - c ./ lambda0 * 1e-9, '.-', 'LineWidth', 2, 'Markersize', 14); 
   hold off 

end 

subplot(211) 
hold on 
plot(lambdadeltas * 1e9, lambdadeltas * 1e9 * dat.roundtrips, 'k--', 'LineWidth', 2);
hold off 
xticks([-5:2:5] * 1e-2)
xlabel('Applied wavelength shift (nm)')
ylabel('Output shift (nm)')
title({sprintf('Output shift. Trips: %d', dat.roundtrips), ... 
   titlestring})
legend('OC: 0.25', 'OC: 0.50', 'OC: 0.75', sprintf('y = %dx', dat.roundtrips), 'Location', 'SouthEast')
set(gca, 'FontSize', 14)

subplot(212) 
hold on 
plot(dat.freqdeltas * 1e-9, dat.freqdeltas * 1e-9 * dat.roundtrips, 'k--', 'LineWidth', 2); 
hold off 
xticks([-5:1:5])
xlabel('Applied frequency shift (GHz)')
ylabel('Output shift (GHz)')
set(gca, 'FontSize', 14)

% ------------------
saveas(gcf, ['data/' titlestring '.fig'])
% save(['data/' titlestring '.mat'], 'freqdeltas', 'peaklambdas', 'roundtrips', 'feedbacks', 'lambda0', 'f0')
% ------------------