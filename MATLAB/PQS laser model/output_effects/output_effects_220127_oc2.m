% output_effects_220121.m
%
% ------------------
% Created: 2022-01-27 11:50pm
% Author: Cory
% Title: Output Effects (Output Coupling)
% Description:
%   Parametric sweep of output coupling to determine the effect on the output of a PQS laser (SMF + WS).
%   Stores the peak wavelength at multiple points during the round trip to create an animation.
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
peaklambdas = zeros(roundtrips, N); 
% lambda0s = zeros(1, N);   
load('data/02__output_effects__202201270725__oc.mat', 'lambda0s');
% ------------------

% ------------------
titlestring = ['02__output_effects__' datestr(now, 'yyyymmddHHMM') '__oc']
% ------------------

close all;

for i = 1:N 
   tic

   initialise;
   E = E0;

   feedback = feedbacks(i);
   freqdelta = FREQDELTA; 

   for mmm = 1:roundtrips 

       %  ------------------
       SMF;
       gain_PQS;
       SMF;
       E = E*sqrt(connector);
       for kk = 1:4
       SMF;
       end
       SA;
       for kk = 1:4
       SMF;
       end
       E = E*sqrt(connector);
       for kk = 1:4
       SMF;
       end
       E = E*sqrt(lWS);
       WaveShaper
       FreqShift
       for kk = 1:4
       SMF;
       end
       E = E*sqrt(connector);
       SMF;
       OC;
       SMF;
       E = E*sqrt(connector);
       %  ------------------

      spect = 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2)); 
      [~, ind] = max(spect);
      peaklambdas(mmm, i) = lambdanm(ind);  % In nm
      %  ------------------    

   end

   % freqdelta = 0;
   % initialise;
   % laserLoop;
   % spect = 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2)); 
   % [~, ind] = max(spect);
   % lambda0s(i) = lambdanm(ind);  % In nm

   toc
   disp(i);

   save(['data/' titlestring '.mat'], 'peaklambdas', 'roundtrips', 'feedback', 'lambda0s')

end

return

% ------------------
figure('color', 'white');
subplot(211)
p1 = plot(feedbacks * 100, peaklambdas(1, :) - lambda0s, '.-', 'MarkerSize', 20); 
xlabel('Energy feedback per trip (%)')
ylabel('Output shift (nm)')
title({sprintf('Output shift. Trip: %d/%d, Applied shift: %.2f GHz', 1, roundtrips, FREQDELTA * 1e-9)})
set(gca, 'FontSize', 14)
ylim([-5, -1])
xlim([10, 90])

subplot(212)
p2 = plot(feedbacks * 100, c ./ peaklambdas(1, :) - c ./ lambda0s, '.-', 'MarkerSize', 20); 
xlabel('Energy feedback per trip (%)')
ylabel('Output shift (GHz)')
set(gca, 'FontSize', 14)
ylim([200, 600])
xlim([10, 90])

for i = 1:roundtrips 

   subplot(211)
   p1.YData = peaklambdas(i, :) - lambda0s;
   title({sprintf('Output shift. Trip: %d/%d, Applied shift: %.2f GHz', i, roundtrips, FREQDELTA * 1e-9)})
   subplot(212)
   p2.YData =  c ./ peaklambdas(i, :) - c ./ lambda0s; 

   drawnow
end 
% ------------------

% ------------------
saveas(gcf, ['data/' titlestring '.fig'])
save(['data/' titlestring '.mat'], 'peaklambdas', 'roundtrips', 'feedback', 'lambda0s')
% ------------------
