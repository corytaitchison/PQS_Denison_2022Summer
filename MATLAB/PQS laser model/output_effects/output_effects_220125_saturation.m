% output_effects_220125_saturation.m
%
% ------------------
% Created: 2022-01-22 12:30
% Author: Cory
% Title: Saturation
% Description:
%     Determine how many round trips are needed for saturation to be reached (the wavelength no longer changes)
% ------------------
% 

% ------------------
% beta coeff for WS, for quadratic dispersion 
beta2_WS = 0;               % beta2  (s^2/m) 
beta3_WS = -0.12e-39;       % beta3 with opposite sign for Waveshaper (s^3/m)
beta4_WS = 2.2e-54;         % beta4 with opposite sign for Waveshape (s^4/m)
% ------------------

% ------------------
freqMethod = "MULT";       % Frequency shift method
freqdelta = 5e9;          % Frequency shift per round trip (Hz)
feedback = 0.1;
SEED = 222;                % Random seed
% ------------------

% ------------------
peaklambdas = []; 
indices = [];
% ------------------

close all;
initialise;
E = E0;

figure(1)
set(gcf, 'color', 'white')
p = plot(0, 0, 'k.-', 'MarkerSize', 14, 'LineWidth', 2);
xlabel('Trip')
ylabel('Peak Wavelength (nm)')
set(gca, 'FontSize', 14)

mmm = 0;
while true
    mmm = mmm + 1;

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
    EnP(mmm) = (1-feedback)*sum(abs(E).^2)*dt; % Energy
    figure(2)
    plotplot;
    SMF;
    E = E*sqrt(connector);
    %  ------------------

    if mmm >= 5
        spect = 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2)); 
        [~, ind] = max(spect);
        peaklambdas = [peaklambdas lambdanm(ind)];  % In nm
        indices = [5:mmm];
        %  ------------------    
        figure(1)
        p.XData = indices;
        p.YData = peaklambdas;
        drawnow
    else 
        disp(mmm);
    end
end

return

titlestring = ['02__output_effects__' datestr(now, 'yyyymmddHHMM') '__saturation']

figure('color', 'white')
plot(indices, c ./ peaklambdas, '.-', 'LineWidth', 2, 'Markersize', 14)
xlabel('Trip')
ylabel('Peak Frequency (GHz)')
title({sprintf('Applied shift: %.2f GHz, OC: %.2f', freqdelta * 1e-9, feedback), titlestring})
saveas(gcf, ['data/' titlestring '.fig'])

save(['data/' titlestring '.mat'], 'peaklambdas', 'feedback', 'indices', 'freqdelta');