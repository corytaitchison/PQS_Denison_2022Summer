% output_effects__220131__grid.m
%
% ------------------
% Created: 2022-01-31 16:30
% Author: Cory
% Title: Output Effects - Grid Sweep
% Description:
%     Sweep a 2D grid of applied frequency shifts and feedback
%     to determine where the critical points are for creating 
%     stable solitons
% ------------------
% 

% ------------------                        % beta coeff for WS, for quadratic dispersion 
beta2_WS = 0;                               % beta2  (s^2/m) 
beta3_WS = -0.12e-39;                       % beta3 with opposite sign for Waveshaper (s^3/m)
beta4_WS = 2.2e-54;                         % beta4 with opposite sign for Waveshape (s^4/m)
% ------------------

% ------------------
freqMethod = "MULT";                        % Frequency shift method
maxroundtrips = 250;                        % Max number of loops for the laser 
SEED = 222;                                 % RNG Seed
N = 11;                                     % Grid size 
% ------------------

% ------------------
feedbacks = linspace(0.3, 0.8, N);          % Percent energy retained each trip
freqdeltas = logspace(8, 10, N);            % Frequency shift (Hz) per trip
peaklambdas = zeros(N, N, ...
    maxroundtrips)*NaN;                     % Peak wavelength (nm)
% ------------------

% ------------------
titlestring = ['02__output_effects__' datestr(now, 'yyyymmddHHMM') '__grid']
% ------------------

close all;

for i = 1:N 

    for j = 1:N 

        tic

        initialise;
        E = E0;

        feedback = feedbacks(i);
        freqdelta = freqdeltas(j); 

        runningcount = 0; 

        for mmm = 1:maxroundtrips 

            %  ------------------ Laser Loop
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
            %  ------------------ /Laser Loop

            %  ------------------
            spect = 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2)); 
            [~, ind] = max(spect);
            %  ------------------    

            peaklambdas(i, j, mmm) = lambdanm(ind);  % In nm

            if mmm > 50
                if std(peaklambdas(i,j, (mmm-20):(mmm-1))) < 1e-2
                    break 
                end
            end

        end

        save(['data/' titlestring '.mat'], 'peaklambdas', 'feedbacks', 'freqdeltas')
        disp([i, j, mmm]);
        toc

    end
end

return 

% --- Zero data ---
load('data/02__output_effects__202202021029__grid.mat')
baselambdas = zeros(N, 1);

for i = 1:N  
    for mmm = 50:MMM 
        if isnan(peaklambdas(i,1,mmm))
            baselambdas(i) = peaklambdas(i, 1, mmm-1);
            break
        end
    end 
end 

return 

% --- Actual data ---
load('data/02__output_effects__202201311758__grid.mat') 
% load('data/02__output_effects__202202011634__grid.mat')
[N, ~, MMM] = size(peaklambdas); 

peaklambdasfixed = peaklambdas;

for i = 1:N 
    for j = 1:N 
        for mmm = 50:MMM 
            if isnan(peaklambdas(i,j,mmm))
                peaklambdasfixed(i, j, mmm:end) = peaklambdas(i, j, mmm-1);
                break
            end
        end 
    end 
end 

return 

close all
figure('color', 'white')
s = surf(log10(freqdeltas), feedbacks, peaklambdasfixed(:, :, 1) - baselambdas(), 'EdgeColor', 'flat')
xlim(log10(freqdeltas([1, end])))
ylim(feedbacks([1, end]))
view(2)
colorbar
colormap jet
% caxis([1550, 1562])
caxis([-5, 0])
ylabel('Feedback')
xlabel('Log10(Frequency Shift)')
title('Wavelength shift (nm); Trip: 1')
set(gca, 'FontSize', 14)

for mmm = 2:MMM 
    s.ZData = peaklambdasfixed(:, :, mmm) - baselambdas;
    title(sprintf('Wavelength shift (nm); Trip: %d', mmm))
    drawnow 
    pause(0.1)
end

return 

% Peak wavelengths

didterminate = isnan(peaklambdas(:, :, end));

peaklambdasfinal = peaklambdasfixed(:, :, end); 
peaklambdasfinal(~didterminate) = NaN;

titlestring = ['02__output_effects__' datestr(now, 'yyyymmddHHMM') '__grid'];

figure('color', 'white')
s = surf(log10(freqdeltas), feedbacks, peaklambdasfinal, 'EdgeColor', 'flat')
xlim(log10(freqdeltas([1, end])))
ylim(feedbacks([1, end]))
view(2)
colorbar
colormap jet
caxis([1555, 1562])
ylabel('Feedback')
xlabel('Log10(Frequency Shift)')
title({'Peak Wavelengths after 250 trips', 'White = did not saturate', titlestring})
set(gca, 'FontSize', 14)

saveas(gcf, ['data/' titlestring '.fig'])

return 

% Peak wavelength shift

didterminate = isnan(peaklambdas(:, :, end));

peaklambdasfinal = peaklambdasfixed(:, :, end) - baselambdas; 
peaklambdasfinal(~didterminate) = NaN;

titlestring = ['02__output_effects__' datestr(now, 'yyyymmddHHMM') '__grid'];

figure('color', 'white')
s = surf(log10(freqdeltas), feedbacks, peaklambdasfinal, 'EdgeColor', 'flat')
xlim(log10(freqdeltas([1, end])))
ylim(feedbacks([1, end]))
view(2)
colorbar
colormap jet
caxis([-3, 0])
ylabel('Feedback')
xlabel('Log10(Frequency Shift, Hz)')
title({'Peak wavelength shift (nm) after 250 trips', 'White = did not saturate', titlestring})
set(gca, 'FontSize', 14)

saveas(gcf, ['data/' titlestring '.fig'])