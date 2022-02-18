% quartic__220208__grid.m
%
% ------------------
% Created: 2022-02-08 23:00
% Author: Cory
% Title: Quartic - Grid Sweep
% Description:
%     Sweep a 2D grid of applied frequency shifts and feedback
%     using the new E0 from quadratic__220208__zero.m
% ------------------
% 

% ------------------                        % beta coeff for WS, for quadratic dispersion 
beta2_WS = 0;                               % beta2  (s^2/m) 
beta3_WS = -0.12e-39;                       % beta3 with opposite sign for Waveshaper (s^3/m)
beta4_WS = 2.2e-54;                         % beta4 with opposite sign for Waveshape (s^4/m)
% ------------------

% ------------------
freqMethod = 'MULT';                        % Frequency shift method
maxroundtrips = 250;                        % Max number of loops for the laser 
SEED = 222;                                 % RNG Seed
N = 11;                                     % Grid size 
% ------------------

load('data/04__quadratic__202202082154__zero.mat')

% ------------------
feedbacks = linspace(0.3, 0.8, N);          % Percent energy retained each trip
freqdeltas = logspace(8, 10, N);            % Frequency shift (Hz) per trip
peaklambdas = zeros(N, N, ...
    maxroundtrips)*NaN;                     % Peak wavelength (nm)
peakamps = peaklambdas;                     % Peak amplitude
% ------------------

% ------------------
titlestring = ['04__quadratic__' datestr(now, 'yyyymmddHHMM') '__grid']
% ------------------

close all;

for i = 1:N 

    for j = 1:N

        tic

        initialise;
        E = E0s(i, :);

        feedback = feedbacks(i);
        freqdelta = freqdeltas(j); 

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
%             EnP(mmm) = (1-feedback)*sum(abs(E).^2)*dt;
%             plotplot
            SMF;
            E = E*sqrt(connector);
            %  ------------------ /Laser Loop

            %  ------------------
            spect = abs(fftshift(ifft(fftshift(E)))).^2; 
            [val, ind] = max(spect);
            %  ------------------    

            peaklambdas(i, j, mmm) = lambdanm(ind);  % In nm
            peakamps(i, j, mmm) = val; 

            if mmm > 50
                if std(peaklambdas(i, j, (mmm-20):(mmm-1))) < 1e-2
                    break 
                end
            end

        end

        save(['data/' titlestring '.mat'], 'peaklambdas', 'peakamps', 'feedbacks', 'freqdeltas')
        disp([i, j, mmm]);
        toc

    end

end

return 


% --- Zero data ---
load('data/04__quadratic__202202082154__zero.mat')
[N, ~, MMM] = size(peaklambdas); 
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
load('data/04__quadratic__202202082343__grid.mat')
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

% Peak wavelength shift

didterminate = isnan(peaklambdas(:, :, end));

peaklambdasfinal = peaklambdasfixed(:, :, end) - baselambdas; 
peaklambdasfinal(~didterminate) = NaN;

titlestring = ['04__quadratic__' datestr(now, 'yyyymmddHHMM') '__grid'];

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
title({'Peak wavelength shift (nm)', 'White = did not saturate', titlestring})
set(gca, 'FontSize', 14)

saveas(gcf, ['data/' titlestring '.fig'])

return 

% Bifurcation diagrams

% trips = [5 10 15 20 25 30];
trips = [4 8 12 16 20 24];

titlestring = ['04__quadratic__' datestr(now, 'yyyymmddHHMM') '__bif'];

newcolours = zeros(length(trips), 3);
newcolours(:, 2) = linspace(0, 0.8, length(trips));

figure('color', 'white')

subplot(211)
% feedbackindex = 8; 
feedbackindex = N;
colororder(newcolours);
semilogx(freqdeltas, reshape(peaklambdasfixed(feedbackindex, :, trips), [N, length(trips)]) - baselambdas(feedbackindex), '.-', 'LineWidth', 2, 'Markersize', 14)
% ylim([-1, -0.002])
xlabel('Frequency shift (Hz)')
ylabel('Peak wavelength shift (nm)')
% legend('5', '10', '15', '20', '25', '30', 'Location', 'southwest')
legend('4', '8', '12', '16', '20', '24', 'Location', 'southwest')
title({'Peak wavelengths at various round trips', titlestring, sprintf('feedback = %.2f', feedbacks(feedbackindex))})
set(gca, 'Fontsize', 14)

subplot(212)
% feedbackindex = 4; 
feedbackindex = 1; 
colororder(newcolours);
semilogx(freqdeltas, reshape(peaklambdasfixed(feedbackindex, :, trips), [N, length(trips)]) - baselambdas(feedbackindex), '.-', 'LineWidth', 2, 'Markersize', 14)
% ylim([-1, -0.002])
xlabel('Frequency shift (Hz)')
ylabel('Peak wavelength shift (nm)')
% legend('5', '10', '15', '20', '25', '30', 'Location', 'southwest')
legend('4', '8', '12', '16', '20', '24', 'Location', 'southwest')
title(sprintf('feedback = %.2f', feedbacks(feedbackindex)))
set(gca, 'Fontsize', 14)

saveas(gcf, ['data/' titlestring '.fig'])

return