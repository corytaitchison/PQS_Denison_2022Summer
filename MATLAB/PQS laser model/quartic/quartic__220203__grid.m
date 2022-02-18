% quartic__220203__grid.m
%
% ------------------
% Created: 2022-02-03 13:00
% Author: Cory
% Title: Quartic Grid
% Description:
%     Plot the evolution of the wavelengths for quartic dispersion
% ------------------
% 

% --- Zero data ---
load('data/03__quartic__202202031125__grid.mat')
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
load('data/03__quartic__202202022257__grid.mat')
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
s = surf(log10(freqdeltas), feedbacks, peaklambdasfixed(:, :, 1) - baselambdas, 'EdgeColor', 'flat')
xlim(log10(freqdeltas([1, end])))
ylim(feedbacks([1, end]))
view(2)
colorbar
colormap jet
% caxis([1550, 1562])
caxis([-1, 0])
ylabel('Feedback')
xlabel('Log10(Frequency Shift, Hz)')
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

titlestring = ['03__quartic__' datestr(now, 'yyyymmddHHMM') '__grid'];

figure('color', 'white')
s = surf(log10(freqdeltas), feedbacks, peaklambdasfinal, 'EdgeColor', 'flat')
xlim(log10(freqdeltas([1, end])))
ylim(feedbacks([1, end]))
view(2)
colorbar
colormap jet
caxis([-1, 0])
ylabel('Feedback')
xlabel('Log10(Frequency Shift, Hz)')
title({'Peak wavelength shift (nm) after 250 trips', 'White = did not saturate', titlestring})
set(gca, 'FontSize', 14)

saveas(gcf, ['data/' titlestring '.fig'])