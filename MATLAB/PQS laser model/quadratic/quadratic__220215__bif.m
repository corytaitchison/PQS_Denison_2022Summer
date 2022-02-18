% quadratic__220215__bif.m
%
% ------------------
% Created: 2022-02-15 14:30
% Author: Cory
% Title: Quadratic Bifurcation diagrams
% Description:
%     Create bifurcation diagrams for the quadratic data
% ------------------
% 

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
whenfinished = zeros(N, N); 

for i = 1:N 
    for j = 1:N 
        for mmm = 50:MMM 
            if isnan(peaklambdas(i,j,mmm))
                peaklambdasfixed(i, j, mmm:end) = peaklambdas(i, j, mmm-1);
                whenfinished(i, j) = mmm; 
                break
            end
        end 
    end 
end 

return 

% Bifurcation diagrams

startindex = min(whenfinished(whenfinished ~= 0)) - 10;
% startindex = 175;
endindex = min(startindex + 100, MMM);

titlestring = ['04__quadratic__' datestr(now, 'yyyymmddHHMM') '__bif'];
figure('color', 'white')

subplot(211)
feedbackindex = 2; 
semilogx(freqdeltas, reshape(peaklambdasfixed(feedbackindex, :, startindex:endindex), [N, endindex - startindex + 1]) - baselambdas(feedbackindex), 'k.', 'MarkerSize', 12)
% ylim([-1, -0.002])
xlabel('Frequency shift (Hz)')
ylabel('Peak wavelength shift (nm)')
title({'Peak wavelengths', titlestring, sprintf('feedback = %.2f', feedbacks(feedbackindex))})
set(gca, 'Fontsize', 14)

subplot(212)
feedbackindex = 3; 
semilogx(freqdeltas, reshape(peaklambdasfixed(feedbackindex, :, startindex:endindex), [N, endindex - startindex + 1]) - baselambdas(feedbackindex), 'k.', 'MarkerSize', 12)
% ylim([-1, -0.002])
xlabel('Frequency shift (Hz)')
ylabel('Peak wavelength shift (nm)')
title({sprintf('feedback = %.2f', feedbacks(feedbackindex))})
set(gca, 'Fontsize', 14)

saveas(gcf, ['data/' titlestring '.fig'])

return
