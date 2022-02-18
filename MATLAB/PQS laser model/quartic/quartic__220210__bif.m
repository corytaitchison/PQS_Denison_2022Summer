% quartic__220210__bif.m
%
% ------------------
% Created: 2022-02-10 1130
% Author: Cory
% Title: Quartic Bifurcation diagrams
% Description:
%     Create bifurcation diagrams for the new quartic data 
% ------------------
% 

% --- Zero data ---
load('data/03__quartic__202202062259__zero.mat')
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
load('data/03__quartic__202202062351__grid.mat')
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
endindex = min(startindex + 100, MMM);

titlestring = ['03__quartic__' datestr(now, 'yyyymmddHHMM') '__bif'];
figure('color', 'white')

subplot(211)
feedbackindex = 1; 
semilogx(freqdeltas, reshape(peaklambdasfixed(feedbackindex, :, startindex:endindex), [N, endindex - startindex + 1]) - baselambdas(feedbackindex), 'k.', 'MarkerSize', 12)
% ylim([-1, -0.002])
xlabel('Frequency shift (Hz)')
ylabel('Peak wavelength shift (nm)')
title({'Peak wavelengths', titlestring, sprintf('feedback = %.2f', feedbacks(feedbackindex))})
set(gca, 'Fontsize', 14)

subplot(212)
feedbackindex = 11; 
semilogx(freqdeltas, reshape(peaklambdasfixed(feedbackindex, :, startindex:endindex), [N, endindex - startindex + 1]) - baselambdas(feedbackindex), 'k.', 'MarkerSize', 12)
% ylim([-1, -0.002])
xlabel('Frequency shift (Hz)')
ylabel('Peak wavelength shift (nm)')
title({sprintf('feedback = %.2f', feedbacks(feedbackindex))})
set(gca, 'Fontsize', 14)

saveas(gcf, ['data/' titlestring '.fig'])

return

% Bifurcation diagrams

startindex = min(whenfinished(whenfinished ~= 0)) - 10;
endindex = min(startindex + 100, MMM);

titlestring = ['03__quartic__' datestr(now, 'yyyymmddHHMM') '__bif'];
figure('color', 'white')

% newcolours = zeros(N, 3);
% newcolours(:, 2) = linspace(0, 0.8, N);
% colororder(newcolours);

for feedbackindex = 1:2:N 
    semilogx(freqdeltas, reshape(peaklambdasfixed(feedbackindex, :, startindex:endindex), [N, endindex - startindex + 1]) - baselambdas(feedbackindex), '.', 'MarkerSize', 12)
    hold on 
end 
hold off 
% ylim([-1, -0.002])
xlabel('Frequency shift (Hz)')
ylabel('Peak wavelength shift (nm)')
title({'Peak wavelengths', titlestring, sprintf('feedback = %.2f', feedbacks(feedbackindex))})
% legend('0.3', '0.35', '0.4', '0.45', '0.5', '0.55', '0.6', '0.65', '0.7', '0.75', '0.8', 'Location', 'southwest')
legend('0.3',  '0.4', '0.5', '0.6', '0.7', '0.8', 'Location', 'southwest')
set(gca, 'Fontsize', 14)

saveas(gcf, ['data/' titlestring '.fig'])
