% quartic__220215__instfreq.m
%
% ------------------
% Created: 2022-02-15 1500
% Author: Cory
% Title: Quartic - Instantaneous Frequency
% Description:
%     Compute the electric field for a given quartic soliton and
%     then plot the instantaneous frequency
% ------------------
% 

% ------------------
% beta coeff for WS
beta2_WS = 21.4e-27;        % beta2 with opposite sign for Waveshaper (s^2/m) 
beta3_WS = -0.12e-39;       % beta3 with opposite sign for Waveshaper (s^3/m)
beta4_WS = -80e-51;         % beta4 for quartic dispersion (s^4/m)
% ------------------

% ------------------
freqMethod = 'MULT';                        % Frequency shift method
maxroundtrips = 250;                        % Max number of loops for the laser 
SEED = 222;                                 % RNG Seed
N = 11;                                     % Grid size 
% ------------------

load('data/03__quartic__202202062259__zero.mat')

% ------------------
feedbacks = 0.5;                            % Percent energy retained each trip
freqdeltas = logspace(8, 10, N);            % Frequency shift (Hz) per trip
peaklambdas = zeros(1, N, ...
    maxroundtrips)*NaN;                     % Peak wavelength (nm)
peakamps = peaklambdas;                     % Peak amplitude
fields = zeros(N, npts);                    % Complex electric field
% ------------------

% ------------------
titlestring = ['03__quartic__' datestr(now, 'yyyymmddHHMM') '__instfreq']
% ------------------

close all;

for i = 1:1 

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
            spect = 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2)); 
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
        
        fields(j, :) = E;

        save(['data/' titlestring '.mat'], 'peaklambdas', 'peakamps', 'feedbacks', 'freqdeltas', 'fields')
        disp([i, j, mmm]);
        toc

    end

end

return 

%  Initial plots

initialise;
load('data/03__quartic__202202151512__instfreq.mat')
titlestring = ['03__quartic__' datestr(now, 'yyyymmddHHMM') '__instfreq']

i = 1; j = 9;
feedback = feedbacks(i);
freqdelta = freqdeltas(j); 

E = fields(j, :);
[~, ind] = max(abs(E).^2);

phi = angle(E);
phi(phi < 0) = phi(phi < 0) + 2 * pi; 
omegas = - diff(phi) ./ diff(t * 1e12); 
omegas(abs(omegas) > 10) = NaN; 

TE = fftshift(fft(fftshift(E)));  

figure('color', 'white')

subplot(211)
plot(t*1e12, abs(E).^2, 'LineWidth', 2); 
hold on 
plot(t*1e12, phi, '--', 'LineWidth', 2); 
plot(t(1:end-1)*1e12, omegas, 'LineWidth', 2); 
hold off 
xlim([t(ind)*1e12 - 3, t(ind)*1e12 + 3])
legend('|E|^2', '\phi', '\omega', 'Location', 'southwest')
xlabel('Time (ps)')
ylabel('Amplitude')
title({sprintf('feedback = %.2f, freq shift = %.2f GHz', feedback, freqdelta * 1e-9), titlestring, 'Field, phase, and instantaneous frequency'})
set(gca, 'FontSize', 14)

subplot(212)
plot(f*1e-9, abs(TE).^2, 'LineWidth', 2); 
xlim([-1000, 1000])
xlabel('Frequency (GHz)')
ylabel('Squared modulus')
title('Fourier transform of E')
set(gca, 'FontSize', 14)

saveas(gcf, ['data/' titlestring '.fig'])


return 

% More feedbacks and better plot 

initialise;
load('data/03__quartic__202202151700__instfreq.mat')
titlestring = ['03__quartic__' datestr(now, 'yyyymmddHHMM') '__instfreq']

[N, ~, TTT] = size(fields)

i = 2; j = 6;
feedback = feedbacks(i);
freqdelta = freqdeltas(j); 

E = reshape(fields(i, j, :), [1, TTT]);
[~, ind] = max(abs(E).^2);
t = t - t(ind); 

phi = angle(E);
phi(phi < 0) = phi(phi < 0) + 2 * pi; 
omegas = - diff(phi) ./ diff(t * 1e12); 
omegas(abs(omegas) > 10) = NaN; 

TE = fftshift(fft(fftshift(E)));  

figure('color', 'white')
yyaxis left
plot(t*1e12, abs(E).^2, 'LineWidth', 2); 
% ylabel('|E|^2')
hold on 
plot(t(1:end-1)*1e12, omegas, 'k-', 'LineWidth', 2); 
plot(-t(1:end-1)*1e12, -omegas, 'k:', 'LineWidth', 3); 
hold off 
ylabel('Amplitude')

yyaxis right 
plot(t*1e12, phi / pi, '--', 'LineWidth', 2); 
ylabel('\phi (units of \pi)')

xlim([- 3, 3])
legend('|E|^2', '\omega_i', '\omega_i, flipped', 'Location', 'south')
xlabel('Time (ps)')
title({'Field and instantaneous frequency', sprintf('feedback = %.2f, freq shift = %.2f GHz', feedback, freqdelta * 1e-9), titlestring})
set(gca, 'FontSize', 14)

saveas(gcf, ['data/' titlestring '__a.fig'])

figure('color', 'white')
plot(f*1e-9, abs(TE).^2, 'LineWidth', 2); 
hold on 
plot(-f*1e-9, abs(TE).^2, ':', 'LineWidth', 3); 
xlim([-1000, 1000])
xlabel('Frequency (GHz)')
ylabel('Squared modulus')
legend('Original', 'Flipped')
title({'Fourier transform of E', sprintf('feedback = %.2f, freq shift = %.2f GHz', feedback, freqdelta * 1e-9), titlestring})
set(gca, 'FontSize', 14)

saveas(gcf, ['data/' titlestring '__b.fig'])

return 

% Animation

initialise;
load('data/03__quartic__202202151700__instfreq.mat')
titlestring = ['03__quartic__' datestr(now, 'yyyymmddHHMM') '__instfreq']

[N, ~, TTT] = size(fields)

i = 5;
feedback = feedbacks(i);

figure('color', 'white', 'Position', [100 100 1200 450])

j = 1;
freqdelta = freqdeltas(j); 
E = reshape(fields(i, j, :), [1, TTT]);
[~, ind] = max(abs(E).^2);
t = t - t(ind); 

phi = angle(E);
omegas = - diff(phi) ./ diff(t * 1e12); 
omegas(abs(omegas) > 10) = NaN; 

TE = fftshift(fft(fftshift(E)));  

subplot(121)
yyaxis left
p1 = plot(t*1e12, abs(E).^2, 'LineWidth', 2); 
hold on 
p2 = plot(t(1:end-1)*1e12, omegas, 'k-', 'LineWidth', 2); 
p3 = plot(-t(1:end-1)*1e12, -omegas, 'k:', 'LineWidth', 3); 
hold off 
ylabel('Amplitude')
ylim([-10, 10])

yyaxis right 
p4 = plot(t*1e12, phi / pi, '--', 'LineWidth', 2); 
ylabel('\phi (units of \pi)')
ylim([-1.2, 1.2])

xlim([- 3, 3])
legend('|E|^2', '\omega_i', '\omega_i, flipped', 'Location', 'south')
xlabel('Time (ps)')
title({sprintf('feedback = %.2f, freq shift = %.2f GHz', feedback, freqdelta * 1e-9), 'Field, phase and instantaneous frequency'})
set(gca, 'FontSize', 14)

subplot(122)
q1 = plot(f*1e-9, abs(TE).^2, 'LineWidth', 2); 
hold on 
q2 = plot(-f*1e-9, abs(TE).^2, ':', 'LineWidth', 3); 
xlim([-1000, 1000])
xlabel('Frequency (GHz)')
ylabel('Squared modulus')
legend('Original', 'Flipped')
title({titlestring, 'Fourier transform of E'})
set(gca, 'FontSize', 14)

pause(1)

for j = 2:N 
    freqdelta = freqdeltas(j); 
    E = reshape(fields(i, j, :), [1, TTT]);
    [~, ind] = max(abs(E).^2);
    t = t - t(ind); 

    phi = angle(E);
    omegas = - diff(phi) ./ diff(t * 1e12); 
    omegas(abs(omegas) > 10) = NaN; 

    TE = fftshift(fft(fftshift(E)));  

    p1.YData = abs(E).^2; 
    p1.XData = t*1e12;
    p2.YData = omegas;
    p2.XData = t(1:end-1)*1e12;
    p3.YData = -omegas;
    p3.XData = -t(1:end-1)*1e12; 

    p4.YData = phi / pi; 
    p4.XData = t*1e12;

    q1.YData = abs(TE).^2;
    q1.XData = f*1e-9;
    q2.YData = abs(TE).^2;
    q2.XData = -f*1e-9;

    subplot(121)
    title({sprintf('feedback = %.2f, freq shift = %.2f GHz', feedback, freqdelta * 1e-9), 'Field, phase and instantaneous frequency'})

    pause(0.4)
    drawnow 
end 

return 