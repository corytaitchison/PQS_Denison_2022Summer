% quartic__220206__grid.m
%
% ------------------
% Created: 2022-02-06 22:30
% Author: Cory
% Title: Quartic Grid
% Description:
%     Compute the initial E0 profiles for different feedbacks
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

% ------------------
% feedbacks = linspace(0.5, 0.8, N);          % Percent energy retained each trip
% freqdeltas = [0];
% freqdeltas = logspace(10, 12, N);           % Frequency shift (Hz) per trip

feedbacks = linspace(0.3, 0.8, N);          % Percent energy retained each trip
freqdeltas = logspace(8, 10, N);            % Frequency shift (Hz) per trip
peaklambdas = zeros(N, 1, ...
    maxroundtrips)*NaN;                     % Peak wavelength (nm)
E0s = [];
% ------------------

% ------------------
titlestring = ['03__quartic__' datestr(now, 'yyyymmddHHMM') '__zero']
% ------------------


% ------------------ Get initial profile 
for i = 1:N 
    tic 
    initialise;
    E = E0;

    feedback = feedbacks(i);
    freqdelta = 0; 

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
        EnP(mmm) = (1-feedback)*sum(abs(E).^2)*dt;
        plotplot
        SMF;
        E = E*sqrt(connector);
        %  ------------------ /Laser Loop

        %  ------------------
        spect = 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2)); 
        [~, ind] = max(spect);
        %  ------------------    

        peaklambdas(i, 1, mmm) = lambdanm(ind);  % In nm

        if mmm > 50
            if std(peaklambdas(i,1, (mmm-20):(mmm-1))) < 1e-2
                break 
            end
        end

    end

    E0s = [E0s; E];

    save(['data/' titlestring '.mat'], 'peaklambdas', 'feedbacks', 'E0s')
    disp([i, mmm]);
    toc

end