% Code that simulate the PQS laser (SMF + WS)
% A. F. J. Runge, March 2019

clear all,
close all,
clc,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval('InputParameters')
addpath('cavity_elements')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of time and frequency grids 

t = (-npts/2:npts/2-1)/npts*tspan; dt = t(2)-t(1);
f = (-npts/2:npts/2-1)/(npts*dt); df = f(2)-f(1);

f0 = c/lambda0;
wrel = 2*pi*f;                      % relative angular frequency
wabs = 2*pi*(f+f0);                 % absolute angular frequency
lambda = 2*pi*c./wabs;              % wavelength array
lambdanm = lambda*1e9;              % wavelength array in nm
lambdanm(lambdanm<0) = nan;         % if neg. freqs => neg. lambdas

V = 0*circshift(-1000*1i*sech(t/2e-12),npts/2);

figure('color', 'white')
set(gca, 'FontSize', 14)
plot(t*1e12,abs(V).^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input pulse definition

randn('state',221)
E0 = sech(t/1000e-15).*randn(1,npts);

% figure(444)
% plot(t,abs(E0).^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cavity specific parameters

roundtrips = 100;
splice = 0.97;                  % power loss in splice 
connector = 0.90;               % power loss in conncector
feedback = 0.5;                 % percent output coupler
lWS = 0.3;                      % insertion loss WS (~5dB)  

Ltot = 21.5;                  % Total cavity length (m), to be manually adjusted

% For PQS mode-locking
Esat = 99e-12;          % Sat energy of gain medium to simulate pump power
q0 = 0.7;             % unsaturated loss of the saturable absorber          A. Chong, H. Renninger, and F. W. Wise, JOSAB 25, 140-148 (2008)
Psat = 100;               % saturation power of the saturable absorber 

% % For PQS mode-locking
% Esat = 400e-12;          % Sat energy of gain medium to simulate pump power
% q0 = 0.7;             % unsaturated loss of the saturable absorber          A. Chong, H. Renninger, and F. W. Wise, JOSAB 25, 140-148 (2008)
% Psat = 500;               % saturation power of the saturable absorber 


% beta coeff for WS
beta2_WS = 21.4e-27;        % beta2 with opposite sign for Waveshaper (s^2/m) 
beta3_WS = -0.12e-39;       % beta3 with opposite sign for Waveshaper (s^3/m)
beta4_WS = -80e-51;         % beta4 for quartic dispersion (s^4/m)

% -------
freqdelta = 1e10;             % Frequency shift per round trip (in Hz)
freqMethod = "MULT";        % Frequency shift method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through roundtrips or until s or p is pressed. s stops, p pauses

% Es = zeros(npts,roundtrips);
% EsFT = zeros(npts,roundtrips);

E = E0;

for mmm = 1:roundtrips
   tic 
   
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
   plotplot;
   SMF;
   
   E = E*sqrt(connector);
   
   toc
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Fitting the laser output stable solution
% addpath('C:\Users\arun2330\Documents\MATLAB\Pure-Quartic Soliton\Article_figures')
% addpath('C:\Users\arun2330\Documents\MATLAB\Pure-Quartic Soliton\Laser WS\after_fixing')

lambdanum = lambdanm;
Epqs = E;
% save('Epqs.mat','lambdanum','Epqs')

load WaveData20190807_229.csv 
lambda_exp = WaveData20190807_229(:,1);
Ift_exp = WaveData20190807_229(:,2);

figure()
plot(lambda_exp,10*log10(Ift_exp./max(Ift_exp)),'color',[0 0 1],'linewidth',1)
hold on
plot(lambdanm, 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2))-6.17,'color',[1 0 0.2],'linewidth',1)
ylim([-38 5])
% xlim([1558 1563])
xlabel('Wavelength (nm)')
ylabel('Intensity (dB)')
legend('Exp','Sim')

figure()
plot(1:roundtrips, EnP*1e12, 'o')
xlabel('Roundtrips')
ylabel('Energy (pJ)')
