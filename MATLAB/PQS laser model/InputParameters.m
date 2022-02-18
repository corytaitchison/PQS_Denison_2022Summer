c = 2.99792458e8;                   % Speed of light (m/s)
hplanck = 6.626e-34;                % Planck's constant
kT = 1.38e-23*300;                  % Boltzmann term kT

npts =  4*4096;  %5*                   % Number of time-frequency points
tspan = 400e-12; % 600                 % s

lambda0 = 1560.35e-9;%28.2                  % Input Pulse Wavelength
% lambdaR = 1075e-9;                    % Raman shift silica - 13THz
% fR = 1*0.18;                          % Raman fraction

% taushock0 = 1*lambda0/(2*pi*c);       % 1/w0 = 0.563e-15 for 1060 nm;
% 
% ramanfile = 'chicom1';              % Raman data filename
% ramannoise = 0;                     % 0 for no spontaneous Raman noise, 1 for spontaneous Raman noise 
% 
randn('state',0);                   % state of random number generator
rand('state',0);                    % state of random number generator
statezero = 0;                      % set state offset value
%------------------------------------------------------------------------
