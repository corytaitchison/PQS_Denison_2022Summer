%% Fibre segment

nstep = 100;
L_fibre = 1;
Nplots = 101;

% Dispersion paramters for SMF28

beta2 = -21.4e-27;        % GVD in s^2/m       Opt. Lett. 36, 112 (2011) "Peregrine soliton generation and breakup in standard telecom fiber" 
beta3 = 0.12e-39;         % TOD in s^3/m       Opt. Lett. 36, 112 (2011) "Peregrine soliton generation and breakup in standard telecom fiber" 
beta4 = -0.0022e-51;       % FOD in s^4/m    

maxorder = 4;

beta = 0;
% beta = beta + beta2/2*wrel.^2 + beta3/6*wrel.^3 + beta4/24*wrel.^4;

for ii = 2:maxorder
    
        eval(['beta = beta + beta' num2str(ii) '*wrel.^' num2str(ii) ...
            '/factorial(' num2str(ii) ');']);
end

% Optical loss

Alpha_l = 0.19e-3;      % Loss SMF28 in dB/m
LdB = Alpha_l*L_fibre;  % Loss in dB
Lss = 10^(LdB/10);      
alpha_l = log(Lss)/L_fibre;   % Loss in 1/m

% Nonlinear parameters for SMF28

n2 = 2.6e-20;               % Nonlinear refractive index silica @ 1550 nm (m^2/W)
MFD = 10.4e-6;              % Mode field diameter Corning SMF28  
Aeff = 0.944*pi/4*MFD^2;    % Effective waveguide area (m^2)
k0 = 2*pi/lambda0;
gamma = (k0*n2/Aeff);       % Nonlinear coeff (W.m)^-1
% gamma = 0;

hSMF = L_fibre/nstep;

betaop_full = exp(1i*hSMF*beta);         
betaop_half = exp(1i*hSMF*beta/2);          

PropagationFibre;

