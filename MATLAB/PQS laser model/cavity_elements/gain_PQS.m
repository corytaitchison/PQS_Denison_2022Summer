%% Gain fibre

LG = 1.5;         % Length of the Gain fiber (m)
nstepG = 300;   % Number of steps 
Nplots = 201;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gain profile / 50 nm Lorentzian profile @ 1550 nm

gainband = 50e-9;           % Gain spectral bandwidth Erbium
dlam = gainband;    
lam_FWHM = lambda0+dlam/2;  % Center the gain band around central wavelength
dlamg = sqrt((lam_FWHM-lambda0)^2/0.5);
domg = 2*pi*c/lambda0^2*dlamg;

gprofile = (1-wrel.^2/domg.^2);
gprofile(gprofile < 0) = 0;

% plot(lambdanm,gprofile,'ro')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dispersion paramters for SMF28

beta2 = -21.4e-27;        % GVD in s^2/m       from Opt. Lett. 36, 112 (2011) "Peregrine soliton generation and breakup in standard telecom fiber" 
beta3 = 0.12e-39;         % TOD in s^3/m       from Opt. Lett. 36, 112 (2011) "Peregrine soliton generation and breakup in standard telecom fiber" 
beta4 = -0.0022e-51;      % FOD in s^4/m    

maxorder = 4;

beta = 0;

for ii = 2:maxorder
    
        eval(['betaG = beta + beta' num2str(ii) '*wrel.^' num2str(ii) ...
            '/factorial(' num2str(ii) ');']);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optical loss

Alpha_l = 0.19e-3;      % Loss SMF28 in dB/m
LdB = Alpha_l*LG;  % Loss in dB
Lss = 10^(LdB/10);      
alpha_l = log(Lss)/LG;   % Loss in 1/m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear parameters for SMF28

n2 = 2.6e-20;               % Nonlinear refractive index silica @ 1550 nm (m^2/W)
MFD = 10.4e-6;              % Mode field diameter Corning SMF28  
Aeff = 0.944*pi/4*MFD^2;    % Effective waveguide area (m^2)
k0 = 2*pi/lambda0;
gamma = (k0*n2/Aeff);       % Nonlinear coeff (W.m)^-1
% gamma = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dispersion operator definition etc.

hG = LG/nstepG;
betaop_fullG = exp(1i*hG*betaG);
betaop_halfG = exp(1i*hG*betaG/2);
  
L = LG;
nstep = nstepG;
h = L/nstep;

betaop_full = betaop_fullG;
betaop_half = betaop_halfG;
% gamma = gammaG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of gain at z = 0

Epulse = sum(abs(E).^2).*dt;

% g0_dB = 30;       % 30 dB small-signal gain
% g0 = 10^(g0_dB/10);
g0 = 3.45;      % Small-signal gain         % Nat. Photonics 4, 307 (2010) "Soliton-similariton fibre laser"
gg = g0/(1+Epulse/Esat);
g = gg*gprofile;

PropagationGain;



