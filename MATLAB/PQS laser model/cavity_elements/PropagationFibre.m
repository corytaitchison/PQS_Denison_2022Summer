

TE = fftshift(ifft(fftshift(E)));
TE = TE.*betaop_half;              % First dispersion step
E  = fftshift(fft(fftshift(TE)));

for mm = 1:nstep;

    % Nonlinear step
    A = E;
    I = abs(A).^2;
    E = A.*exp(1i*gamma*I*hSMF-V);
    
    % Dispersion step
    TE = fftshift(ifft(fftshift(E)));
    
    if (mm<nstep)
        TE = TE.*betaop_full; 
    else 
        TE = TE.*betaop_half; 
    end
    E  = fftshift(fft(fftshift(TE)));
    
    loss = exp(-0.5*alpha_l*hSMF);
    E = E.*loss;
end
