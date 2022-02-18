

TE = fftshift(ifft(fftshift(E)));
beta_WS = wrel.^2*beta2_WS/2 + wrel.^3*beta3_WS/6 + wrel.^4*beta4_WS/24;
Wsop = exp(1i*beta_WS*Ltot);
TE = TE.*Wsop;
E = fftshift(fft(fftshift(TE)));

% figure()
% plot(lambdanm,angle(Wsop))