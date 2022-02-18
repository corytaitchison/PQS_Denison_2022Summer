

% figure
subplot(211)
E2 = abs(E).^2;
plot(t*1e12, E2, 'LineWidth', 2)
xlabel('Time (ps)')
title(['Roundtrip = ' num2str(mmm),'   ','Energy = ' num2str(EnP(mmm)*1e12),'pJ'])
[maxE2, arg] = max(E2); 
tcentre = t(arg)*1e12;
xlim([tcentre-5 tcentre+5])

% xlim([-5 5])

% axis([-10 10 0 1.2])

subplot(212)
plot(lambdanm, 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2)), 'LineWidth', 2)
hold on
xlabel('Wavelength (nm)')
xlim([1530 1590])
% axis([1540 1580 0 1.2])
drawnow
hold off