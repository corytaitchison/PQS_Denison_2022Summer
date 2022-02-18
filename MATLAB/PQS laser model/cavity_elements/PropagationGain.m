% Relevant arrays when plotting while simulating

sel = round(nstep/(Nplots-1));      % plot selection and arrays

TE = fftshift(ifft(fftshift(E)));
TE = TE.*betaop_half.*exp(g*h/4);
E  = fftshift(fft(fftshift(TE)));

for mm = 1:nstep;
   
    Epulse = sum(abs(E).^2).*dt;
    gg = g0/(1+Epulse/Esat);
    g = gg*gprofile;
    
    % Nonlinear step
    A = E;
    I = abs(A).^2;
    E = A.*exp(1i*gamma*I*h-V);
    
    % Dispersion step
    TE = fftshift(ifft(fftshift(E)));
    
    if (mm<nstep)
        TE = TE.*betaop_full.*exp(g*abs(h)/2); 
    else
        TE = TE.*betaop_half.*exp(g*abs(h)/4);
    end
    
    E  = fftshift(fft(fftshift(TE)));
    loss = exp(-0.5*alpha_l*h);
    E = E.*loss;
    
end

% E  = fftshift(E);
% Especout = fftshift(TE);

% Relevant arrays when plotting while simulating --------------------------

% sel = round(nstep/(Nplots-1));      % plot selection and arrays
% 
% E = fftshift(E);
% TE = ifft(E).*opdisphalf.*exp(gcoeff*g*abs(h)/4);      % First dispersion step
% E  = fft(TE);
% 
% nphot0 = sum( abs(TE*dt).^2 * df ./hplanck ./ffw);
% % 
% % for k = 1:nstep,
%         
%    Epulse = sum(abs(E).^2).*dt;
%    Pavg = Epulse*frep;
% %    [a b] = min(abs(Pavg - P_signalIns));
%    [aa bb] = min(abs(zarray(k) - Larray));
% 
%    g = fftshift(gprofile.*gmatrix(b,bb));   
%    
%    U0 = E;
%    TI = ifft(abs(U0).^2);
%    
%    int0 = (1-fR)*abs(U0).^2+fR*fft(chi.*TI);
%    Nd = U0.*int0;
%    U1 = U0 - h/2*gamma*taushock*(Nd(idxp)-Nd(idxm))/(2*dt);
%    TI = ifft(abs(U1).^2);
%    int1 = (1-fR)*abs(U1).^2+fR*fft(chi.*TI);
%    Nd = U1.*int1;
%    U1 = U0-h*1i*gamma*U1.*(int1-int0)-h*gamma*taushock*(Nd(idxp)-Nd(idxm))/(2*dt);
%    
%    E  = U1.*exp(1i*gamma*h*int0);
% 
%    TE = ifft(E);
%         
%    if (k<nstep),
%        TE = TE.*opdispfull.*exp(gcoeff*g*abs(h)/2);
%        
%    else TE = TE.*(opdisphalf).*exp(gcoeff*g*abs(h)/4);
%    end;
% 
% 
%    E  = fft(TE);
%         if saveplot == 1
%          if (rem(k,sel)==0),
%              TEstep = TE./opdisphalf;
%              Estep  = fftshift(fft(TEstep));
%              
%              Up =  [Up; multcoeff*Estep];
%              Ufp = [Ufp fftshift(TEstep)];
% 
%             subplot(211),
%             plot(t,abs(Estep).^2,'k','linewidth',2), hold on
%             hold off
%             xlabel('Time (ps)', 'Fontsize',24),ylabel('Intensity', 'Fontsize',24),title(['z = ' num2str(k*h) 'm'],'Fontsize',24)
%             set(gca,'Fontsize',24)
%             axis([-tspan/2 tspan/2 0 1.5*max(abs(Estep).^2)])
%             Pmax(1+(k/sel)) = max(abs(Estep).^2);
%             set(gca,'YTickLabel',[])
%             drawnow
%             
%            
%             subplot(212),
%             plot(omw,abs(fftshift(ifft(Estep))).^2,'k','linewidth',2)
%             xlabel('Wavelength (nm)','Fontsize',24),ylabel('Spectrum', 'Fontsize',24)
%             set(gca,'YTickLabel',[])
%             set(gca,'Fontsize',24)
%             drawnow       
% 
%         end;
%         end
% end;
%         
% E  = fftshift(E);
% Especout = fftshift(TE);



