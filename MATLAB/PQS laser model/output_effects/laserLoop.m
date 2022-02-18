%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop through roundtrips or until s or p is pressed. s stops, p pauses

% Es = zeros(npts,roundtrips);
% EsFT = zeros(npts,roundtrips);

E = E0;

for mmm = 1:roundtrips
   % tic 
   
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
   % plotplot;
   SMF;
   
   E = E*sqrt(connector);
   
   % toc
end