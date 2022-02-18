% savePlot.m
%
% ------------------
% Created: 2022-01-17 10:52am
% Author: Cory
% Title: Save Plot
% Description:
%   Saves the plot and the data for the plot to a text file
% ------------------

function savePlot(E, t, lambdanm, filename, folder) 
    arguments
        E 
        t 
        lambdanm 
        filename (1,:) char
        folder (1,:) char = "data"
    end
    % Saves the plotplot.m figure and the data for the plot to a text file

    % Calcualte wavelength spectrum
    spec = 10*log10(abs(fftshift(ifft(fftshift(E)))).^2/max(abs(fftshift(ifft(fftshift(E)))).^2));

    display(['Saving to ' folder '/' filename '.mat'])
    % Save data
    save([folder '/' filename '.mat'], 't', 'E', 'lambdanm', 'spec')
    % Save figure
    saveas(gcf, [folder '/' filename '.fig'])

end