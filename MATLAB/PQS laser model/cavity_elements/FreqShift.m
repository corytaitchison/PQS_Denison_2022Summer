% FreqShift.m
%
% ------------------
% Created: 2022-01-17 10:52am
% Author: Cory
% Title: Frequency Shifter
% Description:
%   Cavity element for applying a small frequency shift to the pulse each rounnd trip
% ------------------
% 
% --------
% MODIFIES
% --------
% - E: [1, n] complex float
%       The electric field 
% 
% -------
% GLOBALS
% -------
% - freqMethod: "FFT" or "MULT"
%       Whether to perform Fourier transform shifts, or multiply
%       by exp(-i omega t)
% - freqdelta: float
%       The frequency shift in Hertz to be applied
% - df: float
%       The frequency interval in Hz for the Fourier transform frequencies

if freqMethod == "FFT"

    % Get the spectrum
    TE = fftshift(ifft(fftshift(E)));  
    % Shift the spectrum by a given frequency shift
    shift = round(freqdelta / df); 
    if shift == 0 
        error('Frequency shift is smaller than frequency resolution');
    end
    TE = [zeros(1, shift) TE(1:end-shift)];
    % Recompute the field
    E = fftshift(fft(fftshift(TE)));

elseif freqMethod == "MULT"

    E = E .* exp(-1i * freqdelta * 2 * pi * t);

end 