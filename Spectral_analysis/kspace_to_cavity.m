function [cavity_data, cavity_distance] = kspace_to_cavity(kspace_data,k_even_spacing,padding)
%KSPACE_TO_CAVITY Summary of this function goes here
%   Every column is a vector of k-space values, but the function returnes a
%   matrix where every column is the complex value in time for a given fft
%   bin

%hardcoded value for number of bins, f(spectrometer)
n_lambda = 512;

cavity_data = fft(kspace_data,n_lambda*padding)';

%just the first half
cavity_data = cavity_data(:,1:256*padding);

%from k-space to distance in [um]
space_freq = 1/(k_even_spacing(2)-k_even_spacing(1));
Fs = 10*space_freq/(pi*(n_lambda*padding)); % Sampling frequency [nm]
cavity_distance = Fs*(0:(n_lambda*padding/2-1))/1000; %size of the cavity [um]

end

