function [LIV_estimate] = LIV_calc(intensity_over_time_matrix,n_frames_averaging)
%Estimate logarithmic intensity variance, used as a proxy for cellular
%activity
%
%   intensity matrix format: every column is a frame, every row a cavity
%%
%calculate time average log(variance) for every cavity
time_avg_LI = mean(10*log10(intensity_over_time_matrix),2);

%remove time indipendent part, calculate log intensity variance
LI_dynamic_only = 10*log10(intensity_over_time_matrix)-time_avg_LI;
LIV_estimate = 1 / n_frames_averaging * sum(LI_dynamic_only.^2,2);

%remove symmetric part of FFT
LIV_estimate = LIV_estimate(1:length(time_avg_LI));

end