function [OCDS_estimate] = OCDS_calc(intensity_over_time_matrix,delay_range,acquisition_frequency,fit_delay)
%Estimate OCT correlation decay speed, to infer dynamics of activity in
%spheroids
%
%intensity matrix format: every column is a frame, every row a cavity
%delay range is expressed in [s]
%acquisition frequency is expressend in [Hz]
%%
%calculate autocorrelation f(tau,z)

%autocorr_est = zeros(size(intensity_over_time_matrix,1),delay_range*acquisition_frequency*2-1);

temp = zeros(size(intensity_over_time_matrix,1),1);
for i=1:size(intensity_over_time_matrix,1)
    ACF = xcorr((intensity_over_time_matrix(i,1:end)),'normalized');
    ACF_cut = ACF(round(length(ACF)/2)+1+fit_delay:round(length(ACF)/2)+delay_range+1+fit_delay);
    [p, ~] = polyfit(linspace(0,10,delay_range+1),ACF_cut,1);
    temp(i) = p(1);
    
    %abs(mean(diff(ACF(round(length(ACF)/2):timespan*resampling_freq))));
end

OCDS_estimate = temp;

% 
% for i = 1:size(intensity_over_time_matrix,1)
%     [autocorr_est(i,:), ~] = xcorr(log(intensity_over_time_matrix(i,1:delay_range*acquisition_frequency)),'normalized');
% end
% 
% %calculate OCDS
% OCDS_estimate = abs(mean(diff(autocorr_est(:,delay_range*acquisition_frequency:end)')))/delay_range;
% % 
% figure
% mesh(diff(autocorr_est(:,delay_range*acquisition_frequency:end)))
end