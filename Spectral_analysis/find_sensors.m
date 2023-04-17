function [idx_pressure_sensor, idxs_sample] = find_sensors(cavity_data)
%FIND_SENSORS Summary of this function goes here
%   Detailed explanation goes here

% take abs, log10 transform
cavity_data = log10(abs(cavity_data));
% average out intensity
cavity_data = mean(cavity_data);
padding_est = round(length(cavity_data)/256);

% find maximum to normalize
[max_intensity, idx_max_intensity] = max(cavity_data);
cavity_data_norm = cavity_data/max_intensity;

% take mean + 3 *std value of noise using the right end of the array
treshold = 0.7;%mean(cavity_data_norm(end-50:end)) + 10 * std(cavity_data_norm(end-50:end));

% put everything below to 0
idx_tresholding = cavity_data_norm < treshold;
cavity_data_tresholded = cavity_data_norm;
cavity_data_tresholded(idx_tresholding) = 0;
% DC part is also useless
cavity_data_tresholded(1:5) = 0;

% pressure sensor is the local maximum
idx_pressure_sensor = find(islocalmax(cavity_data_tresholded(8*padding_est:11*padding_est)),1);
idx_pressure_sensor = idx_pressure_sensor + 8*padding_est;

%indices of the samples are nonzero events from left and right
narrowed_search_idx = idx_pressure_sensor*4;
idxs_sample(1) = narrowed_search_idx + ...
    find(cavity_data_tresholded(narrowed_search_idx:end),1,'first');
idxs_sample(2) = find(cavity_data_tresholded,1,'last');

% temporary, just to make sure it didn't mess up
figure
plot(cavity_data_norm)
hold on
plot(idxs_sample(1),cavity_data_tresholded(idxs_sample(1)),'ro')
plot(idxs_sample(2),cavity_data_tresholded(idxs_sample(2)),'ro')
plot(idx_pressure_sensor,cavity_data_tresholded(idx_pressure_sensor),'ro')

end

