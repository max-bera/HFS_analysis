function [pressure_signal,sample_displ_signal] = get_sensor_data(cavity_data,mean_lambda,idx_pressure_cavity, idxs_sample_cavity)
%GET_SENSOR_DATA Summary of this function goes here
%   Detailed explanation goes here

pressure_signal = phase_to_displacement(cavity_data,mean_lambda,idx_pressure_cavity);
sample_displ_signal = phase_to_displacement(cavity_data,mean_lambda,idxs_sample_cavity);

end

