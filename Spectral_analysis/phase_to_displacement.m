function [sensor_signal] = phase_to_displacement(cavity_data,mean_lambda,idx_sensor)
%PHASE_TO_DISPLACEMENT input is matrix of complex values
%   This function calculated the displacement (nm) field of an HFS exp.
%   Important: the displacement is to be intended as the Optical path
%   length. For the pressure sensor, this is effectively the geometrical
%   displacement (as n=1). For the rest of the cavities, correction for
%   refractive index of the sample still need to be applied.

%%
phase = zeros(size(cavity_data));
det_treshold = 2000;
pressure_sens = 76.50; %Pa/nm
RI_sample = 1.4;
K_convolution_size = 4;

for i = 1:length(cavity_data)
    z = cavity_data(i,:);
    threshold =max(abs(z))/det_treshold; %tolerance threshold
    z(abs(z)<threshold) = 0; %mask values below the threshold
    phase(i,:) = angle(z);
end


if length(idx_sensor) == 1
    displacement_field = unwrap(phase,[],1) * mean_lambda / (4*pi);
    sensor_signal = displacement_field(:,idx_sensor)*pressure_sens-...
        displacement_field(1,idx_sensor)*pressure_sens;
else
    K = (1/K_convolution_size^2)*ones(K_convolution_size);
    displacement_field = conv2(unwrap(phase,[],1),K,'same') * mean_lambda / (4*pi);
    sensor_signal = (displacement_field(:,idx_sensor(1):idx_sensor(2))-...
        displacement_field(1,idx_sensor(1):idx_sensor(2)))/RI_sample;
end

end

