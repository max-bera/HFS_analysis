function [displacement_field] = phase_to_displacement(cavity_data,mean_lambda)
%PHASE_TO_DISPLACEMENT input is matrix of complex values
%   This function calculated the displacement (nm) field of an HFS exp.
%   Important: the displacement is to be intended as the Optical path
%   length. For the pressure sensor, this is effectively the geometrical
%   displacement (as n=1). For the rest of the cavities, correct for
%   refractive index of the sample.

%%
displacement_field = unwrap(angle(cavity_data)) * mean_lambda / (4*pi);
end

