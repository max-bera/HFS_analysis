path = 'C:\Users\massimiliano.berardi\Desktop\test_data\HCV29_3\20220324_155747_S5_01.tdms';
path_bg = 'C:\Users\massimiliano.berardi\Desktop\test_data\20220323_112003_background_01.tdms';
resampling_freq = 100; %Hz, max = acq_freq
padding = 5; %n times the original array length
pressure_sens = 76.50; %Pa/nm
pip_rad = 25e-6;

%% load spectrum, transform

[spectrum_lambda, mean_lambda, k_space] = load_spectrum(path,resampling_freq);

% background removal
spectrum_bg_lambda = load_spectrum(path,resampling_freq,'BG');

spectrum_noBG_lambda = spectrum_lambda-spectrum_bg_lambda;

% this needs to be done twice. BG removal messes up the phase demodulation
% k-space transform
[spectrum_k, k_even_spacing] = k_resample_M(spectrum_lambda,k_space);
% k-space to cavity
[cavity_data, cavity_distance] = kspace_to_cavity(spectrum_k,k_even_spacing,padding);

% calculate dOPL field
displacement_field = phase_to_displacement(cavity_data,mean_lambda);
time_array = linspace(0,length(displacement_field)/resampling_freq,length(displacement_field));

%% find sensors

% isolate pressure sensor and sample via intensity analysis
[spectrum_k_noBG, ~] = k_resample_M(spectrum_noBG_lambda,k_space);
[cavity_data_noBG, ~] = kspace_to_cavity(spectrum_k_noBG,k_even_spacing,padding);
%%
[idx_pressure_cavity, idxs_sample_cavity] = find_sensors(cavity_data_noBG);

%% extract pressure

pressure_signal = displacement_field(:,idx_pressure_cavity)*pressure_sens-...
    displacement_field(1,idx_pressure_cavity)*pressure_sens;

%% extract displacement
displ_signal = displacement_field(:,idxs_sample_cavity(1):idxs_sample_cavity(2)+5)/1.45 - ...\
    displacement_field(1,idxs_sample_cavity(1):idxs_sample_cavity(2)+5)/1.45;

figure(10)
subplot(1,2,1)
plot(time_array,pressure_signal,'k')
grid on
xlabel('time [s]')
ylabel('pressure [Pa]')
subplot(1,2,2)
hold off
for i = 1000:50:1000
    plot3(linspace(i,i,size(displ_signal,2))/resampling_freq,linspace(1,cavity_distance(2)*size(displ_signal,2),size(displ_signal,2)),displ_signal(i,:))
    xlabel('time [s]')
    ylabel('distance [\mum]')
    zlabel('displacement [nm]')
    hold on
    grid on
end
zlim([0 10000])