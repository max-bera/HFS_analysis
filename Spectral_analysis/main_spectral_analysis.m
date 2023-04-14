% macos paths
path = '/Users/massimilianoberardi/Desktop/Spheroids_results/CONTROL_3_14/HCV_DAY3/20220324_142758_S1_01.tdms';
path_bg = '/Users/massimilianoberardi/Desktop/Spheroids_results/20220323_112003_background_01.tdms';
% win paths
path = 'C:\Users\massimiliano.berardi\Desktop\test_data\HCV29_3\20220324_150902_S3_01.tdms';
path_bg = 'C:\Users\massimiliano.berardi\Desktop\test_data\20220323_112003_background_01.tdms';

resampling_freq = 500; %Hz, max = acq_freq
padding = 2; %n times the original array length
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

[idx_pressure_cavity, idxs_sample_cavity] = find_sensors(cavity_data_noBG);

opti_cross_sec = (idxs_sample_cavity(2) - idxs_sample_cavity(1)) * cavity_distance(2);
%% extract pressure

pressure_signal = displacement_field(:,idx_pressure_cavity)*pressure_sens-...
    displacement_field(1,idx_pressure_cavity)*pressure_sens;
%make it nice
freq = 60; %cutoff
     Order = 5;
     [b, a] = butter(Order, freq / (resampling_freq * 2), 'low');
     pressure_signal = filter(b,a,pressure_signal);

%% extract displacement
displ_signal = displacement_field(:,idxs_sample_cavity(1):idxs_sample_cavity(2))/1.45 - ...
    displacement_field(1,idxs_sample_cavity(1):idxs_sample_cavity(2))/1.45;
displ_signal = filter(b,a,displ_signal);

K = (1/9)*ones(3);
displ_signal = conv2(displ_signal(1:1:end,:),K,'same');
figure(10)
subplot(1,2,1),
plot(time_array,pressure_signal,'k')
grid on
xlabel('time [s]')
ylabel('pressure [Pa]')
%limit plot to interesting part, i.e. 1.8 times pip_rad
idx_decay = find(cavity_distance/1.45>pip_rad*4*1e6*1.45,1,'first')
subplot(1,2,2)
hold off
j=1;
lprp=0;
ptoP=0;
for i = 3*resampling_freq:resampling_freq/10:6*resampling_freq
    plot(cavity_distance(1:idx_decay)/1.45,displ_signal(i,1:idx_decay),'k.-')
    lprp(j) = mean(displ_signal(i,2:4))*1e-9/pip_rad/1.45*1.33;
    ptoP(j) = mean(displ_signal(i,idx_decay-5:idx_decay))/(mean(displ_signal(i,2:4))/1.45*1.33);
    j = j+1;
    xlabel('distance [\mum]')
    ylabel('displacement [nm]')
    hold on
    grid on
end
%%
figure
plot(lprp,ptoP)