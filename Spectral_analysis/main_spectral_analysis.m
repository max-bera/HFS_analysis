% macos paths
path = '/Users/massimilianoberardi/Desktop/Spheroids_results/CONTROL_3_14/HCV_DAY14/20220330_112835_S5_01.tdms';
path_bg = '/Users/massimilianoberardi/Desktop/Spheroids_results/20220323_112003_background_01.tdms';
% win paths
%path = 'C:\Users\massimiliano.berardi\Desktop\test_data\HCV29_3\20220324_150902_S3_01.tdms';
%path_bg = 'C:\Users\massimiliano.berardi\Desktop\test_data\20220323_112003_background_01.tdms';

resampling_freq = 100; %Hz, max = acq_freq
padding = 2; %n times the original array length
pip_rad = 25e-6;

%% load spectrum, transform

[spectrum_lambda, mean_lambda, k_space] = load_spectrum(path,resampling_freq);

% k-space transform
[spectrum_k, k_even_spacing] = k_resample_M(spectrum_lambda,k_space);
% k-space to cavity
[cavity_data, cavity_distance] = kspace_to_cavity(spectrum_k,k_even_spacing,padding);

%% find sensors

% isolate pressure sensor and sample via intensity analysis
spectrum_noBG_lambda = spectrum_lambda-load_spectrum(path,resampling_freq,'BG');
[spectrum_k_noBG, ~] = k_resample_M(spectrum_noBG_lambda,k_space);
[cavity_data_noBG, ~] = kspace_to_cavity(spectrum_k_noBG,k_even_spacing,padding);

% indices
[idx_pressure_cavity, idxs_sample_cavity] = find_sensors(cavity_data_noBG);

%% extract pressure and displacement in time

[pressure_signal, displ_signal] = get_sensor_data(cavity_data,mean_lambda,idx_pressure_cavity,idxs_sample_cavity);
time_array = linspace(0,length(displ_signal)/resampling_freq,length(displ_signal));

%add a lowpass
freq = 20; %cutoff
Order = 5;
[b, a] = butter(Order, freq / (resampling_freq * 2), 'low');

pressure_signal = filter(b,a,pressure_signal);
displ_signal = filter(b,a,displ_signal);


%%
figure(10)
subplot(1,2,1),
plot(time_array,pressure_signal,'k')
yyaxis right
plot(time_array,displ_signal(:,5),'r')
grid on
xlabel('time [s]')
ylabel('pressure [Pa]')
%limit plot to interesting part, i.e. 1.8 times pipette diameter
idx_decay = find(cavity_distance/1.45>pip_rad*5*1e6,1,'first');
subplot(1,2,2)
hold off
j=1;
lprp=0;
ptoP=0;
for i = 2*resampling_freq:resampling_freq:9*resampling_freq
    plot(cavity_distance(1:idx_decay)/1.45,displ_signal(i,1:idx_decay),'ko-')
    lprp(j) = mean(displ_signal(i,4:6))*1e-9/pip_rad/1.45*1.33;
    ptoP(j) = mean(displ_signal(i,idx_decay-5:idx_decay))/(mean(displ_signal(i,4:6)));
    j = j+1;
    xlabel('distance [\mum]')
    ylabel('displacement [nm]')
    hold on
    grid on
end

%% prepare structure for analysis in ALVA
j=1;
for i=20:20:300
exp.diameter = (idxs_sample_cavity(2)-idxs_sample_cavity(1))*cavity_distance(2)/1.45;
exp.load = i; %Pa
idx_profile = find(pressure_signal<-exp.load,1,'first');
% at times segmentation is too generous, double check based on max displ
[~,idx_max_displacement] = max(displ_signal(idx_profile,:));
if idx_max_displacement == 1
    exp.zdisp = displ_signal(idx_profile,1:idx_decay)';
else
    exp.zdisp = displ_signal(idx_profile,idx_max_displacement-1:idx_decay+idx_max_displacement-2)';
end
exp.zpos = cavity_distance(1:idx_decay)'/1.45;

save(strcat('HCV_3_',num2str(j)), '-struct', 'exp');
j = j+1;
end