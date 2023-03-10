%% deltasens data
linepath = '/Volumes/T7/PhD/KRAKOW/DOCETAXEL/HT1376/';
resampling_freq =500; %Hz, max = acq_freq
padding = 1; %n times the original array length
autocorr_interval = round(2*resampling_freq); %seconds*frames per sec
%% import file list
timepoint_list = dir(linepath);

%% load file
for k = 3:3%size(timepoint_list,1)-2
    dosage_list = dir(strcat(timepoint_list(k+2).folder,'/',timepoint_list(k+2).name,'/'));
    for j = 4:4%size(dosage_list,1)-3
        file_list = dir(strcat(dosage_list(j+2).folder,'/',dosage_list(j+2).name,'/'));
        LIV_estimate = zeros(256,size(file_list,1)-2);
        OCDS_estimate = zeros(256,size(file_list,1)-2);
        intensity_over_time = zeros(256,5000);
        fit_delay = 10;
        autolag_e = zeros(autocorr_interval,size(file_list,1)-2);
        autolag_c = zeros(autocorr_interval,size(file_list,1)-2);
        for i = 1:size(file_list,1)-2
            filename = strcat(file_list(i+2).folder,'/',file_list(i+2).name);
            intensity_over_time = HFS_struct_load(filename,resampling_freq,padding);
            %% LIV calculation
            LIV_estimate(:,i) = LIV_calc(intensity_over_time,1000);
            %% OCDS calculation
            OCDS_estimate(:,i) = OCDS_calc(intensity_over_time,400,resampling_freq,0);
            disp(i);
            %% Autocorr calc
      %      [t1 lags] = autocorr(intensity_over_time(44,1:autocorr_interval),'NumLags',autocorr_interval-1);
      %      [t2 lags] = autocorr(intensity_over_time(45,1:autocorr_interval),'NumLags',autocorr_interval-1);
           % [t3 lags] = autocorr(intensity_over_time(46,1:autocorr_interval),'NumLags',autocorr_interval-1);
        %    autolag_e(:,i) = (t2)';%/3;
         %   [t1 lags] = autocorr(intensity_over_time(56,1:autocorr_interval),'NumLags',autocorr_interval-1);
        %    [t2 lags] = autocorr(intensity_over_time(59,1:autocorr_interval),'NumLags',autocorr_interval-1);
         %   [t3 lags] = autocorr(intensity_over_time(60,1:autocorr_interval),'NumLags',autocorr_interval-1);
         %   autolag_c(:,i) = (t2)';%/3;
        end
    end
end

%%
%figure
%errorbar(mean(HT1376_3D_LIV(1:100,:)'),std(HT1376_3D_LIV(1:100,:)'))
%hold on
%errorbar(mean(HT1376_14D_LIV(1:100,:)'),std(HT1376_14D_LIV(1:100,:)'))
% 
% for i = 1:size(intensity_over_time,1)
%     [autocorr_est(i,:), lags] = xcorr(log(intensity_over_time(i,3001:5500)));
% end
% 
% OCDS_estimate2 = abs(mean(diff(autocorr_est(:,3100:end)')));
% %%
%%
 startm = 1;
  stopm = 23;
 figure(10)
 subplot(2,1,1)
 hold on
 errorbar(mean(LIV_estimate(:,startm:stopm)'),std(LIV_estimate(:,startm:stopm)'))
 subplot(2,1,2)
 hold on
 errorbar(abs(mean(OCDS_estimate(:,startm:stopm)')),std(OCDS_estimate(:,startm:stopm)'))
%%
%prova = cov(intensity_over_time(45,1:end),[intensity_over_time(45,3:end),0,0])./sqrt(var(intensity_over_time(45,1:end))*var(intensity_over_time(45,1:end)))

%%
sample_number = 2;
%%
HSV_spheroid = zeros(256,23,3);
%%
for j = 1:23
    sample_number = j;
HUE_spheroid = zeros(256,1);
SATURATION_spheroid = zeros(256,1);
VALUE_spheroid = zeros(256,1);
for i=1:size(intensity_over_time,1)
    [pxx,w] = pwelch((intensity_over_time(i,1:1000)),25);
    HUE_spheroid(i) = dot(w,log(pxx/vecnorm(pxx,1)));
    SATURATION_spheroid(i) = dot(pxx/vecnorm(pxx,1),(w).^2)-dot(pxx/vecnorm(pxx,1),(w)).^2;
    VALUE_spheroid(i) = (LIV_estimate(i,sample_number));
end
    HUE_spheroid = rescale((HUE_spheroid),0,0.66);
    VALUE_spheroid = rescale(abs(VALUE_spheroid),0.05,0.8);
    SATURATION_spheroid = rescale((1./SATURATION_spheroid),0.0,1);
    HSV_spheroid(:,sample_number,:) = hsv2rgb(cat(3,HUE_spheroid,SATURATION_spheroid,VALUE_spheroid));
end
%%
figure
imagesc(HSV_spheroid(10:120,:,:))

title("TEST- A")
%%
   % HUE_spheroid = OCDS_estimate(:,sample_number);

