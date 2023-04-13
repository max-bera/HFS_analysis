function [spectrum,mean_lambda,k_space] = load_spectrum(path,resampling_freq,BG)
%LOAD_SPECTRUM takes the path of a DeltaSens *.tdms file, and extracts the
%Lambda space resampled data, in matrix form
%%

if nargin<3, BG = 'data'; end

    % import file
    dsens_struct = TDMS_getStruct(path);

    % extract sampling freq
    acq_freq = str2num(dsens_struct.Raw_data.Data.Props.Effective_F_Sps);
    if resampling_freq>acq_freq
        resampling_freq = acq_freq;
    end
    
    % calculate time interval [ms]
    time_ms = length(dsens_struct.Raw_data.Time.data);

    % define length of experiment
    N = length(dsens_struct.Raw_data.LambdaArray.data);
    mean_lambda = mean(dsens_struct.Raw_data.LambdaArray.data);
    acq_interval = round(acq_freq/resampling_freq);

    % if background import, return only a time-avg vector of
    % intensity(lambda)
    if strcmp(BG,'BG')
        bg_envelope = zeros(1,512);
        for i=1:512:length(dsens_struct.Raw_data.Data.data)
            bg_envelope = bg_envelope+double(dsens_struct.Raw_data.Data.data(i:i+511));
        end        
        spectrum = bg_envelope'/(i/512);
    else
        % rewrite raw data in matrix form 
        raw_time_domain = reshape(dsens_struct.Raw_data.Data.data,...
                [N time_ms]);
    
        % downsample as per user input
        downsampling_idx = 1:acq_interval:time_ms;
        spectrum= double(raw_time_domain(:,downsampling_idx));
        k_space = 2*pi./dsens_struct.Raw_data.LambdaArray.data;
    end

end

