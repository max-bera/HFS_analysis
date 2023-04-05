%% import experimental data (formatted as SLINE TIME #MEAS FREQ ES EL)

path = '/Volumes/T7/PhD/EXPERIMENTAL DATA/SPHEROIDS/Aging/Analysis/Mechanical analysis/DMA_HFS_HCV3.csv'

%select line to model
s_line = 'T24';
%select timepoint to model
s_time = '3';
%select the lumped parameter model to use
LP_model = '';


DATA_M = readtable(path);

%% filter by line and day

DATA_M.line = categorical(DATA_M.line)
idx_line = DATA_M.line == s_line;
data_oneline = DATA_M(idx_line,:);

data_oneline.day = categorical(data_oneline.day);
idx_day = data_oneline.day == s_time;
data_DMA_tab = data_oneline(idx_day,:);

%extract number of indipendent samples measured
n_ind_samples = length(unique(data_DMA_tab.measure));

%% clean data and model

%move form table to array, for convenience
data_DMA = table2array(data_DMA_tab(:,4:6));
%remove failed measurements
fail_idx = isnan(data_DMA(:,2)) | isnan(data_DMA(:,3));
rheo_data = data_DMA(~fail_idx,:);
rheo_freq = rheo_data(:,1);
rheo_compl = rheo_data(:,2) + 1i * rheo_data(:,3);


%do the fitting. 2PL require 4 pars, PT 6
x0 = [1000 0.1 10 1 10000 1 1];
lb = [0 0.06 0 0 0 0];
ub = [1e6 1 1e6 1 1e6 1];

options = optimoptions(@lsqnonlin,'Display','off','MaxFunctionEvaluations',3e4,'MaxIteration',1e3,'FunctionTolerance',1e-12,'StepTolerance',5e-12);
[model_val,resnorm,al,exitflag] = lsqnonlin(@(x) fractional_poynting_thompson_objfun(x,rheo_freq,rheo_compl/1000),x0,lb,ub,options);

%show me what's up
figure
visuals_complex_model('PT','SE',rheo_freq,rheo_data(:,2)/1000,rheo_data(:,3)/1000,model_val,3);
