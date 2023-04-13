%% import experimental data (formatted as SLINE TIME #MEAS FREQ ES EL)

path = 'D:\PhD\EXPERIMENTAL DATA\SPHEROIDS\Aging\Analysis\Mechanical analysis\DMA_HFS_complete.csv';
%path = 'C:\Users\massimiliano.berardi\Downloads\AFM_REFORMATTED.csv'

% select method (the syntax between files is slightly different)
s_method = 'HFS';
% select line to model
s_line = 'T24';
% select timepoint to model
s_time = [3,14];

DATA_M = readtable(path);

%%

for j=1:length(s_time)
    %% filter by line and day
    
    DATA_M.line = categorical(DATA_M.line);
    idx_line = DATA_M.line == s_line;
    data_oneline = DATA_M(idx_line,:);
    
    data_oneline.day = categorical(data_oneline.day);
    idx_day = data_oneline.day == num2str(s_time(j));
    data_DMA_tab = data_oneline(idx_day,:);
    
    %% clean data and model
    
    % move form table to array, for convenience
    data_DMA = table2array(data_DMA_tab(:,4:6));
    % remove failed measurements
    fail_idx = isnan(data_DMA(:,2)) | isnan(data_DMA(:,3));
    rheo_data = data_DMA(~fail_idx,:);
    rheo_freq = rheo_data(:,1);
    if strcmp(s_method,'HFS')
        [rheo_stor, rheo_loss] = hydrodynamic_drag_correction_HFS(rheo_freq,rheo_data(:,2)/1000,rheo_data(:,3)/1000);
    else
        rheo_stor = rheo_data(:,2)/1000;
        rheo_loss = rheo_data(:,3)/1000;
    end
    rheo_compl = rheo_stor + 1i * rheo_loss;

    %%
    
    % do the fitting. 2PL require 4 pars, PT 6
    x0 = [17.9 0.1  .023 0.95 1000 0.8];
    lb = [9 .1 0 0.77 0 0];
    ub = [1e6 1 1e6 1 1e6 1];

  %  x0 = [0.1 0 800 1 .20 0 200 1];
 %   lb = [0 0 0 1 0 0 0 1];
%    ub = [1e6 0 1e6 1 1e6 0 1e6 1];


%      x0 = [17.9 0  .023 1 1 1000];
 %     lb = [0 0 0 1 0 1];
  %    ub = [1e6 1 1e6 1 1e6 1];
  %   x0 = [1 0.1 0.1 .8];
   %  lb = [0 0 0 0];
    % ub = [1e6 1 1e6 1];
%     
    options = optimoptions(@lsqnonlin,'Display','off','MaxFunctionEvaluations',3e4,'MaxIteration',1e3,'FunctionTolerance',1e-18,'StepTolerance',5e-16);
    [model_val,resnorm,resid,exitflag(j),~,~,jacobian] = lsqnonlin(@(x) fractional_poynting_thompson_objfun(x,rheo_freq,rheo_compl),x0,lb,ub,options);
    
    %calculate confidence bounds
    CI = nlparci(model_val,resid,'jacobian',jacobian); 
    % show me what's up
    figure(2)
    subplot(1,2,j)
    visuals_complex_model('PT','SE',rheo_freq,rheo_stor,rheo_loss,model_val,s_time(j));
    title(strcat(s_method,' - ',s_line,' - day ',num2str(s_time(j))))
end
CI(:,2)-model_val'
model_val