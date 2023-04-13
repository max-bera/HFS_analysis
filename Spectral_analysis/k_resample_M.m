function [raw_data_resample, k_even_sp] = k_resample_M(wl_data, k_space)
%K_SPACE_RESAMPLE takes Deltasense data and outputs evenly spaced raw data
%in k space 
%   INPUT ARGUMENTS:
%            -deltasense structure;
%            -pointer (int>=1, defines which time point to analyze).

if nargin<2, help K_resample_M; return; end

line = zeros(2,length(wl_data));
wl_data_corr = wl_data;

for j=1:length(wl_data)
    line(1,j) = (wl_data(end,j)-wl_data(1,j))./(k_space(end)-k_space(1));
    line(2,j) =(wl_data(1,j));
    wl_data_corr(:,j) = wl_data(:,j)-polyval(line(:,j),k_space)';
end

[raw_data_resample, k_even_sp] = resample(wl_data_corr,k_space,...
    'Dimension',1);

for j=1:length(wl_data)
    raw_data_resample(:,j) = raw_data_resample(:,j) + ...
        polyval(line(:,j),k_even_sp)';
end

end