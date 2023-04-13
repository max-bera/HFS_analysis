function [Ecomp_mean] = visuals_complex_model(rheomodel,errortype,freq_list,Estor_list,Eloss_list,fit_params,day)
%fractional_poynting_thompson.m fits E* vs frequency with a lumped parameter model
%equivalent to two parallel springpots + a thid spring pot in series.
%
%INPUT (* are mandatory)
%model           *    either 2PL or PT
%errortype            either STD or SE, defaults to STD
%freq_list       *
%E_comp_list     * 
%fit_params      *
%OUTPUT
%visuals. Errorplot + fit overlay
%%

%obtain mean/SE values for plotting
unique_osc = unique(freq_list);
Ecomp_mean = zeros(length(unique_osc),1);
Ecomp_SE = zeros(length(unique_osc),1);


theor_freq = logspace(-2,3);

for i=1:length(unique_osc)
    idx_list = freq_list == unique_osc(i);
    Ecomp_mean(i) = mean(Estor_list(idx_list))...
        + 1i * mean(Eloss_list(idx_list));
    if strcmp(errortype,'SE')
        Ecomp_SE(i) = std(Estor_list(idx_list))/sqrt(length(Estor_list)/length(unique_osc))...
            + 1i * std(Eloss_list(idx_list))/sqrt(length(Estor_list)/length(unique_osc));
    else
        Ecomp_SE(i) = std(Estor_list(idx_list))...
            + 1i * std(Eloss_list(idx_list));
end


no_model_plot = false;
if strcmp('2PL',rheomodel)
    fun = @(x)x(1).*(1i*theor_freq).^x(2)+x(3).*(1i*theor_freq).^x(4);
   % fun = @(x)x(1).*cos(pi/2*x(2)).*(unique_osc).^x(2)+x(3).*cos(pi/2*x(3)).*(unique_osc).^x(4)+1i*(x(1).*sin(pi/2*x(2)).*(unique_osc).^x(2)+x(3).*sin(pi/2*x(3)).*(unique_osc).^x(4));
elseif strcmp('PT',rheomodel)
    fun = @(x)(x(5)*(1i.*theor_freq).^x(6)...
        .*(x(1)*(1i.*theor_freq).^x(2)+x(3)*(1i.*theor_freq).^x(4)))...
        ./(x(5)*(1i.*theor_freq).^x(6)+x(1)*(1i.*theor_freq).^x(2)+x(3)*(1i.*theor_freq).^x(4));
elseif strcmp('FB',rheomodel)
    fun = @(x)(x(1)*(1i.*theor_freq).^x(2).*(x(3)*(1i.*theor_freq).^x(4)))./...
    (x(1)*(1i.*theor_freq).^x(2)+x(3)*(1i.*theor_freq).^x(4))+...
    (x(5)*(1i.*theor_freq).^x(6).*(x(7)*(1i.*theor_freq).^x(8)))./...
    (x(5)*(1i.*theor_freq).^x(6)+x(7)*(1i.*theor_freq).^x(8))
else
    no_model_plot = true;
end

hold off
errorbar(unique_osc,real(Ecomp_mean),real(Ecomp_SE),real(Ecomp_SE),'ok')
hold on
errorbar(unique_osc,imag(Ecomp_mean),imag(Ecomp_SE),imag(Ecomp_SE),'*k')
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log';
ylim([10e-2 10e2])
grid on
xlabel('frequency [Hz]')
ylabel('E^* [Pa]')

if no_model_plot == false
    if day == 14
        storlinepat = '[0.4940 0.1840 0.5560]';
        losslinepat = '[0.3010 0.7450 0.9330]';
        disp_style = '-';
    else 
        storlinepat = '[0.6350 0.0780 0.1840]';
        losslinepat = '[0.4660 0.6740 0.1880]';
        disp_style = '--';
    end

    loglog(theor_freq,real(fun(fit_params)),...
        'Color',storlinepat,'LineStyle',disp_style);

    loglog(theor_freq,imag(fun(fit_params)),...
        'Color',losslinepat,'LineStyle',disp_style);

    xlim([10e-3 10e2])
    ylim([10e-1 10e1])
end

end

