function visuals_complex_model(rheomodel,freq_list,Estor_list,Eloss_list,fit_params)
%fractional_poynting_thompson.m fits E* vs frequency with a lumped parameter model
%equivalent to two parallel springpots + a thid spring pot in series.
%
%INPUT (* are mandatory)
%model           *    either 2PL or PT
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

for i=1:length(unique_osc)
    idx_list = freq_list == unique_osc(i);
    Ecomp_mean(i) = mean(Estor_list(idx_list))...
        + 1i * mean(Eloss_list(idx_list));
    Ecomp_SE(i) = std(Estor_list(idx_list))/sqrt(sum(idx_list))...
        + 1i * std(Eloss_list(idx_list))/sqrt(sum(idx_list));
end



if strcmp('2PL',rheomodel)
    fun = @(x)x(1).*(1i*unique_osc).^x(2)+x(3).*(1i*unique_osc).^x(4);
elseif strcmp('PT',rheomodel)
    fun = @(x)(x(5)*(1i.*unique_osc).^x(6).*(x(1)*(1i.*unique_osc).^x(2)+...
    x(3)*(1i.*unique_osc).^x(4)))./(x(5)*(1i.*unique_osc).^x(6)+...
    x(1)*(1i.*unique_osc).^x(2)+x(3)*(1i.*unique_osc).^x(4));
end

figure(1)
errorbar(unique_osc,real(Ecomp_mean),real(Ecomp_SE),real(Ecomp_SE),'ok')
hold on
errorbar(unique_osc,imag(Ecomp_mean),imag(Ecomp_SE),imag(Ecomp_SE),'*k')
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log';
grid on
xlabel('frequency [Hz]')
ylabel('E^* [Pa]')

loglog(unique_osc,real(fun(fit_params)+Ecomp_mean),'g-');
loglog(unique_osc,imag(fun(fit_params)+Ecomp_mean),'r--');
end

