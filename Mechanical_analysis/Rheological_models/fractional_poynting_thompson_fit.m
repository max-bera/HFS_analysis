function F = fractional_poynting_thompson_fit(params, oscfreq,Ecomp)
%fractional_poynting_thompson.m fits E* vs frequency with a lumped parameter model
%equivalent to two parallel springpots + a thid spring pot in series.
%
%INPUT (* are mandatory)
%osc_range  *    measured frequencies
%Ecomp       *    complex modulus
%
%OUTPUT
%fit function
%%
storagemod = real((params(5)*(1i.*oscfreq).^params(6).*(params(1)*(1i.*oscfreq).^params(2)+...
    params(3)*(1i.*oscfreq).^params(4)))./(params(5)*(1i.*oscfreq).^params(6)+...
    params(1)*(1i.*oscfreq).^params(2)+params(3)*(1i.*oscfreq).^params(4))-Ecomp);

lossmod = imag((params(5)*(1i.*oscfreq).^params(6).*(params(1)*(1i.*oscfreq).^params(2)+...
    params(3)*(1i.*oscfreq).^params(4)))./(params(5)*(1i.*oscfreq).^params(6)+...
    params(1)*(1i.*oscfreq).^params(2)+params(3)*(1i.*oscfreq).^params(4))-Ecomp);
 
F=[storagemod,lossmod];
 
end

