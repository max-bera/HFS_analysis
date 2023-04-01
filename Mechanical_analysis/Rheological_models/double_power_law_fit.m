function F = double_power_law_fit(params, oscfreq, Ecomp)
%double_power_law_fit.m fits E* vs frequency with a lumped parameter model
%equivalent to two parallel springpots. E* = A*(iw)^alpha)+B*(iw)^beta
%
%INPUT (* are mandatory)
%oscfreq  *    measured frequencies
%Ecomp       *    complex modulus
%
%OUTPUT
%fit function
%%
storagemod = real(params(1).*(1i*oscfreq).*params(2)+params(3).*(1i*oscfreq).*(params(4))-Ecomp);
lossmod = imag(params(1).*(1i*oscfreq).*params(2)+params(3).*(1i*oscfreq).*(params(4))-Ecomp);
 
F=[storagemod,lossmod];
 
end