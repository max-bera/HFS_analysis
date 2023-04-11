function F = double_power_law_objfun(x, w, Ecomp)
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
storagemod = real(x(1).*(1i*w).^x(2)+x(3).*(1i*w).^(x(4))-Ecomp);

lossmod = imag(x(1).*(1i*w).^x(2)+x(3).*(1i*w).^(x(4))-Ecomp);

F=[storagemod;lossmod];
end