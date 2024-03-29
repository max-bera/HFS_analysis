function F = fractional_poynting_thompson_objfun(x, w, Ecomp)
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
storagemod = real((x(5)*(1i.*w).^x(6).*(x(1)*(1i.*w).^x(2)+...
    x(3)*(1i.*w).^x(4)))./(x(5)*(1i.*w).^x(6)+...
    x(1)*(1i.*w).^x(2)+x(3)*(1i.*w).^x(4))-Ecomp);

lossmod = imag((x(5)*(1i.*w).^x(6).*(x(1)*(1i.*w).^x(2)+...
    x(3)*(1i.*w).^x(4)))./(x(5)*(1i.*w).^x(6)+...
    x(1)*(1i.*w).^x(2)+x(3)*(1i.*w).^x(4))-Ecomp);
 
F=[storagemod;lossmod];
 
end

