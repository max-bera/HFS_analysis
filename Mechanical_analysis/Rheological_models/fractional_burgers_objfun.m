function F = fractional_burgers_objfun(x, w, Ecomp)
%burgers_objfun.m fits E* vs frequency with a lumped parameter model
%equivalent a SLS + a dashpot in series, but every term is fractional.
%
%INPUT (* are mandatory)
%osc_range  *    measured frequencies
%Ecomp       *    complex modulus
%
%OUTPUT
%fit function
%%
storagemod = real((x(1)*(1i.*w).^x(2).*(x(3)*(1i.*w).^x(4)))./...
    (x(1)*(1i.*w).^x(2)+x(3)*(1i.*w).^x(4))+...
    (x(5)*(1i.*w).^x(6).*(x(7)*(1i.*w).^x(8)))./...
    (x(5)*(1i.*w).^x(6)+x(7)*(1i.*w).^x(8))-Ecomp);

lossmod = imag((x(1)*(1i.*w).^x(2).*(x(3)*(1i.*w).^x(4)))./...
    (x(1)*(1i.*w).^x(2)+x(3)*(1i.*w).^x(4))+...
    (x(5)*(1i.*w).^x(6).*(x(7)*(1i.*w).^x(8)))./...
    (x(5)*(1i.*w).^x(6)+x(7)*(1i.*w).^x(8))-Ecomp);
 
F=[storagemod;lossmod];
 
end

