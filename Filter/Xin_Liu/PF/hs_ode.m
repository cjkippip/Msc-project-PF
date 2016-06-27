function hsState = hs_ode(t,y)
%%%----------------------------------------------------------------
%       
%   Description: This is the ODE functions for Heat Shock Model, including
%                three state variables and ten parameters. 
%   
%   Date: 14/09/2010
%
%   Author: Xin Liu
%
%%%---------------------------------------------------------------------

Dt = y(1);
St = y(2);
Uf = y(3);

%   Parameter declare
kd = 3; ad = 0.015; nt1 = 10; nt2 = 60; a0 = 0.03; as = 3; ks = 0.05;
ku = 0.0254; kt1 = 40; kt2 = 80; kfold = 6000; pt = 2*10^6;

if (t < 200)
    dDtdt = kd * (St/(1+(ks*Dt/(1+ku*Uf)))) - ad * Dt;
    dStdt = nt1 - a0*St - as * (ks*Dt/(1+ku*Uf))/(1 + ( ks*Dt/(1 + ku*Uf)))*St;
    Dufdt = kt1*(pt - Uf) - (kt1 + kfold)*Dt;
else
    dDtdt = kd * (St/(1+(ks*Dt/(1+ku*Uf)))) - ad * Dt;
    dStdt = nt2 - a0*St - as * (ks*Dt/(1+ku*Uf))/(1 + ( ks*Dt/(1 + ku*Uf)))*St;
    Dufdt = kt2*(pt - Uf) - (kt2 + kfold)*Dt;
end

hsState = [dDtdt; dStdt; Dufdt];