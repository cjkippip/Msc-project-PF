function hsState = hs_ode(t,y)
%%%-----------------------------------------------------------------------
%       
%   Description: This is the ODE functions for Heat Shock Model, including
%                three state variables, six fixed values parameters and 
%                six time-dependent parameters
%   
%   
%   Input: Three initial values for state variables have to be input, 
%          they are for Dt(0), St(0) and Uf(0).
%
%   Output: A vector namely hsState, which contains the integration respect
%           to the input time interval.
%
%   Date: 14/09/2010
%
%   Author: Xin Liu
%   
%%%------------------------------------------------------------------------

Dt = y(1);
St = y(2);
Uf = y(3);
%   Parameter declare
kd = 3; ad = 0.015; nt1 = 10; nt2 = 60; a0 = 0.03; as = 3; ks = 0.05;
ku = 0.0254; kt1 = 40; kt2 = 80; pt1 = 6000; pt2 = 2*10^6;

if (t < 200)
    dDtdt  = kd * (St/(1+(ks*Dt/(1+ku*Uf)))) - ad * Dt;
    dStdt  = nt1 - a0*St - as * (ks*Dt/(1+ku*Uf))/(1 + ( ks*Dt/(1 + ku*Uf)))*St;
    Dufdt  = kt1*(pt2 - Uf) - (kt1 + pt1)*Dt;
else
    dDtdt  = kd * (St/(1+(ks*Dt/(1+ku*Uf)))) - ad * Dt;
    dStdt  = nt2 - a0*St - as * (ks*Dt/(1+ku*Uf))/(1 + ( ks*Dt/(1 + ku*Uf)))*St;
    Dufdt  = kt2*(pt2 - Uf) - (kt2 + pt1)*Dt;
end

hsState = [dDtdt; dStdt; Dufdt];