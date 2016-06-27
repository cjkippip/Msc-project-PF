function hsState = hs_odeEKF_Single_Para(t,y)
%%%-----------------------------------------------------------------------
%       
%   Description: This is the ODE functions for Heat Shock Model used for
%                generating the synthetic data. 
%                This model includes three state variables and six 
%                parameters. This version is for single parameter unknown
%                (Ks)
%   
%   Input:       t: Time interval for integration
%                y: Initial values for integration.
%
%   Output:      hsState is the numerical solution for heat shock ODEs, 
%                The first three are state variables. The rest are
%                parameters.
%   
%   Date: 23/04/2011
%
%   Author: Xin Liu
%
%%%-----------------------------------------------------------------------

Dt = y(1);
St = y(2);
Uf = y(3);
ks = y(4);

%   Parameter declare
nt1 = 10; nt2 = 60; as = 3; kd = 3; ad = 0.015; a0 = 0.03;
kt1 = 40; kt2 = 80; kfold = 6000; pt = 2*10^6; ku = 0.0254;

if (t < 200)
    dDtdt = kd * (St/(1+(ks*Dt/(1+ku*Uf)))) - ad * Dt;
    
    dStdt = nt1 - a0*St - as * (ks*Dt/(1+ku*Uf))/...
            (1 + ( ks*Dt/(1 + ku*Uf)))*St;
        
    Dufdt = kt1*(pt - Uf) - (kt1 + kfold)*Dt;
else
    dDtdt = kd * (St/(1+(ks*Dt/(1+ku*Uf)))) - ad * Dt;
    
    dStdt = nt2 - a0*St - as * (ks*Dt/(1+ku*Uf))/...
            (1 + ( ks*Dt/(1 + ku*Uf)))*St;
        
    Dufdt = kt2*(pt - Uf) - (kt2 + kfold)*Dt;
end

hsState = [dDtdt; dStdt; Dufdt;zeros(1,1)];