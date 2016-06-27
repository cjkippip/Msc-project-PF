function hsStateUKF = hs_odeUKF_allPara(t,x,para)
%%%----------------------------------------------------------------
%       
%   Description: This is the ODE functions for Heat Shock Model used for
%                parameter estimation by UKF.
%                This model includes three state variables and six
%                parameters. 
%                This version is for all parameters unknown.
%   
%   Input:       
%                t:     The time for running the experiment
%                x:     The values for underlying states.
%                para:  The values for recent parameters
%   Output:
%                hsStateUKF: The numerical solution for integrating heat
%                shock model, respect to sampling interval.
%
%   Date: 23/04/2011
%
%   Author: Xin Liu
%
%%%---------------------------------------------------------------------

Dt = x(1);
St = x(2);
Uf = x(3);
% 
% para(1,:) = as;
% para(2,:) = kd;
% para(3,:) = ks;
% para(4,:) = ad;
% para(5,:) = a0;
% para(6,:) = ku;

%   Parameter declare
nt1 = 10; nt2 = 60;
kt1 = 40; kt2 = 80; kfold = 6000; pt = 2*10^6;

if (t < 200)
    dDtdt = para(2,:) * (St/(1+(para(3,:)*Dt/(1+para(6,:)*Uf)))) ...
            - para(4,:) * Dt;
        
    dStdt = nt1 - para(5,:)*St - para(1,:) * (para(3,:)*Dt/(1+para(6,:)...
            *Uf))/(1 + ( para(3,:)*Dt/(1 + para(6,:)*Uf)))*St;
        
    Dufdt = kt1*(pt - Uf) - (kt1 + kfold)*Dt;
else
    
    dDtdt = para(2,:) * (St/(1+(para(3,:)*Dt/(1+para(6,:)*Uf))))...
            - para(4,:) * Dt;
        
    dStdt = nt2 - para(5,:)*St - para(1,:) * (para(3,:)*Dt/(1+para(6,:)...
            *Uf))/(1 + ( para(3,:)*Dt/(1 + para(6,:)*Uf)))*St;
        
    Dufdt = kt2*(pt - Uf) - (kt2 + kfold)*Dt;
end

hsStateUKF = [dDtdt; dStdt; Dufdt];