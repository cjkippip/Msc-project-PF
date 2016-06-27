function hsState = hs_odeSmc_All_Para(t,y,paraSmc)
Dt = y(1);
St = y(2);
Uf = y(3);

%   Allocation of parameter space
%   paraSmc(1,:) is Ku
%   paraSmc(2,:) is As
%   paraSmc(3,:) is Ks
%   paraSmc(4,:) is A0
%   paraSmc(5,:) is Kd
%   paraSmc(6,:) is Ad

%   Parameter declare
nt1 = 10; nt2 = 60; 
kt1 = 40; kt2 = 80; kfold = 6000; pt = 2*10^6;

if (t < 200)
    dDtdt = paraSmc(5,:) * (St/(1+(paraSmc(3,:)*Dt/(1+paraSmc(1,:)*Uf)))) - paraSmc(6,:) * Dt;
    dStdt = nt1 - paraSmc(4,:)*St - paraSmc(2,:) * (paraSmc(3,:)*Dt/(1+paraSmc(1,:)*Uf))/(1 + ( paraSmc(3,:)*Dt/(1 + paraSmc(1,:)*Uf)))*St;
    Dufdt = kt1*(pt - Uf) - (kt1 + kfold)*Dt;
else
    dDtdt = paraSmc(5,:) * (St/(1+(paraSmc(3,:)*Dt/(1+paraSmc(1,:)*Uf)))) - paraSmc(6,:) * Dt;
    dStdt = nt2 - paraSmc(4,:)*St - paraSmc(2,:) * (paraSmc(3,:)*Dt/(1+paraSmc(1,:)*Uf))/(1 + ( paraSmc(3,:)*Dt/(1 + paraSmc(1,:)*Uf)))*St;
    Dufdt = kt2*(pt - Uf) - (kt2 + kfold)*Dt;
end

hsState = [dDtdt; dStdt; Dufdt];











