function J = Jacobian_HeatShock_allparams(y)
%%%------------------------------------------------------------------------
%
%   Description: This is evaluation of the current Jocobian of process
%                model.
%
%   Input: 
%           y: The previous posterior of state variables.
%
%   Output:
%           J: is current Jocobian matrix of system.
%%%-----------------------------------------------------------------------
% Recover state values
Dt=y(1);
St=y(2);
Uf=y(3);
alpha_s=y(4);
K_d=y(5);
K_s=y(6);
alpha_d=y(7);
alpha_0=y(8);
K_u=y(9);

% Preallocate Jacobian
Kfold=6000;
J=zeros(9,9);

% Set nonzero elements
J(1,1)=-(K_d*St*K_s+K_d*St*K_s*K_u*Uf+alpha_d+2*alpha_d*K_u*Uf+2*...
    alpha_d*K_s*Dt+alpha_d*K_u^2*Uf^2+2*alpha_d*K_u*Uf*K_s*Dt+alpha_d*...
    K_s^2*Dt^2)/(1+K_u*Uf+K_s*Dt)^2;
J(1,2)=K_d*(1+K_u*Uf)/(1+K_u*Uf+K_s*Dt);
J(1,3)=K_u*Dt*K_d*St*K_s/(1+K_u*Uf+K_s*Dt)^2;
J(1,5)=St*(1+K_u*Uf)/(1+K_u*Uf+K_s*Dt);
J(1,6)=-Dt*(1+K_u*Uf)*K_d*St/(1+K_u*Uf+K_s*Dt)^2;
J(1,7)=-Dt;
J(1,9)=Uf*Dt*K_d*St*K_s/(1+K_u*Uf+K_s*Dt)^2;
 
J(2,1)=-alpha_s*K_s*(1+K_u*Uf)*St/(1+K_u*Uf+K_s*Dt)^2;
J(2,2)=-(alpha_0+alpha_0*K_u*Uf+alpha_0*K_s*Dt+alpha_s*K_s*Dt)/...
    (1+K_u*Uf+K_s*Dt);
J(2,3)=K_u*St*alpha_s*K_s*Dt/(1+K_u*Uf+K_s*Dt)^2;
J(2,4)=-St*K_s*Dt/(1+K_u*Uf+K_s*Dt);
J(2,6)=-alpha_s*Dt*(1+K_u*Uf)*St/(1+K_u*Uf+K_s*Dt)^2;
J(2,8)=-St;
J(2,9)=Uf*St*alpha_s*K_s*Dt/(1+K_u*Uf+K_s*Dt)^2;
 
J(3,1)=-40-Kfold;
J(3,3)=-40;
