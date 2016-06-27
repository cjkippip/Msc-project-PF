function P = DifferentialLyapunov(F,Q,tspan,P0)
%%%------------------------------------------------------------------------
%
%   Description: This is calculation of the current prior covariance,
%                obtaining by calling differential Lyapunov equation.
%
%   Input: 
%           F: is the Jacobian of process model evaluated at the previous
%           posterior state estimate.
%           Q: is the noise covariance in process model.
%           tspan: The time interval for integrating
%           P0: the posterior covariance for previous state
%
%   Output:
%           P: is matrix for the posterior covariance, being formed by
%           integrating the differential Lyapunov equation respect to the
%           fixed time interval.
%%%-----------------------------------------------------------------------
p0=reshape(P0,length(F)^2,1);
opts=odeset('RelTol',1e-2,'AbsTol',1e-4);
[tsim,psim] = ode45(@(t,y)dreof(t,y,F,Q),tspan,p0,opts);

P=cell(length(tsim),1);

for ii=1:length(tsim)
    P{ii}=reshape(psim(ii,:),size(F));
end




