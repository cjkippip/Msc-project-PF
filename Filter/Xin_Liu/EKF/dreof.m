function yp=dreof(t,y,F,Q)
%%%------------------------------------------------------------------------
%
%   Description: This is differential Lyapunov equation.
%
%   Input: 
%           F: is the Jacobian of process model evaluated at the previous
%           posterior state estimate.
%           Q: is the noise covariance in process model.
%           tspan: The time interval for integrating
%           P0: the posterior covariance for previous state
%
%   Output:
%           yp: posterior covariance matrix
%%%-----------------------------------------------------------------------
Y=reshape(y,size(F));
Y=(Y+transpose(Y))./2; % Symmetrize
YP=F*Y+Y*transpose(F)+Q;

yp=reshape(YP,length(F)^2,1);