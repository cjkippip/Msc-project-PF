
function [xhat,Pxx,K] = ukf(xhat,Pxx,y,state_fct,obs_fct,P_obs_noise)
%%%------------------------------------------------------------------------
%
%   Description: This is implementation of UKF.
%
%   Input: 
%           xhat: The current prior state.
%           Pxx : The current prior covariance.
%           y   : The generated synthetic data, which is used in correction
%                 step.
%           state_fct: The process propogation.
%           obs_fct: The observation function.
%           xhat: Observation noise.
%
%   Output:
%           xhat: The current posterior state.
%           Pxx : The current posterior state.
%           K   : The current Kalman gain.
%   
%
%   Author: Xin Liu
%
%   Date:   23/04/2011
%
%%%------------------------------------------------------------------------

global aug_state_dim;   % dimension of the augmented state
global obs_dim;         % dimension of observations

  N=2*aug_state_dim;    % number of sigma points
  
  xsigma = chol( aug_state_dim*Pxx )'; 
  
  Xa = xhat*ones(1,N) + [xsigma, -xsigma];
   
  X = feval(state_fct,Xa);  % No need to input time interval, because the 
                            % state_fct calculate the integration with 
                            % respect to deltaT
  
  xtilde = sum(X')'/N;

  Pxx = zeros(aug_state_dim,aug_state_dim);
  
  for i=1:N;
    Pxx = Pxx +(X(:,i)-xtilde)*(X(:,i)-xtilde)'/N;
  end;
  
  Y=feval(obs_fct,X);
   
  ytilde=sum(Y')'/N;
  
  Pyy = P_obs_noise;
  
  for i=1:N;
    Pyy = Pyy+(Y(:,i)-ytilde)*(Y(:,i)-ytilde)'/N;
  end;
  
  Pxy = zeros(aug_state_dim,obs_dim); 
  
  for i = 1:N;
    Pxy = Pxy + (X(:,i)-xtilde)*(Y(:,i)-ytilde)'/N;
  end;
    
  K = Pxy * inv(Pyy);
  
  xhat = xtilde + K*(y - ytilde);
 
  Pxx =Pxx - K*Pxy';
  
