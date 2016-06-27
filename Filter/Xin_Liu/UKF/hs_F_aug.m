function x_next = hs_F_aug(x)
%%%------------------------------------------------------------------------
%
%   Description: This is process propogation for system. The extended state
%   space is updating by integrating respect to fixed sampling interval.
%   The starting point is the previous posterior states.
%
%   Input: 
%           x: is the extended state space. The first row until par_dim is
%           representing the parameter estimation. The rest are state
%           variables.
%
%   Output:
%           x_next: the matrix is the posterior states.
%   
%
%   Author: Xin Liu
%
%   Date:   23/04/2011
%
%%%------------------------------------------------------------------------

% Augmented state function

global deltaT
global para_dim;
global state_dim;

%   Initialization of parameters
p = x(1:para_dim,:);  

%   Assign memory for state fitering 
x_next = zeros(size(x));

x_next(1:para_dim,:) = p;  % Keep the parameter unchange
options = odeset('RelTol',1e-02);   % the real error tolerance is 1e-02

for i=1:size(x,2)
    %   Obtain solutions from integrating heat shock model, respect to
    %   sampling frequency.
    sol = ode45(@(t,x) hs_odeUKF_allPara(t,x,p(:,i)),[0 deltaT],...
                            x(para_dim+1 : para_dim+state_dim,i),options);
    %   Update state                    
    x_next(para_dim+1 : para_dim+state_dim,i) = deval(sol,deltaT);
end

