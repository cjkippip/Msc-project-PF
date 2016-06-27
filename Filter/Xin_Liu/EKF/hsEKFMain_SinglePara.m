%%%-----------------------------------------------------------------------
%
%   Description: The parameter estimation used EKF. 
%                This implementation is for single parameters (Ks)
%                unknown case. There are three sections in this programme,
%                which are generating synthetic data, filtering and plot.
%
%   Author: Xin Liu
%
%   programme goes through
%
%   Date:   23/04/2011
%
%--------------------------------------------------------------------------
clear all;
clc;
echo off;

startTime = cputime;    % For measuring time
p = 0.05;   % The true value of parameter for parameter Ks

global state_dim;       % Dimension for state variables
global para_dim;        % Dimension for parameter variables
global aug_state_dim;   % Dimension for augmented state space
global deltaT;          % Sampling frequency

para_dim      = length(p);
state_dim     = 3;
obs_dim       = 2;  % the dimension of observation
obs_index     = [1 3];  % the index of observation
aug_state_dim = state_dim + para_dim;   % the dimension of augmented...
                                        % state space
                                        
deltaT = 0.2;   
%%  Generate True State

init = [0 0 0 0.05];    % true values for Dt St Uf...
                                              % as kd ks ad a0 ku. First
                                              % three are states and rest
                                              % of them are parameters.
                                              
options = odeset('RelTol',1e-02);   % the real error tolerance is 1e-02

%   Generate the synthetic data
sol     = ode45(@hs_odeEKF_Single_Para,[0 200],init,options);                                                                
tspan   = linspace(0,200,1000);
xtrue   = deval(sol,tspan);

%   The observation function, only the first and third one can be observed
% h = [1 0 1 0 0 0 0 0 0];

%   This is the observation noise
P_obs_noise = 0.25 * diag(var(xtrue(obs_index,:)'));    

%   Only first and third state is observable
ytrue   = xtrue(obs_index,:);

%   The covariance matrix for observation noise.
R       = diag([(mean(ytrue(1,:)/3)*(5/100))^2 (mean(ytrue(2,:)/3)*(50/100))^2]);
ytrue(1,:) = ytrue(1,:) + sqrt(R(1,1))*randn(1,size(ytrue,2));
ytrue(2,:) = ytrue(2,:) + sqrt(R(2,2))*randn(1,size(ytrue,2));

%%  EKF Algorithm

Q = 0.0000001 * eye(aug_state_dim);    % the process noise 0.0000001

%   The initialization for state variables
xe_init(1:state_dim,:) = [0 0 0]';

%   The initialization for parameter Ks
xe_init(state_dim+1:para_dim+state_dim,:) = 10; 

%   Assign the memory for filtering
xe_ekf = zeros(aug_state_dim,size(ytrue,2));

%   First column is initialization
xe_ekf(:,1) = xe_init;

%   The initial covariance matrix of EKF
Pe_init = diag([25 25 25 0.01]);

%   Assign memory for covariance matrices
Pe = zeros(aug_state_dim,aug_state_dim,size(ytrue,2));

%   First element is the initialization
Pe(:,:,1) = Pe_init;

%   Jacobian matrix for observation model
H = [1 0 0;
     0 0 1];
H = [H zeros(size(H,1),para_dim)];

I  = eye(aug_state_dim);

for t = 2:size(ytrue,2)

    % **1** Compute linearization around previous a posteriori estimate
    Fk = Jacobian_HeatShock_singleParams(xe_ekf(:,t-1));
    
    % **2** Advance time
    solEKF = ode45(@hs_odeEKF_Single_Para,[0 deltaT],xe_ekf(:,t-1),options);
    xp = deval(solEKF,deltaT);
    
    %   prior covariance is up to integrate a differential Lyapunov
    %   equation using the previous posterior covariance as initialization
    Pp = DifferentialLyapunov(Fk,Q,[0 deltaT],Pe(:,:,t-1));
    Pp=Pp{numel(Pp)};
    Pp=(Pp+Pp')/2; % Symmetrize
    
    % **3** Compute optimal gain (try to avoid inverse...)
    Lk = Pp*H'*(H*Pp*H'+R)^(-1);
    
    % **4** Incorporate currente measurement
    xe_ekf(:,t)= xp + Lk*(ytrue(:,t)-H*xp);
    Pe(:,:,t)=(I-Lk*H)*Pp*(I-Lk*H)'+Lk*R*Lk';
    Pe(:,:,t)=(Pe(:,:,t) + Pe(:,:,t)')/2; % Symmetrize
  
    t
    
    
end

finishTime = cputime - startTime

%%  Plot Graphs
figure,
plot(0.05*ones(1,size(ytrue,2)),'b--','LineWidth',2);hold on;
plot(xe_ekf(4,:),'r--','LineWidth',2);hold on;
xlabel('Time','Fontsize',20);
ylabel('K_{s}','Fontsize',20);
legend('True','EKF');
set(gca,'Fontsize',20);



