%%%------------------------------------------------------------------------
%
%   Description: This is the implementation for running UKF on Heat Shock,
%   which refer to Florence's work.
%
%   Author: Xin Liu
%
%   Date:   23/04/2011
%
%%%-----------------------------------------------------------------------
clear all;
clc
echo off;
% === PARAMETER DECLARATION

global para_dim;        % Number of parameters to be estimated
global state_dim;       % Dimension of the state vector
global aug_state_dim;   % Dimension of the augmented state
global obs_index;       % Index of observation vector
global deltaT;          % Sampling time interval
global obs_dim;         % Dimension of observations

% Dimensions
para_dim      = 6;
state_dim     = 3;
obs_dim       = 2;
obs_index     = [1 3];
aug_state_dim = state_dim + para_dim;

T = 200;                % the time duration for integration
deltaT = 0.2;
L = T/deltaT + 1;       % Number of data samples

% True parameter value
p_true = [3 3 0.05 0.015 0.03 0.0254];    % The values for parameters:
                                          % as kd ks ad a0 ku.
                                          
                                          
startTime = cputime;   % Measure time

%% Generate True Data

% Allocate memory space for generated state vector
x_true = zeros(state_dim,L);    

x_true(:,1) = [0 0 0]';   % Initial state values

options = odeset('RelTol',1e-02);   % the real error tolerance is 1e-02

%   The numerical solutions are viewed as synthetic data
sol = ode45(@hs_ode,[0 T],x_true(:,1)',options);

x_true = deval(sol,0:deltaT:T);

%   Covariance matrix for observation noise
P_obs_noise = 0.25 * diag(var(x_true(obs_index,:)'));

y_true = x_true(obs_index,:) + sqrtm(P_obs_noise) * randn(obs_dim,L);   % Noisy data

%%  UKF Filtering
% -------------- Initialization of Parameters -----------------------------
p_init = [4 0 1 1 1 1];  % Initialization of parameter estimation

P_params = 0.5 * diag(p_true);               

P_params = P_params * P_params;      % Covariance matrix of parameters

% Initialization the augmented state

x_init = zeros(aug_state_dim,1);

x_init(1:para_dim,1) = p_init;    % First input the parameter value

% the rest component in x comes from true observation

x_init(para_dim + obs_index,:) = y_true(:,1);  

% =================== Initialization =========================

% Allocate memory for state variables
mu_ukf = zeros(aug_state_dim,L);    

% Initialization of state variables
mu_ukf(:,1) = x_init;   

% Allocate memory for error covariance for first time step 
Pxx_init = zeros(aug_state_dim,aug_state_dim);  

%   Initialization of covariance for parameter estimation
Pxx_init(1:para_dim,1:para_dim) = diag([1e-1 1e-1  1e-1  1e-1  1e-1  1e-1]);

%   Initialization of covariance for state
Pxx_init(para_dim+1:para_dim+state_dim,para_dim+1:...
         para_dim+state_dim) = 0.25*diag(var((x_true)'));

% Allocate memory for error covariance for all time steps 
P_ukf = zeros(aug_state_dim,aug_state_dim,L);   

% First element is assigned to initialization
P_ukf(:,:,1) = Pxx_init;    

% Kalman gain
K = zeros(aug_state_dim,obs_dim,L);     

% =================== Main Loop for UKF ============================

% Heat shock model by associating augmented states
state_fct = 'hs_F_aug';     

% Observation model by associating augmented states
obs_fct   = 'hs_H_aug';

for k = 2:L
    [mu_ukf(:,k),P_ukf(:,:,k)]=ukf(mu_ukf(:,k-1),P_ukf(:,:,k-1),y_true(:,k),...
                                   state_fct,obs_fct,P_obs_noise);
    
    k
end


finishTime = cputime - startTime;
%%  Plot Graph
figure,
plot(3*ones(1,1000),'b--','LineWidth',2);hold on;
plot(mu_ukf(1,:),'r--','LineWidth',2);hold on;
xlabel('Time','Fontsize',20);
ylabel('Alpha_s','Fontsize',20);
set(gca,'Fontsize',20);
legend('True','UKF');
xlim([0 1000])

figure,
plot(3*ones(1,1000),'b--','LineWidth',2);hold on;
plot(mu_ukf(2,:),'r--','LineWidth',2);hold on;
xlabel('Time','Fontsize',20);
ylabel('K_D','Fontsize',20);
set(gca,'Fontsize',20);
legend('True','UKF');
xlim([0 1000])

figure,
plot(0.05*ones(1,1000),'b--','LineWidth',2);hold on;
plot(mu_ukf(3,:),'r--','LineWidth',2);hold on;
xlabel('Time','Fontsize',20);
ylabel('K_S','Fontsize',20);
set(gca,'Fontsize',20);
legend('True','UKF');
xlim([0 1000])

figure,
plot(0.015*ones(1,1000),'b--','LineWidth',2);hold on;
plot(mu_ukf(4,:),'r--','LineWidth',2);hold on;
xlabel('Time','Fontsize',20);
ylabel('Alpha_D','Fontsize',20);
set(gca,'Fontsize',20);
legend('True','UKF');
xlim([0 1000])

figure,
plot(0.03*ones(1,1000),'b--','LineWidth',2);hold on;
plot(mu_ukf(5,:),'r--','LineWidth',2);hold on;
xlabel('Time','Fontsize',20);
ylabel('Alpha_0','Fontsize',20);
set(gca,'Fontsize',20);
legend('True','UKF');
xlim([0 1000])

figure,
plot(0.0254*ones(1,1000),'b--','LineWidth',2);hold on;
plot(mu_ukf(6,:),'r--','LineWidth',2);hold on;
xlabel('Time','Fontsize',20);
ylabel('K_U','Fontsize',20);
set(gca,'Fontsize',20);
legend('True','UKF');
xlim([0 1000])






