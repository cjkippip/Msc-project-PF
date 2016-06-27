% Project Msc 
% Multi-model version
% SIR estimate all parameters, one guess of initial state
%%
ID_mdl='BIOMD0000000002'; % model ID
%{ 
run ext_SBML first
need **libSBML** and **SBMLToolbox-4.1.0**
website: http://sbml.org/Downloads
%}
% run ext_SBML

% if no **libSBML** and **SBMLToolbox-4.1.0**
% run the below statement
addpath (genpath(['./model/',ID_mdl]))
load (['./model/',ID_mdl,'/BIOMD02_para'])
%% Parameter
% Global variables declaration
global N;
global a;
global covParaCoefVect;

% initial states
state_init = x0;

state_dim = length(state_init);         % state dimension
para_dim = length(allParaValue);        % parameter dimension
ext_state_dim = state_dim+para_dim;     % extended dimension
obs_dim = 3;                            % observed dimension
obs_index = [1 3 5];                    % index of observation
N = 500;                                % number of particles {{***AD***}}


% time length and delta_t {{***AD***}}
tspan = 0:0.1:300; % time range
delta_t = tspan(2)-tspan(1); % delta time
T = length(tspan); % length of iteration time

% options set
% opts = odeset('AbsTol',1e-3);
opts = odeset('RelTol',1e-4);

% solve ode
[~,trueState] = ode23tb(@BIOMD0002,tspan,state_init,opts);
trueState = trueState';
% plot(tspan,trueState);

% perturbation of parameter {{***AD***}}
e = 0.96;  
a = (3*e-1)/(2*e);
covParaCoef = 1-a^2;
covParaCoefVect = covParaCoef*ones(para_dim,1);

% Convariance of state {{***AD***}}
covState=100;

% Convariance of measurement {{***AD***}}
covMeasure = 0.25*diag(var(trueState(obs_index,:),0,2));

% Convariance of particles of state {{***AD***}}
covP_State=100;

%% Generate true state, parameter and observation
% true states
for i=1:state_dim
    trueState(i,:)=trueState(i,:) + sqrt(covState)*randn(1,T);
end

% true parameters
truePara = ones(para_dim,T);
for i=1:para_dim
    truePara(i,:) = allParaValue(i)*ones(1,T); 
end

% generate observation data
trueOb=ones(obs_dim,T);
for i=1:obs_dim
    trueOb(i,:)=trueState(obs_index(i),:)+sqrt(covMeasure(i,i))*randn(1,T);
end
%% Figure to show parameter
% Figure 1 for parameter 1-3
figure(1);
subplot(3,1,1);
h1 = plot(zeros(2,2));title('Parameter 1','FontSize',11);
set(h1(1),'XData',1:size(truePara,2),'YData',truePara(1,:),'LineStyle',...
           '--','color','blue','LineWidth',0.75);
xlabel('Time');
ylabel('Parameter 1');

subplot(3,1,2);
h2 = plot(zeros(2,2));title('Parameter 2','FontSize',11);
set(h2(1),'XData',1:size(truePara,2),'YData',truePara(2,:),'LineStyle',...
           '--','color','blue','LineWidth',0.75);
xlabel('Time');
ylabel('Parameter 2');

subplot(3,1,3);
h3 = plot(zeros(2,2));title('Parameter 3','FontSize',11);
set(h3(1),'XData',1:size(truePara,2),'YData',truePara(3,:),'LineStyle',...
           '--','color','blue','LineWidth',0.75);
xlabel('Time');
ylabel('Parameter 3');

% Figure 2 for parameter 4-6
figure(2);
subplot(3,1,1);
h4 = plot(zeros(2,2));title('Parameter 4','FontSize',11);
set(h4(1),'XData',1:size(truePara,2),'YData',truePara(4,:),'LineStyle',...
           '--','color','blue','LineWidth',0.75);
xlabel('Time');
ylabel('Parameter 4');

subplot(3,1,2);
h5 = plot(zeros(2,2));title('Parameter 5','FontSize',11);
set(h5(1),'XData',1:size(truePara,2),'YData',truePara(5,:),'LineStyle',...
           '--','color','blue','LineWidth',0.75);
xlabel('Time');
ylabel('Parameter 5');

subplot(3,1,3);
h6 = plot(zeros(2,2));title('Parameter 6','FontSize',11);
set(h6(1),'XData',1:size(truePara,2),'YData',truePara(6,:),'LineStyle',...
           '--','color','blue','LineWidth',0.75);
xlabel('Time');
ylabel('Parameter 6');
%% SMC
%=======================================================
% initial particles of parameters
P_Para=ones(para_dim,N);

% gaussian {{***AD***}}
mean1=2;
var1=2;
for i=1:para_dim
%     P_Para(i,:)=allParaValue(i)+sqrt(var1)*randn(1,N);
    P_Para(i,:)=2*allParaValue(i)+sqrt(var1)*randn(1,N);
end

% uniform {{***AD***}}
% marginal=10;
% for i=1:para_dim
%     P_Para(i,:)=unifrnd(-marginal,marginal,[1,N]);
% end
%=======================================================

% memory for particles of parameters in all iteration
P_ParaSpace=ones(para_dim,N,T);
P_ParaSpace(:,:,1)=P_Para; % 1st iter

% memory for means of particles(parameters) in all iteration
meanP_Para=ones(para_dim,T);
meanP_Para(:,1)=mean(P_Para,2); % 1st iter

% covariance of particles(parameters) 1st iter
covP_Para=var(P_Para,0,2);

%=======================================================
% initial particles of state
P_State=sqrt(25)*randn(state_dim,N);
%=======================================================

% memory for particles(states) in all iteration
P_StateSpace=ones(state_dim,N,T);
P_StateSpace(:,:,1)=P_State; % 1st iter

% memory for means of particles(states) in all iteration
meanP_State=ones(state_dim,T);
meanP_State(:,1)=mean(P_State,2); % 1st iter

% memory for weights and normal weights of particles in all iteration
P_WeightSpace = ones(1,N,T);
P_Weight_NorSpace = ones(1,N,T);

% weights and normal weights
P_Weight=ones(1,N)/N;
P_Weight_Nor=ones(1,N)/N;

% number of effective particles
Neff=ones(1,T);

% == SMC main loop ==
for t=2:T
    % update states by integrating previous states respect to sampling
    % interval
    
    for n=1:N
        [~,y] = ode23tb(@(t,y) BIOMD0002_SMC(t,y,P_Para(:,n)),...
            [0 delta_t],P_State(:,n),opts);
        y=y';     
%         y = deval(solSmc,delta_t);
        P_State(:,n) = y(:,end);
        fprintf('%d particle %d iteration \n',n,t);
    end
    
    % add noise to underlying state
    for i=1:state_dim
        P_State(i,:) = P_State(i,:) + sqrt(covP_State)*randn(1,N);
    end    
    
    % update parameters (perturbation)
    P_Para=updatePara_all_para(P_Para,covP_Para,meanP_Para(:,t-1),para_dim);
    
    % weight, likelihood      
    for n=1:N
        P_Ob = P_State([1,3,5],n);
        P_Weight(n) = exp( -0.5 * ...
            (trueOb(:,t)-P_Ob)' / ...
            (covMeasure) * ...
            (trueOb(:,t)-P_Ob) ) + eps;      
    end
    P_WeightSpace(:,:,t)=P_Weight;
    
    % normalise weight
    P_Weight_Nor=P_Weight/sum(P_Weight);
    P_Weight_NorSpace(:,:,t)=P_Weight_Nor;
    %% Resampling
    
    Neff(t)=1/sum(P_Weight_Nor.*2);
    Neff_hat=N*0.7; % threshold
    
    if (Neff(t) < Neff_hat)
        % resample state
        for i=1:state_dim
            P_State(i,:)=randsample(P_State(i,:),N,true,P_Weight_Nor);
        end
        % resample parameter
        for i=1:para_dim
            P_Para(i,:)=randsample(P_Para(i,:),N,true,P_Weight_Nor);
        end
    end
    
    
    %% Saving   
    % update covariance of parameters
    covP_Para=var(P_Para,0,2);
    
    % store particles (parameter)
    P_ParaSpace(:,:,t)=P_Para;
    
    % store mean 
    meanP_Para(:,t)=mean(P_Para,2);
    
    % store particles (state)
    P_StateSpace(:,:,t)=P_State;
    
    % store mean
    meanP_State(:,t)=mean(P_State,2);
    
    %% plot parameters every iteration
    set(h1(2),'XData',1:t,'YData',meanP_Para(1,1:t),'LineStyle','-',...
               'Color','red','LineWidth',0.5);
    set(h2(2),'XData',1:t,'YData',meanP_Para(2,1:t),'LineStyle','-',...
               'Color','red','LineWidth',0.5);
    set(h3(2),'XData',1:t,'YData',meanP_Para(3,1:t),'LineStyle','-',...
               'Color','red','LineWidth',0.5);
    set(h4(2),'XData',1:t,'YData',meanP_Para(4,1:t),'LineStyle','-',...
               'Color','red','LineWidth',0.5);
    set(h5(2),'XData',1:t,'YData',meanP_Para(5,1:t),'LineStyle','-',...
               'Color','red','LineWidth',0.5);
    set(h6(2),'XData',1:t,'YData',meanP_Para(6,1:t),'LineStyle','-',...
               'Color','red','LineWidth',0.5);
    drawnow;
end
%% error of parameters
err_Para=abs(truePara(:,1)-meanP_Para(:,T));
disp('error of parameters:')
disp(err_Para)




