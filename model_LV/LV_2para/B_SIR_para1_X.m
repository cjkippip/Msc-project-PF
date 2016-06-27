% Project Msc 
% Lotka-Volterra model
% SIR estimate 2 parameters, one guess of initial state
%% Parameter
clear, clc;

%   Global variables declaration
global N;
global a;
global covParaCoefVect;

state_dim = 2;                      % state dimension
para_dim = 2;                       % parameter dimension
ext_state_dim = state_dim+para_dim; % extended dimension
obs_dim = 2;                        % 1 and 2 state can be observed
N = 1000;                            % number of particles 


e = 0.96;
a = (3*e-1)/(2*e);
covParaCoef = 1-a^2;
covParaCoefVect = [covParaCoef;covParaCoef];

% Convariance of state
covState1 = 150;
covState2 = 50;

% Convariance of measurement
covMeasure = 10000;

% time length and sampling rate
endT = 30;
tspan = linspace(0,endT,301);
delta_t = tspan(2)-tspan(1);
T = length(tspan); % iteration time
%% Generate true state, parameter and observation
state1=10;
state2=10;
alphaStar=1;
betaStar=0.05;
deltaStar=0.02;
gammaStar=0.5;

init=[state1 state2 alphaStar betaStar];

options = odeset('RelTol', 1e-2);
[~,trueExtState] = ode45('LV_2para', tspan, init, options);
trueExtState=trueExtState';

% true state
trueState=trueExtState(1:2,:);
trueState(1,:)=trueState(1,:) + sqrt(covState1)*randn(1,T);  
trueState(2,:)=trueState(2,:) + sqrt(covState1)*randn(1,T);   

% true parameter
truePara=trueExtState(3:end,:);

% observation matrix
H=diag([1 1]);

% generate observation data
trueOb=H*trueState+sqrt(covMeasure)*randn(2,T);

% Figure one for showing Ku and Alphs_s
figure(1);
subplot(2,1,1);
h11 = plot(zeros(2,2));title('Parameter \alpha','FontSize',20);
set(h11(1),'XData',1:size(truePara,2),'YData',truePara(1,:),'LineStyle',...
           '--','color','blue','LineWidth',1);
xlabel('Time');
ylabel('Parameter \alpha');

subplot(2,1,2);
h12 = plot(zeros(2,2));title('Parameter \beta','FontSize',20);
set(h12(1),'XData',1:size(truePara,2),'YData',truePara(2,:),'LineStyle',...
           '--','color','blue','LineWidth',1);
xlabel('Time');
ylabel('Parameter \beta');

%% SMC

% particles of parameters
P_Para=ones(2,N);
P_Para(1,:)=5+sqrt(2)*randn(1,N);
P_Para(2,:)=1+sqrt(2)*randn(1,N);

% memory for particles of parameters in all iteration
P_ParaSpace=ones(2,N,T);
P_ParaSpace(:,:,1)=P_Para; % 1st iter

% memory for means of particles of parameters in all iteration
meanP_Para=ones(2,T);
meanP_Para(:,1)=mean(P_Para,2); % 1st iter

% covariance of particles of parameters 1st iter
covP_Para=var(P_Para,0,2);

% particles of state
P_State=randn(2,N);

% memory for particles of states in all iteration
P_StateSpace=ones(2,N,T);
P_StateSpace(:,:,1)=P_State; % 1st iter

% memory for means of particles of states in all iteration
meanP_State=ones(2,T);
meanP_State(:,1)=mean(P_State,2); % 1st iter

% == SMC main loop ==
for t=2:T
    % generate state by parameter for every particles
    for n=1:N
        solSmc = ode45(@(t,y) LV_2para_SMC(t,y,P_Para(:,n)),...
            [0 delta_t],P_State(:,n),options);
        y = deval(solSmc,delta_t);
        P_State(:,n) = y;
    end
%     for n=1:N
%         y1=P_State(1,n)+delta_t*...
%             (P_Para(1,n)*P_State(1,n)-P_Para(2,n)*P_State(1,n)*P_State(2,n));
%         y2=P_State(2,n)+delta_t*...
%             (0.02*P_State(1,n)*P_State(2,n)-0.5*P_State(2,n));
%         P_State(:,n) = [y1,y2];
%     end
    
    % add noise to underlying state
    P_State(1,:) = P_State(1,:) + sqrt(covState1)*randn(1,N);
    P_State(2,:) = P_State(2,:) + sqrt(covState2)*randn(1,N); 
                 
    % update parameters
    P_Para=updatePara_all_para(P_Para,covP_Para,meanP_Para(:,t-1));
    
    % weight
    P_Weight=ones(1,N);
    for n=1:N
        P_Weight(n)=exp(-0.5*(trueOb(:,t)-H*P_State(:,n))'*...
            inv(covMeasure*eye(2))*...
            (trueOb(:,t)-H*P_State(:,n)))+eps;      
    end
    
    % normalise weight
    norP_Weight=P_Weight/sum(P_Weight);
    
    %% Resampling
    % resample state
    for i=1:2
        P_State(i,:)=randsample(P_State(i,:),N,true,norP_Weight);
    end
    
    % resample parameter
    for i=1:2
        P_Para(i,:)=randsample(P_Para(i,:),N,true,norP_Weight);
    end
    
    %% Saving   
    % update covariance of parameters
    covP_Para=var(P_Para,0,2);
    
    % store particles parameter
    P_ParaSpace(:,:,t)=P_Para;
    
    % store mean 
    meanP_Para(:,t)=mean(P_Para,2);
    
    % store particles state
    P_StateSpace(:,:,t)=P_State;
    
    % store mean
    meanP_State(:,t)=mean(P_State,2);
    
    %% plot
    %   alpha
    set(h11(2),'XData',1:t,'YData',meanP_Para(1,1:t),'LineStyle','-',...
               'Color','red','LineWidth',1);
    %   beta
    set(h12(2),'XData',1:t,'YData',meanP_Para(2,1:t),'LineStyle','-',...
               'Color','red','LineWidth',1);
    drawnow;
end

%% Plot result
para_name={'\alpha','\beta'};

figure(2),clf,
for k=1:state_dim                               
    subplot(state_dim,1,k)
    plot(tspan, trueState(k,:), 'b-', tspan, meanP_State(k,:), 'r-')
    if k==1
        title('PF estimate state','FontSize',15);
    end
    xlabel('Time','FontSize',13,'FontWeight','bold'); 
    ylabel('Number of spicies','FontSize',13,'FontWeight','bold'); 
    legend({'ture','estimation'},...
        'Location','northeast','FontSize',11,'FontWeight','bold'); 
    grid on
    grid minor
end
figure(3),clf,
for k=1:para_dim                              
    subplot(para_dim,1,k)
    plot(tspan, truePara(k,:), 'b-', tspan, meanP_Para(k,:), 'r-')
    if k==1
        title('PF estimate parameter','FontSize',15);
    end
    xlabel('Time','FontSize',13,'FontWeight','bold'); 
    ylabel(['Parameter',para_name(k)],'FontSize',13,'FontWeight','bold'); 
    legend({'ture','estimation'},...
        'Location','northeast','FontSize',11,'FontWeight','bold');  
    grid on
    grid minor
end

















