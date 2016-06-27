% Project Msc 
% Lotka-Volterra model
% SIR estimate 4 parameters, one guess of initial state
%% Parameter
clear, clc;

% Global variables declaration
global N;
global a;
global covParaCoefVect;

state_dim = 2;                      % state dimension
para_dim = 4;                       % parameter dimension
ext_state_dim = state_dim+para_dim; % extended dimension
obs_dim = 2;                        % 1 and 2 state can be observed
N = 500;                           % number of particles {{***AD***}}

e = 0.96; % {{***AD***}}
a = (3*e-1)/(2*e);
covParaCoef = 1-a^2;
% covParaCoefVect = [covParaCoef;covParaCoef;covParaCoef;covParaCoef];
covParaCoefVect = [0.05;0.05;0.05;0.05];

% Convariance of state {{***AD***}}
covState1 = 20;
covState2 = 200;

% Convariance of measurement {{***AD***}}
covMeasure = 1000;

% time length and sampling rate {{***AD***}}
endT = 200;
N_step = 1001;
tspan = linspace(0,endT,N_step);
delta_t = tspan(2)-tspan(1);
T = length(tspan); % iteration time (equal to N_step)
%% Generate true state, parameter and observation
state1=10;
state2=10;
alphaStar=1;
betaStar=0.05;
deltaStar=0.02;
gammaStar=0.5;

init=[state1 state2 alphaStar betaStar deltaStar gammaStar];

options = odeset('RelTol', 1e-2);
sol = ode45(@LV_4para, [0 endT], init, options);
trueExtState = deval(sol,tspan);

% true state
trueState=trueExtState(1:2,:);
trueState(1,:)=trueState(1,:) + sqrt(covState1)*randn(1,T);  
trueState(2,:)=trueState(2,:) + sqrt(covState2)*randn(1,T);   

% true parameter
truePara=trueExtState(3:end,:);

% observation matrix
H=diag([1 1]);

% generate observation data
trueOb=H*trueState+sqrt(covMeasure)*randn(2,T);
%% Figure to show parameter
% Figure 1 for alpha and beta
figure(1);
subplot(2,1,1);
h1 = plot(zeros(2,2));title('Parameter \alpha','FontSize',20);
set(h1(1),'XData',1:size(truePara,2),'YData',truePara(1,:),'LineStyle',...
           '--','color','blue','LineWidth',0.75);
xlabel('Time');
ylabel('Parameter \alpha');

subplot(2,1,2);
h2 = plot(zeros(2,2));title('Parameter \beta','FontSize',20);
set(h2(1),'XData',1:size(truePara,2),'YData',truePara(2,:),'LineStyle',...
           '--','color','blue','LineWidth',0.75);
xlabel('Time');
ylabel('Parameter \beta');

% Figure 2 for delta and gamma
figure(2);
subplot(2,1,1);
h3 = plot(zeros(2,2));title('Parameter \delta','FontSize',20);
set(h3(1),'XData',1:size(truePara,2),'YData',truePara(3,:),'LineStyle',...
           '--','color','blue','LineWidth',0.75);
xlabel('Time');
ylabel('Parameter \delta');

subplot(2,1,2);
h4 = plot(zeros(2,2));title('Parameter \gamma','FontSize',20);
set(h4(1),'XData',1:size(truePara,2),'YData',truePara(4,:),'LineStyle',...
           '--','color','blue','LineWidth',0.75);
xlabel('Time');
ylabel('Parameter \gamma');
%% SMC
%=======================================================
% initial particles of parameters
P_Para=ones(para_dim,N);
% gaussian {{***AD***}}
P_Para(1,:)=1+rand(1,N);
P_Para(2,:)=2+sqrt(1.5)*randn(1,N);
P_Para(3,:)=3+rand(1,N);
P_Para(4,:)=4+sqrt(1.2)*randn(1,N);
% uniform {{***AD***}}
% P_Para(1,:)=unifrnd(-10,10,[1,N]);
% P_Para(2,:)=unifrnd(-100,100,[1,N]);
% P_Para(3,:)=unifrnd(-100,100,[1,N]);
% P_Para(4,:)=unifrnd(-10,10,[1,N]);
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
% initial particles of state {{***AD***}}
P_State=randn(state_dim,N);
%=======================================================

% memory for particles(states) in all iteration
P_StateSpace=ones(state_dim,N,T);
P_StateSpace(:,:,1)=P_State; % 1st iter

% memory for means of particles(states) in all iteration
meanP_State=ones(state_dim,T);
meanP_State(:,1)=mean(P_State,2); % 1st iter

% memory for weights and norweights of particles in all iteration
P_WeightSpace = ones(1,N,N_step);
norP_WeightSpace = ones(1,N,N_step);

% weight
P_Weight=ones(1,N)/N;
norP_Weight=ones(1,N)/N;

% == SMC main loop ==
for t=2:T
    % generate state by parameter for every particles
    for n=1:N
        solSmc = ode45(@(t,y) LV_4para_SMC(t,y,P_Para(:,n)),...
            [0 delta_t],P_State(:,n),options);
        
        y = deval(solSmc,delta_t);
        P_State(:,n) = y;
    end

    % faster one than above
%     for n=1:N
%         y1=P_State(1,n)+delta_t*...
%             (P_Para(1,n)*P_State(1,n) - P_Para(2,n)*P_State(1,n)*P_State(2,n));
%         y2=P_State(2,n)+delta_t*...
%             (P_Para(3,n)*P_State(1,n)*P_State(2,n) - P_Para(4,n)*P_State(2,n));
% %         if ~((y1>-1e12)&&(y1<1e12)&&(y2>-1e12)&&(y2<1e12))
% %             y1=P_State(1,n) + sqrt(covState1)*randn(1);
% %             y2=P_State(2,n) + sqrt(covState1)*randn(1);
% %         end
%         P_State(:,n) = [y1;y2];
%     end
    
    % add noise to underlying state
    P_State(1,:) = P_State(1,:) + sqrt(covState1)*randn(1,N);
    P_State(2,:) = P_State(2,:) + sqrt(covState2)*randn(1,N); 
                 
    % update parameters (perturbation)
    P_Para=updatePara_all_para(P_Para,covP_Para,meanP_Para(:,t-1));
    
    % weight, likelihood   
    for n=1:N
        P_Weight(n)=exp( -0.5*(trueOb(:,t)-H*P_State(:,n))' * ...
            inv(covMeasure*eye(obs_dim)) * ...
            (trueOb(:,t)-H*P_State(:,n)) ) + eps;      
    end
    P_WeightSpace(:,:,t)=P_Weight;
    
    % normalise weight
    norP_Weight=P_Weight/sum(P_Weight);
    norP_WeightSpace(:,:,t)=norP_Weight;
    %% Resampling
    % resample state
    for i=1:state_dim
        P_State(i,:)=randsample(P_State(i,:),N,true,norP_Weight);
    end
    
    % resample parameter
    for i=1:para_dim
        P_Para(i,:)=randsample(P_Para(i,:),N,true,norP_Weight);
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
    %   alpha
    set(h1(2),'XData',1:t,'YData',meanP_Para(1,1:t),'LineStyle','-',...
               'Color','red','LineWidth',0.5);
    %   beta
    set(h2(2),'XData',1:t,'YData',meanP_Para(2,1:t),'LineStyle','-',...
               'Color','red','LineWidth',0.5);
    %   delta
    set(h3(2),'XData',1:t,'YData',meanP_Para(3,1:t),'LineStyle','-',...
               'Color','red','LineWidth',0.5);
    %   gamma
    set(h4(2),'XData',1:t,'YData',meanP_Para(4,1:t),'LineStyle','-',...
               'Color','red','LineWidth',0.5);
    drawnow;
end
%% error of parameters
err_Para=abs(truePara(:,1)-...
    [meanP_Para(1,N_step);...
    meanP_Para(2,N_step);...
    meanP_Para(3,N_step);...
    meanP_Para(4,N_step)]);
disp('error of parameters:')
disp(err_Para)
%% Plot result
para_name={'\alpha','\beta','\delta','\gamma'};

% state
% figure(3),clf,
% for k=1:state_dim                               
%     subplot(state_dim,1,k)
%     plot(tspan, trueState(k,:), 'b-', tspan, meanP_State(k,:), 'r-')
%     if k==1
%         title('PF estimate state','FontSize',15);
%     end
%     xlabel('Time','FontSize',13,'FontWeight','bold'); 
%     ylabel('Number of spicies','FontSize',13,'FontWeight','bold'); 
%     legend({'ture','estimation'},...
%         'Location','eastoutside','FontSize',11,'FontWeight','bold'); 
%     grid on
%     grid minor
% end

% parrameter
figure(4),clf,
for k=1:para_dim                              
    subplot(para_dim,1,k)
    plot(tspan, truePara(k,:), 'b--', tspan, meanP_Para(k,:), 'r-')
    if k==1
        title('PF estimate parameter','FontSize',15);
    end
    xlabel('Time','FontSize',13,'FontWeight','bold'); 
    ylabel(['Parameter',para_name(k)],'FontSize',13,'FontWeight','bold'); 
    legend({'ture','estimation'},...
        'Location','eastoutside','FontSize',11,'FontWeight','bold');  
    grid on
    grid minor
end
%% plot distribution of states and parameters
% state
% figure(5),clf,
% subplot(211)
% hist(P_StateSpace(1,:,N_step),50);
% hold on
% plot(meanP_State(1,N_step),0,'r*');
% 
% subplot(212)
% hist(P_StateSpace(2,:,N_step),50);
% hold on
% plot(meanP_State(2,N_step),0,'r*');

% parameter
figure(6),clf,
subplot(221)
hist(P_ParaSpace(1,:,N_step),50);
hold on
plot(meanP_Para(1,N_step),0,'r*');
plot(alphaStar,0,'ro');
title('true \alpha = 1')

subplot(222)
hist(P_ParaSpace(2,:,N_step),50);
hold on
plot(meanP_Para(2,N_step),0,'r*');
plot(betaStar,0,'ro');
title('true \beta = 0.05')

subplot(223)
hist(P_ParaSpace(3,:,N_step),50);
hold on
plot(meanP_Para(3,N_step),0,'r*');
plot(deltaStar,0,'ro');
title('true \delta = 0.02')

subplot(224)
hist(P_ParaSpace(4,:,N_step),50);
hold on
plot(meanP_Para(4,N_step),0,'r*');
plot(gammaStar,0,'ro');
title('true \gamma = 0.5')
%%










