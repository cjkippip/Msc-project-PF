%%%-----------------------------------------------------------------------
%
%   Description: This is the main file for applying SMC to heat shock model
%   and the initial values of state variables are 0,0 and 0
%   respectively. The time interval for integration is from zero to 200 and
%   it has been seperated into 1000 points, which means the sampling 
%   interval is 0.2 second.
%
%   For SMC part,initially the particle number is chose as 1000, this is
%   the version for all parameters unknown case. The estimations are
%   simultaneously updating when the programme is running. 
%
%   Date: 25/04/2011
%
%   Author: Xin Liu
%%%-----------------------------------------------------------------------

clear all;
clc

%   Global variables
global N;
global a;
global covParaCoefVect;

deltaT = 0.2;   % The sampling time interval
N = 500;        % number of used particles

% Covariance of states
covState1 = 10000;
covState2 = 10;
covState3 = 1250000; 

% Covariance of measurement
covMeasure = 1000000000; 

%   For measuring time cost
startTime = cputime;                     

%   Parameters for SMC algorithm
e = 0.96;
a = (3*e - 1)/(2*e);
covParaCoeff = 1 - a^2;
covParaCoefVect = [0.045;covParaCoeff;0.05;0.045;covParaCoeff;0.045];

%%  Generate the true state

init = [0 0 0];    % the initial values
options = odeset('RelTol',1e-02);   % the real error tolerance is 1e-02

%   Generate synthetic data for states
sol = ode45(@hs_ode,[0 200],init,options);
tspan = linspace(0,200,1000);
trueState = deval(sol,tspan);

% Adding noise for state
trueState(1,:)=trueState(1,:) + sqrt(covState3)*randn(1,size(trueState,2));
trueState(2,:)=trueState(2,:) + sqrt(covState2)*randn(1,size(trueState,2));  
trueState(3,:)=trueState(3,:) + sqrt(covState3)*randn(1,size(trueState,2));  

% The observation matrix, means the second state is unobservable
h = [1 0 1];    
h = diag(h);

%   Generate synthetic data for observations
trueObser = h * trueState + sqrt(covMeasure)* randn(size(h * trueState));

%%  Plot Graph

%   Aim to simultaneously show the results when the programme is going
trueParaKu = 0.0254*ones(1,length(trueState(1,:)));
trueParaAs = 3*ones(1,length(trueState(1,:)));
trueParaKs = 0.05*ones(1,length(trueState(1,:)));
trueParaAd = 0.015*ones(1,length(trueState(1,:)));
trueParaKd = 3*ones(1,length(trueState(1,:)));
trueParaA0 = 0.03*ones(1,length(trueState(1,:)));
set(gcf,'doublebuffer','on');

% Figure one for showing Ku and Alphs_s
figure(1);
subplot(2,1,1);
h21 = plot(zeros(2,2));title('Parameter K_{u}','FontSize',20);
set(h21(1),'XData',1:size(trueState,2),'YData',trueParaKu,'LineStyle',...
           '--','color','blue','LineWidth',1);
xlabel('Time');
ylabel('Parameter K_{u}');

subplot(2,1,2);
h22 = plot(zeros(2,2));title('Parameter \alpha_{s}','FontSize',20);
set(h22(1),'XData',1:size(trueState,2),'YData',trueParaAs,'LineStyle',...
           '--','color','blue','LineWidth',1);
xlabel('Time');
ylabel('Parameter \alpha_{s}');

% Figure two for showing Ks and Alphs_d
figure(2);
subplot(2,1,1);
h31 = plot(zeros(2,2));title('Parameter K_{s}','FontSize',20);
set(h31(1),'XData',1:size(trueState,2),'YData',trueParaKs,'LineStyle',...
           '--','color','blue','LineWidth',1);
xlabel('Time');
ylabel('Parameter K_{s}');

subplot(2,1,2);
h32 = plot(zeros(2,2));title('Parameter \alpha_{d}','FontSize',20);
set(h32(1),'XData',1:size(trueState,2),'YData',trueParaAd,'LineStyle',...
           '--','color','blue','LineWidth',1);
xlabel('Time');
ylabel('Parameter \alpha_{d}');

% Figure three for showing Kd and Alphs_0
figure(3);
subplot(2,1,1);
h41 = plot(zeros(2,2));title('Parameter K_{d}','FontSize',20);
set(h41(1),'XData',1:size(trueState,2),'YData',trueParaKd,'LineStyle',...
           '--','color','blue','LineWidth',1);
xlabel('Time');
ylabel('Parameter K_{d}');

subplot(2,1,2);
h42 = plot(zeros(2,2));title('Parameter \alpha_{0}','FontSize',20);
set(h42(1),'XData',1:size(trueState,2),'YData',trueParaA0,'LineStyle',...
           '--','color','blue','LineWidth',1);
xlabel('Time');
ylabel('Parameter \alpha_{0}');


%%  SMC

paraSmc = zeros(6,N);
paraSmc(1,:) = 1 + rand(1,N);               % Initialization for Ku
paraSmc(2,:) = 4 + sqrt(5.2)*randn(1,N);    % Initialization for As
paraSmc(3,:) = 3 + rand(1,N);               % Initialization for Ks
paraSmc(4,:) = 2 + rand(1,N);               % Initialization for A0
paraSmc(5,:) = 0 + sqrt(1.95)*randn(1,N);      % Initialization for Kd
paraSmc(6,:) = 2 + rand(1,N);               % Initialization for Ad

%   Obtain the variance of the first iteration
covParaSmc = var((paraSmc'))';             

%   Create memory for mean of particles for parameters
realParaSmc = zeros(6,size(trueState,2));

%   Calculate the mean of particles for parameters in the first iteration
realParaSmc(:,1) = mean((paraSmc'));

%   Create memory for saving particles for parameters in all iterations
paraSmcSpace = zeros(6,N,size(trueState,2));

%   Save particles for parameters in the first iterations
paraSmcSpace(:,:,1) = paraSmc;

%   Initilize the particle for states
concenSmc = randn(3,N);

%   Create memory for saving particles for states in all iterations
concenSpace = zeros(3,N,size(trueState,2));

%   Save particles for states in the first iterations
concenSpace(:,:,1) = concenSmc;

%   Create memory for mean of particles for staets
realSmc = zeros(3,size(trueState,2));

%   Save particles for states in the first iterations
realSmc(:,1) = mean((concenSmc'));

%   == Main loop of SMC ==
for t = 2:length(trueState)
    %   Update the states by integrating previous states respect to
    %   sampling interval.
    for n = 1:N
        solSmc = ode45(@(t,y) hs_odeSmc_All_Para(t,y,paraSmc(:,n)),...
                              [0 deltaT],concenSmc(:,n),options);
        y = deval(solSmc,deltaT);
        concenSmc(:,n) = y;
        fprintf('%d particles in %d iteration \n',n,t);
    end
    
    %   Add noise to underlying state
    concenSmc(1,:) = concenSmc(1,:) + sqrt(covState1)*...
                     randn(1,size(concenSmc,2));
    concenSmc(2,:) = concenSmc(2,:) + sqrt(covState2)*...
                     randn(1,size(concenSmc,2));    
    concenSmc(3,:) = concenSmc(3,:) + sqrt(covState3)*...
                     randn(1,size(concenSmc,2));
    
    %   Updating the parameters through artificial dynamics
    paraSmc = updatePara_all_para(paraSmc,covParaSmc,realParaSmc(:,t-1));
    
    %   Create memory for particle weights 
    smcWeight = zeros(N,1);
    
    %   Calculate the weights
    for n = 1:N
        smcWeight(n,1) = exp(-0.5* (trueObser(:,t) - h*concenSmc(:,n))'*...
                         inv(covMeasure*eye(size(trueState,1)))*...
                         (trueObser(:,t) - h*concenSmc(:,n))) + eps;
    end
    
    %   Normalized weight
    norWeights = (1/sum(smcWeight))*smcWeight;
    
    %   Resample states
    for c = 1:3
        concenSmc(c,:) = randsample(concenSmc(c,:),N,true,norWeights);
    end
    
    %   Resample parameters
    for k = 1:6
        paraSmc(k,:) = randsample(paraSmc(k,:),N,true,norWeights);
    end        
    
    %   Update the variance vector based on underlying parameters
    covParaSmc = var((paraSmc'))'; 
    
    %   Calculate the mean of underlying parameters    
    realParaSmc(:,t) = mean((paraSmc'));
    
    %   Store the values of parameters
    paraSmcSpace(:,:,t) = paraSmc;
    
    %   Calculate the mean of underlying states 
    realSmc(:,t) = mean((concenSmc'));
    
    %   Store the values of states
    concenSpace(:,:,t) = concenSmc;
    
    %   ===   simultaneously plot ====== %
    %   Ku
    set(h21(2),'XData',1:t,'YData',realParaSmc(1,1:t),'LineStyle','-',...
               'Color','red','LineWidth',1);
    %   As
    set(h22(2),'XData',1:t,'YData',realParaSmc(2,1:t),'LineStyle','-',...
               'Color','red','LineWidth',1);
    %   Ks
    set(h31(2),'XData',1:t,'YData',realParaSmc(3,1:t),'LineStyle','-',...
               'Color','red','LineWidth',1);
    %   Ad
    set(h32(2),'XData',1:t,'YData',realParaSmc(4,1:t),'LineStyle','-',...
               'Color','red','LineWidth',1);
    %   Kd
    set(h41(2),'XData',1:t,'YData',realParaSmc(5,1:t),'LineStyle','-',...
               'Color','red','LineWidth',1);
    %   A0
    set(h42(2),'XData',1:t,'YData',realParaSmc(6,1:t),'LineStyle','-',...
               'Color','red','LineWidth',1);
    drawnow;
    
    t
end

%   calculate the cpu time
finishTime = cputime - startTime 

%%  Plot figures
figure,
plot(trueParaAd,'b--','LineWidth',2);hold on;
plot(realParaSmc(4,:),'r--','LineWidth',2);
ylabel('alpha_{d}','Fontsize',20);
xlabel('Time','Fontsize',20);
legend('True','PF');
set(gca,'Fontsize',20);

figure,
plot(trueParaKs,'b--','LineWidth',2);hold on;
plot(realParaSmc(3,:),'r--','LineWidth',2);
ylabel('K_{s}','Fontsize',20);
xlabel('Time','Fontsize',20);
legend('True','PF');
set(gca,'Fontsize',20);

figure,
plot(trueParaKu,'b--','LineWidth',2);hold on;
plot(realParaSmc(1,:),'r--','LineWidth',2);
ylabel('K_{u}','Fontsize',20);
xlabel('Time','Fontsize',20);
legend('True','PF');
set(gca,'Fontsize',20);

figure,plot(trueParaAs,'b--','LineWidth',2);hold on;
plot(realParaSmc(2,:),'r--','LineWidth',2);
ylabel('alpha_{s}','Fontsize',20);
xlabel('Time','Fontsize',20);
legend('True','PF');
set(gca,'Fontsize',20);

figure,plot(trueParaA0,'b--','LineWidth',2);hold on;
plot(realParaSmc(6,:),'r--','LineWidth',2);
ylabel('alpha_{0}','Fontsize',20);
xlabel('Time','Fontsize',20);
legend('True','PF');
set(gca,'Fontsize',20);

figure,plot(trueParaKd,'b--','LineWidth',2);hold on;
plot(realParaSmc(5,:),'r--','LineWidth',2);
ylabel('K_{d}','Fontsize',20);
xlabel('Time','Fontsize',20);
legend('True','PF');
set(gca,'Fontsize',20);


