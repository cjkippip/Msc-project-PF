% Project Msc 
% Lotka-Volterra model
clear, clc;

state_dim=2; % state dimension
para_dim=2; % parameter dimension
ext_state_dim=state_dim+para_dim; % extended dimension
obs_dim=2; % 1 and 2 state can be observed
%% generate data by ODEs solver
state1=10;
state2=10;
alphaStar=1;
betaStar=0.05;
deltaStar=0.02;
gammaStar=0.5;

% state1=10;
% state2=10;
% alphaStar=2;
% betaStar=0.6;
% deltaStar=0.4;
% gammaStar=2.2;

init=[state1 state2 alphaStar betaStar];

endT=30;
t=linspace(0,endT,300);
delta_t=t(2)-t(1);
T=length(t);

options = odeset('RelTol', 1e-3);
[~,trueState] = ode45('LV_2para', t, init, options);

figure(1),clf,
plot(t,trueState(:,1),'b',t,trueState(:,2),'r','LineWidth',1.5);
xlabel('Time','FontSize',13,'FontWeight','bold'); 
ylabel('Number of spicies','FontSize',13,'FontWeight','bold'); 
legend({'prey', 'predators'},...
'Location','northeast','FontSize',11,'FontWeight','bold');
grid on
grid minor
%%
figure(2),clf,
plot3(trueState(:,1),trueState(:,2),t,'LineWidth',1.5);
xlabel('Prey','FontSize',13,'FontWeight','bold'); 
ylabel('Predators','FontSize',13,'FontWeight','bold'); 
zlabel('Time','FontSize',13,'FontWeight','bold'); 
grid on
grid minor
%%
T=length(t);
numGrid=50;
alpha=linspace(alphaStar*3/5,alphaStar*7/5,numGrid);
beta=linspace(betaStar*3/5,betaStar*7/5,numGrid);

[meshA,meshB] = meshgrid(alpha,beta);
change=ones(numGrid,numGrid);

sigma1=std(trueState(:,1));
sigma2=std(trueState(:,2));
allEstiState=cell(numGrid,numGrid);

%%
for i=1:numGrid
    for j =1:numGrid
        init2=[state1 state2 alpha(j) beta(i)];
        options = odeset('RelTol', 1e-3);
        [~,estiState] = ode45('LV_2para',t,init2,options);
        allEstiState{i,j}=estiState;     
        risd1=sumsqr(estiState(:,1)-trueState(:,1))/sigma1^2;
        risd2=sumsqr(estiState(:,2)-trueState(:,2))/sigma2^2;
        change(i,j)=(risd1+risd2)/(T*2*state_dim);
    end
end
%% draw contour
figure(3),clf,
% contour(log(meshA),log(meshB),change);
contour(meshA,meshB,change);
hold on
% plot(log(alphaStar),log(betaStar),'o','MarkerSize',9,...
%     'MarkerEdgeColor','k','MarkerFaceColor','r');
plot(alphaStar,betaStar,'o','MarkerSize',9,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');
xlabel('alpha','FontSize',13,'FontWeight','bold'); 
ylabel('beta','FontSize',13,'FontWeight','bold'); 
grid on
grid minor
%% 3D
figure(4),clf,
% surf(log(meshA),log(meshB),change);
surf(meshA,meshB,change);
hold on
% plot3(log(alphaStar),log(betaStar),0,'o','MarkerSize',9,...
%     'MarkerEdgeColor','k','MarkerFaceColor','r');
plot3(alphaStar,betaStar,0,'o','MarkerSize',9,...
    'MarkerEdgeColor','k','MarkerFaceColor','r');
xlabel('alpha','FontSize',13,'FontWeight','bold'); 
ylabel('beta','FontSize',13,'FontWeight','bold'); 
zlabel('Behavior Change','FontSize',13,'FontWeight','bold'); 
grid on
grid minor










