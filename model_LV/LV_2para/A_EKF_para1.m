% Project Msc 
% Lotka-Volterra model
% EKF estimate 2 parameters, one guess of initial state
%%
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

endT=50;
t=linspace(0,endT,500);
delta_t=t(2)-t(1);
T=length(t);

options = odeset('RelTol', 1e-4);
[~,trueState] = ode45('LV_2para', t, init, options);

figure(1),clf,
plot(t,trueState(:,1),'b',t,trueState(:,2),'r','LineWidth',1.5);
xlabel('Time','FontSize',13,'FontWeight','bold'); 
ylabel('Number of spicies','FontSize',13,'FontWeight','bold'); 
legend({'prey', 'predators'},...
'Location','northeast','FontSize',11,'FontWeight','bold');
grid on
grid minor
%% EKF estimate state and parameter alpha and beta

% covariance of process 
q1=1e-6; q2=1e-6; q3=1e-6; q4=1e-6;
Q=diag([q1 q2 q3 q4]); 

% true states and parameters
aV = trueState'; 
aV(1,:) = aV(1,:)+sqrt(q1)*randn(1,T);
aV(2,:) = aV(2,:)+sqrt(q2)*randn(1,T);
% aV(3,:) = aV(3,:)+q3*randn(1,T);
% aV(4,:) = aV(4,:)+q4*randn(1,T);

% covariance of measurement 
varX1=var(aV(1,:));
varX2=var(aV(2,:));
% r1=0.25*varX1;
% r2=0.25*varX2;
r1=3;
r2=3;
R=diag([r1 r2]); 

% state equations
f=@(x)[x(1)+delta_t*(x(1)*x(3)-x(4)*x(1)*x(2));...
        x(2)+delta_t*(deltaStar*x(1)*x(2)-gammaStar*x(2));...
        x(3);...
        x(4)]; 
% measurement equation
h=@(x)[x(1);x(2)]; 

% estmation
xV = zeros(ext_state_dim,T); 
%*******change initial estimation**************
x = [15;5;3;1]; % initial estmation 
%**********************************************
xV(:,1) = x; % save estimate

% measurments
zV=zeros(para_dim,T); 
z=h(aV(:,1)) + [sqrt(r1)*randn;sqrt(r2)*randn]; % measurments
zV(:,1)=z; % measurments

% initial covariance matrix
P=diag([25 25 4 1]);
PM=ones(ext_state_dim,ext_state_dim,T);
PM(:,:,1)=P;

for k=2:T
    z = h(aV(:,k)) + [sqrt(r1)*randn;sqrt(r2)*randn]; % measurments
    zV(:,k) = z; % save measurment
    [x, P] = ekf(f,x,P,h,z,Q,R); % ekf 
    PM(:,:,k) = P;
    xV(:,k) = x; % save estimate
end
%%
figure(2),clf,
for k=1:2      
    subplot(2,1,k)    
    plot(t, aV(k,:), 'b-', t, xV(k,:), 'r-')
    if k==1
        title('EKF estimate states','FontSize',15);
    end
    xlabel('Time','FontSize',13,'FontWeight','bold'); 
    ylabel('Number of spicies','FontSize',13,'FontWeight','bold'); 
    legend({'ture','estimation'},...
        'Location','northeast','FontSize',11,'FontWeight','bold'); 
    grid on
    grid minor
end
figure(3),clf,
for k=1:2                              
    subplot(2,1,k)
    plot(t, aV(k+2,:), 'b-', t, xV(k+2,:), 'r-')
    if k==1
        title('EKF estimate parameter','FontSize',15);
    end
    xlabel('Time','FontSize',13,'FontWeight','bold'); 
    ylabel(['Parameter',num2str(k)],'FontSize',13,'FontWeight','bold'); 
    legend({'ture','estimation'},...
        'Location','northeast','FontSize',11,'FontWeight','bold');  
    grid on
    grid minor
end
%%







