% Project Msc 
% Lotka-Volterra model
% RPF estimate 2 parameters, one guess of initial state
%%
clear, clc;

state_dim=2; % state dimension
para_dim=2; % parameter dimension
ext_state_dim=state_dim+para_dim; % extended dimension
obs_dim=2; % 1 and 2 state can be observed
N_P = 200; % number of particles 
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
%% RPF Filter estimate state and parameter alpha and beta
% covariance of process 
q1=1e-6; q2=1e-6; q3=1e-6; q4=1e-6;
Q=diag([q1 q2 q3 q4]); 

% true states and parameters
aV = trueState'; 
aV(1,:) = aV(1,:)+sqrt(q1)*randn(1,T);
aV(2,:) = aV(2,:)+sqrt(q2)*randn(1,T);

% covariance of measurement 
varX1=var(aV(1,:));
varX2=var(aV(2,:));
% r1=0.25*varX1;
% r2=0.25*varX2;
r1=3;
r2=3;
R=diag([r1 r2]); 

% measurments
zV=zeros(para_dim,T); 
for i=1:T
    z=h(aV(:,i)) + [sqrt(r1)*randn;sqrt(r2)*randn]; % measurments
    zV(:,i)=z; % measurments vector
end

% state equations
f=@(x)[x(1)+delta_t*(x(1)*x(3)-x(4)*x(1)*x(2));...
        x(2)+delta_t*(deltaStar*x(1)*x(2)-gammaStar*x(2));...
        x(3);...
        x(4)]; 
% measurement equation
h=@(x)[x(1);x(2)]; 

% estimate
xV = zeros(4,T);

P=ones(4,N_P);
PWeight = zeros(1,N_P); 
for i = 1 : N_P
    P(:, i) = [10+10*rand; 10+10*rand; 1+1*rand; 0.05+0.02*rand];
    z_update=h(P(:, i))+[sqrt(r1)*randn;sqrt(r2)*randn];  
    dist = norm(z_update-zV(:, 1));     
    PWeight(i) = (1/sqrt(2*pi*r1))*exp(-(dist)^2/(2 * r1)); 
end
PWeight = PWeight./sum(PWeight);

for i=1:N_P
    xV(:,1) = xV(:,1)+ P(:,i)*PWeight(i); % save estimate
end
%%


















