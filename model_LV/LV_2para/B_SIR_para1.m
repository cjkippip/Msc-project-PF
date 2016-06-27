% Project Msc 
% Lotka-Volterra model
% SIR estimate 2 parameters, one guess of initial state
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
t=linspace(0,endT,501);
delta_t=t(2)-t(1);
T=length(t);

options = odeset('RelTol', 1e-3);
[~,trueExtState] = ode45('LV_2para', t, init, options);

figure(1),clf,
plot(t,trueExtState(:,1),'b',t,trueExtState(:,2),'r','LineWidth',1.5);
xlabel('Time','FontSize',13,'FontWeight','bold'); 
ylabel('Number of spicies','FontSize',13,'FontWeight','bold'); 
legend({'prey', 'predators'},...
'Location','northeast','FontSize',11,'FontWeight','bold');
grid on
grid minor
%% SIR Filter estimate state and parameter alpha and beta
% covariance of process 
q1=1e-6; q2=1e-6; q3=1e-6; q4=1e-6;
Q=diag([q1 q2 q3 q4]); 

% true states and parameters
aV = trueExtState'; 
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
    
% state equations
f=@(x)[x(1)+delta_t*(x(1)*x(3)-x(4)*x(1)*x(2));...
        x(2)+delta_t*(deltaStar*x(1)*x(2)-gammaStar*x(2));...
        x(3);...
        x(4)]; 
% measurement equation
h=@(x)[x(1);x(2)]; 

% estimate
xV = zeros(4,T);

% measurments
zV=zeros(para_dim,T); 
z=h(aV(:,1)) + [sqrt(r1)*randn;sqrt(r2)*randn]; % measurments
zV(:,1)=z; % measurments vector

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
for k = 2 : T
    zV(:,k) = h(aV(:,k)) + [sqrt(r1)*randn;sqrt(r2)*randn];    
    for i = 1 : N_P
        P(:, i)=f(P(:, i))+[sqrt(q1)*randn;sqrt(q2)*randn;...
            sqrt(q3)*randn;sqrt(q4)*randn];
%         P(:, i)=f(P(:, i));

        z_update=h(P(:, i))+[sqrt(r1)*randn;sqrt(r2)*randn];  
%         z_update=h(P(:, i)); 
        
        dist = norm(z_update-zV(:,k));     
        PWeight(i) = (1/sqrt(2*pi*r1))*exp(-(dist)^2/(2 * r1)); 
    end
    PWeight = PWeight./sum(PWeight);
    
%     aaCum=cumsum(w);
%     aaindex=find(rand <= aaCum,1);
%     disp('breakpoint');
    
%     resamling
%     for i = 1 : N_P  
%         P(:,i) = P(:,find(rand <= cumsum(w),1));      
%     end 
    Neff(k) = 1/sum(PWeight.^2);
    Nt = N_P*0.5;
    if Neff(k)<Nt      
        edges = min([0 cumsum(PWeight)],1); % protect against accumulated round-off
        edges(end) = 1;                 % get the upper edge exact
        u1 = rand/N_P;
        [~, idx] = histc(u1:1/N_P:1, edges);
        P = P(:,idx);                    % extract new particles
        PWeight = repmat(1/N_P, 1, N_P); 
        
%         for i = 1 : N_P
%             wmax = 2 * max(PWeight) * rand; 
%             index = randi(N_P, 1);
%             while(wmax > PWeight(index))
%                 wmax = wmax - PWeight(index);
%                 index = index + 1;
%                 if index > N_P
%                     index = 1;
%                 end          
%             end
%             P(:, i) = P(:, index);     
%         end

    end
    
%     for i = 1:N_P;
%         xV(:, k) = xV(:, k) + PWeight(i)*P(:,i);
%     end 
    
    xV(:, k) = mean(P,2);
    
end
%%
figure(2),clf,
for k=1:state_dim                               
    subplot(state_dim,1,k)
    plot(t, aV(k,:), 'b-', t, xV(k,:), 'r-')
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
    plot(t, aV(k+2,:), 'b-', t, xV(k+2,:), 'r-')
    if k==1
        title('PF estimate parameter','FontSize',15);
    end
    xlabel('Time','FontSize',13,'FontWeight','bold'); 
    ylabel(['Parameter',num2str(k)],'FontSize',13,'FontWeight','bold'); 
    legend({'ture','estimation'},...
        'Location','northeast','FontSize',11,'FontWeight','bold');  
    grid on
    grid minor
end
%%





