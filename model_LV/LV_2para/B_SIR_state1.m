% Project Msc 
% Lotka-Volterra model
% SIR estimate 2 states
%%
clear, clc;

state_dim=2; % state dimention
para_dim=2; % parameter dimention
ext_state_dim=state_dim+para_dim; % extend state dimention
delta_t=0.1; % delta time
obs_dim=2; % observatin dimention
N_P = 200; % number of particles 
%%
init=[10 10]; % LV model initial state
endT=50;

options = odeset('RelTol', 1e-4, 'NonNegative', [1 2],...
    'Refine',1,'MaxStep ',0.1);
[t,trueState] = ode45('LV', [0 endT], init, options);
% plot LV model
figure(1),clf,
plot(t,trueState(:,1),'b',t,trueState(:,2),'r','LineWidth',1.5);
xlabel('Time','FontSize',13,'FontWeight','bold'); 
ylabel('Number of spicies','FontSize',13,'FontWeight','bold'); 
legend({'prey', 'predators'},...
'Location','northeast','FontSize',11,'FontWeight','bold');
grid on
grid minor
%%
f=@(x)[x(1)+delta_t*(1*x(1)-0.05*x(1)*x(2));...
        x(2)+delta_t*(0.02*x(1)*x(2)-0.5*x(2))]; % state equations
h=@(x)[x(1);x(2)];  % measurement equation

T=length(t); 
q=1; %std of process 
r=3; %std of measurement

Q=diag([q^2 q^2]); % covariance of process 
R=diag([r^2 r^2]); % covariance of measurement 

aV = trueState'; % actuality
aV(1,:) = aV(1,:)+q*randn(1,T);
aV(2,:) = aV(2,:)+q*randn(1,T);
%%
xV = zeros(2,T); % estmation
% x_init=[h(aV(:,1)) + r*randn(2,1);2;0.1];
x_init=[30;30];
x = x_init; % initial estmation 
 
zV = zeros(2,T); % measurments
z = h(aV(:,1)) + r*randn(2,1); % measurments
zV(:,1) = z; % measurments

P=ones(2,N_P);
w = zeros(1,N_P); 
for i = 1 : N_P
    P(:, i) = [10+5*rand; 10+5*rand];
    dist = norm(P(:, i)-zV(:, 1));     
    w(i) = (1 / sqrt(r) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / r); 
end

xV(:,1) = mean(P,2); % save estimate

%%

for k = 2 : T
    zV(:,k) = h(aV(:,k)) + r*randn(2,1);    
    for i = 1 : N_P
        P(:, i) = f(P(:, i)) + q*randn(2,1);
        dist = norm(P(:, i)-zV(:,k));     
        w(i) = (1 / sqrt(r) / sqrt(2 * pi)) * exp(-(dist)^2 / 2 / r); 
    end
    wsum = sum(w);
    w = w ./ wsum;
    
%     % resamling
%     for i = 1 : N_P  
%         P(:,i) = P(:,find(rand <= cumsum(w),1));      
%     end 
    
    for i = 1 : N_P
        wmax = 2 * max(w) * rand; 
        index = randi(N_P, 1);
        while(wmax > w(index))
            wmax = wmax - w(index);
            index = index + 1;
            if index > N_P
                index = 1;
            end          
        end
        P(:, i) = P(:, index);     
    end
    
    xV(:, k) = sum(P, 2) / N_P;    
end


%%
figure(2),clf,
for k=1:2                               
    subplot(2,1,k)
    plot(t, aV(k,:), 'b-', t, xV(k,:), 'r-')
    if k==1
        title('PF estimate states','FontSize',15);
    end
    xlabel('Time','FontSize',13,'FontWeight','bold'); 
    ylabel('Number of spicies','FontSize',13,'FontWeight','bold'); 
    legend({'reality','estimation'},...
        'Location','northeast','FontSize',11,'FontWeight','bold'); 
    grid on
    grid minor
end







