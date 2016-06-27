%% LV
% % Set our preferences for ode45
% % The default relative tolerance is 1e-3.
% % To set our output to non-negative, we provide an array listing
% %   each population that we want to constrain.  Since this example
% %   has two populations, we pass the array [1 2]
% options = odeset('RelTol', 1e-4, 'NonNegative', [1 2]);
% % Use ode45 to solve our ODE
% % Place the time points in a vector 't'
% % Place the solution in a vector 'x'
% [t,x] = ode45('lotka_volterra', [0 20], [10 10], options);
% plot(t,x);
% legend('prey', 'predators');
%% SIR
% options = odeset('RelTol', 1e-4, 'NonNegative', [1 2 3]);
% % Use ode45 to solve our ODE
% % Place the time points in a vector 't'
% % Place the solution in a vector 'x'
% [t,x] = ode45('sir', [0 10], [1000 1 0], options);
% plot(t,x);
% legend('S', 'I', 'R');
%% rigid
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
% [T,Y] = ode45(@rigid,[0 12],[0 1 1],options);
% plot(T,Y(:,1),'-',T,Y(:,2),'-.',T,Y(:,3),'.')
%% sin(x) and diff 1,2
% x=0:0.01:(4*pi);
% y=sin(x);
% yx=diff(y)*100;
% yxx=diff(yx)*100;
% 
% figure(1),clf,
% plot(x,y,'r');
% hold on
% plot(x(2:end),yx,'b');
% plot(x(3:end),yxx,'g');
% legend('sin(x)','yx','yxx');
% grid on
% grid minor
% hold off
%% y=1/(1+exp(-x))
% x=-5:0.01:5;
% y=1./(1+exp(-x));
% yx=diff(y)*100;
% yxx=diff(yx)*100;
% 
% figure(1),clf,
% plot(x,y,'r');
% hold on
% plot(x(2:end),yx,'b');
% plot(x(3:end),yxx,'g');
% legend('y','yx','yxx');
% grid on
% grid minor
% hold off
%% y=x^3
% x=-2:0.01:2;
% y=x.^3;
% yx=diff(y)*100;
% yxx=diff(yx)*100;
% 
% figure(1),clf,
% plot(x,y,'r');
% hold on
% plot(x(2:end),yx,'b');
% plot(x(3:end),yxx,'g');
% legend('y','yx','yxx');
% grid on
% grid minor
% hold off
%% normal distribution
% m=0;
% C=1;
% x=-5:0.01:5;
% y=1/(sqrt(2*pi)*C)*exp(-1/2*(x-m).^2/C^2);
% yx=diff(y)*100;
% yxx=diff(yx)*100;
% 
% figure(1),clf,
% plot(x,y,'r');
% hold on
% plot(x(2:end),yx,'b');
% plot(x(3:end),yxx,'g');
% legend('y','yx','yxx');
% grid on
% grid minor
% hold off
%% cholesky decomposition
% x=[2 -2;-2 5];
% [R,p] = chol(x);

% mu=[0;0];
% sigma=[1 0;0 1];
% aa=[3 2;2 3];
% aa1=sqrt(aa);
% aa2=chol(aa);
% bb1=mvnrnd(mu,sigma,1000)*aa1;
% bb2=mvnrnd(mu,sigma,1000)*aa2;
% cc=mvnrnd(mu,aa,1000);
% 
% figure(6),clf,
% scatter(bb1(:,1),bb1(:,2));
% figure(7),clf,
% scatter(bb2(:,1),bb2(:,2));
% figure(8),clf,
% scatter(cc(:,1),cc(:,2));
%% Poisson Distribution
% x=0:10;
% lambda=2;
% y=exp(-lambda).*lambda.^x./factorial(x);
% z=cumsum(y);
%% EKF example
% n=2; %number of state
% delta_t=0.1;
% t=0:delta_t:5;
% N=length(t); 
% g=5;
% 
% q=0.9; %std of process 
% r=3; %std of measurement
% Q=q^2*eye(n); % covariance of process 
% R=r^2*eye(n); % covariance of measurement 
% 
% f=@(x)[x(1)+delta_t*x(2)+delta_t^2*g/2;x(2)+delta_t*g];
% h=@(x)[x(1);x(2)];  % measurement equation
% 
% s=[0;0]; % initial state actual
% x=s+q*randn(n,1); % initial state with noise
% P = zeros(n); % initial state covraiance
% 
% xV = zeros(n,N); % estmation
% sV = zeros(n,N); % actuality
% zV = zeros(n,N); % measurments
% for k=1:N
%   z = h(s) + r*randn(n,1); % measurments
%   sV(:,k)= s; % save actual state
%   zV(:,k)  = z; % save measurment
%   [x, P] = ekf(f,x,P,h,z,Q,R); % ekf 
%   xV(:,k) = x; % save estimate
%   s = f(s) + q*randn(n,1); % update process 
% end
% figure(2),clf,
% for k=1:n                               
%     subplot(n,1,k)
%     plot(t, sV(k,:), 'b-', t, xV(k,:), 'r--')
%     legend({'reality','estimation'},...
%     'Location','northwest','FontSize',11,'FontWeight','bold'); 
%     grid on  
%     grid minor
% end
%% UKF example
% n=3;      %number of state
% q=0.1;    %std of process 
% r=0.1;    %std of measurement
% Q=q^2*eye(n); % covariance of process
% R=r^2;        % covariance of measurement  
% f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
% h=@(x)x(1);                               % measurement equation
% s=[0;0;1];                                % initial state
% x=s+q*randn(3,1); %initial state          % initial state with noise
% P = eye(n);                               % initial state covraiance
% N=20;                                     % total dynamic steps
% xV = zeros(n,N);          %estmate        % allocate memory
% sV = zeros(n,N);          %actual
% zV = zeros(1,N);
% for k=1:N
%   z = h(s) + r*randn;                     % measurments
%   sV(:,k)= s;                             % save actual state
%   zV(k)  = z;                             % save measurment
%   [x, P] = ukf(f,x,P,h,z,Q,R);            % ekf 
%   xV(:,k) = x;                            % save estimate
%   s = f(s) + q*randn(3,1);                % update process 
% end
% for k=1:3                                 % plot results
%   subplot(3,1,k)
%   plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
% end
%%
% n = 10;
% aa1 = randsample(n,5,true);
% aa2 = randsample(6:n,5,true);
%%
% aa1=[2 1];
% aa2=size(aa1,1);
%%
% init=[10 10];
% endT=50;
% 
% options = odeset('RelTol', 1e-4);
% sol = ode45(@LV, [0 endT], init, options);
% t = linspace(0,50,1000);
% trueState = deval(sol,t);
% 
% figure(1),clf,
% plot(t,trueState(1,:),'b',t,trueState(2,:),'r','LineWidth',1.5);
% xlabel('Time','FontSize',13,'FontWeight','bold'); 
% ylabel('Number of spicies','FontSize',13,'FontWeight','bold'); 
% legend({'prey', 'predators'},...
% 'Location','northeast','FontSize',11,'FontWeight','bold');
% grid on
% grid minor
%% cell {}
% aa=cell(3,1);
% aa{1}=1;
% aa{2}='aeg';
% aa{3}=1.35;
% 
% bb1=aa(1);
% bb2=aa(2);
% cc1=aa{3};
% cc2=aa{2};
%%
% guess_para1=[1,4,7,10,13,16,19,22,25,28];
% guess_para2=[4,8,12,16,20,24,28,32,36,40];
% guess_para1_a=randperm(guess_para1,10);
% guess_para2_a=randperm(guess_para2,10);
%%















