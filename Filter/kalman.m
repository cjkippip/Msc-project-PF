clc  
clear
  
% ��ʼ������  
delta_t=0.1;   %����ʱ��  
t=0:delta_t:5;  
N = length(t); % ���еĳ���  
sz = [2,N];    % �ź��迪�ٵ��ڴ�ռ��С  2��*N��  2:Ϊ״̬������ά��n  
g=5;          %���ٶ�ֵ  
x=1/2*g*t.^2;      %ʵ����ʵλ��  
z = x + sqrt(10).*randn(1,N); % ����ʱ�������������  
  
Q =[0 0;0 9e-1]; %���轨����ģ��  ��������������ٶ��� ��СΪn*n���� n=״̬������ά��  
R = 10;    % λ�ò���������ƣ����Ըı���������ͬЧ��  m*m      m=z(i)��ά��  
  
A=[1 delta_t;0 1];  % n*n  
B=[1/2*delta_t^2;delta_t];  
H=[1,0];            % m*n  
  
n=size(Q);  %nΪһ��1*2������  QΪ����  
m=size(R);  
  
% ����ռ�  
xhat=zeros(sz);       % x�ĺ������  
P=zeros(n);           % ���鷽�����  n*n  
xhatminus=zeros(sz);  % x���������  
Pminus=zeros(n);      % n*n  
K=zeros(n(1),m(1));   % Kalman����  n*m  
I=eye(n);  
  
% ���Ƶĳ�ʼֵ��ΪĬ�ϵ�0����P=[0 0;0 0],xhat=0  
for k = 9:N           %���賵���Ѿ��˶�9��delta_T�ˣ����ǲſ�ʼ����  
    % ʱ����¹���  
    xhatminus(:,k) = A*xhat(:,k-1)+B*g;  
    Pminus= A*P*A'+Q;  
      
    % �������¹���  
    K = Pminus*H'*inv( H*Pminus*H'+R );  
    xhat(:,k) = xhatminus(:,k)+K*(z(k)-H*xhatminus(:,k));  
    P = (I-K*H)*Pminus;  
end  
%%  
figure(1),clf,  
plot(t,z,'LineWidth',1.5);  
hold on  
plot(t,xhat(1,:),'r-','LineWidth',1.5)  
plot(t,x(1,:),'g-','LineWidth',1.5);  
legend({'measurement(with noise)', 'estimation', 'reality'},...
    'Location','northwest','FontSize',11,'FontWeight','bold');  
xlabel('Iteration(t=1:0.1:5)','FontSize',13,'FontWeight','bold'); 
ylabel('Distancce','FontSize',13,'FontWeight','bold'); 
hold off
grid on
grid minor
%%
clc  
clear  
  
% ��ʼ������  
delta_t=0.1;  
t=0:delta_t:5;  
g=5;%���ٶ�ֵ  
n_iter = length(t); % ���еĳ���  
sz = [n_iter, 1]; % �ź��迪�ٵ��ڴ�ռ��С  
x=1/2*g*t.^2;  
x=x';  
z = x + sqrt(10).*randn(sz); % ����ʱ�������������  
  
Q = 5; % ���̼�����������     
         %ע��Qֵ�øı�  ��������2������Ч�����Աȿ�Ч��ʱ���޸Ĵ��벻Ҫ�ı�z��ֵ  
R = 10; % ����������ƣ����Ըı���������ͬЧ��  
   
% ����ռ�  
xhat=zeros(sz);      % x�ĺ������  
P=zeros(sz);         % ���鷽�����  
xhatminus=zeros(sz); % x���������  
Pminus=zeros(sz);    % ���鷽�����  
K=zeros(sz);         % Kalman����  
   
% ���Ƶĳ�ʼֵ  
xhat(1) = 0.0;  
P = 1.0;  
for k = 2:n_iter   %  
    % ʱ����¹���  
    xhatminus(k) = xhat(k-1);  
    Pminus(k) = P(k-1)+Q;  
      
    % �������¹���  
    K(k) = Pminus(k)/( Pminus(k)+R );  
    xhat(k) = xhatminus(k)+K(k)*(z(k)-xhatminus(k));  
    P(k) = (1-K(k))*Pminus(k);  
end  
   
figure(2),clf,  
plot(t,z,'LineWidth',1.5);  
hold on  
plot(t,xhat,'r-','LineWidth',1.5)  
plot(t,x,'g-','LineWidth',1.5);  
legend({'measurement(with noise)', 'estimation', 'reality'},...
    'Location','northwest','FontSize',11,'FontWeight','bold');    
xlabel('Iteration(t=1:0.1:5)','FontSize',13,'FontWeight','bold'); 
ylabel('Distancce','FontSize',13,'FontWeight','bold');  
hold off
grid on
grid minor
