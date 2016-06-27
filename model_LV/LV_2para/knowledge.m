%% diffenential, derivative, Jacobian, Hessian



%% sin(x) and diff 1,2
x=0:0.01:(4*pi);
y=sin(x);
yx=diff(y)*100;
yxx=diff(yx)*100;

figure(1),clf,
plot(x,y,'r');
hold on
plot(x(2:end),yx,'b');
plot(x(3:end),yxx,'g');
legend('sin(x)','yx','yxx');
grid on
grid minor
hold off
%% y=1/(1+exp(-x))
x=-5:0.01:5;
y=1./(1+exp(-x));
yx=diff(y)*100;
yxx=diff(yx)*100;

figure(1),clf,
plot(x,y,'r');
hold on
plot(x(2:end),yx,'b');
plot(x(3:end),yxx,'g');
legend('y','yx','yxx');
grid on
grid minor
hold off
%% y=x^3
x=-5:0.01:5;
y=x.^4;
yx=diff(y)*100;
yxx=diff(yx)*100;

figure(1),clf,
plot(x,y,'r');
hold on
plot(x(2:end),yx,'b');
plot(x(3:end),yxx,'g');
legend('y','yx','yxx');
grid on
grid minor
hold off
%% normal distribution
m=0;
C=1;
x=-5:0.01:5;
y=1/(sqrt(2*pi)*C)*exp(-1/2*(x-m).^2/C^2);
yx=diff(y)*100;
yxx=diff(yx)*100;

figure(1),clf,
plot(x,y,'r');
hold on
plot(x(2:end),yx,'b');
plot(x(3:end),yxx,'g');
legend('y','yx','yxx');
grid on
grid minor
hold off
%% cholesky decomposition
x=[2 -2;-2 5];
[R,p] = chol(x);

mu=[0;0];
sigma=[1 0;0 1];
C=[3 2;2 3];
C1=sqrt(C);
C2=chol(C);
aa1=mvnrnd(mu,sigma,1000)*C1;
aa2=mvnrnd(mu,sigma,1000)*C2;
aa3=mvnrnd(mu,C,1000);

figure(6),clf,
scatter(aa1(:,1),aa1(:,2));
figure(7),clf,
scatter(aa2(:,1),aa2(:,2));
figure(8),clf,
scatter(aa3(:,1),aa3(:,2));
%% Poisson Distribution
x=0:10;
lambda=5;
y=exp(-lambda).*lambda.^x./factorial(x);
z=cumsum(y);
figure(10),clf,
plot(x,y);
%% generate gaussian distribution
% a=[2 1;1 2];
% [b,c]=eig(a);
mu=[0,0];
sigma=[2 1.5;1.5 2];
[eigVec,lambda]=eig(sigma);
pointGaus=mvnrnd(mu,sigma,500);
figure(1),clf,
scatter(pointGaus(:,1),pointGaus(:,2));
grid on
grid minor
%% Hessian
syms x y 
f = x^2*y + 2*y*x;
hess=hessian(f,[x,y]);
jaco=jacobian(f,[x,y]);

% bb=double(aa2(1,1));
%% Jacobian
% syms x y z
% jacobian([x*y*z, y^2, x + z], [x, y, z])
% jacobian([x*y*z, y^2, x + z], [x; y; z])
% syms x y z
% jacobian(2*x + 3*y + 4*z, [x, y, z])
% gradient(2*x + 3*y + 4*z, [x, y, z])
% syms x y
% jacobian([x^2*y, x*sin(y)], x)
% diff([x^2*y, x*sin(y)], x)

delta_t=0.1;
f=@(x)[x(1)+delta_t*(x(1)*x(3)-x(4)*x(1)*x(2));...
            x(2)+delta_t*(0.02*x(1)*x(2)-0.5*x(2));...
            x(3);...
            x(4)];
x=[30;30;27;19];
z=f(x);
n=numel(x);% elements' number of x
m=numel(z);
A=zeros(m,n);
h=n*eps;
for k=1:n
    x1=x;
    x1(k)=x1(k)+h*1i;
    A(:,k)=imag(f(x1))/h;
end
%% Jacobian
delta_t=0.1;
f=@(x)[x(1)+delta_t*(x(1)*x(3)-x(4)*x(1)*x(2));...
            x(2)+delta_t*(0.02*x(1)*x(2)-0.5*x(2));...
            x(3);...
            x(4)];

x=[30;30;27;19];
syms x1 x2 x3 x4
xSym=[x1;x2;x3;x4];
F=f(xSym);
grad=jacobian(F,xSym);
double(grad(x))
% Y=subs(grad,xSym,x);
% aa=double(Y);





