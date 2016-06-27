% MATLAB code for implementing a simple importance sampling scheme
% Ryan Guerra
% Nov 17, 2008
clear;
mu = 6; % gaussian mean parameter
sigma = sqrt(1.5); % gaussian standard deviation parameter
k=1.65; % chi DOF parameter
c=2; % multiplicative constant to make p < cq
num_iterations = 100; % number of Monte Carlo runs to execute
p = @(x) x^(k-1)*exp(-x^2/2); % target distribution
q = @(x) c/sqrt(2*pi*sigma^2)... % sampling distribution
*exp(-(mu-x)^2/(2*sigma^2));
f = @(x) 2*sin(pi/1.5*x); % function we want to take expectation of

for iter=1:num_iterations
    samp_size = 1000; % number of initial samples to draw from q
    X = normrnd(mu,sigma,[samp_size 1]); % generate vector of samples from q
    for i=1:length(X)
        if X(i)>=0 % we will discard all sampled X that are not in the support of p
            W(i) = p(X(i))/q(X(i)); % calculate importance weight for sample i
            I(i) = W(i)*f(X(i)); % weigh our function
        else
            W(i)=0; % these values will be ignored
            I(i)=0;
            samp_size = samp_size - 1;
        end
    end
    I_hat = sum(I)/sum(W); % perform summation in (5)
    [tmp tmp2 non_zero_weights] = find(W); % remove all discarded importance weights
    variance = var(non_zero_weights); % calculate the variance of the importance weights
    eff_samp_size = samp_size/(1 + variance); % calculate the effective sample size
    format short g % store the results of this run
    results(iter,:) = [I_hat variance samp_size eff_samp_size];
end

% print compiled results
totals = [mean(results(:,1)) mean(results(:,2)) mean(results(:,3)) mean(results(:,4));...
var(results(:,1)) var(results(:,2)) var(results(:,3)) var(results(:,4))];
%%
% figure(1),clf,
% ezplot(p,'b');
% hold on
% ezplot(q,'g');
% ezplot(f,'r');
% hold off
%%
% xx1=linspace(-6,-6,500);
% yy1=ones(500,1);
% yy2=ones(500,1);
% yy3=ones(500,1);
% for i=1:500
%     yy1(i)=p(xx1(i));
%     yy2(i)=q(xx1(i));
%     yy3(i)=f(xx1(i));
% end
% figure(2),clf,
% plot(xx1.yy1,'b');
% hold on
% plot(xx1.yy2,'g');
% plot(xx1.yy3,'r');
% hold off




