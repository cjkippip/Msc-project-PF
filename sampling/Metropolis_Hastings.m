n = 250000;
x = zeros(n, 1);
x(1) = 0.5;
for i = 1: n-1
    x_c = normrnd(x(i), 0.05);% q(x_c|x(i))
%     ratio=(normpdf(x_c)*normpdf(x(i),x_c,0.05))/(normpdf(x(i))*normpdf(x_c,x(i),0.05));
    ratio=normpdf(x_c)/normpdf(x(i));
    if rand < min(1, ratio)
        x(i+1) = x_c;% accept
    else
        x(i+1) = x(i);% reject
    end
end
%%
m=mean(x);
v=var(x);
s=std(x);
%%
figure(1),clf,
hist(x,50);
% hold on
% plot();




