% lotka_volterra.m
function lvState = LV_4para(t, state)
x=state(1);
y=state(2);

alpha=state(3);
beta=state(4);
delta=state(5);
gamma=state(6);

% alpha=1;
% beta=0.05;
% delta=0.02;
% gamma=0.5;

% alpha=2;
% beta=0.6;
% delta=0.4;
% gamma=2.2;

% alpha=3;
% beta=2;
% delta=0.5;
% gamma=2;

dxdt = alpha * x - beta * x * y;
dydt = delta * x * y - gamma * y;

lvState=[dxdt;dydt;zeros(4,1)];
end

  






