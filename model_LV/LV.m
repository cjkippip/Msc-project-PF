% lotka_volterra2.m
function lvState = LV(t, state)

x=state(1);
y=state(2);

alpha=1;
beta=0.05;
delta=0.02;
gamma=0.5;

dxdt = alpha * x - beta * x * y;
dydt = delta * x * y - gamma * y;
lvState=[dxdt;dydt];
end