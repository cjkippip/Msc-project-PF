function lvState = LV_2para_SMC(t, state, P_Para)
x=state(1);
y=state(2);


% alpha=1;
% beta=0.05;
delta=0.02;
gamma=0.5;

dxdt = P_Para(1) * x - P_Para(2) * x * y;
dydt = delta * x * y - gamma * y;

lvState=[dxdt;dydt];

end