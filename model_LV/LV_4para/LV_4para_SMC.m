function lvState = LV_4para_SMC(t, state, P_Para)
x=state(1);
y=state(2);


% alpha=1;
% beta=0.05;
% delta=0.02;
% gamma=0.5;

dxdt = P_Para(1) * x - P_Para(2) * x * y;
dydt = P_Para(3) * x * y - P_Para(4) * y;

lvState=[dxdt;dydt];

end