state_dim = 2;                      % state dimension
para_dim = 4;                       % parameter dimension
ext_state_dim = state_dim+para_dim; % extended dimension
obs_dim = 2;                        % 1 and 2 state can be observed
N = 2000;                            % number of particles 


e = 0.88;
a = (3*e-1)/(2*e);
covParaCoef = 1-a^2;
covParaCoefVect = [covParaCoef;covParaCoef;covParaCoef;covParaCoef];

% Convariance of state
covState1 = 100;
covState2 = 50;

% Convariance of measurement
covMeasure = 1000;

% time length and sampling rate
endT = 30;
tspan = linspace(0,endT,301);
delta_t = tspan(2)-tspan(1);
T = length(tspan); % iteration time

%% Generate true state, parameter and observation
state1=10;
state2=10;
alphaStar=1;
betaStar=0.05;
deltaStar=0.02;
gammaStar=0.5;

% particles of parameters
P_Para=ones(para_dim,N);
P_Para(1,:)=3+sqrt(3)*randn(1,N);
P_Para(2,:)=6+sqrt(2)*randn(1,N);
P_Para(3,:)=6+sqrt(2)*randn(1,N);
P_Para(4,:)=3+sqrt(2)*randn(1,N);