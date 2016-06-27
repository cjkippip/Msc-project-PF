deltaT = 0.2;   % The sampling time interval
N = 500;        % number of used particles

% Covariance of states
covState1 = 10000;
covState2 = 10;
covState3 = 1250000; 

% Covariance of measurement
covMeasure = 1000000000;                    

%   Parameters for SMC algorithm
e = 0.96;
a = (3*e - 1)/(2*e);
covParaCoeff = 1 - a^2;
covParaCoefVect = [0.045;covParaCoeff;0.05;0.045;covParaCoeff;0.045];