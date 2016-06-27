function parameter = updatePara_all_para(P_Para,covP_Para,meanP_Para)
%%%-----------------------------------------------------------------------
%       
%   Description: This is for updating the parameter. Refer to the paper, we
%                update the parameters through the constructed artificial 
%                dynamics. 
%
%   Input:
%                 para:        Parameters are going to be updated.
%                 covParaSmc:  Variance vector of particles for parameters
%                 realParaSmc: The mean of particles for parameters.
%   
%   Output:
%                 parameter:   Matrix contains updated values for
%                              parameters.
%
%   Date: 25/04/2011
%
%   Author: Xin Liu
%
%%%------------------------------------------------------------------------

%   Global variables declaration
global N;
global a;
global covParaCoefVect;

%   Update for alpha with (1-a)
P_Para(1,:) = a*P_Para(1,:) + (1-a) * meanP_Para(1,:) * ones(1,N) +...
    sqrt(covParaCoefVect(1,:) * covP_Para(1,:)) * randn(1,N);

%   Update for alpha
% P_Para(1,:) = a*P_Para(1,:) +...
%     sqrt(covParaCoefVect(1,:) * covP_Para(1,:)) * randn(1,N);

%   Update for beta with (1-a)
P_Para(2,:) = a*P_Para(2,:) + (1-a) * meanP_Para(2,:) * ones(1,N) +...
    sqrt(covParaCoefVect(2,:) * covP_Para(2,:)) * randn(1,N);

%   Update for beta
% P_Para(2,:) = a*P_Para(2,:) +...
%     sqrt(covParaCoefVect(2,:) * covP_Para(2,:)) * randn(1,N);

parameter = P_Para;        






