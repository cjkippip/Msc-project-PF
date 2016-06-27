function parameter = updatePara_all_para(P_Para,covP_Para,meanP_Para,para_dim)
%%%-----------------------------------------------------------------------
%       
%   Description: This is for updating the parameter. Refer to the paper, we
%                update the parameters through the constructed artificial 
%                dynamics. 
%
%   Input:
%                 P_Para:      Parameters are going to be updated.
%                 covP_Para:   Variance vector of particles for parameters
%                 meanP_Para:  The mean of particles for parameters.
%                 para_dim:    number of parameters
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
%% parameter 1-para_dim
for i=1:para_dim
    P_Para(i,:) = a*P_Para(i,:) + (1-a) * meanP_Para(i,:) * ones(1,N) +...
        sqrt(covParaCoefVect(i) * covP_Para(i,:)) * randn(1,N);
end
%%
% % parameter 1
% P_Para(1,:) = a*P_Para(1,:) + (1-a) * meanP_Para(1,:) * ones(1,N) +...
%     sqrt(covParaCoefVect(1) * covP_Para(1,:)) * randn(1,N);
% 
% % parameter 2
% P_Para(2,:) = a*P_Para(2,:) + (1-a) * meanP_Para(2,:) * ones(1,N) +...
%     sqrt(covParaCoefVect(2) * covP_Para(2,:)) * randn(1,N);
% 
% % parameter 3
% P_Para(3,:) = a*P_Para(3,:) + (1-a) * meanP_Para(3,:) * ones(1,N) +...
%     sqrt(covParaCoefVect(3) * covP_Para(3,:)) * randn(1,N);
% 
% % parameter 4
% P_Para(4,:) = a*P_Para(4,:) + (1-a) * meanP_Para(4,:) * ones(1,N) +...
%     sqrt(covParaCoefVect(4) * covP_Para(4,:)) * randn(1,N);
% 
% % parameter 5
% P_Para(5,:) = a*P_Para(5,:) + (1-a) * meanP_Para(5,:) * ones(1,N) +...
%     sqrt(covParaCoefVect(5) * covP_Para(5,:)) * randn(1,N);
% 
% % parameter 6
% P_Para(6,:) = a*P_Para(6,:) + (1-a) * meanP_Para(6,:) * ones(1,N) +...
%     sqrt(covParaCoefVect(6) * covP_Para(6,:)) * randn(1,N);
% 
% % parameter 7
% P_Para(7,:) = a*P_Para(7,:) + (1-a) * meanP_Para(7,:) * ones(1,N) +...
%     sqrt(covParaCoefVect(7) * covP_Para(7,:)) * randn(1,N);
% 
% % parameter 8
% P_Para(8,:) = a*P_Para(8,:) + (1-a) * meanP_Para(8,:) * ones(1,N) +...
%     sqrt(covParaCoefVect(8) * covP_Para(8,:)) * randn(1,N);
% 
% % parameter 9
% P_Para(4,:) = a*P_Para(9,:) + (1-a) * meanP_Para(9,:) * ones(1,N) +...
%     sqrt(covParaCoefVect(9) * covP_Para(9,:)) * randn(1,N);
% 
% % parameter 10
% P_Para(10,:) = a*P_Para(10,:) + (1-a) * meanP_Para(10,:) * ones(1,N) +...
%     sqrt(covParaCoefVect(10) * covP_Para(10,:)) * randn(1,N);

parameter = P_Para;        






