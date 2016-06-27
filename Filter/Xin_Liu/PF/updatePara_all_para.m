function parameter = updatePara_all_para(para,covParaSmc,realParaSmc)
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

%   Update for Ku
para(1,:) = a*para(1,:) + sqrt(covParaCoefVect(1,:) * covParaSmc(1,:))* ...
            randn(1,N);

%   Update for As
para(2,:) = a*para(2,:) + (1-a) * realParaSmc(2,:) * ones(1,N)...
            + sqrt(covParaCoefVect(2,:) * covParaSmc(2,:)) * randn(1,N);

%   Update for Ks
para(3,:) = a*para(3,:) + sqrt(covParaCoefVect(3,:) * covParaSmc(3,:))* ...
            randn(1,N);
        
%   Update for A0
para(4,:) = a*para(4,:) + sqrt(covParaCoefVect(4,:) * covParaSmc(4,:))* ...
            randn(1,N);
        
%   Update for Kd
para(5,:) = a*para(5,:) + (1-a) * realParaSmc(5,:) * ones(1,N)...
            + sqrt(covParaCoefVect(5,:) * covParaSmc(5,:)) * randn(1,N);

%   Update for Ad
para(6,:) = a*para(6,:) + sqrt(covParaCoefVect(6,:) * covParaSmc(6,:))* ...
            randn(1,N);

parameter = para;        






