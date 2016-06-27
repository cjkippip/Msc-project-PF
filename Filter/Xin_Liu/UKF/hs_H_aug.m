function y = hs_H_aug(x)
%%%------------------------------------------------------------------------
%
%   Description: This is observation function for system. Only the first 
%   and the third one can be observed. This is implemented by setting the
%   obs_index which is 1 and 3. The dimension of input x is 5 and para_dim
%   is 2, thus para_dim+obs_index is the third and fifth element in input
%   variable x.
%
%   Input: 
%           x: is the extended state space. The first row until par_dim is
%           representing the parameter estimation. The rest are state
%           variables.
%
%   Output:
%           y: the matrix is for observation
%   
%
%   Author: Xin Liu
%
%   Date:   08/10/2010
%
%%%------------------------------------------------------------------------
global obs_index;
global para_dim;

%   The first and third element can be observed
y = x(para_dim+obs_index,:);