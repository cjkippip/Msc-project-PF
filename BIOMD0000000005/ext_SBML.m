%% extract resource from SBML file
% ID of biomodel
ID_mdl='BIOMD0000000005';

% add path
addpath (genpath('/Users/CJ/Documents/MATLAB/libSBML'))
addpath (genpath('/Users/CJ/Documents/MATLAB/SBMLToolbox-4.1.0'))
addpath (genpath(['./model/',ID_mdl]))

% translate SBML to Matlab struct
SBMLModel=TranslateSBML(['./model/',ID_mdl,'/',ID_mdl,'.xml']);

% get all parameters
[allParaName, allParaValue] = GetGlobalParameters(SBMLModel);
allParaName=allParaName';
allParaValue=allParaValue';



