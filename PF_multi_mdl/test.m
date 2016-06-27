

clc, clear;
% ID of biomodel
ID_mdl='BIOMD0000000002';
% add path
addpath (genpath('/Users/CJ/Documents/MATLAB/libSBML'))
addpath (genpath('/Users/CJ/Documents/MATLAB/SBMLToolbox-4.1.0'))
addpath (genpath(['./model/',ID_mdl]))
% translate SBML to Matlab struct
SBMLModel=TranslateSBML(['./model/',ID_mdl,'/',ID_mdl,'.xml']);
% get all parameters
[allParaName1, allParaValue1] = GetAllParameters(SBMLModel);
% get all parameters unique
[allParaName2, allParaValue2] = GetAllParametersUnique(SBMLModel);
% get global parameters unique
[allParaName3, allParaValue3] = GetGlobalParameters(SBMLModel);
% Get Parameter Assignment Rules
[allParaName4, allParaValue4] = GetParameterAssignmentRules(SBMLModel);
% Get Parameter Algebraic Rules
[allParaName5, allParaValue5] = GetParameterAlgebraicRules(SBMLModel);
% get all species
[allSpeciesName, allSpeciesValue] = GetSpecies(SBMLModel);

%% determine species role in reaction
SR=DetermineSpeciesRoleInReaction(SBMLModel.species(1), SBMLModel.reaction(1));
%% compartment
[compartName, compartValue] = GetCompartments(SBMLModel);
%% compartment type
% compartmentType = GetCompartmentTypes(mdl_NEGF);
%% get global parameters
[globalParaName, globalParaValue] = GetGlobalParameters(SBMLModel);
%% Algebraic
[Algebraic, algebraicRules] = GetParameterAlgebraicRules(SBMLModel);
%% Assignment
[Assignment, assignmentRules] = GetParameterAssignmentRules(SBMLModel);
%% parameter from reaction
[paraName1, paraValue1] = GetParameterFromReaction(SBMLModel.reaction(9));
%%












