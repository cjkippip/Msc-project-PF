function xdot=BIOMD0033_SMC(t,x,P_Para)
% End Matlab code

% Start Octave code
%function xdot=f(x,t)
% End Octave code

% Compartment: id = cell, name = cell, constant
	compartment_cell=1.0;
    
% Species:   id = RasGapActive, name = RasGapActive, constant	
const_species_RasGapActive=120000.0;

% Species:   id = RapGapActive, name = RapGapActive, constant	
const_species_RapGapActive=120000.0;

% Species:   id = PP2AActive, name = PP2AActive, constant	
const_species_PP2AActive=120000.0;

% Species:   id = Raf1PPtase, name = Raf1PPtase, constant	
const_species_Raf1PPtase=120000.0;

% Reaction: id = EGFBindingReaction, name = EGF binding
	reaction_EGFBindingReaction=compartment_cell*P_Para(1)*x(1)*x(3);

% Reaction: id = EGFUnbindingReaction, name = EFG unbinding
	reaction_EGFUnbindingReaction=compartment_cell*P_Para(2)*x(4);

% Reaction: id = NGFBindingReaction, name = NGF binding
	reaction_NGFBindingReaction=P_Para(3)*x(2)*x(5)*compartment_cell;

% Reaction: id = NGFUnbindingReaction, name = NGF unbinding
	reaction_NGFUnbindingReaction=P_Para(4)*x(6)*compartment_cell;

% Reaction: id = SosActivationByEGFReaction, name = SOS activation by EGF
	reaction_SosActivationByEGFReaction=compartment_cell*P_Para(5)*x(4)*x(7)/(x(7)+P_Para(6));

% Reaction: id = SosActivationByNGFReaction, name = SOS activation by NGF
	reaction_SosActivationByNGFReaction=compartment_cell*P_Para(7)*x(6)*x(7)/(x(7)+P_Para(8));

% Reaction: id = SosDeactivationReaction, name = SOS deactivation
	reaction_SosDeactivationReaction=compartment_cell*P_Para(9)*x(10)*x(8)/(x(8)+P_Para(10));

% Reaction: id = RasActivationReaction, name = Ras activation
	reaction_RasActivationReaction=compartment_cell*P_Para(11)*x(8)*x(11)/(x(11)+P_Para(12));

% Reaction: id = RasDeactivationReaction, name = Ras deactivation
	reaction_RasDeactivationReaction=compartment_cell*P_Para(13)*const_species_RasGapActive*x(12)/(x(12)+P_Para(14));

% Reaction: id = Raf1ByRasActivationReaction, name = Raf1 activation by Ras
	reaction_Raf1ByRasActivationReaction=compartment_cell*P_Para(15)*x(12)*x(13)/(x(13)+P_Para(16));

% Reaction: id = MekbyRaf1ActivationReaction, name = Mek activation by Raf1
	reaction_MekbyRaf1ActivationReaction=compartment_cell*P_Para(17)*x(14)*x(17)/(x(17)+P_Para(18));

% Reaction: id = MekbyBRafActivationReaction, name = Mek activation by B-Raf
	reaction_MekbyBRafActivationReaction=compartment_cell*P_Para(19)*x(16)*x(17)/(x(17)+P_Para(20));

% Reaction: id = ErkActivationReaction, name = Erk activation
	reaction_ErkActivationReaction=compartment_cell*P_Para(23)*x(18)*x(19)/(x(19)+P_Para(24));

% Reaction: id = MekDeactivationReaction, name = Mek deactivation
	reaction_MekDeactivationReaction=compartment_cell*P_Para(21)*const_species_PP2AActive*x(18)/(x(18)+P_Para(22));

% Reaction: id = ErkDeactivationReaction, name = Erk deactivation
	reaction_ErkDeactivationReaction=compartment_cell*P_Para(25)*const_species_PP2AActive*x(20)/(x(20)+P_Para(26));

% Reaction: id = Raf1byPPtaseDeactivationReaction, name = Raf1 deactivation by PPase
	reaction_Raf1byPPtaseDeactivationReaction=compartment_cell*P_Para(45)*const_species_Raf1PPtase*x(14)/(x(14)+P_Para(46));

% Reaction: id = P90RskActivationReaction, name = P90Rsk activation
	reaction_P90RskActivationReaction=compartment_cell*P_Para(27)*x(20)*x(9)/(x(9)+P_Para(28));

% Reaction: id = PI3KbyEGFRActivationReaction, name = PI3K activation by EGFR
	reaction_PI3KbyEGFRActivationReaction=compartment_cell*P_Para(29)*x(4)*x(21)/(x(21)+P_Para(30));

% Reaction: id = PI3KbyRasActivationReaction, name = PI3K activation by Ras
	reaction_PI3KbyRasActivationReaction=compartment_cell*P_Para(31)*x(12)*x(21)/(x(21)+P_Para(32));

% Reaction: id = AktActivationReaction, name = Akt activation
	reaction_AktActivationReaction=compartment_cell*P_Para(33)*x(22)*x(23)/(x(23)+P_Para(34));

% Reaction: id = Raf1ByAktDeactivationReaction, name = Raf1 deactivation by Akt
	reaction_Raf1ByAktDeactivationReaction=compartment_cell*P_Para(35)*x(24)*x(14)/(x(14)+P_Para(36));

% Reaction: id = C3GActivationReaction, name = C3G activation
	reaction_C3GActivationReaction=compartment_cell*P_Para(37)*x(6)*x(25)/(x(25)+P_Para(38));

% Reaction: id = Rap1ActivationReaction, name = Rap1 activation
	reaction_Rap1ActivationReaction=compartment_cell*P_Para(39)*x(26)*x(27)/(x(27)+P_Para(40));

% Reaction: id = Rap1DeactivationReaction, name = Rap1 deactivation
	reaction_Rap1DeactivationReaction=compartment_cell*P_Para(41)*const_species_RapGapActive*x(28)/(x(28)+P_Para(42));

% Reaction: id = BRafByRap1ActivationReaction, name = BRaf activation by Rap1
	reaction_BRafByRap1ActivationReaction=compartment_cell*P_Para(43)*x(28)*x(15)/(x(15)+P_Para(44));

% Reaction: id = BRafbyPPtaseDeactivationReaction, name = BRaf deactivation by PPase
	reaction_BRafbyPPtaseDeactivationReaction=compartment_cell*P_Para(47)*const_species_Raf1PPtase*x(16)/(x(16)+P_Para(48));



	xdot=zeros(28,1);
	
% Species:   id = EGF, name = EGF, affected by kineticLaw
	xdot(1) = (1/(compartment_cell))*((-1.0 * reaction_EGFBindingReaction) + ( 1.0 * reaction_EGFUnbindingReaction));
	
% Species:   id = NGF, name = NGF, affected by kineticLaw
	xdot(2) = (1/(compartment_cell))*((-1.0 * reaction_NGFBindingReaction) + ( 1.0 * reaction_NGFUnbindingReaction));
	
% Species:   id = freeEGFReceptor, name = freeEGFReceptor, affected by kineticLaw
	xdot(3) = (1/(compartment_cell))*((-1.0 * reaction_EGFBindingReaction) + ( 1.0 * reaction_EGFUnbindingReaction));
	
% Species:   id = boundEGFReceptor, name = boundEGFReceptor, affected by kineticLaw
	xdot(4) = (1/(compartment_cell))*(( 1.0 * reaction_EGFBindingReaction) + (-1.0 * reaction_EGFUnbindingReaction));
	
% Species:   id = freeNGFReceptor, name = freeNGFReceptor, affected by kineticLaw
	xdot(5) = (1/(compartment_cell))*((-1.0 * reaction_NGFBindingReaction) + ( 1.0 * reaction_NGFUnbindingReaction));
	
% Species:   id = boundNGFReceptor, name = boundNGFReceptor, affected by kineticLaw
	xdot(6) = (1/(compartment_cell))*(( 1.0 * reaction_NGFBindingReaction) + (-1.0 * reaction_NGFUnbindingReaction));
	
% Species:   id = SosInactive, name = SosInactive, affected by kineticLaw
	xdot(7) = (1/(compartment_cell))*((-1.0 * reaction_SosActivationByEGFReaction) + (-1.0 * reaction_SosActivationByNGFReaction) + ( 1.0 * reaction_SosDeactivationReaction));
	
% Species:   id = SosActive, name = SosActive, affected by kineticLaw
	xdot(8) = (1/(compartment_cell))*(( 1.0 * reaction_SosActivationByEGFReaction) + ( 1.0 * reaction_SosActivationByNGFReaction) + (-1.0 * reaction_SosDeactivationReaction));
	
% Species:   id = P90RskInactive, name = P90RskInactive, affected by kineticLaw
	xdot(9) = (1/(compartment_cell))*((-1.0 * reaction_P90RskActivationReaction));
	
% Species:   id = P90RskActive, name = P90RskActive, affected by kineticLaw
	xdot(10) = (1/(compartment_cell))*(( 1.0 * reaction_P90RskActivationReaction));
	
% Species:   id = RasInactive, name = RasInactive, affected by kineticLaw
	xdot(11) = (1/(compartment_cell))*((-1.0 * reaction_RasActivationReaction) + ( 1.0 * reaction_RasDeactivationReaction));
	
% Species:   id = RasActive, name = RasActive, affected by kineticLaw
	xdot(12) = (1/(compartment_cell))*(( 1.0 * reaction_RasActivationReaction) + (-1.0 * reaction_RasDeactivationReaction));
	
% Species:   id = Raf1Inactive, name = Raf1Inactive, affected by kineticLaw
	xdot(13) = (1/(compartment_cell))*((-1.0 * reaction_Raf1ByRasActivationReaction) + ( 1.0 * reaction_Raf1byPPtaseDeactivationReaction) + ( 1.0 * reaction_Raf1ByAktDeactivationReaction));
	
% Species:   id = Raf1Active, name = Raf1Active, affected by kineticLaw
	xdot(14) = (1/(compartment_cell))*(( 1.0 * reaction_Raf1ByRasActivationReaction) + (-1.0 * reaction_Raf1byPPtaseDeactivationReaction) + (-1.0 * reaction_Raf1ByAktDeactivationReaction));
	
% Species:   id = BRafInactive, name = BRafInactive, affected by kineticLaw
	xdot(15) = (1/(compartment_cell))*((-1.0 * reaction_BRafByRap1ActivationReaction) + ( 1.0 * reaction_BRafbyPPtaseDeactivationReaction));
	
% Species:   id = BRafActive, name = BRafActive, affected by kineticLaw
	xdot(16) = (1/(compartment_cell))*(( 1.0 * reaction_BRafByRap1ActivationReaction) + (-1.0 * reaction_BRafbyPPtaseDeactivationReaction));
	
% Species:   id = MekInactive, name = MekInactive, affected by kineticLaw
	xdot(17) = (1/(compartment_cell))*((-1.0 * reaction_MekbyRaf1ActivationReaction) + (-1.0 * reaction_MekbyBRafActivationReaction) + ( 1.0 * reaction_MekDeactivationReaction));
	
% Species:   id = MekActive, name = MekActive, affected by kineticLaw
	xdot(18) = (1/(compartment_cell))*(( 1.0 * reaction_MekbyRaf1ActivationReaction) + ( 1.0 * reaction_MekbyBRafActivationReaction) + (-1.0 * reaction_MekDeactivationReaction));
	
% Species:   id = ErkInactive, name = ErkInactive, affected by kineticLaw
	xdot(19) = (1/(compartment_cell))*((-1.0 * reaction_ErkActivationReaction) + ( 1.0 * reaction_ErkDeactivationReaction));
	
% Species:   id = ErkActive, name = ErkActive, affected by kineticLaw
	xdot(20) = (1/(compartment_cell))*(( 1.0 * reaction_ErkActivationReaction) + (-1.0 * reaction_ErkDeactivationReaction));
	
% Species:   id = PI3KInactive, name = PI3KInactive, affected by kineticLaw
	xdot(21) = (1/(compartment_cell))*((-1.0 * reaction_PI3KbyEGFRActivationReaction) + (-1.0 * reaction_PI3KbyRasActivationReaction));
	
% Species:   id = PI3KActive, name = PI3KActive, affected by kineticLaw
	xdot(22) = (1/(compartment_cell))*(( 1.0 * reaction_PI3KbyEGFRActivationReaction) + ( 1.0 * reaction_PI3KbyRasActivationReaction));
	
% Species:   id = AktInactive, name = AktInactive, affected by kineticLaw
	xdot(23) = (1/(compartment_cell))*((-1.0 * reaction_AktActivationReaction));
	
% Species:   id = AktActive, name = AktActive, affected by kineticLaw
	xdot(24) = (1/(compartment_cell))*(( 1.0 * reaction_AktActivationReaction));
	
% Species:   id = C3GInactive, name = C3GInactive, affected by kineticLaw
	xdot(25) = (1/(compartment_cell))*((-1.0 * reaction_C3GActivationReaction));
	
% Species:   id = C3GActive, name = C3GActive, affected by kineticLaw
	xdot(26) = (1/(compartment_cell))*(( 1.0 * reaction_C3GActivationReaction));
	
% Species:   id = Rap1Inactive, name = Rap1Inactive, affected by kineticLaw
	xdot(27) = (1/(compartment_cell))*((-1.0 * reaction_Rap1ActivationReaction) + ( 1.0 * reaction_Rap1DeactivationReaction));
	
% Species:   id = Rap1Active, name = Rap1Active, affected by kineticLaw
	xdot(28) = (1/(compartment_cell))*(( 1.0 * reaction_Rap1ActivationReaction) + (-1.0 * reaction_Rap1DeactivationReaction));
end


