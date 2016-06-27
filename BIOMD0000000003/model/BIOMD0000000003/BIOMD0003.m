function xdot=BIOMD0003(t,x)
% End Matlab code

% Start Octave code
%function xdot=f(x,t)
% End Octave code

% Compartment: id = cell, name = cell, constant
	compartment_cell=1.0;
% Parameter:   id =  V1, name = V1
% Parameter:   id =  V3, name = V3
% Parameter:   id =  VM1, name = VM1
	global_par_VM1=3.0;
% Parameter:   id =  VM3, name = VM3
	global_par_VM3=1.0;
% Parameter:   id =  Kc, name = Kc
	global_par_Kc=0.5;
% assignmentRule: variable = V1
	global_par_V1=x(1)*global_par_VM1*(x(1)+global_par_Kc)^(-1);
% assignmentRule: variable = V3
	global_par_V3=x(2)*global_par_VM3;

% Reaction: id = reaction1, name = creation of cyclin	% Local Parameter:   id =  vi, name = vi
	reaction_reaction1_vi=0.025;

	reaction_reaction1=compartment_cell*reaction_reaction1_vi;

% Reaction: id = reaction2, name = default degradation of cyclin	% Local Parameter:   id =  kd, name = kd
	reaction_reaction2_kd=0.01;

	reaction_reaction2=x(1)*compartment_cell*reaction_reaction2_kd;

% Reaction: id = reaction3, name = cdc2 kinase triggered degration of cyclin	% Local Parameter:   id =  vd, name = vd
	reaction_reaction3_vd=0.25;
	% Local Parameter:   id =  Kd, name = Kd
	reaction_reaction3_Kd=0.02;

	reaction_reaction3=x(1)*compartment_cell*reaction_reaction3_vd*x(3)*(x(1)+reaction_reaction3_Kd)^(-1);

% Reaction: id = reaction4, name = activation of cdc2 kinase	% Local Parameter:   id =  K1, name = K1
	reaction_reaction4_K1=0.005;

	reaction_reaction4=compartment_cell*(1- 1*x(2))*global_par_V1*(reaction_reaction4_K1- 1*x(2)+1)^(-1);

% Reaction: id = reaction5, name = deactivation of cdc2 kinase	% Local Parameter:   id =  V2, name = V2
	reaction_reaction5_V2=1.5;
	% Local Parameter:   id =  K2, name = K2
	reaction_reaction5_K2=0.005;

	reaction_reaction5=compartment_cell*x(2)*reaction_reaction5_V2*(reaction_reaction5_K2+x(2))^(-1);

% Reaction: id = reaction6, name = activation of cyclin protease	% Local Parameter:   id =  K3, name = K3
	reaction_reaction6_K3=0.005;

	reaction_reaction6=compartment_cell*global_par_V3*(1- 1*x(3))*(reaction_reaction6_K3- 1*x(3)+1)^(-1);

% Reaction: id = reaction7, name = deactivation of cyclin protease	% Local Parameter:   id =  K4, name = K4
	reaction_reaction7_K4=0.005;
	% Local Parameter:   id =  V4, name = V4
	reaction_reaction7_V4=0.5;

	reaction_reaction7=compartment_cell*reaction_reaction7_V4*x(3)*(reaction_reaction7_K4+x(3))^(-1);

	xdot=zeros(3,1);
	
% Species:   id = C, name = Cyclin, affected by kineticLaw
	xdot(1) = (1/(compartment_cell))*(( 1.0 * reaction_reaction1) + (-1.0 * reaction_reaction2) + (-1.0 * reaction_reaction3));
	
% Species:   id = M, name = CDC-2 Kinase, affected by kineticLaw
	xdot(2) = (1/(compartment_cell))*(( 1.0 * reaction_reaction4) + (-1.0 * reaction_reaction5));
	
% Species:   id = X, name = Cyclin Protease, affected by kineticLaw
	xdot(3) = (1/(compartment_cell))*(( 1.0 * reaction_reaction6) + (-1.0 * reaction_reaction7));
end