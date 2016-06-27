function xdot=BIOMD0002_SMC(t,x,P_Para)
% End Matlab code

% Start Octave code
%function xdot=f(x,t)
% End Octave code

% Compartment: id = comp1, name = compartment1, constant
	compartment_comp1=1.0E-16;
% % Parameter:   id =  kf_0, name = kf_0
% 	P_Para(1)=3.0E8;
% % Parameter:   id =  kr_0, name = kr_0
% 	P_Para(2)=8000.0;
% % Parameter:   id =  kf_1, name = kf_1
% 	P_Para(3)=1.5E8;
% % Parameter:   id =  kr_1, name = kr_1
% 	P_Para(4)=16000.0;
% % Parameter:   id =  kf_2, name = kf_2
% 	P_Para(5)=30000.0;
% % Parameter:   id =  kr_2, name = kr_2
% 	P_Para(6)=700.0;
% % Parameter:   id =  kf_3, name = kf_3
% 	P_Para(7)=3.0E8;
% % Parameter:   id =  kr_3, name = kr_3
% 	P_Para(8)=8.64;
% % Parameter:   id =  kf_4, name = kf_4
% 	P_Para(9)=1.5E8;
% % Parameter:   id =  kr_4, name = kr_4
% 	P_Para(10)=17.28;
% % Parameter:   id =  kf_5, name = kf_5
% 	P_Para(11)=0.54;
% % Parameter:   id =  kr_5, name = kr_5
% 	P_Para(12)=10800.0;
% % Parameter:   id =  kf_6, name = kf_6
% 	P_Para(13)=130.0;
% % Parameter:   id =  kr_6, name = kr_6
% 	P_Para(14)=2740.0;
% % Parameter:   id =  kf_7, name = kf_7
% 	P_Para(15)=3.0E8;
% % Parameter:   id =  kr_7, name = kr_7
% 	P_Para(16)=4.0;
% % Parameter:   id =  kf_8, name = kf_8
% 	P_Para(17)=1.5E8;
% % Parameter:   id =  kr_8, name = kr_8
% 	P_Para(18)=8.0;
% % Parameter:   id =  kf_9, name = kf_9
% 	P_Para(19)=19.7;
% % Parameter:   id =  kr_9, name = kr_9
% 	P_Para(20)=3.74;
% % Parameter:   id =  kf_10, name = kf_10
% 	P_Para(21)=19.85;
% % Parameter:   id =  kr_10, name = kr_10
% 	P_Para(22)=1.74;
% % Parameter:   id =  kf_11, name = kf_11
% 	P_Para(23)=20.0;
% % Parameter:   id =  kr_11, name = kr_11
% 	P_Para(24)=0.81;
% % Parameter:   id =  kf_12, name = kf_12
% 	P_Para(25)=3.0E8;
% % Parameter:   id =  kr_12, name = kr_12
% 	P_Para(26)=4.0;
% % Parameter:   id =  kf_13, name = kf_13
% 	P_Para(27)=1.5E8;
% % Parameter:   id =  kr_13, name = kr_13
% 	P_Para(28)=8.0;
% % Parameter:   id =  kf_14, name = kf_14
% 	P_Para(29)=0.05;
% % Parameter:   id =  kr_14, name = kr_14
% 	P_Para(30)=0.0012;
% % Parameter:   id =  kf_15, name = kf_15
% 	P_Para(31)=0.05;
% % Parameter:   id =  kr_15, name = kr_15
% 	P_Para(32)=0.0012;
% % Parameter:   id =  kf_16, name = kf_16
% 	P_Para(33)=0.05;
% % Parameter:   id =  kr_16, name = kr_16
% 	P_Para(34)=0.0012;

% Reaction: id = React0, name = React0
	reaction_React0=compartment_comp1*(P_Para(1)*x(6)*x(13)-P_Para(2)*x(5));

% Reaction: id = React1, name = React1
	reaction_React1=compartment_comp1*(P_Para(3)*x(5)*x(13)-P_Para(4)*x(1));

% Reaction: id = React2, name = React2
	reaction_React2=compartment_comp1*(P_Para(5)*x(1)-P_Para(6)*x(12));

% Reaction: id = React3, name = React3
	reaction_React3=compartment_comp1*(P_Para(7)*x(4)*x(13)-P_Para(8)*x(3));

% Reaction: id = React4, name = React4
	reaction_React4=compartment_comp1*(P_Para(9)*x(3)*x(13)-P_Para(10)*x(12));

% Reaction: id = React5, name = React5
	reaction_React5=compartment_comp1*(P_Para(11)*x(6)-P_Para(12)*x(4));

% Reaction: id = React6, name = React6
	reaction_React6=compartment_comp1*(P_Para(13)*x(5)-P_Para(14)*x(3));

% Reaction: id = React7, name = React7
	reaction_React7=compartment_comp1*(P_Para(15)*x(11)*x(13)-P_Para(16)*x(2));

% Reaction: id = React8, name = React8
	reaction_React8=compartment_comp1*(P_Para(17)*x(2)*x(13)-P_Para(18)*x(9));

% Reaction: id = React9, name = React9
	reaction_React9=compartment_comp1*(P_Para(19)*x(4)-P_Para(20)*x(11));

% Reaction: id = React10, name = React10
	reaction_React10=compartment_comp1*(P_Para(21)*x(3)-P_Para(22)*x(2));

% Reaction: id = React11, name = React11
	reaction_React11=compartment_comp1*(P_Para(23)*x(12)-P_Para(24)*x(9));

% Reaction: id = React12, name = React12
	reaction_React12=compartment_comp1*(P_Para(25)*x(8)*x(13)-P_Para(26)*x(10));

% Reaction: id = React13, name = React13
	reaction_React13=compartment_comp1*(P_Para(27)*x(10)*x(13)-P_Para(28)*x(7));

% Reaction: id = React14, name = React14
	reaction_React14=compartment_comp1*(P_Para(29)*x(11)-P_Para(30)*x(8));

% Reaction: id = React15, name = React15
	reaction_React15=compartment_comp1*(P_Para(31)*x(2)-P_Para(32)*x(10));

% Reaction: id = React16, name = React16
	reaction_React16=compartment_comp1*(P_Para(33)*x(9)-P_Para(34)*x(7));

	xdot=zeros(13,1);
	
% Species:   id = BLL, name = BasalACh2, affected by kineticLaw
	xdot(1) = (1/(compartment_comp1))*(( 1.0 * reaction_React1) + (-1.0 * reaction_React2));
	
% Species:   id = IL, name = IntermediateACh, affected by kineticLaw
	xdot(2) = (1/(compartment_comp1))*(( 1.0 * reaction_React7) + (-1.0 * reaction_React8) + ( 1.0 * reaction_React10) + (-1.0 * reaction_React15));
	
% Species:   id = AL, name = ActiveACh, affected by kineticLaw
	xdot(3) = (1/(compartment_comp1))*(( 1.0 * reaction_React3) + (-1.0 * reaction_React4) + ( 1.0 * reaction_React6) + (-1.0 * reaction_React10));
	
% Species:   id = A, name = Active, affected by kineticLaw
	xdot(4) = (1/(compartment_comp1))*((-1.0 * reaction_React3) + ( 1.0 * reaction_React5) + (-1.0 * reaction_React9));
	
% Species:   id = BL, name = BasalACh, affected by kineticLaw
	xdot(5) = (1/(compartment_comp1))*(( 1.0 * reaction_React0) + (-1.0 * reaction_React1) + (-1.0 * reaction_React6));
	
% Species:   id = B, name = Basal, affected by kineticLaw
	xdot(6) = (1/(compartment_comp1))*((-1.0 * reaction_React0) + (-1.0 * reaction_React5));
	
% Species:   id = DLL, name = DesensitisedACh2, affected by kineticLaw
	xdot(7) = (1/(compartment_comp1))*(( 1.0 * reaction_React13) + ( 1.0 * reaction_React16));
	
% Species:   id = D, name = Desensitised, affected by kineticLaw
	xdot(8) = (1/(compartment_comp1))*((-1.0 * reaction_React12) + ( 1.0 * reaction_React14));
	
% Species:   id = ILL, name = IntermediateACh2, affected by kineticLaw
	xdot(9) = (1/(compartment_comp1))*(( 1.0 * reaction_React8) + ( 1.0 * reaction_React11) + (-1.0 * reaction_React16));
	
% Species:   id = DL, name = DesensitisedACh, affected by kineticLaw
	xdot(10) = (1/(compartment_comp1))*(( 1.0 * reaction_React12) + (-1.0 * reaction_React13) + ( 1.0 * reaction_React15));
	
% Species:   id = I, name = Intermediate, affected by kineticLaw
	xdot(11) = (1/(compartment_comp1))*((-1.0 * reaction_React7) + ( 1.0 * reaction_React9) + (-1.0 * reaction_React14));
	
% Species:   id = ALL, name = ActiveACh2, affected by kineticLaw
	xdot(12) = (1/(compartment_comp1))*(( 1.0 * reaction_React2) + ( 1.0 * reaction_React4) + (-1.0 * reaction_React11));
	
% Species:   id = L, name = ACh, affected by kineticLaw
	xdot(13) = (1/(compartment_comp1))*((-1.0 * reaction_React0) + (-1.0 * reaction_React1) + (-1.0 * reaction_React3) + (-1.0 * reaction_React4) + (-1.0 * reaction_React7) + (-1.0 * reaction_React8) + (-1.0 * reaction_React12) + (-1.0 * reaction_React13));
end