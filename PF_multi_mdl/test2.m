
%% normrnd
a1=normrnd(2,2,[2,1000]);
a2=mean(a1(1,:));
a3=mean(a1(2,:));
figure(20),clf,
hist(a1(1,:),50);
%% mvnrnd
mu = [2,3];
sigma = [1,1.5;1.5,3];
% rng default  % For reproducibility
r = mvnrnd(mu,sigma,100);
figure(20),clf,
plot(r(:,1),r(:,2),'+')
%%

% Parameter:   id =  krbEGF, name = krbEGF
	P_Para(1)=2.18503E-5;
% Parameter:   id =  kruEGF, name = kruEGF
	P_Para(2)=0.0121008;
% Parameter:   id =  krbNGF, name = krbNGF
	P_Para(3)=1.38209E-7;
% Parameter:   id =  kruNGF, name = kruNGF
	P_Para(4)=0.00723811;
% Parameter:   id =  kEGF, name = kEGF
	P_Para(5)=694.731;
% Parameter:   id =  KmEGF, name = KmEGF
	P_Para(6)=6086070.0;
% Parameter:   id =  kNGF, name = kNGF
	P_Para(7)=389.428;
% Parameter:   id =  KmNGF, name = KmNGF
	P_Para(8)=2112.66;
% Parameter:   id =  kdSos, name = kdSos
	P_Para(9)=1611.97;
% Parameter:   id =  KmdSos, name = KmdSos
	P_Para(10)=896896.0;
% Parameter:   id =  kSos, name = kSos
	P_Para(11)=32.344;
% Parameter:   id =  KmSos, name = KmSos
	P_Para(12)=35954.3;
% Parameter:   id =  kRasGap, name = kRasGap
	P_Para(13)=1509.36;
% Parameter:   id =  KmRasGap, name = KmRasGap
	P_Para(14)=1432410.0;
% Parameter:   id =  kRasToRaf1, name = kRasToRaf1
	P_Para(15)=0.884096;
% Parameter:   id =  KmRasToRaf1, name = KmRasToRaf1
	P_Para(16)=62464.6;
% Parameter:   id =  kpRaf1, name = kpRaf1
	P_Para(17)=185.759;
% Parameter:   id =  KmpRaf1, name = KmpRaf1
	P_Para(18)=4768350.0;
% Parameter:   id =  kpBRaf, name = kpBRaf
	P_Para(19)=125.089;
% Parameter:   id =  KmpBRaf, name = KmpBRaf
	P_Para(20)=157948.0;
% Parameter:   id =  kdMek, name = kdMek
	P_Para(21)=2.83243;
% Parameter:   id =  KmdMek, name = KmdMek
	P_Para(22)=518753.0;
% Parameter:   id =  kpMekCytoplasmic, name = kpMekCytoplasmic
	P_Para(23)=9.85367;
% Parameter:   id =  KmpMekCytoplasmic, name = KmpMekCytoplasmic
	P_Para(24)=1007340.0;
% Parameter:   id =  kdErk, name = kdErk
	P_Para(25)=8.8912;
% Parameter:   id =  KmdErk, name = KmdErk
	P_Para(26)=3496490.0;
% Parameter:   id =  kpP90Rsk, name = kpP90Rsk
	P_Para(27)=0.0213697;
% Parameter:   id =  KmpP90Rsk, name = KmpP90Rsk
	P_Para(28)=763523.0;
% Parameter:   id =  kPI3K, name = kPI3K
	P_Para(29)=10.6737;
% Parameter:   id =  KmPI3K, name = KmPI3K
	P_Para(30)=184912.0;
% Parameter:   id =  kPI3KRas, name = kPI3KRas
	P_Para(31)=0.0771067;
% Parameter:   id =  KmPI3KRas, name = KmPI3KRas
	P_Para(32)=272056.0;
% Parameter:   id =  kAkt, name = kAkt
	P_Para(33)=0.0566279;
% Parameter:   id =  KmAkt, name = KmAkt
	P_Para(34)=653951.0;
% Parameter:   id =  kdRaf1ByAkt, name = kdRaf1ByAkt
	P_Para(35)=15.1212;
% Parameter:   id =  KmRaf1ByAkt, name = KmRaf1ByAkt
	P_Para(36)=119355.0;
% Parameter:   id =  kC3GNGF, name = kC3GNGF
	P_Para(37)=146.912;
% Parameter:   id =  KmC3GNGF, name = KmC3GNGF
	P_Para(38)=12876.2;
% Parameter:   id =  kC3G, name = kC3G
	P_Para(39)=1.40145;
% Parameter:   id =  KmC3G, name = KmC3G
	P_Para(40)=10965.6;
% Parameter:   id =  kRapGap, name = kRapGap
	P_Para(41)=27.265;
% Parameter:   id =  KmRapGap, name = KmRapGap
	P_Para(42)=295990.0;
% Parameter:   id =  kRap1ToBRaf, name = kRap1ToBRaf
	P_Para(43)=2.20995;
% Parameter:   id =  KmRap1ToBRaf, name = KmRap1ToBRaf
	P_Para(44)=1025460.0;
% Parameter:   id =  kdRaf1, name = kdRaf1
	P_Para(45)=0.126329;
% Parameter:   id =  KmdRaf1, name = KmdRaf1
	P_Para(46)=1061.71;
% Parameter:   id =  kdBRaf, name = kdBRaf
	P_Para(47)=441.287;
% Parameter:   id =  KmdBRaf, name = KmdBRaf
	P_Para(48)=1.08795E7;
%%
    
    
    
    
    
    