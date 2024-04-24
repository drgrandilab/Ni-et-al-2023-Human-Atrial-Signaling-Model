%% LTCC Mode 1
%% Input: Caj, Vm, najLCC, kjLCC, ISO_PKA;

function [ydot, LTCC_out] = LTCC_mode_1(state, Caj, Vm, najLCC, kjLCC, ISO_PKA, Frdy, FoRT,  Cao, Nao, Ko)


ISO_Shift = ISO_PKA * 3;
flag_7_state = 1;
Zca = 2;
% cAF = 0;
% ICa_scale = 1;
% aff = 1;
gammaCai = 0.0341;
gammaCao = 0.341;
Zk = 1;
gammaKi = 0.75;
gammaKo = 0.75;
Zna = 1;
gammaNai = 0.75;
gammaNao = 0.75;

% LTCC Current - Fixed Parameters
cpt = 3.75e-3;
% cat = 7.617e-3;
% s1o = 1.2 * 1.4 * 1.3 * 5.0 * 0.0182688;

k1o = 1.2 * 1.4 * 1.3 * 5.0 * 0.024168 ;




% k2o = 0.000103615;

% sp0 = 1.5;
sp1 = 3;
sp2 = 27;
sp3 = 3;
% sp4 = 4;
% sp5 = 7.1;
% sp6 = 15.6;
% sp7 = 10;
% sp8 = 4954;
% sp9 = 78.0329;
% sp10 = 0.1;
aR2 = 1;

sR2 = -9.5;
pR2 = 1.0 / 6.5;
% aT2 = 1;
% sT2 = -1000;
% pT2 = 0.100;
aR1 = 0.09091;
sR1 = -1000;
pR1 = 0.100;
aT1 = 0.30303;
sT1 = -1000;
pT1 = 0.100;
% aRv2 = 0.9;
% sRv2 = -29;
% pRv2 = 0.135;
% aTv2 = 500;
% sTv2 = -25;
% pTv2 = 0.050;
aRv1 = 0.85;
sRv1 = -27;  % -180;
pRv1 = 0.090;
% aTv1 = 270;
% sTv1 = -180;
% pTv1 = 0.100;
% aTrev1 = 205.12;
% sTrev1 = -65;
% pTrev1 = 0.100;
% aTrev2 = 7e8;
% sTrev2 = 60;
% pTrev2 = 0.130;
ICa_speed = 1;



fcp = 0.2 / (1 + (cpt / Caj)^1)  + 0.3 / (1 + (cpt / Caj)^2) + 0.5 / (1 + (cpt / Caj)^3);
R2 = aR2 / (1 + exp(-(Vm + ISO_Shift - sR2) * pR2));

T2 = 0.6 * 1.5 * (0.59 + (0.8 * exp(0.052 * (Vm + 13.0 + ISO_Shift))) / (1 + exp(0.132 * (Vm + 13.0 + ISO_Shift))));

PT = 1 - (1 / (1 + exp(-(Vm + ISO_Shift + sp2) / sp3)));
R1 = aR1 / (1 + exp(-(Vm + ISO_Shift - sR1) * pR1));
T1 = aT1 / (1 + exp(-(Vm + ISO_Shift - sT1) * pT1));

Rv1 = aRv1 / (1 + exp(-(Vm + ISO_Shift - sRv1) * pRv1));

alphaLCC = ICa_speed * R2 / T2;
betaLCC = ICa_speed * (1 - R2) / T2;
r1 = ICa_speed * R1 / T1;
r2 = ICa_speed * (1 - R1) / T1;


SLOPE = 8;
tau_k1k2 = 1.0 / (k1o * fcp);


	% dependent on PKA, Findlay, I. (2002). β-Adrenergic stimulation modulates Ca2+- and voltage-dependent inactivation of L-type Ca2+ channel currents in guinea-pig ventricular myocytes. 
	% J Physiol 541, 741–751. https://doi.org/10.1113/jphysiol.2002.019737.
	% Fig. 3 B vs A, ISO enhances CDI

Rv1_k1k2 = (0.75 +  ISO_PKA * 0.10) / (1 + exp((Vm + 25 - 14  + 5 * fcp + ISO_Shift) / SLOPE)) + (1 - (0.75 +  ISO_PKA * 0.10));

k1 = flag_7_state * ICa_speed * (1 - Rv1_k1k2) / tau_k1k2;
k2 = flag_7_state * ICa_speed * Rv1_k1k2  / tau_k1k2 ;

k3 = ICa_speed * PT / sp1;


% Recovery from inactivation was slowed by both depolarization and high [Ca2+]i: 
% Altamirano, J., and Bers, D.M. (2007). Effect of intracellular Ca2+ and action potential duration on L-type Ca2+ channel inactivation and recovery from inactivation in rabbit cardiac myocytes. 
% American Journal of Physiology - Heart and Circulatory Physiology 293, H563–H573. https://doi.org/10.1152/ajpheart.00469.2006.


Tv1_K5  = 20 + 0.45 * (1 - fcp) * 100 / (1 + exp((Vm + 40 + ISO_Shift) / -5)) ... 
         +  (0.8 + 3 * fcp) * ( 150 * exp(-(Vm + 45 + ISO_Shift) * (Vm + 45 + ISO_Shift) / 50) + 30.0 / (1 + exp((Vm + 40 + ISO_Shift) / 5))) ;

fcp_k5 =  fcp;
Is_Vss = (1.0 / (1 + exp( (Vm  + 25 - 14 + 5 * fcp + ISO_Shift) / SLOPE )));
k5 =  Is_Vss / Tv1_K5 / (fcp_k5 * 1.0 + 1.0);
k6 = (1 - Is_Vss) / Tv1_K5 * (fcp_k5 * 1 + 1.0) * 1;

k4 = k3 * (alphaLCC / betaLCC) * (k1 / k2) * (k5 / k6);

Tv1 =  56 + 3000 * exp(-(Vm + 50 + ISO_Shift) * (Vm + 50 + ISO_Shift) / 150) + 200 * 1.0 / (1 + exp((Vm + 60 + ISO_Shift) / -3)); 

Is_Vss_V = (1 / (1 + exp( (Vm + 27.0 + ISO_Shift) / SLOPE )) + 0);

k1p = ICa_speed * Rv1 / Tv1;
k2p = ICa_speed * (1 - Rv1) / Tv1;

Tv1_K5p = 56 + 1700 * exp(-(Vm + 37 + ISO_Shift) * (Vm + 37 + ISO_Shift) / 300) + 200 * 1.0 / (1 + exp((Vm + 60 + ISO_Shift) / -3));
k5p = Is_Vss_V / Tv1_K5p ;
k6p =  (1 - Is_Vss_V) / Tv1_K5p;

k3p = k3;
k4p = k3p * k5p * alphaLCC * k1p / (k2p * betaLCC * k6p);
s1p = k1p;

s2p = s1p * k2p * r1 / (r2 * k1p); 
% s2 = 0.01;

s1 = k1;
s2 = k2;
k13 = r2;
k14 =  k13 * k2 * r1 * s1 / ( s2 * r2 * k1);

Tv1_K9  = (50 + 30 * exp(-(Vm + 55 + ISO_Shift) * (Vm + 55 + ISO_Shift) / 150));
Rv1_k1k2 = 1 / (1 + exp((Vm + 27 + 5+10 - 14 + 5 * fcp + ISO_Shift) / SLOPE));

k7 = flag_7_state * ICa_speed * (1 - Rv1_k1k2) / Tv1_K9 ; % (50); % tau_k1k2;  % /(1+0.1*fcp)
k8 = flag_7_state * ICa_speed * Rv1_k1k2 / Tv1_K9 ;% (50); %  tau_k1k2 / (1);  % /(1+0.1*fcp)

Is_Vss = Rv1_k1k2;
k9 = 6 * Is_Vss / Tv1_K9 ;  % change to Caj dependent caused alternans /(1.0 + 3.0*fcp)
k10 = 6 * (1 - Is_Vss) / Tv1_K9;

k12 = k3;
k11 = k12 * (k9 * k4 * k7) / (k8 * k3 * k10 );
fcp_k15 = fcp;
k15_factor = 1;
tau_k15 = 150;
k15 = k15_factor * (1 - fcp_k15) / tau_k15;
k15p = k15_factor * fcp_k15 / tau_k15;
% tau_k16 = 150;

k16 = k15;
k16p = k15p;



Pc2_LCCj_m1   = state(1);
Pc1_LCCj_m1   = state(2);
Pi1Ca_LCCj_m1 = state(3);
Pi2Ca_LCCj_m1 = state(4);
Pi1Ba_LCCj_m1 = state(5);
Pi2Ba_LCCj_m1 = state(6);
Pi3Ca_LCCj_m1 = state(7);  
Pi4Ca_LCCj_m1 = state(8);   
PiOCa_LCCj_m1 = state(9);

Po_LCCj_m1 =  1 - Pc2_LCCj_m1 - Pc1_LCCj_m1 - Pi1Ca_LCCj_m1 - Pi2Ca_LCCj_m1 - Pi1Ba_LCCj_m1 - Pi2Ba_LCCj_m1 - Pi3Ca_LCCj_m1 - Pi4Ca_LCCj_m1 - PiOCa_LCCj_m1;

dPc2_LCCj_m1 = betaLCC * Pc1_LCCj_m1 + k5 * Pi2Ca_LCCj_m1 + k5p * Pi2Ba_LCCj_m1 - (k6 + k6p + alphaLCC) * Pc2_LCCj_m1;          %  C2_m1j
dPc1_LCCj_m1 = alphaLCC * Pc2_LCCj_m1 + k2 * Pi1Ca_LCCj_m1 + k2p * Pi1Ba_LCCj_m1 + r2 * Po_LCCj_m1 - (r1 + betaLCC + k1 + k1p) * Pc1_LCCj_m1; %  C1_m1j
dPi1Ca_LCCj_m1 = k1 * Pc1_LCCj_m1 + k4 * Pi2Ca_LCCj_m1 + k13 * PiOCa_LCCj_m1 + k8 * Pi3Ca_LCCj_m1 + k16p * Pi1Ba_LCCj_m1 - (k2 + k3 + k14 + k7 + k16) * Pi1Ca_LCCj_m1;             %  I1Ca_m1j
dPi2Ca_LCCj_m1 = k3 * Pi1Ca_LCCj_m1 + k6 * Pc2_LCCj_m1 + k9 * Pi4Ca_LCCj_m1 + k15p * Pi2Ba_LCCj_m1 - (k4 + k5 + k10  + k15) * Pi2Ca_LCCj_m1; %  I2Ca_m1j
dPi1Ba_LCCj_m1 = k1p * Pc1_LCCj_m1 + k4p * Pi2Ba_LCCj_m1 + s1p * Po_LCCj_m1 + k16 * Pi1Ca_LCCj_m1 - (k2p + k3p + s2p + k16p) * Pi1Ba_LCCj_m1;        %  I1Ba_m1j
dPi2Ba_LCCj_m1 = k3p * Pi1Ba_LCCj_m1 + k6p * Pc2_LCCj_m1 + k15 * Pi2Ca_LCCj_m1 - (k5p + k4p + k15p) * Pi2Ba_LCCj_m1;                                 %  I2Ba_m1j
dPi3Ca_LCCj_m1 = k7 * Pi1Ca_LCCj_m1 + k11 * Pi4Ca_LCCj_m1  - (k12 + k8 ) * Pi3Ca_LCCj_m1;
dPi4Ca_LCCj_m1 = k12 * Pi3Ca_LCCj_m1 + k10 * Pi2Ca_LCCj_m1 - (k11 + k9) * Pi4Ca_LCCj_m1;
dPiOCa_LCCj_m1 = s1 * Po_LCCj_m1 + k14 * Pi1Ca_LCCj_m1 - (s2 + k13) * PiOCa_LCCj_m1;

ydot = zeros(10, 1);
ydot(1) = dPc2_LCCj_m1;
ydot(2) = dPc1_LCCj_m1;
ydot(3) = dPi1Ca_LCCj_m1;
ydot(4) = dPi2Ca_LCCj_m1;
ydot(5) = dPi1Ba_LCCj_m1;
ydot(6) = dPi2Ba_LCCj_m1;
ydot(7) = dPi3Ca_LCCj_m1;
ydot(8) = dPi4Ca_LCCj_m1;
ydot(9) = dPiOCa_LCCj_m1;


Pna = 1.05* 0.65  * (0.70 - 0.20) * 100 * 0.675e-9; % [cm/s]
Pk = 1.05* 0.65  * (0.70 - 0.20) * 100 * 12.15e-9; % [cm/s]
Pca = 1.05* 0.65  * (0.70 - 0.20) * 100 * 24.3e-6; % [cm/s] 10*0.45*5.4e-6

ibarca_jm1 = Zca * Zca  * Pca * Frdy * FoRT * Vm / (exp(Zca * Vm * FoRT) - 1) * (gammaCai * Caj * exp(Zca * Vm * FoRT) - gammaCao * Cao);
LTCC_out.I_Ca_junc_m1 = (ibarca_jm1 * Po_LCCj_m1);

ibarna_jm1 = Zna * Zna * Pna * Frdy * FoRT * Vm / (exp(Zna * Vm * FoRT) - 1) * (gammaNai * najLCC * exp(Zna * Vm * FoRT) - gammaNao * Nao);
LTCC_out.I_Na_junc_m1 = (ibarna_jm1 * Po_LCCj_m1);

ibark_jm1 = Zk * Zk * Pk * Frdy * FoRT * Vm / (exp(Zk * Vm * FoRT) - 1) * (gammaKi * kjLCC * exp(Zk * Vm * FoRT) - gammaKo * Ko);
LTCC_out.I_K_junc_m1 = (ibark_jm1 * Po_LCCj_m1);


end


