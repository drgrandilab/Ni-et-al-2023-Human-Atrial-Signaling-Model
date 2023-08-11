// #include "HAM_Cell.hpp"
#include <cmath>

#include "LTCC_Markov.hpp"



double LTCC_Markov::update_states_v9_mode_2(double Caj, double Vm, double najLCC, double kjLCC, double ISO_PKA) {
	double ISO_Shift = ISO_PKA * 3;
	// double ISO_Shift = 0;
	// Vm = Vm + ISO_PKA * 3;
	int flag_7_state = 1;
	// Caj=2e-3;
	double Zca = 2;
	int cAF = 0;
	double ICa_scale = 1;
	double aff = 1;
	double gammaCai = 0.0341;
	double gammaCao = 0.341;
	double Zk = 1;
	double gammaKi = 0.75;
	double gammaKo = 0.75;
	double Zna = 1;
	double gammaNai = 0.75;
	double gammaNao = 0.75;

	// LTCC Current - Fixed Parameters
	double cpt = 3.75e-3 ;// [mM]
	double cat = 7.617e-3;    // [mM]
	double s1o = 1.2*1.4*1.3 * 5.0 * 0.0182688; // [1/ms]  
	double k1o = 1.2*1.4*1.3 * 5.0 * 0.024168; // [1/ms]  
	double k2o = 0.000103615; // [1/ms]
	double sp0 = 1.5;
	double sp1 = 3;           // [ms]
	// sp2 = 40;          // [mV]
	double sp2 = 27;          // [mV]
	double sp3 = 3;           // [mV]
	double sp4 = 4;           // [mV]
	// sp5 = 11.32;       // [mV]
	double sp5 = 7.1 ;// based on Li et al. 1997. haibo //11.32;       // [mV]
	double sp6 = 15.6;        // [mV]
	double sp7 = 10;          // [ms]
	double sp8 = 4954;        // [ms]
	double sp9 = 78.0329;     // [ms]
	double sp10 = 0.1;        // [ms]
	double aR2 = 1;
	// sR2 = -2;          // [mV]

	double sR2 = -9.5;       // [mV] // left shift in steady-state activation
	double pR2 = 1.0 / 6.5;//6.5;//6.5;//6.5; //0.145;       // [1/mV]  // was 6.0
	double aT2 = 1;           // [1/ms]
	double sT2 = -1000;       // [mV]
	double pT2 = 0.100;       // [1/mV]
	double aR1 = 0.09091;
	double sR1 = -1000;       // [mV]
	double pR1 = 0.100;       // [1/mV]
	double aT1 = 0.30303;     // [1/ms]
	double sT1 = -1000;       // [mV]
	double pT1 = 0.100;       // [1/mV]
	double aRv2 = 0.9;
	double sRv2 = -29;        // [mV]
	double pRv2 = 0.135;      // [1/mV]
	double aTv2 = 500;        // [1/ms]
	double sTv2 = -25;        // [mV]
	double pTv2 = 0.050;      // [1/mV]
	double aRv1 = 0.85;
	double sRv1 = -27;//-180;       // [mV]
	double pRv1 = 0.090;      // [1/mV]
	double aTv1 = 270;        // [1/ms]
	double sTv1 = -180;       // [mV]
	double pTv1 = 0.100;      // [1/mV]
	double aTrev1 = 205.12;   // [1/ms]
	double sTrev1 = -65;      // [mV]
	double pTrev1 = 0.100;    // [1/mV]
	double aTrev2 = 7e8;      // [1/ms]
	double sTrev2 = 60;       // [mV]
	double pTrev2 = 0.130;    // [1/mV]
	int ICa_speed = 1;
	
	double fcp = 0.2 / (1 + pow(cpt / Caj, 1))  + 0.3 / (1 + pow(cpt / Caj, 2.0)) +0.5 / (1 + pow(cpt / Caj, 3)); 
	double R2 = aR2 / (1 + exp(-(Vm+ISO_Shift - sR2) * pR2));

	double T2 = 0.6 * 1.5 * (0.59 + (0.8 * exp(0.052 * (Vm + 13.0 + ISO_Shift))) / (1 + exp(0.132 * (Vm + 13.0 + ISO_Shift))));

	double PT = 1 - (1 / (1 + exp(-(Vm+ISO_Shift + sp2) / sp3)));
	double R1 = aR1 / (1 + exp(-(Vm+ISO_Shift - sR1) * pR1));
	double T1 = aT1 / (1 + exp(-(Vm+ISO_Shift - sT1) * pT1));
	double Rv1 = aRv1 / (1 + exp(-(Vm +ISO_Shift- sRv1) * pRv1));


	//  Transition Rates (20 rates)
	double alphaLCC = ICa_speed * R2 / T2;
	double betaLCC = ICa_speed * (1 - R2) / T2;
	r1 = ICa_speed * R1 / T1;
	r2 = ICa_speed * (1 - R1) / T1;
	
	// mode 2, // 10:58:32, Tue, 08-December-2020, By Haibo
	// T1 = 1.0;  // was 0.3030
	// R1 = 0.1;
	r1 = ICa_speed * 3. * R1 / T1;
	r2 = ICa_speed * (1 - 3.0 * R1) / T1; //  * (1.0 + 1.5 / (1 + exp((Vm + 50) / 1)));



	double SLOPE = 8;//7.1;
	//  k1 = flag_7_state*ICa_speed*k1o*fcp;
	double tau_k1k2 = 1.0 / (k1o * fcp );
	double Rv1_k1k2 = (0.75 +  ISO_PKA * 0.10) / (1 + exp((Vm + 25 - 14  + 5 * fcp + ISO_Shift) / SLOPE)) + (1 - (0.75 + ISO_PKA * 0.10));


	k1 = flag_7_state * ICa_speed * (1 - Rv1_k1k2) / tau_k1k2;
	k2 = flag_7_state * ICa_speed * Rv1_k1k2  / tau_k1k2 ;

	k3 = ICa_speed * PT / sp1;

	double Tv1_K5  = 20 + 0.45*(1 - fcp) * 100 / (1 + exp((Vm + 40 + ISO_Shift) / -5))
	                 +  (0.8 + 3 * fcp) * ( 150 * exp(-(Vm + 45 + ISO_Shift) * (Vm + 45 + ISO_Shift) / 50) + 30.0 / (1 + exp((Vm + 40 + ISO_Shift) / 5))) ; 
	
	double fcp_k5 =  fcp;
	double Is_Vss = (1.0 / (1 + exp( (Vm  + 25 - 14 + 5 * fcp + ISO_Shift) / SLOPE )));
	k5     =  Is_Vss / Tv1_K5 / (fcp_k5 * 1.0 + 1.0);


	k6      = (1 - Is_Vss) / Tv1_K5 * (fcp_k5 * 1 + 1.0) * 1;


	k4 = k3 * (alphaLCC / betaLCC) * (k1 / k2) * (k5 / k6); //  REV

	double Tv1  =  56 + 3000 * exp(-(Vm + 50 + ISO_Shift) * (Vm + 50 + ISO_Shift) / 150) + 200 * 1.0 / (1 + exp((Vm + 60 + ISO_Shift) / -3)); 

	double Is_Vss_V = (1 / (1 + exp( (Vm + 27.0 + ISO_Shift) / SLOPE )) + 0);


	k1p = ICa_speed * Rv1 / Tv1;
	k2p = ICa_speed * (1 - Rv1) / Tv1;
	
	double Tv1_K5p = 56 + 1700 * exp(-(Vm + 37 + ISO_Shift) * (Vm + 37 + ISO_Shift) / 300) + 200 * 1.0 / (1 + exp((Vm + 60 + ISO_Shift) / -3));
	k5p = Is_Vss_V / Tv1_K5p ;
	k6p =  (1 - Is_Vss_V) / Tv1_K5p;
	
	k3p = k3;
	k4p = k3p * k5p * alphaLCC * k1p / (k2p * betaLCC * k6p);
	s1p = k1p;
	s2p = s1p * k2p * r1 / (r2 * k1p); 

	s1 = flag_7_state * ICa_speed * k1o * fcp ;
	s2 = 0.01;

	
	s1 = k1;
	s2 = k2;

	k13 = r2;
	k14 =  k13 * k2 * r1 * s1 / ( s2 * r2 * k1);

	//  adding new states here
	double Tv1_K9  = (50 + 30 * exp(-(Vm + 55 + ISO_Shift) * (Vm + 55 + ISO_Shift) / 150));
	Rv1_k1k2 = 1 / (1 + exp((Vm + 27 + 5+10 - 14 + 5 * fcp + ISO_Shift) / SLOPE));

	k7 = flag_7_state * ICa_speed * (1 - Rv1_k1k2) / Tv1_K9 ;
	k8 = flag_7_state * ICa_speed * Rv1_k1k2 / Tv1_K9 ;


	Is_Vss = Rv1_k1k2;
	k9 = 6 * Is_Vss / Tv1_K9 ;  
	k10 = 6 * (1 - Is_Vss) / Tv1_K9;



	k12 = k3;
	k11 = k12 * (k9 * k4 * k7) / (k8 * k3 * k10 );
	double fcp_k15 = fcp;
	double k15_factor = 1;
	double tau_k15 = 150;
	double k15 = k15_factor * (1 - fcp_k15) / tau_k15;
	double k15p = k15_factor * fcp_k15 / tau_k15;
	double tau_k16 = 150;

	double k16 = k15;
	double k16p = k15p;

	// global ICaL_inac_tau ICaL_rates
	// ICaL_inac_tau(1,tStep) = 1.0 /(s1+s2);
	// ICaL_inac_tau(2,tStep) = 1.0 /(s1p+s2p);
	// ICaL_inac_tau(3,tStep) = 1.0 /(k1+k2);
	// ICaL_inac_tau(4,tStep) = 1.0 /(k1p+k2p);
	// ICaL_inac_tau(5,tStep) = 1.0 /(s1p+s2p);

	// ICaL_rates(1,tStep) = s1;
	// ICaL_rates(2,tStep) = s2;
	// ICaL_rates(3,tStep) = s1p;
	// ICaL_rates(4,tStep) = s2p;
	// ICaL_rates(5,tStep) = k1;
	// ICaL_rates(6,tStep) = k2;
	// ICaL_rates(7,tStep) = k1p;
	// ICaL_rates(8,tStep) = k2p;
	// ICaL_rates(9,tStep) = k3;
	// ICaL_rates(10,tStep) = k4;
	// ICaL_rates(11,tStep) = k6;
	// ICaL_rates(12,tStep) = k5;
	// ICaL_rates(13,tStep) = r1;
	// ICaL_rates(14,tStep) = r2;
	double Pc2_LCCj_m1   = state[0];
	double Pc1_LCCj_m1   = state[1];
	double Pi1Ca_LCCj_m1 = state[2];
	double Pi2Ca_LCCj_m1 = state[3];
	double Pi1Ba_LCCj_m1 = state[4];
	double Pi2Ba_LCCj_m1 = state[5];
	double Pi3Ca_LCCj_m1 = state[6];  
	double Pi4Ca_LCCj_m1 = state[7];   
	double PiOCa_LCCj_m1 = state[8];   


	double Po_LCCj_m1 =  1 - Pc2_LCCj_m1 - Pc1_LCCj_m1 - Pi1Ca_LCCj_m1 - Pi2Ca_LCCj_m1 - Pi1Ba_LCCj_m1 - Pi2Ba_LCCj_m1 - Pi3Ca_LCCj_m1 - Pi4Ca_LCCj_m1 - PiOCa_LCCj_m1;

	//  State transitions for mode-1 junctional LCCs
	double dPc2_LCCj_m1 = betaLCC * Pc1_LCCj_m1 + k5 * Pi2Ca_LCCj_m1 + k5p * Pi2Ba_LCCj_m1 - (k6 + k6p + alphaLCC) * Pc2_LCCj_m1;          //  C2_m1j
	double dPc1_LCCj_m1 = alphaLCC * Pc2_LCCj_m1 + k2 * Pi1Ca_LCCj_m1 + k2p * Pi1Ba_LCCj_m1 + r2 * Po_LCCj_m1 - (r1 + betaLCC + k1 + k1p) * Pc1_LCCj_m1; //  C1_m1j	
	double dPi1Ca_LCCj_m1 = k1 * Pc1_LCCj_m1 + k4 * Pi2Ca_LCCj_m1 + k13 * PiOCa_LCCj_m1 + k8 * Pi3Ca_LCCj_m1 + k16p * Pi1Ba_LCCj_m1 - (k2 + k3 + k14 + k7 + k16) * Pi1Ca_LCCj_m1;             //  I1Ca_m1j
	double dPi2Ca_LCCj_m1 = k3 * Pi1Ca_LCCj_m1 + k6 * Pc2_LCCj_m1 + k9 * Pi4Ca_LCCj_m1 + k15p * Pi2Ba_LCCj_m1 - (k4 + k5 + k10  + k15) * Pi2Ca_LCCj_m1; //  I2Ca_m1j
	double dPi1Ba_LCCj_m1 = k1p * Pc1_LCCj_m1 + k4p * Pi2Ba_LCCj_m1 + s1p * Po_LCCj_m1 + k16 * Pi1Ca_LCCj_m1 - (k2p + k3p + s2p + k16p) * Pi1Ba_LCCj_m1;        //  I1Ba_m1j
	double dPi2Ba_LCCj_m1 = k3p * Pi1Ba_LCCj_m1 + k6p * Pc2_LCCj_m1 + k15 * Pi2Ca_LCCj_m1 - (k5p + k4p + k15p) * Pi2Ba_LCCj_m1;                                 //  I2Ba_m1j
	double dPi3Ca_LCCj_m1 = k7 * Pi1Ca_LCCj_m1 + k11 * Pi4Ca_LCCj_m1  - (k12 + k8 ) * Pi3Ca_LCCj_m1;
	// double dPi4Ca_LCCj_m1 = k3 * Pi3Ca_LCCj_m1 + k10 * Pi2Ca_LCCj_m1 - (k4 + k9) * Pi4Ca_LCCj_m1;
	double dPi4Ca_LCCj_m1 = k12 * Pi3Ca_LCCj_m1 + k10 * Pi2Ca_LCCj_m1 - (k11 + k9) * Pi4Ca_LCCj_m1;
	double dPiOCa_LCCj_m1 = s1 * Po_LCCj_m1 + k14 * Pi1Ca_LCCj_m1 - (s2 + k13) * PiOCa_LCCj_m1;



	ydot[0] = dPc2_LCCj_m1;
	ydot[1] = dPc1_LCCj_m1;
	ydot[2] = dPi1Ca_LCCj_m1;
	ydot[3] = dPi2Ca_LCCj_m1;
	ydot[4] = dPi1Ba_LCCj_m1;
	ydot[5] = dPi2Ba_LCCj_m1;
	ydot[6] = dPi3Ca_LCCj_m1;
	ydot[7] = dPi4Ca_LCCj_m1;
	ydot[8] = dPiOCa_LCCj_m1;
	double ibarca_jm1 = Zca * Zca  * Pca * Frdy * FoRT * Vm / (exp(Zca * Vm * FoRT) - 1) * (gammaCai * Caj * exp(Zca * Vm * FoRT) - gammaCao * Cao);
	I_Ca_junc_m1 = Fjunc_CaL * (ibarca_jm1 * Po_LCCj_m1) ;

	double ibarna_jm1 = Zna * Zna * Pna * Frdy * FoRT * Vm / (exp(Zna * Vm * FoRT) - 1) * (gammaNai * najLCC * exp(Zna * Vm * FoRT) - gammaNao * Nao);
	I_Na_junc_m1 = Fjunc_CaL * (ibarna_jm1 * Po_LCCj_m1) ;

	double ibark_jm1 = Zk * Zk * Pk * Frdy * FoRT * Vm / (exp(Zk * Vm * FoRT) - 1) * (gammaKi * kjLCC * exp(Zk * Vm * FoRT) - gammaKo * Ko);
	I_K_junc_m1 = Fjunc_CaL * (ibark_jm1 * Po_LCCj_m1) ;
}

