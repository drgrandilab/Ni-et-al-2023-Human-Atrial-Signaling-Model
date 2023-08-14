#include <math.h>
#include "HAM_ECC.hpp"

// update to include signalling-dependent effects
// 14:27:49, Tue, 05-March-2019, By Haibo

#pragma optimize("", off)
void HAM_ECC::HAM_ECC_update_ODE(double t, signalling_para & para) {
	// This file, which describes the human atrial excitation-contraction
	// coupling, was built upon the code developed by the Grandi et al.
	// Reference: Grandi E, Pandit SV, Voigt N, Workman AJ, Dobrev D, Jalife J,
	// Bers DM. (2011). Human Atrial Action Potential and Ca2+ Model: Sinus
	// Rhythm and Chronic Atrial Fibrillation. Circ Res. 2011 Oct 14;109(9):
	// 1055-66. Epub 2011 Sep 15.

	// MOD1 is for K1 (MOD1: GNa, GNaL, IClbkg) - anything but 3
	// MOD2 is for K2
	// MOD3 is for K3 - 3
	int MOD_ind = 3; // <- now implemented!

	// Optimization coeffcients (MOD_ind = 3)
	//parameters_min = ones(1,25);
	// parameters_min_opt15 = [0.999577536811449,0.980825078550174,1.09603000184601,0.924655894922876,1.16521268185567,1.81090408596575,0.367634769993851,0.357030670487448,1.78217332117361,0.979209273540801,1.60170059166390,1,0.334590150555472,1,1.12891867706223,1.28106061951874,1.05775442816844,1.14616657989970,0.701149131592228,0.973434224380471,0.704636234456055,0.710072856800229,1.28704962772999,1.35326407799692,1.43817543727274];
	// parameters_min_opt16 = parameters_min_opt15.*[1,1,1,1,0.931388465326278,1.24121222425510,0.697723678281066,0.959203070619648,0.698887199522057,1.21790637146654,0.996742520272117,1,1.32714104391583,1,1,1,1,1,1,1,1,1,1,1.12963540145110,0.901842104741936];
	std::vector<double> parameters_min_opt17 {0.999763767766896, 1.02704368132844, 0.904232756929961, 0.912089674296620, 1.06595749235083, 1.00817013225314, 0.561056138848061, 0.354046530214812, 1.35707583730889, 2.84563198875882, 1.52172787019912, 1, 0.569127790357487, 1, 1.06023886956084, 1.28931577624649, 1.02096444680901, 0.956205050080743, 0.814370961024797, 0.854520304154750, 0.701949192397789, 0.833439913899504, 1.29893009653586, 1.42202828590400, 0.951213166214512};
	//parameters_min = parameters_min_opt15;
	//parameters_min = parameters_min_opt16;
	// parameters_min = parameters_min_opt17;

	//// State variables
	// 1       2       3       4       5       6       7       8      9       10
	// m       h       j       d       f       fcaBj   fcaBsl  xkur   ykur    xtof
	// 11      12      13      14      15      16      17      18     19      20
	// ytof    xkr     xks     RyRr    RyRo    RyRi    NaBj    NaBsl  TnCL    TnCHc
	// 21      22      23      24      25      26      27      28     29      30
	// TnCHm   CaM     Myoc    Myom    SRB     SLLj    SLLsl   SLHj   SLHsl   Csqnb
	// 31      32      33      34      35      36      37      38     39      40
	// Ca_sr   Naj     Nasl    Nai     Ki      Caj     Casl    Cai    Vm      m(late)
	// 41      42
	// h(late) xkp2

	// ydot = zeros(size(y));



	//// cell_par Parameters
	//epi        = 1; // EPI or ENDO?
	int RA       = 0; // Right ATRIUM
	int AF       = 0; // p[5]; // AF
	AF = para.AF;
	int ISO      = 0;//p[6]; // ISO (boolean)
	double CCh   = 0;//p(7); // [uM]
	double ACh   = 0;//p(8); // [uM]
	int K2P_cond = 1;//p(9);
	int SK_cond  = 1;//p(10);
	int SK_shift = 0;// p(11);
	// int Ca_clamp = 0;//p(12);
	// int Na_clamp = 0; //p(13);
	// Scaling factors for sensitivity analysis
	// SA_par = p(14:end);
	// if length(SA_par) == 23
	//     SA_par = [SA_par 1 1];
	// end
	// Optimization
	// SA_par = SA_par.*parameters_min;


	std::vector<double> SA_par(parameters_min_opt17.begin(), parameters_min_opt17.end());  // set SA_Par parameters

	// Nernst Potentials
	double ena_junc = (1.0 / FoRT) * log(Nao / y[32 - 1]); // [mV]
	double ena_sl   = (1.0 / FoRT) * log(Nao / y[33 - 1]); // [mV]
	double ek       = (1.0 / FoRT) * log(Ko / y[35 - 1]);	 // [mV]
	double eca_junc = (1.0 / FoRT / 2) * log(Cao / y[35]); // [mV2
	double eca_sl   = (1.0 / FoRT / 2) * log(Cao / y[37 - 1]); // [mV]
	double ecl      = (1.0 / FoRT) * log(Cli / Clo);      // [mV]

	//// I_Na
	// // Grandi et al, 2011
	// GNa = SA_par(1)*23*(1-0.1*AF);  // [mS/uF] Grandi
	//
	// mss = 1 / ((1 + exp( -(56.86 + y[39-1]) / 9.03 ))^2);
	// taum = 0.1292 * exp(-((y[39-1]+45.79)/15.54)^2) + 0.06487 * exp(-((y[39-1]-4.823)/51.12)^2);
	//
	// ah = (y[39-1] >= -40) * (0)...
	//    + (y[39-1] < -40) * (0.057 * exp( -(y[39-1] + 80) / 6.8 ));
	// bh = (y[39-1] >= -40) * (0.77 / (0.13*(1 + exp( -(y[39-1] + 10.66) / 11.1 )))) ...
	//    + (y[39-1] < -40) * ((2.7 * exp( 0.079 * y[39-1]) + 3.1e5 * exp(0.3485 * y[39-1])));
	// tauh = 1 / (ah + bh);
	// hss = 1 / ((1 + exp( (y[39-1] + 71.55)/7.43 ))^2);
	//
	// aj = (y[39-1] >= -40) * (0) ...
	//     +(y[39-1] < -40) * (((-2.5428 * 10^4*exp(0.2444*y[39-1]) - 6.948e-6 * exp(-0.04391*y[39-1])) * (y[39-1] + 37.78)) / ...
	//                      (1 + exp( 0.311 * (y[39-1] + 79.23) )));
	// bj = (y[39-1] >= -40) * ((0.6 * exp( 0.057 * y[39-1])) / (1 + exp( -0.1 * (y[39-1] + 32) ))) ...
	//    + (y[39-1] < -40) * ((0.02424 * exp( -0.01052 * y[39-1] )) / (1 + exp( -0.1378 * (y[39-1] + 40.14) )));
	// tauj = 1 / (aj + bj);
	// jss = 1 / ((1 + exp( (y[39-1] + 71.55)/7.43 ))^2);
	//
	// ydot(1) = (mss - y[0]) / taum;
	// ydot(2) = (hss - y[2-1]) / tauh;
	// ydot(3) = (jss - y[2]) / tauj;



	// para.kCKII_Nav_change has not been implemented yet. // 16:53:03, Tue, 05-March-2019, By Haibo
	// Courtemanche et al, 1998
	//GNa = SA_par(1)*23*(1-0.1*AF);  // [mS/uF] Grandi
	//GNa = SA_par(1)*7.8*(1-0.1*AF);  // [mS/uF] Courtemanche
	double Ina_h_shift = /*0.5 * */(3.25 - 2 * 3.25 * 1 / (1 + exp(-(para.NaV_CKp - 0.12) / 0.1))); // h gate dependent on INa camKII

	if (para.No_CaMKII_NaV) {
		Ina_h_shift = 0.3239;  //based on CaKII p of NaV at baseline 1 Hz // NaV_CKp = 0.1 /*0.5 * */(3.25 - 2 * 3.25 * 1 / (1 + exp(-(0.1 - 0.12) / 0.1)));
	}
	// double Ina_h_shift = /*0.5 * */(5 - 2 * 5 * 1 / (1 + exp(-(para.NaV_CKp - 0.12) / 0.1))); // h gate dependent on INa camKII
	double GNa = (0.25 * para.kPKA_Ina_increase + 1) * (1 - 0.1 * AF) * SA_par[0] * 9/*9*/ ; // [mS/uF] Courtemanche (adjusted for dV/dt max) MOD1

	double am = (y[39 - 1] == -47.13) ? 3.2 : (0.32 * (y[39 - 1] + 47.13) / (1 - exp( -0.1 * (y[39 - 1] + 47.13))));

	double bm = 0.08 * exp(-y[39 - 1] / 11.0);

	double Vm = y[38] - Ina_h_shift;

	double ah = (Vm >= -40) ? 0 : 0.135 * exp( -(Vm + 80) / 6.8 );

	double bh = (Vm >= -40) ? (1.0 / (0.13 * (1 + exp( -(Vm + 10.66) / 11.1 )))) : (3.56 * exp( 0.079 * Vm) + 3.1e5 * exp(0.35 * Vm));

	double aj = (Vm >= -40) ? 0 : ( (-127140 * exp(0.2444 * Vm) - 3.474e-5 * exp(-0.04391 * Vm)) * (Vm + 37.78)) / (1 + exp( 0.311 * (Vm + 79.23) ) );
	double bj = (Vm >= -40) ? ((0.3 * exp(-2.535e-7 * Vm)) / (1 + exp( -0.1 * (Vm + 32) ))) : (0.1212 * exp( -0.01052 * Vm)) / (1 + exp( -0.1378 * (Vm + 40.14) ));

	ydot[0] = am * (1 - y[0]) - bm * y[0];
	ydot[1] = ah * (1 - y[1]) - bh * y[1];
	ydot[2] = aj * (1 - y[2]) - bj * y[2];

	I_Na_junc = para.INa_Scale * Fjunc * GNa * y[0] * y[0] * y[0] * y[1] * y[2] * (y[39 - 1] - ena_junc);
	I_Na_sl = para.INa_Scale * Fsl * GNa * y[0] * y[0] * y[0] * y[1] * y[2] * (y[39 - 1] - ena_sl);
	I_Na = I_Na_junc + I_Na_sl; //#ok<NASGU>

	//// Late I_Na
	//GNaL = SA_par(2)*0.0025*AF; // [mS/uF]
	// (current)  //(1.0 + 0.5 / (1 + exp((para.kCKII_Nav_change - 0.3) / -0.05)))
	double GNaL = para.INaL_Scale * ( (14.0 / 11.0) / (1 + exp(-(para.NaV_CKp - 0.12) / 0.1)) ) * SA_par[1] * (1.2 * 0.003 + 0.0015 * AF) * (1 - 0.1 * AF) * SA_par[0] * 9; // [mS/uF] (adjusted) MOD1  // increased GNaL here // 18:29:43, Thu, 07-February-2019, By Haibo

	am = 0.32 * (y[39 - 1] + 47.13) / (1 - exp(-0.1 * (y[39 - 1] + 47.13)));
	bm = 0.08 * exp(-y[39 - 1] / 11);
	double inf = 1 / (1 + exp((y[39 - 1] + 91) / 6.1));
	double tau = 200;
	ydot[39] = am * (1 - y[39]) - bm * y[39];
	ydot[40] = (inf - y[40]) / tau;

	I_NaL_junc =  Fjunc * GNaL * y[39] * y[39] * y[39] * y[40] * (y[38] - ena_junc);
	I_NaL_sl =  Fsl * GNaL * y[39] * y[39] * y[39] * y[40] * (y[38] - ena_sl);
	I_NaL = I_NaL_junc + I_NaL_sl; //#ok<NASGU>

	//// I_nabk: Na Background Current
	double GNaB = para.INab_Scale * /*SA_par[2] * */0.597e-3; // [mS/uF]

	I_nabk_junc = Fjunc * GNaB * (y[39 - 1] - ena_junc);
	I_nabk_sl = Fsl * GNaB * (y[39 - 1] - ena_sl);
	I_nabk = I_nabk_junc + I_nabk_sl; //#ok<NASGU>

	//// I_nak: Na/K Pump Current
	double IbarNaK = /*0.9 **/ para.INaK_Scale * 1.0 * /*SA_par[3] **/ 1.26; // [uA/uF]
	//IbarNaK = SA_par(4)*1.25*(1.26); // NEW K2P +25// in Schmidt et al 2015
	double KmNaip = 11; // 11 original; 10 - from CRN model; 11.0;/** (1 - 0.25 * ISO)*/; // [mM]
	KmNaip =  KmNaip - para.KmNaip_PKA;
	double KmKo = 1.5; // [mM]
	//Q10NaK = 1.63;
	//Q10KmNai = 1.39;

	double sigma = (exp(Nao / 67.3) - 1) / 7.0;
	double fnak = 1 / (1 + 0.1245 * exp(-0.1 * y[39 - 1] * FoRT) + 0.0365 * sigma * exp(-y[39 - 1] * FoRT));

	I_nak_junc = Fjunc * IbarNaK * fnak * Ko / (1 + pow((KmNaip / y[32 - 1]), 4)) / (Ko + KmKo);
	I_nak_sl = Fsl * IbarNaK * fnak * Ko / (1 + pow((KmNaip / y[33 - 1]), 4)) / (Ko + KmKo);
	I_nak = I_nak_junc + I_nak_sl;

	//// I_to: Transient Outward K Current
	double GtoFast = /*para.kCKII_Itof_G **//* 0.8 * */ 1.2 * 0.75 * para.Itof_Scale * 1.2 * (1 + (0.6 - 1) * para.kPKA_Itof) * SA_par[4] * (1 - 0.7 * AF) * 0.165; // 40% decrease w/ 100 nM ISO (Gonzales De La Fuente et al 2013)
	//GtoFast = SA_par(5)*0.4*(1-0.7*AF)*0.165; // NEW K2P -60// in Schmidt et al 2015


	Vm = y[38] - para.kCKII_Itof_vshift;
	double xtoss = 1 / ( 1 + exp( -(Vm + 1) / 11 ) );
	double tauxtof = 3.5 * exp(-((Vm / 30) * (Vm / 30))) + 1.5;

	Vm = y[38];// - para.kCKII_Itof_vshift;  // no shift in ytoss
	double ytoss = 1 / ( 1 + exp( (Vm + 40.5) / 11.5) ) ;
	double tauytof = (para.kCKII_Itof_tau - 1.0) * 20 / (1.0 + exp( (Vm + 40) / -5.0)) + // slows inactivation
	                 25.635 * exp(-(((Vm + 52.45) / 15.8827) * ((Vm + 52.45) / 15.8827))) + 24.14
	                 - (para.kCKII_Itof_tau - 1.0) * 20 / (1.0 + exp( (Vm + 50) / 5.0)); //14.14    // harsens recovery
	// better to scale 24.14 as well? // 13:34:33, Thu, 07-March-2019, By Haibo
	// Tessier et al. 1999, CaMKII slows inactivation
	// Wagner et al. 2009, CaMKII accerates recovery
	ydot[9] = (xtoss - y[9]) / tauxtof;
	ydot[10] = (ytoss - y[10]) / tauytof;

	I_tof = GtoFast * y[9] * y[10] * (y[39 - 1] - ek);
	I_to = I_tof;

	//// I_kr: Rapidly Activating K Current
	double gkr = para.IKr_Scale * (1 + para.dGkr_PKA) * SA_par[5] * 0.95 * 0.035 * sqrt(Ko / 5.4);

	inf = 1 / (1 + exp(-(y[39 - 1] + 10 + para.Vkr_PKA) / 5.0));
	tau = 550 / (1 + exp((-22 - y[39 - 1]) / 9)) * 6 / (1 + exp((y[39 - 1] - (-11.0)) / 9.0)) + 230 / (1 + exp((y[39 - 1] - (-40)) / 20));
	ydot[11] = (inf - y[11]) / tau;
	double rkr = 1 / (1 + exp((y[39 - 1] + 74.0) / 24.0));
	// double rkr = 1 / (1 + exp((y[39 - 1] + 15.0) / 24.0)); // from CRN model
	I_kr = 0.9 * gkr * y[11] * rkr * (y[39 - 1] - ek);


	for (int i = 0; i < 7; ++i)
	{
		IKr_markov.state[i] = y[85 + i];
	}
	IKr_markov.update_states(y[39 - 1], 0.0, ek, para);

	for (int i = 0; i < 7; ++i)
	{
		ydot[85 + i] = IKr_markov.ydot[i];
	}

	I_kr = /*0.95 **/ 1.1 * 0.9 * para.IKr_Scale * 0.35 * IKr_markov.IKr; // times 0.35 to match IKr magnitude in the original Grandi model.
	//// I_ks: Slowly Activating K Current
	//  old version


	double pNaK = 0.01833;
	double eks = (1 / FoRT) * log((Ko + pNaK * Nao) / (y[35 - 1] + pNaK * y[34 - 1]));
#ifdef OLTIKS
	double gks_junc = 1 * SA_par[6] * (1 + 1 * AF + 0.6 * para.kPKA_Iks) * 0.0035;
	double gks_sl = 1 * SA_par[6] * (1 + 1 * AF + 0.6 * para.kPKA_Iks) * 0.0035;

	inf = 1 / (1 + exp(-(y[39 - 1] + 40 * ISO + 3.8) / 14.25));
	tau = 990.1 / (1 + exp(-(y[39 - 1] + 40 * ISO + 2.436) / 14.12));
	ydot[12] = (inf - y[12]) / tau;

	I_ks_junc = para.IKs_Scale * Fjunc * gks_junc * y[12] * y[12] * (y[39 - 1] - eks);
	I_ks_sl = para.IKs_Scale * Fsl * gks_sl * y[12] * y[12] * (y[39 - 1] - eks);
#endif


	// new version from Bartos et al. 2017
	double gks_factor = /*0.95 **/ 1.5 * 8 * 0.0035;
	double P_g_0 = gks_factor * (0.2 + 0.2 * para.kPKA_Iks); // 0.2->0.4 w/ ISO
	double P_g_max = gks_factor * (0.8 + 0.5 * para.kPKA_Iks); // 0.8->1.3 w/ ISO
	double P_vh_0 = -1 - 10 * para.kPKA_Iks; // -1->-11 w/ ISO
	double P_vh_max = -12 - 9 * para.kPKA_Iks; // -12->-21 w/ ISO
	double P_tau_0 = 26 + 9 * para.kPKA_Iks; // 160->260 w/ ISO
	double P_tau_max = 40 + 4 * para.kPKA_Iks; // 300->370 w/ ISO
	double caks_junc = y[36 - 1];
	double caks_sl = y[37 - 1]; // normal simulation
	double gks_junc = P_g_0 + (P_g_max - P_g_0) / (1 + pow(150e-6 / caks_junc, 1.3)); // Regulated by PKA
	double gks_sl = P_g_0 + (P_g_max - P_g_0) / (1 + pow(150e-6 / caks_sl, 1.3)); // Regulated by PKA
	double VsXs_Ca_junc = P_vh_0 + (P_vh_max - P_vh_0) / (1 + pow(350e-6 / caks_junc, 4.0)); // Regulated by PKA
	double xsss_junc = 1.0 / (1 + exp(-(y[38] - VsXs_Ca_junc) / 25.0));
	double VsTs_Ca_junc = P_tau_0 + (P_tau_max - P_tau_0) / (1 + pow(150e-6 / caks_junc, 3)); // Regulated by PKA
	double tauxs_junc = 2 * (50 + (50 + 350 * exp(-(pow((y[38] + 30), 2)) / 4000.0)) * 1.0 / (1 + exp(-(y[38] + VsTs_Ca_junc) / 10.0)));
	double VsXs_Ca_sl = P_vh_0 + (P_vh_max - P_vh_0) / (1 + pow(350e-6 / caks_sl, 4.0) ); // Regulated by PKA
	double xsss_sl = 1 / (1 + exp(-(y[38] - VsXs_Ca_sl) / 25.0));
	double VsTs_Ca_sl = P_tau_0 + (P_tau_max - P_tau_0) / (1 + pow(150e-6 / caks_sl, 3)); // Regulated by PKA
	double tauxs_sl = 2 * (50 + (50 + 350 * exp(-(pow((y[38] + 30) , 2)) / 4000.0)) * 1 / (1 + exp(-(y[38] + VsTs_Ca_sl) / 10.0)));

	// note that we are using the states orginally for LTCC, which has been replaced with Markov model
	ydot[3] = (xsss_junc - y[3]) / tauxs_junc;
	ydot[4] = (xsss_sl - y[4]) / tauxs_sl;
	// eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)));
	I_ks_junc =  para.IKs_Scale * Fjunc * gks_junc * y[3] * y[3] * (y[38] - eks);
	I_ks_sl = para.IKs_Scale *  Fsl * gks_sl * y[4] * y[4] * (y[38] - eks);


	I_ks = I_ks_junc + I_ks_sl;

	//// I_kur: Ultra-Rapid Delayed Rectifier Outward K Current
	double Gkur =   1.2 * para.IKur_Scale * para.kCKII_Ikur * (1.0 + 2 * para.kPKA_Ikur) * 1.8 * SA_par[7] * (1 - 0.5 * AF) * 0.045 * (1 + 0.2 * RA); //1.5* by haibo
	//Gkur = SA_par(8)*(1+0.2*AF)*(1-0.5*AF)*(1+2*ISO)*0.045*(1+0.2*RA); // NEW K2P +20// cAF in Schmidt et al 2015

	double xkurss = 1 / ( 1 + exp( (y[39 - 1] + 6) / -8.6) );
	double tauxkur = 9.0 / (1 + exp((y[39 - 1] + 5) / 12.0)) + 0.5;
	double ykurss = 1 / ( 1 + exp( (y[39 - 1] + 7.5) / 10) );
	double tauykur = 590.0 / (1 + exp((y[39 - 1] + 60) / 10)) + 3050.0;

	ydot[7] = (xkurss - y[7]) / tauxkur;
	ydot[8] = (ykurss - y[8]) / tauykur;
	I_kur = Gkur * y[7] * y[8] * (y[39 - 1] - ek);

	// NEw IKur model here; // 16:37:09, Mon, 25-February-2019, By Haibo
	double CNZ_gkur = 1.05 * para.IKur_Scale * para.kCKII_Ikur * (1.0 + /*0.2*/1 * para.kPKA_Ikur) * (1 - 0.5 * AF) * 1.07 * 1.2 * 0.006398 * 0.85;
	double K_Q10 = 3.0;//3.5308257834747638;
	//y[42]
	// inf = ((het.IKur_ac1_mult * 1.0) / (1 + exp((V + 17.6684 + het.IKur_ac1_shift) / (-5.75418 * het.IKur_ac1_grad))) ) * ((het.IKur_ac2_mult * 1.0) / (1 + exp((V + 8.4153 + het.IKur_ac2_shift) / (-11.51037561 * het.IKur_ac2_grad)))) + het.IKur_ac_add;

	// inf = 1.0 / (1 + exp(-(y[38] + 5.52) / (8.6))); // change to simple version, haibo.  removed IKur mutations here.
	inf = 1.0 / (1 + exp(-(y[38] + 6) / (8.6))); // change to simple version, haibo.  removed IKur mutations here.
	tau = (45.67 / (1.0 + exp((y[38] + 9.23) / 10.53)) + 4.27) * ( 0.55 / (1.0 + exp((y[38] + 60.87) / -22.88)) + 0.03); // incorporate deactivation Thu 20 Jul 2017 23:24:01 BST Haibo
	tau = tau / K_Q10;


	// y[42] = inf + (y[42] - inf) * exp(-(dt) / tau);
	ydot[42] = (inf - y[42]) / tau;

	// y[43]
	inf = 1.0 / (1.0 + exp( (y[38] + 7.5) / 6.67)); // Thu 20 Jul 2017 23:24:21 BST
	tau = 1300.0 / (exp((y[38] - 133.0) / 97.0) * 5.04 + 0.2 * exp((y[38] + 10.9233071) / (-12.0))) + 370.0; // at 37deg  Thu 20 Jul 2017 23:24:35 BST haibo Feng et al. 1998
	// y[43] = inf + (y[43] - inf) * exp(-(dt) / tau);
	ydot[43] = (inf - y[43]) / tau;


	tau = 1400.0 / (exp((y[38] + 0.1) / 8.86) + 1.0) + 478.0 / (exp((y[38] - 41.0) / -3.0) + 1.0) + 6500.0; // // at 37deg  Thu 20 Jul 2017 23:24:35 BST haibo Feng et al. 1998
	// y[44] = inf + (y[44] - inf) * exp(-(dt) / tau);
	ydot[44] = (inf - y[44]) / tau;

	double IKur = CNZ_gkur * (4.5128 + 1.899769 / (1.0 + exp((y[39 - 1] - 20.5232) / (-8.26597)))) * y[42] * (0.65 * y[43] + 0.35 * y[44]) * (y[39 - 1] - ek);

	// I_kur = IKur;  //

	//Gki = SA_par(9)*(1+0.45*(1-AF)+0.70*AF)*(1+1*AF)*(2.1*0.0525)*sqrt(Ko/5.4); // NEW K2P +45// nSR/70// cAF in Schmidt et al 2015


	// double aki_sl = (0.1 + 0.9 / (1 + (y[33 - 1] / 7) * (y[33 - 1] / 7.0))) * 1 / (1 + exp(0.2385 * (y[39 - 1] - ek - 59.215)));
	// double aki_j = (0.1 + 0.9 / (1 + (y[32 - 1] / 7) * (y[32 - 1] / 7.0))) * 1 / (1 + exp(0.2385 * (y[39 - 1] - ek - 59.215)));

	// Na-dependence (Schmidt et al 2015)
	// aki_sl = 1/(1+(y[33-1]/7)^2);
	// aki_j = 1/(1+(y[32-1]/7)^2);

	//// I_ki: Time-Independent K Current
	double Gki = /*0.95 * */0.75 * SA_par[8] * (1 + 1 * AF) * (2.1 * 0.0525) * sqrt(Ko / 5.4);


	double fracIK1_avail = 1 + (0.55 - 1) * para.kPKA_IK1; // 45 % decrease w / 100 nM ISO (Gonzales De La Fuente et al 2013)
	// increased effect (-55 % ) in cAF - to be included
	Gki = para.kCKII_IK1_G * fracIK1_avail * Gki;
	// Na-dependence (Voigt-Heijman et al 2013)
	double aki_sl = (0.15 + 0.85 / (1 + (y[33 - 1] / 7.0) * (y[33 - 1] / 7.0))) * 1 / (1 + exp(0.2385 * (y[39 - 1] - ek - 59.215)));
	double aki_j = (0.15 + 0.85 / (1 + (y[32 - 1] / 7.0) * (y[32 - 1] / 7.0))) * 1 / (1 + exp(0.2385 * (y[39 - 1] - ek - 59.215)));
	double bki = (0.49124 * exp(0.08032 * (y[39 - 1] + 5.476 - ek)) + exp(0.06175 * (y[39 - 1] - ek - 594.31))) / (1 + exp(-0.5143 * (y[39 - 1] - ek + 4.753)));
	double kiss_sl = aki_sl / (aki_sl + bki);
	double kiss_j = aki_j / (aki_j + bki);

	I_ki_sl = para.IK1_Scale * Fsl * Gki * kiss_sl * (y[39 - 1] - ek);
	I_ki_j = para.IK1_Scale * Fjunc * Gki * kiss_j * (y[39 - 1] - ek);


	// replace with CRN IK1, adding Nai-dependence
	I_ki_j = para.IK1_Scale * 1 / 1.31 * 0.65 * (0.15 + 0.85 / (1 + (y[32 - 1] / 10.0) * (y[32 - 1] / 10.0))) * Fjunc * Gki * (y[39 - 1] - ek) / (1 + exp(0.07 * (y[39 - 1] - (ek + 6.94)))); // CRN IK1
	I_ki_sl = para.IK1_Scale * 1 / 1.31 * 0.65 * (0.15 + 0.85 / (1 + (y[33 - 1] / 10.0) * (y[33 - 1] / 10.0))) * Fsl * Gki * (y[39 - 1] - ek) / (1 + exp(0.07 * (y[39 - 1] - (ek + 6.94)))); // CRN IK1
	// if(y[38] > -50) I_ki += 0.1;

	// I_ki = 0.3*Gki * (y[39 - 1] - ek) / (1 + exp(1.5 * (y[39 - 1] - ek + 3.6 )*FoRT));
	I_ki = I_ki_j + I_ki_sl;
	//// I_kp: Plateau K current
	double gkp = /*10**/ 0.95 * SA_par[9] * 0.002;

	double kp_kp = 1 / (1 + exp(7.488 - y[39 - 1] / 5.98));
	I_kp = para.IKp_Scale * gkp * kp_kp * (y[39 - 1] - ek);

	//// I_k2p: K2P3.1 K current (from Schmidt et al 2015)
	//gk2p = SA_par(11)*K2P_cond*(0.005+0.014*AF); // wrong
	//gk2p = SA_par(11)*K2P_cond*0.0050*(1+(0.0145/0.005-1)*AF); // 2.9-fold increase in AF
	double gk2p = /*1.1 **/ 0.9 * para.IK2p_Scale /** 1.2*/ * 0.95 * SA_par[10] * K2P_cond * 0.0050 * (1 + (SA_par[23] * 0.0145 / 0.005 - 1) * AF);

	inf = 0.2 + 0.8 / (1 + exp(-(y[39 - 1] - 10 + 15 * AF) / 14));
	tau = 2 + 40 / (1 + exp((y[39 - 1] + 25) * (y[39 - 1] + 25) / 80));
	ydot[41] = (inf - y[41]) / tau;

	I_k2p = gk2p * y[41] * (y[39 - 1] - ek); // NEW K2P





	//// I_k,ach: Muscarinic-Receptor-Activated K Current
	// Grandi (ACh)
	// I_KAch = 1/(1+(0.03/Ach)^2.1)*(0.08+0.4./(1+exp((y[39-1]+91)/12))).*(y[39-1]-ek);

	// Schmidt (CCh)
	// fkach = 1.5/(1+pow((9/y[33-1]),4)); // Na-dependence (NaSL)
	// rkach = 0.055+0.4/(1+exp((y[39-1]-ek+9.53)/17.18));
	// if AF == 0
	//     Gkach = 0.10*(1+fkach)/0.9411765;
	// else
	//     Gkach = 0.05/0.9411765;
	// end
	// I_kach = Gkach*sqrt(Ko/5.4)*(CCh/(CCh+0.125))*rkach*(y[39-1]-ek);
	double Gkach_sl = 0, Gkach_j = 0, d_kach = 0, r_kach = 0;
	double fkach_sl = 0, fkach_j = 0;
	if (CCh == 0 and ACh == 0) {// none
		Gkach_sl = SA_par[11] * 0;
		Gkach_j = SA_par[11] * 0;
		d_kach = 0;
		r_kach = 0;
	} else if (CCh > 0 and ACh == 0) { // CCh
		// Na-dependence
		fkach_sl = 1.5 / (1 + pow((9 / y[33 - 1]), 4)); // Na-dependence sl
		fkach_j = 1.5 / (1 + pow((9 / y[32 - 1]), 4)); // Na-dependence j
		if (AF == 0) {
			Gkach_sl = SA_par[11] * (RA * 1 + (1 - RA) * 0.3) * 0.10 * (1 + fkach_sl) / 0.9411765 * sqrt(Ko / 5.4);
			Gkach_j = SA_par[11] * (RA * 1 + (1 - RA) * 0.3) * 0.10 * (1 + fkach_j) / 0.9411765 * sqrt(Ko / 5.4);
		} else {
			Gkach_sl = SA_par[11] * (RA * 0.5 + (1 - RA) * 0.3) * 0.10 / 0.9411765 * sqrt(Ko / 5.4);
			Gkach_j = SA_par[11] * (RA * 0.5 + (1 - RA) * 0.3) * 0.10 / 0.9411765 * sqrt(Ko / 5.4);
		}
		// V-dependence
		r_kach = 0.055 + 0.4 / (1 + exp((y[39 - 1] - ek + 9.53) / 17.18));
		// Dose-dependence
		d_kach = CCh / (CCh + 0.125);
	} else if (CCh == 0 and ACh > 0) { // ACh
		// Na-dependence
		fkach_sl = 1.5 / (1 + pow((9 / y[33 - 1]), 4)); // Na-dependence sl
		fkach_j = 1.5 / (1 + pow((9 / y[32 - 1]), 4)); // Na-dependence j
		if (AF == 0) {
			Gkach_sl = SA_par[11] * (RA * 1 + (1 - RA) * 0.3) * 5 * 0.10 * (1 + fkach_sl) / 0.9411765 * sqrt(Ko / 5.4);
			Gkach_j = SA_par[11] * (RA * 1 + (1 - RA) * 0.3) * 5 * 0.10 * (1 + fkach_j) / 0.9411765 * sqrt(Ko / 5.4);
		} else {
			Gkach_sl = SA_par[11] * (RA * 0.5 + (1 - RA) * 0.3) * 5 * 0.10 / 0.9411765 * sqrt(Ko / 5.4);
			Gkach_j = SA_par[11] * (RA * 0.5 + (1 - RA) * 0.3) * 5 * 0.10 / 0.9411765 * sqrt(Ko / 5.4);
		}
		// V-dependence
		r_kach = 0.08 + 0.4 / (1 + exp((y[39 - 1] + 91) / 12));
		// Dose-dependence
		d_kach = 1 / (1 + pow(0.03 / ACh, 2.1));
	}

	I_kach_sl = para.IKach_Scale * Fsl * Gkach_sl * d_kach * r_kach * (y[39 - 1] - ek);
	I_kach_j = para.IKach_Scale * Fjunc * Gkach_j * d_kach * r_kach * (y[39 - 1] - ek);

	I_kach = I_kach_sl + I_kach_j;

	//// I_sk: Small-Conductance Ca-Activated K Current
	double par_sk[] = {0.0506381114404388, 0.273335569451572, 2.96381060498817, 0.199981221802789, 0.279328126521496, -86.9289059836381, 0.00636311816933264, 5.22915055145375};

	//gsk = SA_par(13)*SK_cond*par_sk(1)*(1+(par_sk(8)-1)*AF);
	double gsk = /*1.05 * */0.95 * para.ISK_Scale  * 0.85 * SA_par[12] * SK_cond * par_sk[0] * (1 + (SA_par[24] * par_sk[7] - 1) * AF);

	double kdsk = SA_par[13] * (pow(10, (SK_shift - 3.45))); // kdsk = (10^(ISK_shift-3.3)); // (mM)
	double gsk_ca_junc = 1.0 / (1 + exp((log10(kdsk) - log10(y[35])) / 0.3));
	double gsk_ca_sl = 1.0 / (1 + exp((log10(kdsk) - log10(y[36])) / 0.3));
	//     // Ca clamp for SK channels
	//     gsk_ca = 500e-6;
	//     gsk_ca_junc = 1/(1+ exp((log10(kdsk)-log10(gsk_ca))/0.3));
	//     gsk_ca_sl = 1/(1+ exp((log10(kdsk)-log10(gsk_ca))/0.3));

	double gsk_vm = par_sk[1] / (1 + exp((y[39 - 1] - ek + par_sk[2]) * par_sk[3])) + par_sk[4] / (1 + exp((-(y[39 - 1] - ek + par_sk[5]) * par_sk[6])));

	I_sk_junc = Fjunc * gsk * gsk_ca_junc * gsk_vm * (y[39 - 1] - ek);
	I_sk_sl = Fsl * gsk * gsk_ca_sl * gsk_vm * (y[39 - 1] - ek);
	I_sk = I_sk_junc + I_sk_sl;

	//// I_ClCa and I_Clbk: Ca-activated and Background Cl Currents
	double GClCa = 1.178030 * SA_par[14] * 0.0548;   // [mS/uF]
	double KdClCa = (1 + (0.704 - 1) * para.kPKA_IClCa) * 100e-3;       // [mM]
	//GClB = SA_par(16)*9e-3;        // [mS/uF]
	double GClB = /*SA_par[15] * 0.5 * */9e-3;    // [mS/uF] MOD1

	I_ClCa_junc = para.IClCa_Scale * Fjunc * GClCa / (1 + KdClCa / y[35]) * (y[39 - 1] - ecl);
	I_ClCa_sl = para.IClCa_Scale * Fsl * GClCa / (1 + KdClCa / y[37 - 1]) * (y[39 - 1] - ecl);
	I_ClCa = I_ClCa_junc + I_ClCa_sl;

	// I_Clbk = GClB * (y[39 - 1] - ecl);


	I_Clbk = /*0.5 **/ /*0.95**/para.IClb_Scale * 1.05 * 0.19e-3 * (y[39 - 1] - ecl) / (1 - 0.94 * exp(2.5e-4 * (y[39 - 1] - ecl))); 	// Li et al. An outwardly rectifying anionic background current in atrial myocytes from the human heart

	double GClCFTR = 0; // 4.9e-3*ISO;     // [mS/uF]
	I_ClCFTR = GClCFTR * (y[39 - 1] - ecl);

	//// I_Ca: L-type Ca Current

#ifdef OLDLTCC
	double pNa = 0 * SA_par[16] * (1 + 0.5 * ISO) * (1 - 0.5 * AF) * 0.75e-8; // [cm/sec]
	double pCa = 0 * SA_par[16] * (1 + 0.5 * ISO) * (1 - 0.5 * AF) * 2.7e-4; // [cm/sec]
	double pK = 0 * SA_par[16] * (1 + 0.5 * ISO) * (1 - 0.5 * AF) * 1.35e-7; // [cm/sec]
	double Q10CaL = 1.8;

	double dss = 1.0 / (1 + exp(-(y[38] + 3 * ISO + 9) / 6.0));
	//dss = 1/(1+exp(-(y[39-1]+3*ISO+9)/(6*(1-AF)+7*AF))); // NEW K2P modifications in Schmidt et al 2015
	double taud = dss * (1 - exp(-(y[38] + 3 * ISO + 9) / 6)) / (0.035 * (y[38] + 3 * ISO + 9));
	double fss = 1 / (1 + exp((y[38] + 3 * ISO + 30) / 7)) + 0.2 / (1 + exp((50 - y[38] - 3 * ISO) / 20.0));
	double tauf = 1 / (0.0197 * exp( -(0.0337 * (y[38] + 3 * ISO + 25)) * (0.0337 * (y[38] + 3 * ISO + 25)) ) + 0.02);
	ydot[3] = (dss - y[3]) / taud;
	ydot[4] = (fss - y[4]) / tauf;
	ydot[5] = 1.7 * y[35] * (1 - y[5]) - 1 * 11.9e-3 * y[5]; // fCa_junc
	ydot[6] = 1.7 * y[36] * (1 - y[6]) - 1 * 11.9e-3 * y[6]; // fCa_sl
	//fcaCaMSL = 0.1/(1+(0.01/y[37-1]));
	//fcaCaj = 0.1/(1+(0.01/y[35]));
	double fcaCaMSL = 0;
	double fcaCaj = 0;
	double ibarca_j = pCa * 4 * (y[39 - 1] * Frdy * FoRT) * (0.341 * y[35] * exp(2 * y[39 - 1] * FoRT) - 0.341 * Cao) / (exp(2 * y[39 - 1] * FoRT) - 1);
	double ibarca_sl = pCa * 4 * (y[39 - 1] * Frdy * FoRT) * (0.341 * y[37 - 1] * exp(2 * y[39 - 1] * FoRT) - 0.341 * Cao) / (exp(2 * y[39 - 1] * FoRT) - 1);
	double ibark = pK * (y[39 - 1] * Frdy * FoRT) * (0.75 * y[35 - 1] * exp(y[39 - 1] * FoRT) - 0.75 * Ko) / (exp(y[39 - 1] * FoRT) - 1);
	double ibarna_j = pNa * (y[39 - 1] * Frdy * FoRT) * (0.75 * y[32 - 1] * exp(y[39 - 1] * FoRT) - 0.75 * Nao) / (exp(y[39 - 1] * FoRT) - 1);
	double ibarna_sl = pNa * (y[39 - 1] * Frdy * FoRT) * (0.75 * y[33 - 1] * exp(y[39 - 1] * FoRT) - 0.75 * Nao) / (exp(y[39 - 1] * FoRT) - 1);

	I_Ca_junc = para.ICaL_Scale * (Fjunc_CaL * ibarca_j * y[3] * y[4] * ((1 - y[5]) + fcaCaj)) * 0.45; // *Q10CaL^Qpow was removed.
	I_Ca_sl = para.ICaL_Scale * (Fsl_CaL * ibarca_sl * y[3] * y[4] * ((1 - y[6]) + fcaCaMSL)) * 0.45; // *Q10CaL^Qpow was removed.
	I_CaK = para.ICaL_Scale * (ibark * y[3] * y[4] * (Fjunc_CaL * (fcaCaj + (1 - y[5])) + Fsl_CaL * (fcaCaMSL + (1 - y[6])))) * 0.45; // *Q10CaL^Qpow was removed.
	I_CaNa_junc = para.ICaL_Scale * (Fjunc_CaL * ibarna_j * y[3] * y[4] * ((1 - y[5]) + fcaCaj)) * 0.45; // *Q10CaL^Qpow was removed.
	I_CaNa_sl = para.ICaL_Scale * (Fsl_CaL * ibarna_sl * y[3] * y[4] * ((1 - y[6]) + fcaCaMSL)) * .45; // *Q10CaL^Qpow was removed.
#endif

	// adding two LTCCs m1 and m2 for junc and sl compartments
	for (int i = 0; i < 10; ++i)
	{
		LTCC.state[i] = y[45 + i];
		LTCC_sl.state[i] = y[45 + 10 + i];
		LTCC_m2.state[i] = y[45 + 20 + i];
		LTCC_sl_m2.state[i] = y[45 + 10 + 20 + i];
	}

	/*LTCC.update_states_v7(y[35], y[38], y[31], y[34], para.K_PKA_LTCC);
	LTCC_sl.update_states_v7(y[36], y[38], y[32], y[34], para.K_PKA_LTCC);
	LTCC_m2.update_states_v7_mode_2(y[35], y[38], y[31], y[34], para.K_PKA_LTCC);
	LTCC_sl_m2.update_states_v7_mode_2(y[36], y[38], y[32], y[34], para.K_PKA_LTCC);*/


	/*LTCC.update_states_v7_BPS_2020(y[35], y[38], y[31], y[34], para.K_PKA_LTCC);
	LTCC_sl.update_states_v7_BPS_2020(y[36], y[38], y[32], y[34], para.K_PKA_LTCC);
	LTCC_m2.update_states_v7_BPS_2020_mode_2(y[35], y[38], y[31], y[34], para.K_PKA_LTCC);
	LTCC_sl_m2.update_states_v7_BPS_2020_mode_2(y[36], y[38], y[32], y[34], para.K_PKA_LTCC);*/

	/*LTCC.update_states_v8(y[35], y[38], y[31], y[34], para.K_PKA_LTCC);
	LTCC_sl.update_states_v8(y[36], y[38], y[32], y[34], para.K_PKA_LTCC);
	LTCC_m2.update_states_v8_mode_2(y[35], y[38], y[31], y[34], para.K_PKA_LTCC);
	LTCC_sl_m2.update_states_v8_mode_2(y[36], y[38], y[32], y[34], para.K_PKA_LTCC);
	*/

	LTCC.update_states_v9(y[35], y[38], y[31], y[34], para.K_PKA_LTCC);
	LTCC_sl.update_states_v9(y[36], y[38], y[32], y[34], para.K_PKA_LTCC);
	LTCC_m2.update_states_v9_mode_2(y[35], y[38], y[31], y[34], para.K_PKA_LTCC);
	LTCC_sl_m2.update_states_v9_mode_2(y[36], y[38], y[32], y[34], para.K_PKA_LTCC);

	for (int i = 0; i < 10; ++i)
	{
		ydot[45 + i] = LTCC.ydot[i];
		ydot[45 + 10 + i] = LTCC_sl.ydot[i];
		ydot[45 + 20 + i] = LTCC_m2.ydot[i];
		ydot[45 + 20 + 10 + i] = LTCC_sl_m2.ydot[i];
	}


	// 0.1 and 0.9 are fractional LTCC in sl and junctional spaces, respectively
	double scale = para.ICaL_Scale * 1.02 * 1.209 * 1.03 * 0.85 * 0.9 * 1.037031 * 1.05; //1.3 * 0.93;
	I_Ca_junc = scale * (1 - 0.3 * AF) * para.G_LTCC * (0.9 * ( (1 - para.LTCC_junc_mode2) * LTCC.I_Ca_junc_m1  + para.LTCC_junc_mode2 * LTCC_m2.I_Ca_junc_m1) );
	I_CaNa_junc = scale * (1 - 0.3 * AF) * para.G_LTCC * (0.9 * ( (1 - para.LTCC_junc_mode2) * LTCC.I_Na_junc_m1  + para.LTCC_junc_mode2 * LTCC_m2.I_Na_junc_m1) );
	I_Ca_sl =  scale * (1 - 0.3 * AF) * para.G_LTCC * (0.1 * ( (1 - para.LTCC_sl_mode2) * LTCC_sl.I_Ca_junc_m1  + para.LTCC_sl_mode2 * LTCC_sl_m2.I_Ca_junc_m1) );
	I_CaNa_sl =  scale * (1 - 0.3 * AF) * para.G_LTCC * (0.1 * ( (1 - para.LTCC_sl_mode2) * LTCC_sl.I_Na_junc_m1  + para.LTCC_sl_mode2 * LTCC_sl_m2.I_Na_junc_m1) );
	I_CaK =  scale * (1 - 0.3 * AF) * para.G_LTCC * (0.9 * ( (1 - para.LTCC_junc_mode2) * LTCC.I_K_junc_m1  + para.LTCC_junc_mode2 * LTCC_m2.I_K_junc_m1)
	         + 0.1 * ( (1 - para.LTCC_sl_mode2) * LTCC_sl.I_K_junc_m1  + para.LTCC_sl_mode2 * LTCC_sl_m2.I_K_junc_m1) );

	LTCC_mode_1_out = scale * (1 - 0.3 * AF) * para.G_LTCC * (
	                      0.9 * ( (1 - para.LTCC_junc_mode2) * LTCC.I_Ca_junc_m1 )
	                      + 0.1 * ( (1 - para.LTCC_sl_mode2) * LTCC_sl.I_Ca_junc_m1 )
	                  );
	LTCC_mode_2_out = scale * (1 - 0.3 * AF) * para.G_LTCC * (
	                      0.9 * ( para.LTCC_junc_mode2 * LTCC_m2.I_Ca_junc_m1 )
	                      + 0.1 * ( para.LTCC_sl_mode2 * LTCC_sl_m2.I_Ca_junc_m1 )
	                  );

	// I_Ca_junc = 0.9 * LTCC_m2.I_Ca_junc_m1;
	// I_CaNa_junc = 0.9 * LTCC_m2.I_Na_junc_m1;

	// I_Ca_sl = 0.1 * LTCC_sl_m2.I_Ca_junc_m1;
	// I_CaNa_sl = 0.1 * LTCC_sl_m2.I_Na_junc_m1;

	// I_CaK = 0.9 * LTCC_m2.I_K_junc_m1 + 0.1 * LTCC_sl_m2.I_K_junc_m1;

	I_Ca = I_Ca_junc + I_Ca_sl;
	I_CaNa = I_CaNa_junc + I_CaNa_sl;
	I_Catot = I_Ca + I_CaK + I_CaNa; //#ok<NASGU>
	//// I_cabk: Ca Background Current
	double GCaB = /*0.85**/para.ICab_Scale * /*1.073250 **/ /*SA_par[17] * */6.0643e-4; // [uA/uF]
	//GCaB = SA_par[17]*6.0643e-4*(1+1*AF); // [uA/uF] MOD2
	if (MOD_ind == 3) {
		GCaB = /*0.85**//*SA_par[17] **/ para.ICab_Scale * 6.0643e-4 * (1 + 0.5 * AF); // [uA/uF] MOD3
	}

	I_cabk_junc = Fjunc * GCaB * (y[39 - 1] - eca_junc);
	I_cabk_sl = Fsl * GCaB * (y[39 - 1] - eca_sl);
	I_cabk = I_cabk_junc + I_cabk_sl; //#ok<NASGU>

	//// I_pca: Sarcolemmal Ca Pump Current
	double IbarSLCaP = para.ICap_Scale * SA_par[18] * 0.0471; // [uA/uF]
	//IbarSLCaP = SA_par[18]*0.0471*4; // [uA/uF] MOD2
	if (MOD_ind == 3) {
		IbarSLCaP = para.ICap_Scale * 1.0 * 1.022257 * 0.9 * SA_par[18] * 0.0471 * 2; // [uA/uF] MOD3
	}

	double KmPCa = 0.5e-3;     // [mM]
	double Q10SLCaP = 2.35;    // [none]

	I_pca_junc = Fjunc * IbarSLCaP * pow(y[35], 1.6) / (pow(KmPCa, 1.6) + pow(y[35], 1.6)); // Q10SLCaP^Qpow* was removed
	I_pca_sl = Fsl * IbarSLCaP * pow(y[36], 1.6) / (pow(KmPCa, 1.6) + pow(y[36], 1.6)); // Q10SLCaP^Qpow* was removed
	I_pca = I_pca_junc + I_pca_sl;

	//// I_ncx: Na/Ca Exchanger flux
	double IbarNCX =  1 * SA_par[19] * (1 + 0.4 * AF) * 3.15; // [uA/uF]
	//IbarNCX = SA_par[19]*(1+1*AF)*3.15; // [uA/uF] MOD2
	if (MOD_ind == 3) {
		IbarNCX =   /*0.95 **/1.1 * 1.05 * 1.2 * SA_par[19] * (1 + 0.8 * AF) * 3.15; // [uA/uF] MOD3
	}

	double KmCai = 3.59e-3;    // [mM]
	double KmCao = 1.3;        // [mM]
	double KmNai = 12.29;      // [mM]
	double KmNao = 87.5;       // [mM]
	double ksat = 0.27;        // [none]
	double nu = 0.35;          // [none]
	double Kdact = 0.384e-3;   // [mM]
	double Q10NCX = 1.57;      // [none]

	double Ka_junc = 1 / (1 + (Kdact / y[35]) * (Kdact / y[35]));
	double Ka_sl = 1 / (1 + (Kdact / y[37 - 1]) * (Kdact / y[37 - 1]));
	double s1_junc = exp(nu * y[39 - 1] * FoRT) * y[32 - 1] * y[32 - 1] * y[32 - 1] * Cao;
	double s1_sl = exp(nu * y[39 - 1] * FoRT) * y[33 - 1] * y[33 - 1] * y[33 - 1] * Cao;
	double Nao_3 = Nao * Nao * Nao;
	double Naj_3 = y[32 - 1] * y[32 - 1] * y[32 - 1];
	double NaSL_3 = y[33 - 1] * y[33 - 1] * y[33 - 1];
	double s2_junc = exp((nu - 1) * y[39 - 1] * FoRT) * Nao_3 * y[35];
	double s3_junc = KmCai * Nao_3 * (1 + (y[32 - 1] / KmNai) * (y[32 - 1] / KmNai) * (y[32 - 1] / KmNai)) + KmNao * KmNao * KmNao * y[35] * (1 + y[35] / KmCai) + KmCao * Naj_3 + Naj_3 * Cao + Nao_3 * y[35];
	double s2_sl = exp((nu - 1) * y[39 - 1] * FoRT) * Nao_3 * y[37 - 1];
	double s3_sl = KmCai * Nao_3 * (1 + (y[33 - 1] / KmNai) * (y[33 - 1] / KmNai) * (y[33 - 1] / KmNai)) + KmNao * KmNao * KmNao * y[37 - 1] * (1 + y[37 - 1] / KmCai) + KmCao * NaSL_3 + NaSL_3 * Cao + Nao_3 * y[37 - 1];

	I_ncx_junc = para.INCX_Scale * Fjunc * IbarNCX * Ka_junc * (s1_junc - s2_junc) / s3_junc / (1 + ksat * exp((nu - 1) * y[39 - 1] * FoRT)); //Q10NCX^Qpow*removed
	I_ncx_sl = para.INCX_Scale * Fsl * IbarNCX * Ka_sl * (s1_sl - s2_sl) / s3_sl / (1 + ksat * exp((nu - 1) * y[39 - 1] * FoRT)); //Q10NCX^Qpow*removed
	I_ncx = I_ncx_junc + I_ncx_sl;

	//// SR fluxes: Calcium Uptake, Release, and Leak
	double Vmax_SRCaP = 1.5 * SA_par[20] * 5.3114e-3; // [mM/msec] (286 umol/L cytosol/sec)
	//Vmax_SRCaP = SA_par[20]*5.3114e-3*(1-0.5*AF); // [mM/msec] MOD2
	if (MOD_ind == 3) {
		Vmax_SRCaP = /*1.2 **/ 1.1 * 1.5 * SA_par[20] * 5.3114e-3 * (1 - 0.25 * AF); // [mM/msec] MOD3
	}

	double Q10SRCaP = 2.6;          // [none]
	double Kmf = para.PLB_kmf_Scale * (2.5 /*- 1.25 * ISO*/) * 0.246e-3; // [mM]
	// std::cout << para.PLB_kmf_Scale << std::endl;
	double Kmr = 1.7;               // [mM]L cytosol
	double hillSRCaP = 1.787;       // [mM]
	double ks = /*SA_par[22 - 1] **/ 1.273543 * 1.1 * 25; // [1/ms]
	double koCa = 10 + 20 * AF /*+ 10 * ISO * (1 - AF)*/; // [mM^-2 1/ms]  // ISO will be modeled by PKA
	double kom = 0.06;              // [1/ms]
	double kiCa = 0.5;              // [1/mM/ms]
	double kim = 0.005;             // [1/ms]
	// double kim = 0.01;             // [1/ms]  // 15:07:36, Thu, 07-February-2019, By Haibo
	double ec50SR = 0.45;           // [mM]
	double MaxSR = 15;              // [mM]
	double MinSR = 1;               // [mM]

	// leak regulation by CaMKII from Bartos et al. 2017
	double kleak = (1.0 / 3.0 + 10 * para.RyR_CKp / 3.0) * SA_par[23 - 1] * (1.0 + 0.25 * AF) * /*1.134138 **/ 5.348e-6;
	// double kleak = (2.0 / 3.0 + 5 * para.RyR_CKp / 3.0) * SA_par[23 - 1] * (1.0 + 0.25 * AF) * 1.134138 * 5.348e-6;

	double kCaSR = MaxSR - (MaxSR - MinSR) / (1 + pow((ec50SR / y[31 - 1]), 2.5));
	double koSRCa = para.RyR_koSRCa_Scale * koCa / kCaSR;
	double kiSRCa = kiCa * kCaSR;

	double KoSR_Ca_2 = koSRCa * y[35] * y[35];//+1e-6;  // original
	// double KoSR_Ca_2 = koSRCa * 1e-4 * 1.0 / (1.0 + pow(0.03/y[35], 2.0));
	// double KoSR_Ca_2 = koSRCa * y[35] * y[35];
	double RI = 1 - y[13] - y[14] - y[15];
	ydot[13] = (kim * RI - kiSRCa * y[35] * y[13]) - (KoSR_Ca_2 * y[13] - kom * y[14]); // R
	ydot[14] = (KoSR_Ca_2 * y[13] - kom * y[14]) - (kiSRCa * y[35] * y[14] - kim * y[15]); // O
	ydot[15] = (kiSRCa * y[35] * y[14] - kim * y[15]) - (kom * y[15] - KoSR_Ca_2 * RI); // 0*
	J_SRCarel = para.Jrel_Scale * ks * y[14] * (y[31 - 1] - y[35]); // [mM/ms]
	J_serca = para.Jserca_Scale * Vmax_SRCaP * (pow(y[37] / Kmf, hillSRCaP) - pow(y[30] / Kmr, hillSRCaP)) / (1 + pow(y[37] / Kmf, hillSRCaP) + pow(y[30] / Kmr, hillSRCaP)); // [mM/ms]   Q10SRCaP^Qpow*
	J_SRleak = /*2 **/ para.Jleak_Scale * kleak * (y[30] - y[35]); // [mM/ms]

	//// Sodium and Calcium Buffering
	// koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
	double BufferScaling_cytosol = /*1.069467 **/ 1 * para.Cytosol_Buffer_Scale;

	double Bmax_Naj = 7.561;       // [mM] // Na buffering
	double Bmax_Nasl = 1.65;       // [mM]
	double koff_na = 1e-3;         // [1/ms]
	double kon_na = 0.1e-3;        // [1/mM/ms]
	double Bmax_TnClow = 70e-3 * BufferScaling_cytosol;  // [mM]                      // TnC low affinity
	double koff_tncl = /*(1 + 0.5 * ISO)*/ para.fPKA_TnI * 19.6e-3; // [1/ms]
	double kon_tncl = 32.7;        // [1/mM/ms]
	double Bmax_TnChigh = 140e-3 * BufferScaling_cytosol; // [mM]                      // TnC high affinity
	double koff_tnchca = 0.032e-3; // [1/ms]
	double kon_tnchca = 2.37;      // [1/mM/ms]
	double koff_tnchmg = 3.33e-3;  // [1/ms]
	double kon_tnchmg = 3e-3;      // [1/mM/ms]
	double Bmax_CaM = 24e-3 * BufferScaling_cytosol;     // [mM] **? about setting to 0 in c-code**   // CaM buffering
	double koff_cam = 238e-3;      // [1/ms]
	double kon_cam = 34;           // [1/mM/ms]
	double Bmax_myosin = 140e-3 * BufferScaling_cytosol; // [mM]                      // Myosin buffering
	double koff_myoca = 0.46e-3;   // [1/ms]
	double kon_myoca = 13.8;       // [1/mM/ms]
	double koff_myomg = 0.057e-3;  // [1/ms]
	double kon_myomg = 0.0157;     // [1/mM/ms]
	double Bmax_SR = 19 * .9e-3 * BufferScaling_cytosol; // [mM] (Bers text says 47e-3) 19e-3
	double koff_sr = 60e-3;        // [1/ms]
	double kon_sr = 100;           // [1/mM/ms]

	double BufferScaling = /*1.111454 **/ 1 * para.Cleft_Buffer_Scale; // scale both cleft and sub-membrane
	double Bmax_SLlowsl = BufferScaling * 37.4e-3 * Vmyo / Vsl;    // [mM]   	// SL buffering
	double Bmax_SLlowj = BufferScaling * 4.6e-3 * Vmyo / Vjunc * 0.1; // [mM]
	double koff_sll = 1300e-3;     // [1/ms]
	double kon_sll = 100;          // [1/mM/ms]
	double Bmax_SLhighsl = BufferScaling * 13.4e-3 * Vmyo / Vsl;   // [mM]
	double Bmax_SLhighj = BufferScaling * 1.65e-3 * Vmyo / Vjunc * 0.1; // [mM]
	double koff_slh = 30e-3;       // [1/ms]
	double kon_slh = 100;          // [1/mM/ms]
	double Bmax_Csqn = 140e-3 * Vmyo / Vsr * para.SR_Buffer_Scale /** 1.064127*/;        // [mM]      // Csqn buffering
	double koff_csqn = 65;         // [1/ms]
	double kon_csqn = 100;         // [1/mM/ms]


	// Junctional and SL Na Buffers
	ydot[16] = kon_na * y[32 - 1] * (Bmax_Naj - y[16]) - koff_na * y[16]; // NaBj      [mM/ms]
	ydot[17] = kon_na * y[33 - 1] * (Bmax_Nasl - y[17]) - koff_na * y[17]; // NaBsl     [mM/ms]

	// Cytosolic Ca Buffers
	ydot[18] = kon_tncl * y[37] * (Bmax_TnClow - y[18]) - koff_tncl * y[18];  // TnCL      [mM/ms]
	ydot[19] = kon_tnchca * y[37] * (Bmax_TnChigh - y[19] - y[20]) - koff_tnchca * y[19]; // TnCHc     [mM/ms]
	ydot[20] = kon_tnchmg * Mgi * (Bmax_TnChigh - y[19] - y[20]) - koff_tnchmg * y[20]; // TnCHm     [mM/ms]
	ydot[21] = 0.0; // kon_cam * y[37] * (Bmax_CaM - y[21]) - koff_cam * y[21];       // CaM       [mM/ms]  // With signaling, this part is taken care of by CaM modules
	ydot[22] = kon_myoca * y[37] * (Bmax_myosin - y[22] - y[23]) - koff_myoca * y[22]; // Myosin_ca [mM/ms]
	ydot[23] = kon_myomg * Mgi * (Bmax_myosin - y[22] - y[23]) - koff_myomg * y[23]; // Myosin_mg [mM/ms]
	ydot[24] = kon_sr * y[37] * (Bmax_SR - y[24]) - koff_sr * y[24];          // SRB       [mM/ms]
	J_CaB_cytosol = ydot[18] + ydot[19] + ydot[21] + ydot[22] + ydot[24];

	// Junctional and SL Ca Buffers
	ydot[25] = kon_sll * y[35] * (Bmax_SLlowj - y[25]) - koff_sll * y[25]; // SLLj      [mM/ms]
	ydot[26] = kon_sll * y[37 - 1] * (Bmax_SLlowsl - y[26]) - koff_sll * y[26]; // SLLsl     [mM/ms]
	ydot[27] = kon_slh * y[35] * (Bmax_SLhighj - y[27]) - koff_slh * y[27]; // SLHj      [mM/ms]
	ydot[28] = kon_slh * y[37 - 1] * (Bmax_SLhighsl - y[28]) - koff_slh * y[28]; // SLHsl     [mM/ms]
	J_CaB_junction = ydot[25] + ydot[27];
	J_CaB_sl = ydot[26] + ydot[28];

	//// Ion concentrations
	// SR Ca Concentrations
	ydot[29] = kon_csqn * y[30] * (Bmax_Csqn - y[29]) - koff_csqn * y[29]; // Csqn      [mM/ms]
	ydot[30] = J_serca - (J_SRleak * Vmyo / Vsr + J_SRCarel) - ydot[29]; // Ca_sr     [mM/ms]

	// Sodium Concentrations
	I_Na_tot_junc = I_Na_junc + I_nabk_junc + 3 * I_ncx_junc + 3 * I_nak_junc + I_CaNa_junc + I_NaL_junc; // [uA/uF]
	I_Na_tot_sl = I_Na_sl + I_nabk_sl + 3 * I_ncx_sl + 3 * I_nak_sl + I_CaNa_sl + I_NaL_sl; // [uA/uF]

	ydot[31] = -I_Na_tot_junc * Cmem / (Vjunc * Frdy) + J_na_juncsl / Vjunc * (y[33 - 1] - y[32 - 1]) - ydot[16];
	ydot[32] = -I_Na_tot_sl * Cmem / (Vsl * Frdy) + J_na_juncsl / Vsl * (y[32 - 1] - y[33 - 1]) + J_na_slmyo / Vsl * (y[33] - y[33 - 1]) - ydot[17];
	ydot[33] = (J_na_slmyo / Vmyo * (y[33 - 1] - y[33])) * (1 - para.Na_clamp); // [mM/msec]

	// Potassium Concentration
	I_K_tot = I_to + I_kr + I_ks + I_ki - 2 * I_nak + I_CaK + I_kp + I_kur + I_kach + I_k2p + I_sk; // [uA/uF] //SVP: added IKur
	ydot[34] = 0; // -I_K_tot*Cmem/(Vmyo*Frdy);         // [mM/msec]

	// Calcium Concentrations
	I_Ca_tot_junc = I_Ca_junc + I_cabk_junc + I_pca_junc - 2 * I_ncx_junc;           // [uA/uF]
	I_Ca_tot_sl = I_Ca_sl + I_cabk_sl + I_pca_sl - 2 * I_ncx_sl;    // [uA/uF]
	ydot[35] = -I_Ca_tot_junc * Cmem / (Vjunc * 2 * Frdy) + J_ca_juncsl / Vjunc * (y[36] - y[35]) - J_CaB_junction + (J_SRCarel) * Vsr / Vjunc + J_SRleak * Vmyo / Vjunc; // Ca_j
	ydot[36] = -I_Ca_tot_sl * Cmem / (Vsl * 2 * Frdy) + J_ca_juncsl / Vsl * (y[35] - y[36]) + J_ca_slmyo / Vsl * (y[37] - y[36]) - J_CaB_sl; // Ca_sl
	ydot[37] = -J_serca * Vsr / Vmyo - J_CaB_cytosol + J_ca_slmyo / Vmyo * (y[36] - y[37]); // [mM/msec]


	if ( para.Ca_clamp == 1 || para.Ca_clamp == 2) { // Ca_clamp, BAPTA
		ydot[35] = 0; ydot[36] = 0; ydot[37] = 0;
	} else if (para.Ca_clamp == 3) { // EGTA
		ydot[37] = 0;
	}



	// //// Simulation type
	// switch p(1)
	//     case 0          // no stimulation
	//         I_app = 0;
	//     case 1          // pace w/ current injection at rate 'rate' (Hz)
	//         rate = p(2); // Hz
	//         period = 1000/rate; // ms
	//         if mod(t,period) <= 5
	//             I_app = 12.5;
	//         else
	//             I_app = 0.0;
	//         end
	//     case 2      // ERP
	//         rate = p(2); // Hz
	//         rec_interval = p(3);
	//         if t <= 5
	//             I_app = 12.5;
	//         elseif t > 5 && t <= rec_interval
	//             I_app = 0.0;
	//         elseif t > rec_interval && t <= rec_interval+5
	//             if rate == 0.5 && AF == 0
	//                 DTE = 12.5*0.125;
	//             else
	//                 DTE = 12.5*0.2;
	//             end
	//             I_app = 2*DTE;
	//         else
	//             I_app = 0.0;
	//         end
	//     case 3      // Voltage step
	//         rate = p(2); // Hz
	//         period = 1000/rate; // ms
	//         step_duration = p(3);
	//         V_test = p(4);
	//         V_hold1 = -80; T_hold1 = 5;
	//         V_hold2 = V_test; T_hold2 = step_duration;
	//         if mod(t,period) <= T_hold1 //#ok<ALIGN>
	//             V_clamp = V_hold1;
	//         elseif mod(t,period) > T_hold1 && mod(t,period) <= T_hold1+T_hold2
	//             V_clamp = V_hold2;
	//         else
	//             V_clamp = V_hold1;
	//         end
	//         R_clamp = 0.02;
	//         I_app = (V_clamp-y[39-1])/R_clamp;
	// end

	//// Membrane Potential
	I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;
	I_Cl_tot = I_ClCa + I_Clbk + I_ClCFTR;
	I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
	I_tot = I_Na_tot + I_Cl_tot + I_Ca_tot + I_K_tot; // [uA/uF]
	ydot[38] = -(I_tot - I_app);
	double vmax = ydot[38];


	dV = ydot[38];
	CaSR = y[30];
	Caj = y[35];
	Casl = y[36];
	Cai = y[37];
	Nai = y[33];
	Naj = y[32];
	Nasl = y[31];



	// dV = ydot[38];
	//// Output adjustment depending on the function call
	// if (nargin == 3)
	//     output = ydot;
	// elseif (nargin == 4) && strcmp(runType,'ydot')
	//     output = ydot;
	// elseif (nargin == 4) && strcmp(runType,'rates')
	//     output = r;
	// elseif (nargin == 4) && strcmp(runType,'currents')
	//     //currents = [I_Na I_nabk I_nak I_kr I_ks I_kp I_tos I_tof I_ki I_ClCa I_Clbk I_Catot I_ncx I_pca I_cabk J_serca*Vmyo/Vsr];
	//     //currents = [I_Na I_tof I_tos I_kr I_ks I_ClCa I_Catot J_SRCarel*Vsr/Vmyo J_SRleak RI I_ncx];
	//     //currents = [I_Na I_Catot I_ncx I_nak I_kr I_ks I_tof I_tos I_ki];
	//     currents = [vmax J_serca*Vsr/Vmyo -I_ncx*2*Cmem/(2*Frdy*Vmyo) I_pca*Cmem/(2*Frdy*Vmyo)];
	//     output = currents;
	// end
}


void HAM_ECC::initialiser() {
	y[0]      = 0.0040891;
	y[1]      = 0.9449988;
	y[2]      = 0.9622242;
	y[3]      = 0.0000084;
	y[4]      = 0.9994133;
	y[5]      = 0.0515998;
	y[6]      = 0.0362242;
	y[7]      = 0.0002035;
	y[8]      = 0.9567656;
	y[9]      = 0.0008244;
	y[10]     = 0.9662879;
	y[11]     = 0.0025541;
	y[12]     = 0.0050455;
	y[13]     = 0.7883706;
	y[14]     = 0.0000028;
	y[15]     = 0.0000007;
	y[16]     = 3.7422141;
	y[17]     = 0.8166018;
	y[18]     = 0.0197666;
	y[19]     = 0.1290834;
	y[20]     = 0.0051135;
	y[21]     = 0.0007777;
	y[22]     = 0.0044295;
	y[23]     = 0.1350567;
	y[24]     = 0.0048088;
	y[25]     = 0.0156659;
	y[26]     = 0.0237294;
	y[27]     = 0.1109636;
	y[28]     = 0.2020127;
	y[29]     = 1.1559742;
	y[30]     = 0.5203395;
	y[31]     = 9.7950634;
	y[32]     = 9.7950461;
	y[33]     = 9.7952094;
	y[34]     = 120.0000000;
	y[35]     = 0.0003778;
	y[36]     = 0.0002588;
	y[37]     = 0.0002343;
	y[38]     = -79.1012084;
	y[39]     = 0.0040891;
	y[40]     = 0.1040307;
	y[41]     = 0.2013751;
	// y[42]     = 1;   // LTCC 7 state - c2
	// y[43]     = 0;
	// y[44]     = 0;
	// y[45]     = 0;
	// y[46]     = 0;
	// y[47]     = 0;
	// y[48]     = 0;
	// y[49]     = 0;
	// y[42 + 8] = 1; // LTCC 7 state - c2
	// y[43 + 8] = 0;
	// y[44 + 8] = 0;
	// y[45 + 8] = 0;
	// y[46 + 8] = 0;
	// y[47 + 8] = 0;
	// y[48 + 8] = 0;
	// y[49 + 8] = 0;

	y[0]     = 0.004702304605;  //0.00399229398;// 0.00445640936;// 0.004080635368;
	y[1]     = 0.9520451134;  //0.9467007117;// 0.9384696031;// 0.9451285505;
	y[2]     = 0.9677505672;  //0.9636109046;// 0.9568180854;// 0.9622833438;
	y[3]     = 0.05503093515;  //0.05460303312;// 0.08131435121;// 8.413437798e-06;
	y[4]     = 0.04592216039;  //0.04688133821;// 0.06740567656;// 0.999413341;
	y[5]     = 0.05762842186;  //0.05762842186;// 0.05762842186;// 0.05762842186;
	y[6]     = 0.03955197377;  //0.03955197377;// 0.03955197377;// 0.03955197377;
	y[7]     = 0.0001934925286;  //0.0002000412669;// 0.0002163330235;// 0.0002032347954;
	y[8]     = 0.9497013683;  //0.9464942804;// 0.9115429599;// 0.9372571103;
	y[9]     = 0.0008245759368;  //0.0008051324621;// 0.0007663098908;// 0.0008234263476;
	y[10]    = 0.966471121;  //0.9671718763;// 0.9686641242;// 0.9663079588;
	y[11]    = 0.002133688834;  //0.002772230208;// 0.01610936656;// 0.003206591541;
	y[12]    = 0.005041386401;  //0.005041386401;// 0.005041386401;// 0.005041386401;
	y[13]    = 0.7690079415;  //0.7937742919;// 0.7916477886;// 0.7834945442;
	y[14]    = 2.954611248e-06;  //3.082711693e-06;// 6.497008409e-06;// 3.83997118e-06;
	y[15]    = 8.87481572e-07;  //8.008856376e-07;// 1.709867405e-06;// 1.06108702e-06;
	y[16]    = 3.898609958;  //3.58731009;// 3.365451293;// 3.800037507;
	y[17]    = 0.8506999403;  //0.782753324;// 0.7343334096;// 0.8292064491;
	y[18]    = 0.01994760791;  //0.02002277168;// 0.01223669363;// 0.02106415581;
	y[19]    = 0.1373865351;  //0.1287392858;// 0.1305759548;// 0.1308372858;
	y[20]    = 0.005777523224;  //0.005288726577;// 0.004318842605;// 0.00425980829;
	y[21]    = 0.0008473481432;  //0.0008473481432;// 0.0008473481432;// 0.0008473481432;
	y[22]    = 0.00445073063;  //0.004298289724;// 0.00489516227;// 0.005291039036;
	y[23]    = 0.1447239502;  //0.1351923896;// 0.1345562962;// 0.1341852704;
	y[24]    = 0.004856136738;  //0.004873380054;// 0.00383671956;// 0.005122289046;
	y[25]    = 0.01862086428;  //0.0161194394;// 0.01684246011;// 0.01752575393;
	y[26]    = 0.02553080338;  //0.02401221582;// 0.02005913571;// 0.02586842913;
	y[27]    = 0.1271215045;  //0.1124443249;// 0.1146257456;// 0.116599258;
	y[28]    = 0.2205195367;  //0.2033109223;// 0.1838000268;// 0.2116176743;
	y[29]    = 1.208931059;  //1.190765872;// 1.222020357;// 1.203914851;
	y[30]    = 0.5044331601;  //0.5492331281;// 0.576433667;// 0.5605281892;
	y[31]    = 10.6396871;  //9.76 + 2; //0.5 + 1.0 + 9.023312821; // 8.012923452;// 9;// 10.10011682;
	y[32]    = 10.63910564;  //9.76 + 2; //0.5 + 1.0 + 9.022782104; // 8.012877676;// 9;// 10.09974519;
	y[33]    = 10.63929166;  //9.76 + 2; //0.5 + 1.0 + 9.022928347; // 8.01319175;// 9;// 10.09990903;
	y[34]    = 140;  //140;//120;//120;// 120;// 120;
	y[35]    = 0.0004049693237;  //0.0003890206036;// 0.0004070258914;// 0.0004240838797;
	y[36]    = 0.0002503917695;  //0.0002619613163;// 0.0002181008294;// 0.0002826454379;
	y[37]    = 0.0002165906556;  //0.0002387248512;// 0.0001730660506;// 0.0002559818175;
	y[38]    = -79.54711902233;  //-79.24833461;// -78.57403977;// -79.11428795;
	y[39]    = 0.004702304605;  //0.00399229398;// 0.00445640936;// 0.004080635368;
	y[40]    = 0.1222806279;  //0.1216252951;// 0.1069020548;// 0.1179555119;
	y[41]    = 0.2013359986;  //0.2013608187;// 0.2014278254;// 0.2013739322;
	y[42]    = 0.0001832041026;  //0.0001891268306;// 0.0002045421467;// 0;
	y[43]    = 0.995757261;  //0.9954683003;// 0.9899943626;// 1;
	y[44]    = 0.9597793584;  //0.9566556703;// 0.9227884593;// 1;


	y[0]     = 0.003597949497;//
	y[1]     = 0.9581027507;//
	y[2]     = 0.972373143;//
	y[3]     = 0.0512017168;//
	y[4]     = 0.0441441933;//
	y[5]     = 0.05762842186;//
	y[6]     = 0.03955197377;//
	y[7]     = 0.0001857903489;//
	y[8]     = 0.9477756353;//
	y[9]     = 0.0008035046359;//
	y[10]    = 0.9684326908;//
	y[11]    = 0.002648425288;//
	y[12]    = 0.005041386401;//
	y[13]    = 0.7896888079;//
	y[14]    = 2.183994596e-06;//
	y[15]    = 5.81637301e-07;//
	y[16]    = 3.765452905;//
	y[17]    = 0.8216421143;//
	y[18]    = 0.01777753152;//
	y[19]    = 0.1274721386;//
	y[20]    = 0.005877212311;//
	y[21]    = 0.0008473481432;//
	y[22]    = 0.003812209607;//
	y[23]    = 0.1356752147;//
	y[24]    = 0.00432988697;//
	y[25]    = 0.01462086751;//
	y[26]    = 0.02105071921;//
	y[27]    = 0.1074860636;//
	y[28]    = 0.1888016889;//
	y[29]    = 1.124895027;//
	y[30]    = 0.4956813367;//
	y[31]    = 9.914946512;//
	y[32]    = 9.914708496;//
	y[33]    = 9.914902042;//
	y[34]    = 140;//
	y[35]    = 0.0003518859336;//
	y[36]    = 0.0002290903651;//
	y[37]    = 0.0002031905044;//
	y[38]    = -79.88436365;//
	y[39]    = 0.003597949497;//
	y[40]    = 0.1300923401;//
	y[41]    = 0.2013004815;//
	y[42]    = 0.0001857283065;//
	y[43]    = 0.9954379337;//
	y[44]    = 0.9576758318;//


	// y[45]    = 0.9895726818;
	// y[43 + 3]    = 1.690746702e-05;
	// y[44 + 3]    = 4.240266403e-10;
	// y[45 + 3]    = 0.0003644860657;
	// y[46 + 3]    = 0.006519349275;
	// y[47 + 3]    = 0.003526009982;
	// y[48 + 3]    = 2.134724476e-15;
	// y[49 + 3]    = 1.299456512e-07;
	// y[45 + 8] = 0.9871462467;
	// y[43 + 3 + 8] = 1.897481412e-05;
	// y[44 + 3 + 8] = 4.090713291e-10;
	// y[45 + 3 + 8] = 0.0003628823657;
	// y[46 + 3 + 8] = 0.008093821615;
	// y[47 + 3 + 8] = 0.004377438036;
	// y[48 + 3 + 8] = 2.125774633e-15;
	// y[49 + 3 + 8] = 1.293105713e-07;

	// states for LTCC and LTCC_sl
	y[45] = 1;
	y[45 + 10] = 1;
	for (int i = 1; i < 10; ++i)
	{
		y[45 + i] = 0;
		y[45 + 10 + i] = 0;
	}

	// states for LTCC_m2 and LTCC_sl_m2
	y[45 + 20] = 1;
	y[45 + 10 + 20] = 1;
	for (int i = 1; i < 10; ++i)
	{
		y[45 + 20 + i] = 0;
		y[45 + 10 + 20 + i] = 0;
	}


	y[45 + 40 + 0] = 1;
	y[45 + 40 + 1] = 0;
	y[45 + 40 + 2] = 0;
	y[45 + 40 + 3] = 0;
	y[45 + 40 + 4] = 0;
	y[45 + 40 + 5] = 0;
	y[45 + 40 + 6] = 0;


	V = y[38];
}



void HAM_ECC::print_to_file(double t, std::ofstream & output_file) {


	output_file <<  std::setprecision(7)
	            << t << " "  << std::setprecision(5)// 1
	            << V << " "
	            << dV << " "
	            << I_app << " " << std::setprecision(7)
	            << CaSR << " "  // 5
	            << Caj << " " // 6
	            << Casl << " "
	            << Cai << " "  << std::setprecision(5)
	            << I_Na << " "  // 9
	            << I_NaL << " "  // 10
	            << I_nabk << " "
	            << I_nak << " "
	            << I_to << " "
	            << I_kur << " "  // 14
	            << I_k2p << " " //15
	            << I_ks  << " "
	            << I_ki << " "
	            << I_kp << " "
	            << I_kach << " "  // 19
	            << I_sk << " "   //20
	            << I_ClCa << " "
	            << I_Clbk << " "
	            << I_ClCFTR << " "
	            << I_Ca << " "  // 24
	            << I_cabk << " "  // 25
	            << I_pca << " "
	            << I_ncx << " "
	            << J_SRCarel << " "
	            << J_serca << " "  // 29
	            << I_Ca_tot_junc << " " // 30
	            << I_Ca_tot_sl << " " // 31
	            << J_CaB_junction << " " //32
	            << Nai << " "  // 33
	            << Naj << " "  // 34
	            << Nasl << " "  // 35
	            << I_kr << " "  // 36
	            << I_Catot << " " // 37
	            << J_SRleak << " " // 38
	            << LTCC_mode_1_out << " "  // 39
	            << LTCC_mode_2_out << " "  // 40
	            << std::endl;
}



void HAM_ECC::print_to_file_Vm_Ca(double t, std::ofstream & output_file) {


	output_file <<  std::setprecision(7)
	            << t << " "  << std::setprecision(5)// 1
	            << V << " "
	            << dV << " "
	            << I_app << " " << std::setprecision(5)
	            << CaSR << " "  // 5
	            << 1000*Caj << " " // 6  // change unit from mM to uM
	            << 1000*Casl << " "
	            << 1000*Cai << " "  << std::setprecision(4)
	            << I_Na << " "  // 9
	            // << I_NaL << " "  // 10
	            // << I_nabk << " "
	            // << I_nak << " "
	            // << I_to << " "
	            // << I_kur << " "  // 14
	            // << I_k2p << " " //15
	            // << I_ks  << " "
	            // << I_ki << " "
	            // << I_kp << " "
	            // << I_kach << " "  // 19
	            // << I_sk << " "   //20
	            // << I_ClCa << " "
	            // << I_Clbk << " "
	            // << I_ClCFTR << " "
	            << I_Ca << " "  // 24
	            // << I_cabk << " "  // 25
	            // << I_pca << " "
	            << I_ncx << " "
	            << J_SRCarel << " "
	            // << J_serca << " "  // 29
	            // << I_Ca_tot_junc << " " // 30
	            // << I_Ca_tot_sl << " " // 31
	            // << J_CaB_junction << " " //32
	            << Nai/*(Nai + Naj + Nasl)/3.0*/ << " "  // 33
	            // << Naj << " "  // 34
	            // << Nasl << " "  // 35
	            // << I_kr << " "  // 36
	            // << I_Catot << " " // 37
	            // << J_SRleak << " " // 38
	            // << LTCC_mode_1_out << " "  // 39
	            // << LTCC_mode_2_out << " "  // 40
	            << std::endl;
}




void HAM_ECC::print_to_file_Vm_only(double t, std::ofstream & output_file) {


	output_file <<  std::setprecision(7)
	            << t << " "  // 1
	            << V << " "
	            // << dV << " "
	            // << I_app << " "
	            // << CaSR << " "  // 5
	            // << Caj << " " // 6
	            // << Casl << " "
	            // << Cai << " "  << std::setprecision(3)
	            // << I_Na << " "  // 9
	            // << I_NaL << " "  // 10
	            // << I_nabk << " "
	            // << I_nak << " "
	            // << I_to << " "
	            // << I_kur << " "  // 14
	            // << I_k2p << " " //15
	            // << I_ks  << " "
	            // << I_ki << " "
	            // << I_kp << " "
	            // << I_kach << " "  // 19
	            // << I_sk << " "   //20
	            // << I_ClCa << " "
	            // << I_Clbk << " "
	            // << I_ClCFTR << " "
	            // << I_Ca << " "  // 24
	            // << I_cabk << " "  // 25
	            // << I_pca << " "
	            // << I_ncx << " "
	            // << J_SRCarel << " "
	            // << J_serca << " "  // 29
	            // << I_Ca_tot_junc << " " // 30
	            // << I_Ca_tot_sl << " " // 31
	            // << J_CaB_junction << " " //3
	            // << Nai << " "  // 33
	            // << Naj << " "  // 34
	            // << Nasl << " "  // 35
	            // << I_kr << " "  // 36
	            // << I_Catot << " " // 37
	            // << J_SRleak << " " // 38
	            << std::endl;
}
