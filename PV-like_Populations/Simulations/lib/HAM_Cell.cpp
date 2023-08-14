

#include <cmath>

#include "HAM_Cell.hpp"

#include  <iomanip>


int human_atrial_ECC_ODE(double t, double *y, double *ydot, Cell * cell_par)
{
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
	//epi = 1; // EPI or ENDO?
	int RA = 0; // Right ATRIUM
	int AF = 0; // p[5]; // AF
	int ISO = 0;//p[6]; // ISO (boolean)
	double CCh = 0;//p(7); // [uM]
	double ACh = 0;//p(8); // [uM]
	int K2P_cond = 1;//p(9);
	int SK_cond = 1;//p(10);
	int SK_shift = 0;// p(11);
	int Ca_clamp = 0;//p(12);
	int Na_clamp = 0; //p(13);
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

	// Courtemanche et al, 1998
	//GNa = SA_par(1)*23*(1-0.1*AF);  // [mS/uF] Grandi
	//GNa = SA_par(1)*7.8*(1-0.1*AF);  // [mS/uF] Courtemanche
	double GNa = SA_par[0] * 9 * (1 - 0.1 * AF); // [mS/uF] Courtemanche (adjusted for dV/dt max) MOD1

	double am = (y[39 - 1] == -47.13) ? 3.2 : (0.32 * (y[39 - 1] + 47.13) / (1 - exp( -0.1 * (y[39 - 1] + 47.13))));

	double bm = 0.08 * exp(-y[39 - 1] / 11.0);

	double ah = (y[39 - 1] >= -40) ? 0 : 0.135 * exp( -(y[39 - 1] + 80) / 6.8 );

	double bh = (y[39 - 1] >= -40) ? (1.0 / (0.13 * (1 + exp( -(y[39 - 1] + 10.66) / 11.1 )))) : (3.56 * exp( 0.079 * y[39 - 1]) + 3.1e5 * exp(0.35 * y[39 - 1]));

	double aj = (y[39 - 1] >= -40) ? 0 : ( (-127140 * exp(0.2444 * y[39 - 1]) - 3.474e-5 * exp(-0.04391 * y[39 - 1])) * (y[39 - 1] + 37.78)) / (1 + exp( 0.311 * (y[39 - 1] + 79.23) ) );
	double bj = (y[39 - 1] >= -40) ? ((0.3 * exp(-2.535e-7 * y[39 - 1])) / (1 + exp( -0.1 * (y[39 - 1] + 32) ))) : (0.1212 * exp( -0.01052 * y[39 - 1] )) / (1 + exp( -0.1378 * (y[39 - 1] + 40.14) ));

	ydot[0] = am * (1 - y[0]) - bm * y[0];
	ydot[1] = ah * (1 - y[1]) - bh * y[1];
	ydot[2] = aj * (1 - y[2]) - bj * y[2];

	cell_par->I_Na_junc = Fjunc * GNa * y[0] * y[0] * y[0] * y[1] * y[2] * (y[39 - 1] - ena_junc);
	cell_par->I_Na_sl = Fsl * GNa * y[0] * y[0] * y[0] * y[1] * y[2] * (y[39 - 1] - ena_sl);
	cell_par->I_Na = cell_par->I_Na_junc + cell_par->I_Na_sl; //#ok<NASGU>

	//// Late I_Na
	//GNaL = SA_par(2)*0.0025*AF; // [mS/uF]
	double GNaL = SA_par[1] * (0.002 + 0.0015 * AF) * GNa; // [mS/uF] (adjusted) MOD1  // increased GNaL here // 18:29:43, Thu, 07-February-2019, By Haibo

	am = 0.32 * (y[39 - 1] + 47.13) / (1 - exp(-0.1 * (y[39 - 1] + 47.13)));
	bm = 0.08 * exp(-y[39 - 1] / 11);
	double inf = 1 / (1 + exp((y[39 - 1] + 91) / 6.1));
	double tau = 200;
	ydot[39] = am * (1 - y[39]) - bm * y[39];
	ydot[40] = (inf - y[40]) / tau;

	cell_par->I_NaL_junc = Fjunc * GNaL * y[39] * y[39] * y[39] * y[40] * (y[38] - ena_junc);
	cell_par->I_NaL_sl = Fsl * GNaL * y[39] * y[39] * y[39] * y[40] * (y[38] - ena_sl);
	cell_par->I_NaL = cell_par->I_NaL_junc + cell_par->I_NaL_sl; //#ok<NASGU>

	//// I_nabk: Na Background Current
	double GNaB = SA_par[2] * 0.597e-3; // [mS/uF]

	cell_par->I_nabk_junc = Fjunc * GNaB * (y[39 - 1] - ena_junc);
	cell_par->I_nabk_sl = Fsl * GNaB * (y[39 - 1] - ena_sl);
	cell_par->I_nabk = cell_par->I_nabk_junc + cell_par->I_nabk_sl; //#ok<NASGU>

	//// I_nak: Na/K Pump Current
	double IbarNaK = SA_par[3] * 1.26; // [uA/uF]
	//IbarNaK = SA_par(4)*1.25*(1.26); // NEW K2P +25// in Schmidt et al 2015
	double KmNaip = 11 * (1 - 0.25 * ISO); // [mM]
	double KmKo = 1.5; // [mM]
	//Q10NaK = 1.63;
	//Q10KmNai = 1.39;

	double sigma = (exp(Nao / 67.3) - 1) / 7.0;
	double fnak = 1 / (1 + 0.1245 * exp(-0.1 * y[39 - 1] * FoRT) + 0.0365 * sigma * exp(-y[39 - 1] * FoRT));

	cell_par->I_nak_junc = Fjunc * IbarNaK * fnak * Ko / (1 + pow((KmNaip / y[32 - 1]), 4)) / (Ko + KmKo);
	cell_par->I_nak_sl = Fsl * IbarNaK * fnak * Ko / (1 + pow((KmNaip / y[33 - 1]), 4)) / (Ko + KmKo);
	cell_par->I_nak = cell_par->I_nak_junc + cell_par->I_nak_sl;

	//// I_to: Transient Outward K Current
	double GtoFast = 1.0 * SA_par[4] * (1 - 0.7 * AF) * 0.165;
	//GtoFast = SA_par(5)*0.4*(1-0.7*AF)*0.165; // NEW K2P -60// in Schmidt et al 2015

	double xtoss = 1 / ( 1 + exp( -(y[39 - 1] + 1) / 11 ) );
	double tauxtof = 3.5 * exp(-((y[39 - 1] / 30) * (y[39 - 1] / 30))) + 1.5;
	double ytoss = 1 / ( 1 + exp( (y[39 - 1] + 40.5) / 11.5) ) ;
	double tauytof = 25.635 * exp(-(((y[39 - 1] + 52.45) / 15.8827) * ((y[39 - 1] + 52.45) / 15.8827))) + 24.14; //14.14

	ydot[9] = (xtoss - y[9]) / tauxtof;
	ydot[10] = (ytoss - y[10]) / tauytof;

	cell_par->I_tof = GtoFast * y[9] * y[10] * (y[39 - 1] - ek);
	cell_par->I_to = cell_par->I_tof;

	//// I_kr: Rapidly Activating K Current
	double gkr = SA_par[5] * 0.035 * sqrt(Ko / 5.4);

	inf = 1 / (1 + exp(-(y[39 - 1] + 10) / 5));
	tau = 550 / (1 + exp((-22 - y[39 - 1]) / 9)) * 6 / (1 + exp((y[39 - 1] - (-11)) / 9)) + 230 / (1 + exp((y[39 - 1] - (-40)) / 20));
	ydot[11] = (inf - y[11]) / tau;
	double rkr = 1 / (1 + exp((y[39 - 1] + 74.0) / 24.0));
	cell_par->I_kr = gkr * y[11] * rkr * (y[39 - 1] - ek);

	//// I_ks: Slowly Activating K Current
	double gks_junc = 2*SA_par[6] * (1 + 1 * AF + 2 * ISO) * 0.0035;
	double gks_sl = 2*SA_par[6] * (1 + 1 * AF + 2 * ISO) * 0.0035;
	double pNaK = 0.01833;
	double eks = (1 / FoRT) * log((Ko + pNaK * Nao) / (y[35 - 1] + pNaK * y[34 - 1]));

	inf = 1 / (1 + exp(-(y[39 - 1] + 40 * ISO + 3.8) / 14.25));
	tau = 990.1 / (1 + exp(-(y[39 - 1] + 40 * ISO + 2.436) / 14.12));
	ydot[12] = (inf - y[12]) / tau;

	cell_par->I_ks_junc = Fjunc * gks_junc * y[12] * y[12] * (y[39 - 1] - eks);
	cell_par->I_ks_sl = Fsl * gks_sl * y[12] * y[12] * (y[39 - 1] - eks);
	cell_par->I_ks = cell_par->I_ks_junc + cell_par->I_ks_sl;

	//// I_kur: Ultra-Rapid Delayed Rectifier Outward K Current
	double Gkur = 1.8 * SA_par[7] * (1 - 0.5 * AF) * (1 + 2 * ISO) * 0.045 * (1 + 0.2 * RA); //1.5* by haibo
	//Gkur = SA_par(8)*(1+0.2*AF)*(1-0.5*AF)*(1+2*ISO)*0.045*(1+0.2*RA); // NEW K2P +20// cAF in Schmidt et al 2015

	double xkurss = 1 / ( 1 + exp( (y[39 - 1] + 6) / -8.6) );
	double tauxkur = 9.0 / (1 + exp((y[39 - 1] + 5) / 12.0)) + 0.5;
	double ykurss = 1 / ( 1 + exp( (y[39 - 1] + 7.5) / 10) );
	double tauykur = 590.0 / (1 + exp((y[39 - 1] + 60) / 10)) + 3050.0;

	ydot[7] = (xkurss - y[7]) / tauxkur;
	ydot[8] = (ykurss - y[8]) / tauykur;
	cell_par->I_kur = Gkur * y[7] * y[8] * (y[39 - 1] - ek);

	// NEw IKur model here; // 16:37:09, Mon, 25-February-2019, By Haibo
	double const CNZ_gkur =1.* 0.006398 * 0.9;
	// double y[42] = y[42];
	// double y[43] = y[43];
	// double y[44] = y[44];
	// IKur = het.GKur * Cm * CNZ_gkur * (4.5128 + 1.899769 / (1.0 + exp((V - 20.5232) / (-8.26597)))) * 0.5 * y[42] * y[43] * (V - Ek);
	//Drug block....
	// double d_IKur_BO = het.drug_IKur_KO * exp(-het.drug_IKur_ZKO * V * F / R / T) * het.drug_IKur_concen * y[43] * y[42] * (1 - IKur_BO - IKur_BC) - het.drug_IKur_LO * IKur_BO * exp(-het.drug_IKur_ZLO * V * F / R / T);
	// double d_IKur_BC = het.drug_IKur_KC * exp(-het.drug_IKur_ZKC * V * F / R / T) * het.drug_IKur_concen * y[43] * (1 - y[42]) * (1 - IKur_BO - IKur_BC) - het.drug_IKur_LC * IKur_BC * exp(-het.drug_IKur_ZLC * V * F / R / T);
	// IKur_BO += dt * d_IKur_BO;
	// IKur_BC += dt * d_IKur_BC;
	// IKur = GKur * IKur;

	/*inf = 1 / (1 + exp(-(V + IKur_ac1_shift + 6) / (8.6 * IKur_ac1_grad)));
	tau = (0.009 / (1.0 + exp((V + 5) / 12)) + 0.0005);
	isusr = inf + (isusr - inf) * exp(-(dt / 1000) / tau);

	inf = 1 / (1 + exp((V + IKur_inac_shift + 7.5) / (10 * IKur_inac_grad)));
	tau = (0.59 / (1 + exp((V + 60.0) / 10)) + 3.05);
	isuss = inf + (isuss - inf) * exp(-(dt / 1000) / tau);*/

	double K_Q10 = 3.0;//3.5308257834747638;
	//y[42]
	// inf = ((het.IKur_ac1_mult * 1.0) / (1 + exp((V + 17.6684 + het.IKur_ac1_shift) / (-5.75418 * het.IKur_ac1_grad))) ) * ((het.IKur_ac2_mult * 1.0) / (1 + exp((V + 8.4153 + het.IKur_ac2_shift) / (-11.51037561 * het.IKur_ac2_grad)))) + het.IKur_ac_add;


	inf = 1.0 / (1 + exp(-(y[38] + 5.52) / (8.6))); // change to simple version, haibo.  removed IKur mutations here.
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

	cell_par->I_kur = IKur;
	//// I_ki: Time-Independent K Current
	double Gki = SA_par[8] * (1 + 1 * AF) * (2.1 * 0.0525) * sqrt(Ko / 5.4);
	//Gki = SA_par(9)*(1+0.45*(1-AF)+0.70*AF)*(1+1*AF)*(2.1*0.0525)*sqrt(Ko/5.4); // NEW K2P +45// nSR/70// cAF in Schmidt et al 2015

	// Na-dependence (Voigt-Heijman et al 2013)
	double aki_sl = (0.1 + 0.9 / (1 + (y[33 - 1] / 7) * (y[33 - 1] / 7.0))) * 1 / (1 + exp(0.2385 * (y[39 - 1] - ek - 59.215)));
	double aki_j = (0.1 + 0.9 / (1 + (y[32 - 1] / 7) * (y[32 - 1] / 7.0))) * 1 / (1 + exp(0.2385 * (y[39 - 1] - ek - 59.215)));
	// Na-dependence (Schmidt et al 2015)
	// aki_sl = 1/(1+(y[33-1]/7)^2);
	// aki_j = 1/(1+(y[32-1]/7)^2);

	double bki = (0.49124 * exp(0.08032 * (y[39 - 1] + 5.476 - ek)) + exp(0.06175 * (y[39 - 1] - ek - 594.31))) / (1 + exp(-0.5143 * (y[39 - 1] - ek + 4.753)));
	double kiss_sl = aki_sl / (aki_sl + bki);
	double kiss_j = aki_j / (aki_j + bki);

	cell_par->I_ki_sl = Fsl * Gki * kiss_sl * (y[39 - 1] - ek);
	cell_par->I_ki_j = Fjunc * Gki * kiss_j * (y[39 - 1] - ek);
	cell_par->I_ki = cell_par->I_ki_j + cell_par->I_ki_sl;

	//// I_kp: Plateau K current
	double gkp = SA_par[9] * 0.002;

	double kp_kp = 1 / (1 + exp(7.488 - y[39 - 1] / 5.98));
	cell_par->I_kp = gkp * kp_kp * (y[39 - 1] - ek);

	//// I_k2p: K2P3.1 K current (from Schmidt et al 2015)
	//gk2p = SA_par(11)*K2P_cond*(0.005+0.014*AF); // wrong
	//gk2p = SA_par(11)*K2P_cond*0.0050*(1+(0.0145/0.005-1)*AF); // 2.9-fold increase in AF
	double gk2p = SA_par[10] * K2P_cond * 0.0050 * (1 + (SA_par[23] * 0.0145 / 0.005 - 1) * AF);

	inf = 0.2 + 0.8 / (1 + exp(-(y[39 - 1] - 10 + 15 * AF) / 14));
	tau = 2 + 40 / (1 + exp((y[39 - 1] + 25) * (y[39 - 1] + 25) / 80));
	ydot[41] = (inf - y[41]) / tau;

	cell_par->I_k2p = gk2p * y[41] * (y[39 - 1] - ek); // NEW K2P





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

	cell_par->I_kach_sl = Fsl * Gkach_sl * d_kach * r_kach * (y[39 - 1] - ek);
	cell_par->I_kach_j = Fjunc * Gkach_j * d_kach * r_kach * (y[39 - 1] - ek);

	cell_par->I_kach = cell_par->I_kach_sl + cell_par->I_kach_j;

	//// I_sk: Small-Conductance Ca-Activated K Current
	double par_sk[] = {0.0506381114404388, 0.273335569451572, 2.96381060498817, 0.199981221802789, 0.279328126521496, -86.9289059836381, 0.00636311816933264, 5.22915055145375};

	//gsk = SA_par(13)*SK_cond*par_sk(1)*(1+(par_sk(8)-1)*AF);
	double gsk = SA_par[12] * SK_cond * par_sk[0] * (1 + (SA_par[24] * par_sk[7] - 1) * AF);

	double kdsk = SA_par[13] * (pow(10, (SK_shift - 3.45))); // kdsk = (10^(ISK_shift-3.3)); // (mM)
	double gsk_ca_junc = 1.0 / (1 + exp((log10(kdsk) - log10(y[35])) / 0.3));
	double gsk_ca_sl = 1.0 / (1 + exp((log10(kdsk) - log10(y[36])) / 0.3));
	//     // Ca clamp for SK channels
	//     gsk_ca = 500e-6;
	//     gsk_ca_junc = 1/(1+ exp((log10(kdsk)-log10(gsk_ca))/0.3));
	//     gsk_ca_sl = 1/(1+ exp((log10(kdsk)-log10(gsk_ca))/0.3));

	double gsk_vm = par_sk[1] / (1 + exp((y[39 - 1] - ek + par_sk[2]) * par_sk[3])) + par_sk[4] / (1 + exp((-(y[39 - 1] - ek + par_sk[5]) * par_sk[6])));

	cell_par->I_sk_junc = Fjunc * gsk * gsk_ca_junc * gsk_vm * (y[39 - 1] - ek);
	cell_par->I_sk_sl = Fsl * gsk * gsk_ca_sl * gsk_vm * (y[39 - 1] - ek);
	cell_par->I_sk = cell_par->I_sk_junc + cell_par->I_sk_sl;

	//// I_ClCa and I_Clbk: Ca-activated and Background Cl Currents
	double GClCa = SA_par[14] * 0.0548;   // [mS/uF]
	double KdClCa = 100e-3;               // [mM]
	//GClB = SA_par(16)*9e-3;        // [mS/uF]
	double GClB = SA_par[15] * 0.5 * 9e-3;    // [mS/uF] MOD1
	double GClCFTR = 0; // 4.9e-3*ISO;     // [mS/uF]

	cell_par->I_ClCa_junc = Fjunc * GClCa / (1 + KdClCa / y[35]) * (y[39 - 1] - ecl);
	cell_par->I_ClCa_sl = Fsl * GClCa / (1 + KdClCa / y[37 - 1]) * (y[39 - 1] - ecl);
	cell_par->I_ClCa = cell_par->I_ClCa_junc + cell_par->I_ClCa_sl;

	cell_par->I_Clbk = GClB * (y[39 - 1] - ecl);

	cell_par->I_ClCFTR = GClCFTR * (y[39 - 1] - ecl);

	//// I_Ca: L-type Ca Current
	double pNa = 0 * SA_par[16] * (1 + 0.5 * ISO) * (1 - 0.5 * AF) * 0.75e-8; // [cm/sec]
	double pCa = 0 * SA_par[16] * (1 + 0.5 * ISO) * (1 - 0.5 * AF) * 2.7e-4; // [cm/sec]
	double pK = 0 * SA_par[16] * (1 + 0.5 * ISO) * (1 - 0.5 * AF) * 1.35e-7; // [cm/sec]

	//pNa = 1.3*SA_par(17)*(1+0.5*ISO)*(1-0.5*AF)*0.75e-8; // [cm/sec] // MOD2
	//pCa = 1.3*SA_par(17)*(1+0.5*ISO)*(1-0.5*AF)*2.7e-4; // [cm/sec] // MOD2
	//pK = 1.3*SA_par(17)*(1+0.5*ISO)*(1-0.5*AF)*1.35e-7; // [cm/sec] // MOD2

	//pNa = SA_par(17)*(1+0.5*ISO)*(1-0.5*AF)*0.75e-8*0.65; // [cm/sec] // NEW K2P -35// in Schmidt et al 2015
	//pCa = SA_par(17)*(1+0.5*ISO)*(1-0.5*AF)*2.7e-4*0.65; // [cm/sec] // NEW K2P -35// in Schmidt et al 2015
	//pK = SA_par(17)*(1+0.5*ISO)*(1-0.5*AF)*1.35e-7*0.65; // [cm/sec] // NEW K2P -35// in Schmidt et al 2015
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

	cell_par->I_Ca_junc = (Fjunc_CaL * ibarca_j * y[3] * y[4] * ((1 - y[5]) + fcaCaj)) * 0.45; // *Q10CaL^Qpow was removed.
	cell_par->I_Ca_sl = (Fsl_CaL * ibarca_sl * y[3] * y[4] * ((1 - y[6]) + fcaCaMSL)) * 0.45; // *Q10CaL^Qpow was removed.
	cell_par->I_CaK = (ibark * y[3] * y[4] * (Fjunc_CaL * (fcaCaj + (1 - y[5])) + Fsl_CaL * (fcaCaMSL + (1 - y[6])))) * 0.45; // *Q10CaL^Qpow was removed.
	cell_par->I_CaNa_junc = (Fjunc_CaL * ibarna_j * y[3] * y[4] * ((1 - y[5]) + fcaCaj)) * 0.45; // *Q10CaL^Qpow was removed.
	cell_par->I_CaNa_sl = (Fsl_CaL * ibarna_sl * y[3] * y[4] * ((1 - y[6]) + fcaCaMSL)) * .45; // *Q10CaL^Qpow was removed.


	for (int i = 0; i < 10; ++i)
	{
		cell_par->LTCC.state[i] = y[42 + 3 + i];
	}


	cell_par->LTCC.update_states_v4_plus(y[35], y[38], y[31], y[34]);


	for (int i = 0; i < 10; ++i)
	{
		ydot[42 + 3 + i] = cell_par->LTCC.ydot[i];
	}


	for (int i = 0; i < 10; ++i)
	{
		cell_par->LTCC_sl.state[i] = y[42 + 3 + 10 + i];
	}


	cell_par->LTCC_sl.update_states_v4_plus(y[36], y[38], y[32], y[34]);


	for (int i = 0; i < 10; ++i)
	{
		ydot[42 + 3 + 10 + i] = cell_par->LTCC_sl.ydot[i];
	}

	cell_par->I_Ca_junc = 0.9 * cell_par->LTCC.I_Ca_junc_m1;
	cell_par->I_CaNa_junc = 0.9 * cell_par->LTCC.I_Na_junc_m1;

	cell_par->I_Ca_sl = 0.1 * cell_par->LTCC_sl.I_Ca_junc_m1;
	cell_par->I_CaNa_sl = 0.1 * cell_par->LTCC_sl.I_Na_junc_m1;

	cell_par->I_CaK = 0.9 * cell_par->LTCC.I_K_junc_m1 + 0.1 * cell_par->LTCC_sl.I_K_junc_m1;

	cell_par->I_Ca = cell_par->I_Ca_junc + cell_par->I_Ca_sl;
	cell_par->I_CaNa = cell_par->I_CaNa_junc + cell_par->I_CaNa_sl;
	cell_par->I_Catot = cell_par->I_Ca + cell_par->I_CaK + cell_par->I_CaNa; //#ok<NASGU>
	//// I_cabk: Ca Background Current
	double GCaB = SA_par[17] * 6.0643e-4; // [uA/uF]
	//GCaB = SA_par[17]*6.0643e-4*(1+1*AF); // [uA/uF] MOD2
	if (MOD_ind == 3) {
		GCaB = SA_par[17] * 6.0643e-4 * (1 + 0.5 * AF); // [uA/uF] MOD3
	}

	cell_par->I_cabk_junc = Fjunc * GCaB * (y[39 - 1] - eca_junc);
	cell_par->I_cabk_sl = Fsl * GCaB * (y[39 - 1] - eca_sl);
	cell_par->I_cabk = cell_par->I_cabk_junc + cell_par->I_cabk_sl; //#ok<NASGU>

	//// I_pca: Sarcolemmal Ca Pump Current
	double IbarSLCaP = SA_par[18] * 0.0471; // [uA/uF]
	//IbarSLCaP = SA_par[18]*0.0471*4; // [uA/uF] MOD2
	if (MOD_ind == 3) {
		IbarSLCaP = SA_par[18] * 0.0471 * 2; // [uA/uF] MOD3
	}

	double KmPCa = 0.5e-3;     // [mM]
	double Q10SLCaP = 2.35;    // [none]

	cell_par->I_pca_junc = Fjunc * IbarSLCaP * pow(y[35], 1.6) / (pow(KmPCa, 1.6) + pow(y[35], 1.6)); // Q10SLCaP^Qpow* was removed
	cell_par->I_pca_sl = Fsl * IbarSLCaP * pow(y[37 - 1], 1.6) / (pow(KmPCa, 1.6) + pow(y[37 - 1], 1.6)); // Q10SLCaP^Qpow* was removed
	cell_par->I_pca = cell_par->I_pca_junc + cell_par->I_pca_sl;

	//// I_ncx: Na/Ca Exchanger flux
	double IbarNCX = 1.3 * SA_par[19] * (1 + 0.4 * AF) * 3.15; // [uA/uF]
	//IbarNCX = SA_par[19]*(1+1*AF)*3.15; // [uA/uF] MOD2
	if (MOD_ind == 3) {
		IbarNCX = 1. * SA_par[19] * (1 + 0.8 * AF) * 3.15; // [uA/uF] MOD3
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

	cell_par->I_ncx_junc = Fjunc * IbarNCX * Ka_junc * (s1_junc - s2_junc) / s3_junc / (1 + ksat * exp((nu - 1) * y[39 - 1] * FoRT)); //Q10NCX^Qpow*removed
	cell_par->I_ncx_sl = Fsl * IbarNCX * Ka_sl * (s1_sl - s2_sl) / s3_sl / (1 + ksat * exp((nu - 1) * y[39 - 1] * FoRT)); //Q10NCX^Qpow*removed
	cell_par->I_ncx = cell_par->I_ncx_junc + cell_par->I_ncx_sl;

	//// SR fluxes: Calcium Uptake, Release, and Leak
	double Vmax_SRCaP = SA_par[20] * 5.3114e-3; // [mM/msec] (286 umol/L cytosol/sec)
	//Vmax_SRCaP = SA_par[20]*5.3114e-3*(1-0.5*AF); // [mM/msec] MOD2
	if (MOD_ind == 3) {
		Vmax_SRCaP = SA_par[20] * 5.3114e-3 * (1 - 0.25 * AF); // [mM/msec] MOD3
	}

	double Q10SRCaP = 2.6;          // [none]
	double Kmf = (2.5 - 1.25 * ISO) * 0.246e-3; // [mM]
	double Kmr = 1.7;               // [mM]L cytosol
	double hillSRCaP = 1.787;       // [mM]
	double ks = SA_par[22 - 1] * 25;  // [1/ms]
	double koCa = 10 + 20 * AF + 10 * ISO * (1 - AF); // [mM^-2 1/ms]
	double kom = 0.06;              // [1/ms]
	double kiCa = 0.5;              // [1/mM/ms]
	double kim = 0.005;             // [1/ms]
	// double kim = 0.01;             // [1/ms]  // 15:07:36, Thu, 07-February-2019, By Haibo
	double ec50SR = 0.45;           // [mM]
	double MaxSR = 15;              // [mM]
	double MinSR = 1;               // [mM]
	double kleak = SA_par[23 - 1] * (1.0 + 0.25 * AF) * 5.348e-6;

	double kCaSR = MaxSR - (MaxSR - MinSR) / (1 + pow((ec50SR / y[31 - 1]), 2.5));
	double koSRCa = koCa / kCaSR;
	double kiSRCa = kiCa * kCaSR;

	double KoSR_Ca_2 = koSRCa * y[35] * y[35];//+1e-6;  // original
	// double KoSR_Ca_2 = koSRCa * 1e-4 * 1.0 / (1.0 + pow(0.03/y[35], 2.0));
	// double KoSR_Ca_2 = koSRCa * y[35] * y[35];
	double RI = 1 - y[13] - y[14] - y[15];
	ydot[13] = (kim * RI - kiSRCa * y[35] * y[13]) - (KoSR_Ca_2 * y[13] - kom * y[14]); // R
	ydot[14] = (KoSR_Ca_2 * y[13] - kom * y[14]) - (kiSRCa * y[35] * y[14] - kim * y[15]); // O
	ydot[15] = (kiSRCa * y[35] * y[14] - kim * y[15]) - (kom * y[15] - KoSR_Ca_2 * RI); // 0*
	cell_par->J_SRCarel =  ks * y[14] * (y[31 - 1] - y[35]); // [mM/ms]

	cell_par->J_serca = Vmax_SRCaP * (pow(y[37] / Kmf, hillSRCaP) - pow(y[30] / Kmr, hillSRCaP)) / (1 + pow(y[37] / Kmf, hillSRCaP) + pow(y[30] / Kmr, hillSRCaP)); // [mM/ms]   Q10SRCaP^Qpow*

	cell_par->J_SRleak =  kleak * (y[30] - y[35]); // [mM/ms]

	//// Sodium and Calcium Buffering
	// koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
	double Bmax_Naj = 7.561;       // [mM] // Na buffering
	double Bmax_Nasl = 1.65;       // [mM]
	double koff_na = 1e-3;         // [1/ms]
	double kon_na = 0.1e-3;        // [1/mM/ms]
	double Bmax_TnClow = 70e-3;    // [mM]                      // TnC low affinity
	double koff_tncl = (1 + 0.5 * ISO) * 19.6e-3; // [1/ms]
	double kon_tncl = 32.7;        // [1/mM/ms]
	double Bmax_TnChigh = 140e-3;  // [mM]                      // TnC high affinity
	double koff_tnchca = 0.032e-3; // [1/ms]
	double kon_tnchca = 2.37;      // [1/mM/ms]
	double koff_tnchmg = 3.33e-3;  // [1/ms]
	double kon_tnchmg = 3e-3;      // [1/mM/ms]
	double Bmax_CaM = 24e-3;       // [mM] **? about setting to 0 in c-code**   // CaM buffering
	double koff_cam = 238e-3;      // [1/ms]
	double kon_cam = 34;           // [1/mM/ms]
	double Bmax_myosin = 140e-3;   // [mM]                      // Myosin buffering
	double koff_myoca = 0.46e-3;   // [1/ms]
	double kon_myoca = 13.8;       // [1/mM/ms]
	double koff_myomg = 0.057e-3;  // [1/ms]
	double kon_myomg = 0.0157;     // [1/mM/ms]
	double Bmax_SR = 19 * .9e-3;   // [mM] (Bers text says 47e-3) 19e-3
	double koff_sr = 60e-3;        // [1/ms]
	double kon_sr = 100;           // [1/mM/ms]

	double BufferScaling = 1;
	double Bmax_SLlowsl = BufferScaling * 37.4e-3 * Vmyo / Vsl;    // [mM]   	// SL buffering
	double Bmax_SLlowj = BufferScaling * 4.6e-3 * Vmyo / Vjunc * 0.1; // [mM]
	double koff_sll = 1300e-3;     // [1/ms]
	double kon_sll = 100;          // [1/mM/ms]
	double Bmax_SLhighsl = BufferScaling * 13.4e-3 * Vmyo / Vsl;   // [mM]
	double Bmax_SLhighj = BufferScaling * 1.65e-3 * Vmyo / Vjunc * 0.1; // [mM]
	double koff_slh = 30e-3;       // [1/ms]
	double kon_slh = 100;          // [1/mM/ms]
	double Bmax_Csqn = 140e-3 * Vmyo / Vsr;        // [mM]      // Csqn buffering
	double koff_csqn = 65;         // [1/ms]
	double kon_csqn = 100;         // [1/mM/ms]


	// OK, here we add a 7-state markov model.	 // 14:14:27, Thu, 07-February-2019, By Haibo


	// Junctional and SL Na Buffers
	ydot[16] = kon_na * y[32 - 1] * (Bmax_Naj - y[16]) - koff_na * y[16]; // NaBj      [mM/ms]
	ydot[17] = kon_na * y[33 - 1] * (Bmax_Nasl - y[17]) - koff_na * y[17]; // NaBsl     [mM/ms]

	// Cytosolic Ca Buffers
	ydot[18] = kon_tncl * y[37] * (Bmax_TnClow - y[18]) - koff_tncl * y[18];  // TnCL      [mM/ms]
	ydot[19] = kon_tnchca * y[37] * (Bmax_TnChigh - y[19] - y[20]) - koff_tnchca * y[19]; // TnCHc     [mM/ms]
	ydot[20] = kon_tnchmg * Mgi * (Bmax_TnChigh - y[19] - y[20]) - koff_tnchmg * y[20]; // TnCHm     [mM/ms]
	ydot[21] = kon_cam * y[37] * (Bmax_CaM - y[21]) - koff_cam * y[21];       // CaM       [mM/ms]
	ydot[22] = kon_myoca * y[37] * (Bmax_myosin - y[22] - y[23]) - koff_myoca * y[22]; // Myosin_ca [mM/ms]
	ydot[23] = kon_myomg * Mgi * (Bmax_myosin - y[22] - y[23]) - koff_myomg * y[23]; // Myosin_mg [mM/ms]
	ydot[24] = kon_sr * y[37] * (Bmax_SR - y[24]) - koff_sr * y[24];          // SRB       [mM/ms]
	cell_par->J_CaB_cytosol = ydot[18] + ydot[19] + ydot[21] + ydot[22] + ydot[24];

	// Junctional and SL Ca Buffers
	ydot[25] = kon_sll * y[35] * (Bmax_SLlowj - y[25]) - koff_sll * y[25]; // SLLj      [mM/ms]
	ydot[26] = kon_sll * y[37 - 1] * (Bmax_SLlowsl - y[26]) - koff_sll * y[26]; // SLLsl     [mM/ms]
	ydot[27] = kon_slh * y[35] * (Bmax_SLhighj - y[27]) - koff_slh * y[27]; // SLHj      [mM/ms]
	ydot[28] = kon_slh * y[37 - 1] * (Bmax_SLhighsl - y[28]) - koff_slh * y[28]; // SLHsl     [mM/ms]
	cell_par->J_CaB_junction = ydot[25] + ydot[27];
	cell_par->J_CaB_sl = ydot[26] + ydot[28];

	//// Ion concentrations
	// SR Ca Concentrations
	ydot[29] = kon_csqn * y[30] * (Bmax_Csqn - y[29]) - koff_csqn * y[29]; // Csqn      [mM/ms]
	ydot[30] = cell_par->J_serca - (cell_par->J_SRleak * Vmyo / Vsr + cell_par->J_SRCarel) - ydot[29]; // Ca_sr     [mM/ms]

	// Sodium Concentrations
	cell_par->I_Na_tot_junc = cell_par->I_Na_junc + cell_par->I_nabk_junc + 3 * cell_par->I_ncx_junc + 3 * cell_par->I_nak_junc + cell_par->I_CaNa_junc + cell_par->I_NaL_junc; // [uA/uF]
	cell_par->I_Na_tot_sl = cell_par->I_Na_sl + cell_par->I_nabk_sl + 3 * cell_par->I_ncx_sl + 3 * cell_par->I_nak_sl + cell_par->I_CaNa_sl + cell_par->I_NaL_sl; // [uA/uF]

	ydot[31] = -cell_par->I_Na_tot_junc * Cmem / (Vjunc * Frdy) + J_na_juncsl / Vjunc * (y[33 - 1] - y[32 - 1]) - ydot[16];
	ydot[32] = -cell_par->I_Na_tot_sl * Cmem / (Vsl * Frdy) + J_na_juncsl / Vsl * (y[32 - 1] - y[33 - 1]) + J_na_slmyo / Vsl * (y[33] - y[33 - 1]) - ydot[17];
	ydot[33] = (J_na_slmyo / Vmyo * (y[33 - 1] - y[33])) * (1 - Na_clamp); // [mM/msec]

	// Potassium Concentration
	cell_par->I_K_tot = cell_par->I_to + cell_par->I_kr + cell_par->I_ks + cell_par->I_ki - 2 * cell_par->I_nak + cell_par->I_CaK + cell_par->I_kp + cell_par->I_kur + cell_par->I_kach + cell_par->I_k2p + cell_par->I_sk; // [uA/uF] //SVP: added IKur
	ydot[34] = 0; // -cell_par->I_K_tot*Cmem/(Vmyo*Frdy);         // [mM/msec]

	// Calcium Concentrations
	cell_par->I_Ca_tot_junc = cell_par->I_Ca_junc + cell_par->I_cabk_junc + cell_par->I_pca_junc - 2 * cell_par->I_ncx_junc;           // [uA/uF]
	cell_par->I_Ca_tot_sl = cell_par->I_Ca_sl + cell_par->I_cabk_sl + cell_par->I_pca_sl - 2 * cell_par->I_ncx_sl;    // [uA/uF]
	ydot[35] = -cell_par->I_Ca_tot_junc * Cmem / (Vjunc * 2 * Frdy) + J_ca_juncsl / Vjunc * (y[36] - y[35]) - cell_par->J_CaB_junction + (cell_par->J_SRCarel) * Vsr / Vjunc + cell_par->J_SRleak * Vmyo / Vjunc; // Ca_j
	ydot[36] = -cell_par->I_Ca_tot_sl * Cmem / (Vsl * 2 * Frdy) + J_ca_juncsl / Vsl * (y[35] - y[36]) + J_ca_slmyo / Vsl * (y[37] - y[36]) - cell_par->J_CaB_sl; // Ca_sl
	ydot[37] = -cell_par->J_serca * Vsr / Vmyo - cell_par->J_CaB_cytosol + J_ca_slmyo / Vmyo * (y[36] - y[37]); // [mM/msec]



	if ( Ca_clamp == 1 || Ca_clamp == 2) { // Ca_clamp, BAPTA
		ydot[35] = 0; ydot[36] = 0; ydot[37] = 0;
	} else if (Ca_clamp == 3) { // EGTA
		ydot[37] = 0;
	}

	// //// Simulation type
	// switch p(1)
	//     case 0          // no stimulation
	//         cell_par->I_app = 0;
	//     case 1          // pace w/ current injection at rate 'rate' (Hz)
	//         rate = p(2); // Hz
	//         period = 1000/rate; // ms
	//         if mod(t,period) <= 5
	//             cell_par->I_app = 12.5;
	//         else
	//             cell_par->I_app = 0.0;
	//         end
	//     case 2      // ERP
	//         rate = p(2); // Hz
	//         rec_interval = p(3);
	//         if t <= 5
	//             cell_par->I_app = 12.5;
	//         elseif t > 5 && t <= rec_interval
	//             cell_par->I_app = 0.0;
	//         elseif t > rec_interval && t <= rec_interval+5
	//             if rate == 0.5 && AF == 0
	//                 DTE = 12.5*0.125;
	//             else
	//                 DTE = 12.5*0.2;
	//             end
	//             cell_par->I_app = 2*DTE;
	//         else
	//             cell_par->I_app = 0.0;
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
	//         cell_par->I_app = (V_clamp-y[39-1])/R_clamp;
	// end

	//// Membrane Potential
	cell_par->I_Na_tot = cell_par->I_Na_tot_junc + cell_par->I_Na_tot_sl;
	cell_par->I_Cl_tot = cell_par->I_ClCa + cell_par->I_Clbk + cell_par->I_ClCFTR;
	cell_par->I_Ca_tot = cell_par->I_Ca_tot_junc + cell_par->I_Ca_tot_sl;
	cell_par->I_tot = cell_par->I_Na_tot + cell_par->I_Cl_tot + cell_par->I_Ca_tot + cell_par->I_K_tot; // [uA/uF]
	ydot[38] = -(cell_par->I_tot - cell_par->I_app);
	double vmax = ydot[38];


	cell_par->dV = ydot[38];
	cell_par->CaSR = y[30];
	cell_par->Caj = y[35];
	cell_par->Casl = y[36];
	cell_par->Cai = y[37];
	cell_par->Nai = y[33];
	cell_par->Naj = y[32];
	cell_par->Nasl = y[31];



	// cell_par->dV = ydot[38];
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


int human_atrial_ECC_initialiser(double *y) {
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

	y[0]     = 0.004080635368;
	y[1]     = 0.9451285505;
	y[2]     = 0.9622833438;
	y[3]     = 8.413437798e-06;
	y[4]     = 0.999413341;
	y[5]     = 0.05762842186;
	y[6]     = 0.03955197377;
	y[7]     = 0.0002032347954;
	y[8]     = 0.9372571103;
	y[9]     = 0.0008234263476;
	y[10]    = 0.9663079588;
	y[11]    = 0.003206591541;
	y[12]    = 0.005041386401;
	y[13]    = 0.7834945442;
	y[14]    = 3.83997118e-06;
	y[15]    = 1.06108702e-06;
	y[16]    = 3.800037507;
	y[17]    = 0.8292064491;
	y[18]    = 0.02106415581;
	y[19]    = 0.1308372858;
	y[20]    = 0.00425980829;
	y[21]    = 0.0008473481432;
	y[22]    = 0.005291039036;
	y[23]    = 0.1341852704;
	y[24]    = 0.005122289046;
	y[25]    = 0.01752575393;
	y[26]    = 0.02586842913;
	y[27]    = 0.116599258;
	y[28]    = 0.2116176743;
	y[29]    = 1.203914851;
	y[30]    = 0.5605281892;
	y[31]    = 10.10011682;
	y[32]    = 10.09974519;
	y[33]    = 10.09990903;
	y[34]    = 120;
	y[35]    = 0.0004240838797;
	y[36]    = 0.0002826454379;
	y[37]    = 0.0002559818175;
	y[38]    = -79.11428795;
	y[39]    = 0.004080635368;
	y[40]    = 0.1179555119;
	y[41]    = 0.2013739322;
	y[42]    = 0;
	y[43]    = 1;
	y[44]    = 1;
	// y[42 + 3]    = 0.9895726818;
	// y[43 + 3]    = 1.690746702e-05;
	// y[44 + 3]    = 4.240266403e-10;
	// y[45 + 3]    = 0.0003644860657;
	// y[46 + 3]    = 0.006519349275;
	// y[47 + 3]    = 0.003526009982;
	// y[48 + 3]    = 2.134724476e-15;
	// y[49 + 3]    = 1.299456512e-07;
	// y[42 + 3 + 8] = 0.9871462467;
	// y[43 + 3 + 8] = 1.897481412e-05;
	// y[44 + 3 + 8] = 4.090713291e-10;
	// y[45 + 3 + 8] = 0.0003628823657;
	// y[46 + 3 + 8] = 0.008093821615;
	// y[47 + 3 + 8] = 0.004377438036;
	// y[48 + 3 + 8] = 2.125774633e-15;
	// y[49 + 3 + 8] = 1.293105713e-07;
	y[42 + 3] = 1;
	y[42 + 3 + 10] = 1;
	for (int i = 1; i < 10; ++i)
	{
		y[42 + 3 + i] = 0;
		y[42 + 3 + 10 + i] = 0;
	}
// 0


}



void Cell::print_to_file(double t, std::ofstream & output_file) {


	output_file <<  std::setprecision(9)
	            << t << " "  // 1
	            << V << " "
	            << dV << " "
	            << I_app << " "
	            << CaSR << " "  // 5
	            << Caj << " " // 6
	            << Casl << " "
	            << Cai << " "
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
	            << std::endl;
}