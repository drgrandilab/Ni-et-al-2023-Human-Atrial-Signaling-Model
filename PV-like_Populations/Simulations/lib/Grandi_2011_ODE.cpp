#include "Grandi_2011_ODE.h"
#include <cmath>

void Grandi_2011_ODE_Initialise(double *state, int celltype) {
	state[0] = 8.72509677797499e-5;
	state[1 ] = Y1;
	state[2 ] = Y2;
	state[3 ] = Y3;
	state[4 ] = Y4;
	state[5 ] = Y5;
	state[6 ] = Y6;
	state[7 ] = Y7;
	state[8 ] = Y8;
	state[9 ] = Y9;
	state[10] = Y10;
	state[11] = Y11;
	state[12] = Y12;
	state[13] = Y13;
	state[14] = Y14;
	state[15] = Y15;
	state[16] = Y16;
	state[17] = Y17;
	state[18] = Y18;
	state[19] = Y19;
	state[20] = Y20;
	state[21] = Y21;
	state[22] = Y22;
	state[23] = Y23;
	state[24] = Y24;
	state[25] = Y25;
	state[26] = Y26;
	state[27] = Y27;
	state[28] = Y28;
	state[29] = Y29;
	state[30] = Y30;
	state[31] = Y31;
	state[32] = Y32;
	state[33] = Y33;
	state[34] = Y34;
	state[35] = Y35;
	state[36] = Y36;
	state[37] = Y37;
	state[38] = Y38;
	state[39] = Y39;
	state[40] = Y40;
// 41
}


void Grandi_2011_ODE(double t, double *Y, double * dY, void *usr_data) {
	SingleCellPara *Data = (SingleCellPara*) usr_data;

	Data->V = Y[38];

	Grandi_2011_ODE( t,  Y,   dY, Data);
	dY[38] =Data->dV;
}


void Grandi_2011_ODE(double t, double *Y, double * dY, SingleCellPara *Data)
{



	// int i;

	// for (i = 0; i < NEQ; i++)
	//     Y[i] = Ith(y, i + 1);

	double I_Ca_tot_junc;   // uA_per_uF (in Ca_Concentrations)
	double I_Ca_tot_sl;   // uA_per_uF (in Ca_Concentrations)
	double J_CaB_cytosol;   // mM_per_msec (in Cytosolic_Ca_Buffers)
	double I_cabk;   // uA_per_uF (in I_CaBK)
	double I_cabk_junc;   // uA_per_uF (in I_CaBK)
	double I_cabk_sl;   // uA_per_uF (in I_CaBK)
	double I_Ca;   // uA_per_uF (in I_Ca)
	double I_CaK;   // uA_per_uF (in I_Ca)
	double I_CaNa;   // uA_per_uF (in I_Ca)
	double I_CaNa_junc;   // uA_per_uF (in I_Ca)
	double I_CaNa_sl;   // uA_per_uF (in I_Ca)
	double I_Ca_junc;   // uA_per_uF (in I_Ca)
	double I_Ca_sl;   // uA_per_uF (in I_Ca)
	double I_Catot;   // uA_per_uF (in I_Ca)
	double dss;   // dimensionless (in I_Ca)
	double fcaCaMSL;   // dimensionless (in I_Ca)
	double fcaCaj;   // dimensionless (in I_Ca)
	double fss;   // dimensionless (in I_Ca)
	double ibarca_j;   // uA_per_uF (in I_Ca)
	double ibarca_sl;   // uA_per_uF (in I_Ca)
	double ibark;   // uA_per_uF (in I_Ca)
	double ibarna_j;   // uA_per_uF (in I_Ca)
	double ibarna_sl;   // uA_per_uF (in I_Ca)
	double taud;   // msec (in I_Ca)
	double tauf;   // msec (in I_Ca)
	double I_ClCa;   // uA_per_uF (in I_ClCa)
	double I_ClCa_junc;   // uA_per_uF (in I_ClCa)
	double I_ClCa_sl;   // uA_per_uF (in I_ClCa)
	double I_Clbk;   // uA_per_uF (in I_ClCa)
	double I_ki;   // uA_per_uF (in I_Ki)
	double aki;   // dimensionless (in I_Ki)
	double bki;   // dimensionless (in I_Ki)
	double kiss;   // dimensionless (in I_Ki)
	double I_kp;   // uA_per_uF (in I_Kp)
	double I_kp_junc;   // uA_per_uF (in I_Kp)
	double I_kp_sl;   // uA_per_uF (in I_Kp)
	double kp_kp;   // dimensionless (in I_Kp)
	double I_kr;   // uA_per_uF (in I_Kr)
	double gkr;   // mS_per_uF (in I_Kr)
	double rkr;   // dimensionless (in I_Kr)
	double tauxr;   // msec (in I_Kr)
	double xrss;   // dimensionless (in I_Kr)
	double I_ks;   // uA_per_uF (in I_Ks)
	double I_ks_junc;   // uA_per_uF (in I_Ks)
	double I_ks_sl;   // uA_per_uF (in I_Ks)
	double eks;   // mV (in I_Ks)
	double gks_junc;   // mS_per_uF (in I_Ks)
	double gks_sl;   // mS_per_uF (in I_Ks)
	double tauxs;   // msec (in I_Ks)
	double xsss;   // dimensionless (in I_Ks)
	double I_ncx;   // uA_per_uF (in I_NCX)
	double I_ncx_junc;   // uA_per_uF (in I_NCX)
	double I_ncx_sl;   // uA_per_uF (in I_NCX)
	double Ka_junc;   // dimensionless (in I_NCX)
	double Ka_sl;   // dimensionless (in I_NCX)
	double s1_junc;   // mM4 (in I_NCX)
	double s1_sl;   // mM4 (in I_NCX)
	double s2_junc;   // mM4 (in I_NCX)
	double s2_sl;   // mM4 (in I_NCX)
	double s3_junc;   // mM4 (in I_NCX)
	double s3_sl;   // mM4 (in I_NCX)
	double I_nabk;   // uA_per_uF (in I_NaBK)
	double I_nabk_junc;   // uA_per_uF (in I_NaBK)
	double I_nabk_sl;   // uA_per_uF (in I_NaBK)
	double I_nak;   // uA_per_uF (in I_NaK)
	double I_nak_junc;   // uA_per_uF (in I_NaK)
	double I_nak_sl;   // uA_per_uF (in I_NaK)
	double fnak;   // dimensionless (in I_NaK)
	double sigma;   // dimensionless (in I_NaK)
	double I_Na;   // uA_per_uF (in I_Na)
	double I_Na_junc;   // uA_per_uF (in I_Na)
	double I_Na_sl;   // uA_per_uF (in I_Na)
	double ah;   // dimensionless (in I_Na)
	double aj;   // dimensionless (in I_Na)
	double bh;   // dimensionless (in I_Na)
	double bj;   // dimensionless (in I_Na)
	double hss;   // dimensionless (in I_Na)
	double jss;   // dimensionless (in I_Na)
	double mss;   // dimensionless (in I_Na)
	double tauh;   // msec (in I_Na)
	double tauj;   // msec (in I_Na)
	double taum;   // msec (in I_Na)
	double I_pca;   // uA_per_uF (in I_PCa)
	double I_pca_junc;   // uA_per_uF (in I_PCa)
	double I_pca_sl;   // uA_per_uF (in I_PCa)
	double GtoFast;   // mS_per_uF (in I_to)
	double GtoSlow;   // mS_per_uF (in I_to)
	double I_to;   // uA_per_uF (in I_to)
	double I_tof;   // uA_per_uF (in I_to)
	double I_tos;   // uA_per_uF (in I_to)
	double tauxtof;   // msec (in I_to)
	double tauxtos;   // msec (in I_to)
	double tauytof;   // msec (in I_to)
	double tauytos;   // msec (in I_to)
	double xtoss;   // dimensionless (in I_to)
	double ytoss;   // dimensionless (in I_to)
	double J_CaB_junction;   // mM_per_msec (in Junctional_and_SL_Ca_Buffers)
	double J_CaB_sl;   // mM_per_msec (in Junctional_and_SL_Ca_Buffers)
	double I_K_tot;   // uA_per_uF (in K_Concentration)
	double dNa_Bj_dt;   // mM_per_msec (in Na_Buffers)
	double dNa_Bsl_dt;   // mM_per_msec (in Na_Buffers)
	double I_Na_tot_junc;   // uA_per_uF (in Na_Concentrations)
	double I_Na_tot_junc2;   // uA_per_uF (in Na_Concentrations)
	double I_Na_tot_sl;   // uA_per_uF (in Na_Concentrations)
	double I_Na_tot_sl2;   // uA_per_uF (in Na_Concentrations)
	double J_SRCarel;   // mM_per_msec (in SR_Fluxes)
	double J_SRleak;   // mM_per_msec (in SR_Fluxes)
	double J_serca;   // mM_per_msec (in SR_Fluxes)
	double MaxSR;   // dimensionless (in SR_Fluxes)
	double MinSR;   // dimensionless (in SR_Fluxes)
	double RI;   // mM (in SR_Fluxes)
	double kCaSR;   // dimensionless (in SR_Fluxes)
	double kiSRCa;   // per_mM_per_msec (in SR_Fluxes)
	double koSRCa;   // per_mM2_per_msec (in SR_Fluxes)
	double I_Ca_tot;   // uA_per_uF (in membrane_potential)
	double I_Cl_tot;   // uA_per_uF (in membrane_potential)
	double I_Na_tot;   // uA_per_uF (in membrane_potential)
	double I_tot;   // uA_per_uF (in membrane_potential)
	double Bmax_CaM;   // mM (in parameters)
	double Bmax_Csqn;   // mM (in parameters)
	double Bmax_Naj;   // mM (in parameters)
	double Bmax_Nasl;   // mM (in parameters)
	double Bmax_SLhighj;   // mM (in parameters)
	double Bmax_SLhighsl;   // mM (in parameters)
	double Bmax_SLlowj;   // mM (in parameters)
	double Bmax_SLlowsl;   // mM (in parameters)
	double Bmax_SR;   // mM (in parameters)
	double Bmax_TnChigh;   // mM (in parameters)
	double Bmax_TnClow;   // mM (in parameters)
	double Bmax_myosin;   // mM (in parameters)
	double Cao;   // mM (in parameters)
	double Cli;   // mM (in parameters)
	double Clo;   // mM (in parameters)
	double Cmem;   // farad (in parameters)
	double DcaJuncSL;   // cm2_per_sec (in parameters)
	double DcaSLcyto;   // cm2_per_sec (in parameters)
	double DnaJuncSL;   // cm2_per_sec (in parameters)
	double DnaSLcyto;   // cm2_per_sec (in parameters)
	double Fjunc;   // dimensionless (in parameters)
	double Fjunc_CaL;   // dimensionless (in parameters)
	double FoRT;   // per_mV (in parameters)
	double Frdy;   // coulomb_per_mole (in parameters)
	double Fsl;   // dimensionless (in parameters)
	double Fsl_CaL;   // dimensionless (in parameters)
	double GCaB;   // mS_per_uF (in parameters)
	double GClB;   // mS_per_uF (in parameters)
	double GClCa;   // mS_per_uF (in parameters)
	double GNa;   // mS_per_uF (in parameters)
	double GNaB;   // mS_per_uF (in parameters)
	double IbarNCX;   // uA_per_uF (in parameters)
	double IbarNaK;   // uA_per_uF (in parameters)
	double IbarSLCaP;   // uA_per_uF (in parameters)
	double J_ca_juncsl;   // liters_per_msec (in parameters)
	double J_ca_slmyo;   // liters_per_msec (in parameters)
	double J_na_juncsl;   // liters_per_msec (in parameters)
	double J_na_slmyo;   // liters_per_msec (in parameters)
	double KdClCa;   // mM (in parameters)
	double Kdact;   // mM (in parameters)
	double KmCai;   // mM (in parameters)
	double KmCao;   // mM (in parameters)
	double KmKo;   // mM (in parameters)
	double KmNai;   // mM (in parameters)
	double KmNaip;   // mM (in parameters)
	double KmNao;   // mM (in parameters)
	double KmPCa;   // mM (in parameters)
	double Kmf;   // mM (in parameters)
	double Kmr;   // mM (in parameters)
	double Ko;   // mM (in parameters)
	double Mgi;   // mM (in parameters)
	double Nao;   // mM (in parameters)
	double Q10CaL;   // dimensionless (in parameters)
	double Q10KmNai;   // dimensionless (in parameters)
	double Q10NCX;   // dimensionless (in parameters)
	double Q10NaK;   // dimensionless (in parameters)
	double Q10SLCaP;   // dimensionless (in parameters)
	double Q10SRCaP;   // dimensionless (in parameters)
	double Qpow;   // dimensionless (in parameters)
	double R;   // joule_per_kelvin_per_kilomole (in parameters)
	double SAjunc;   // um2 (in parameters)
	double SAsl;   // um2 (in parameters)
	double Temp;   // kelvin (in parameters)
	double Vcell;   // liter (in parameters)
	double Vjunc;   // liter (in parameters)
	double Vmax_SRCaP;   // mM_per_msec (in parameters)
	double Vmyo;   // liter (in parameters)
	double Vsl;   // liter (in parameters)
	double Vsr;   // liter (in parameters)
	double cellLength;   // um (in parameters)
	double cellRadius;   // um (in parameters)
	double distJuncSL;   // um (in parameters)
	double distSLcyto;   // um (in parameters)
	double ec50SR;   // mM (in parameters)
	double eca_junc;   // mV (in parameters)
	double eca_sl;   // mV (in parameters)
	double ecl;   // mV (in parameters)
	double ek;   // mV (in parameters)
	double ena_junc;   // mV (in parameters)
	double ena_sl;   // mV (in parameters)
	// double epi;   // dimensionless (in parameters)
	double gkp;   // mS_per_uF (in parameters)
	double hillSRCaP;   // dimensionless (in parameters)
	double junctionLength;   // um (in parameters)
	double junctionRadius;   // um (in parameters)
	double kiCa;   // per_mM_per_msec (in parameters)
	double kim;   // per_msec (in parameters)
	double koCa;   // per_mM2_per_msec (in parameters)
	double koff_cam;   // per_msec (in parameters)
	double koff_csqn;   // per_msec (in parameters)
	double koff_myoca;   // per_msec (in parameters)
	double koff_myomg;   // per_msec (in parameters)
	double koff_na;   // per_msec (in parameters)
	double koff_slh;   // per_msec (in parameters)
	double koff_sll;   // per_msec (in parameters)
	double koff_sr;   // per_msec (in parameters)
	double koff_tnchca;   // per_msec (in parameters)
	double koff_tnchmg;   // per_msec (in parameters)
	double koff_tncl;   // per_msec (in parameters)
	double kom;   // per_msec (in parameters)
	double kon_cam;   // per_mM_per_msec (in parameters)
	double kon_csqn;   // per_mM_per_msec (in parameters)
	double kon_myoca;   // per_mM_per_msec (in parameters)
	double kon_myomg;   // per_mM_per_msec (in parameters)
	double kon_na;   // per_mM_per_msec (in parameters)
	double kon_slh;   // per_mM_per_msec (in parameters)
	double kon_sll;   // per_mM_per_msec (in parameters)
	double kon_sr;   // per_mM_per_msec (in parameters)
	double kon_tnchca;   // per_mM_per_msec (in parameters)
	double kon_tnchmg;   // per_mM_per_msec (in parameters)
	double kon_tncl;   // per_mM_per_msec (in parameters)
	double ks;   // per_msec (in parameters)
	double ksat;   // dimensionless (in parameters)
	double nu;   // dimensionless (in parameters)
	double pCa;   // cm_per_sec (in parameters)
	double pK;   // cm_per_sec (in parameters)
	double pNa;   // cm_per_sec (in parameters)
	double pNaK;   // dimensionless (in parameters)

	double iksmultiply, ik1multiply, ina_c;
	double Gkur, g_Ikur, tau_ac_Ikur, tau_inac_Ikur, ac_Ikur_inf, inac_Ikur_inf;  // Ikur

	/* Late INa, only in AF condition*/
	double I_NaL_junc;
	double I_NaL_sl;
	double I_NaL;
	double GNaL;
	double aml, bml, hlinf, tauhl;
	double Ikur;


	double GClCFTR = 0; //%4.9e-3*ISO;     % [mS/uF]
	double I_ClCFTR ;//= GClCFTR*(y(39)-ecl);
	double ISO, AF, RA;
	ISO = 0.0;
	AF = 0;
	RA = 0;

	double Vm = Data->V;


	//-------------------------------------------------------------------------
	// Computed variables
	//-------------------------------------------------------------------------
	// Enviromental states

	Frdy = 96485.0;
	R = 8314.0;
	Temp = 310.0;
	FoRT = Frdy / (R * Temp);
	Cmem = 1.10e-10;  // was 1.10e-10 here by haibo
	// cell geometry
	cellRadius = 10.25;
	cellLength = 100.0;
	Vcell = 3.14159265358979 * pow(cellRadius, 2.0) * cellLength * 1.0e-15;   // units in liter L
	Vjunc = 0.0539 * 0.01 * Vcell;
	Vmyo = 0.65 * Vcell;
	Vsr = 0.035 * Vcell; // Haibo, change 29/05/2013
	Vsl = 0.02 * Vcell;
	// the values in the non-C code. by haibo
	DcaJuncSL = 1.64e-6;       // in cm
	DcaSLcyto = 1.22e-6;
	DnaJuncSL = 1.09e-5;
	DnaSLcyto = 1.79e-5;
	// comment out by haibo, 29. 05. 2013
	/*DcaJuncSL = 1.2205e-5;
	DcaSLcyto = 2.8914e-6;
	DnaJuncSL = 2.7121e-7;
	DnaSLcyto = 1.2722e-6;
	*/
	junctionLength = 160.0e-3;
	junctionRadius = 15.0e-3;
	distSLcyto = 0.45;
	distJuncSL = 0.5;
	SAjunc = 20150.0 * 3.14159265358979 * 2.0 * junctionLength * junctionRadius; // Haibo
	SAsl = 3.14159265358979 * 2.0 * cellRadius * cellLength; // Haibo

	// fractional parameters:

	Fjunc = 0.11;
	Fjunc_CaL = 0.9;
	Fsl_CaL = 1.0 - Fjunc_CaL;
	Fsl = 1.0 - Fjunc;

	// fixed ion concerntrations.
	Cli = 15;   // Intracellular Cl  [mM]
	Clo = 150;  // Extracellular Cl  [mM]
	Ko  = 5.4;   // Extracellular K   [mM]
	Nao = 140;  // Extracellular Na  [mM]
	Cao = 1.8;  // Extracellular Ca  [mM]
	Mgi = 1.0;    // Intracellular Mg  [mM]


	// Na transport parameters
	GNa      = 23.0 * (1.0 - 0.1 * AF); // [mS/uF]
	GNaB     = 0.597e-3;    // [mS/uF]
	IbarNaK  = 1.26;     // [uA/uF]    // 1.8*0.7
	KmNaip   = 11.0 * (1.0 - 0.25 * ISO);         // [mM]11 by haibo, from matlab code
	KmKo     = 1.5;         // [mM]1.5
	Q10NaK   = 1.63;
	Q10KmNai = 1.39;

	// K current parameters:
	pNaK = 0.01833;
	gkp  = 0.002;

	// Cl current parameters
	GClCa = 0.5 * 0.109625;
	KdClCa = 100.0e-3;
	GClB = 1.0 * 9.0e-3;


	// I_Ca parameters
	/*pNa = 0.75e-8;     //   [cm/sec]  was    pNa = 3.375e-9; here, haibo
	pCa = 2.7e-4;      //  [cm/sec]
	pK = 1.35e-7;      //   [cm/sec]    was    pK = 6.075e-8; by haibo*/
	pNa = (1.0 + 0.5 * ISO) * (1.0 - 0.5 * AF) * 0.75e-8;   // [cm/sec]
	pCa = (1.0 + 0.5 * ISO) * (1.0 - 0.5 * AF) * 2.7e-4;   // [cm/sec]
	pK  = (1.0 + 0.5 * ISO) * (1.0 - 0.5 * AF) * 1.35e-7;    // [cm/sec]
	Q10CaL = 1.8;      //

	IbarNCX = (1.0 + 0.4 * AF) * 3.15;  //  % [uA/uF]4.5 in ventricle
	KmCai = 3.59e-3;    // [mM]
	KmCao = 1.3;        // [mM]
	KmNai = 12.29;      // [mM]
	KmNao = 87.5;       // [mM]
	ksat = 0.27;        // [none]   was 0.32 here, by haibo
	nu = 0.35;          // [none]   // from SK code, not in atrial paper. Haibo.
	Kdact = 0.384e-3;   // [mM] was 0.225 e-3 haibo
	Q10NCX = 1.57;      // [none]
	IbarSLCaP = 0.0471; //  was 0.0673
	KmPCa = 0.5e-3;     // [mM]
	GCaB = 6.0643e-4;   // [uA/uF]  was    GCaB = 5.513e-4; haibo
	Q10SLCaP = 2.35;    // [none]


	// SR SR_Fluxes
	Q10SRCaP = 2.6;
	Vmax_SRCaP = 5.3114e-3;  // [mM/msec] (286 umol/L cytosol/sec)
	Kmf = (2.5 - 1.25 * ISO) * 0.246e-3;        // [mM] was    Kmf = 2.0*0.246e-3; // Haibo
	Kmr = 1.7;               // [mM]L cytosol
	hillSRCaP = 1.787;       // [mM]
	ks = 25.0;                 // [1/ms]
	koCa = 10.0 + 20.0 * AF + 10.0 * ISO * (1.0 - AF);             // [mM^-2 1/ms]   %default 10   modified 20
	kom = 0.06;              // [1/ms]
	kiCa = 0.5;              // [1/mM/ms]
	kim = 0.005;             // [1/ms]
	ec50SR = 0.45;           // [mM]

	// % Buffering parameters
	// % koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
	Bmax_Naj = 7.561;
	Bmax_Nasl = 1.65;
	koff_na = 1.0e-3;
	kon_na = 0.1e-3;
	Bmax_TnClow = 70.0e-3;
	koff_tncl =  (1.0 + 0.5 * ISO) * 19.6e-3;  // according to the matlab code
	kon_tncl = 32.7;
	Bmax_TnChigh = 140.0e-3;

	kon_sll = 100.0;
	Bmax_SLlowj = 4.6e-3 * Vmyo / Vjunc * 0.1;
	koff_sll = 1300.0e-3;
	kon_slh = 100.0;
	Bmax_SLhighj = 1.65e-3 * Vmyo / Vjunc * 0.1;
	koff_slh = 30.0e-3;

	Bmax_SLlowsl = 37.4e-3 * Vmyo / Vsl;
	Bmax_SLhighsl = 13.4e-3 * Vmyo / Vsl;
	Vmax_SRCaP = 5.3114e-3;
	kon_tnchca = 2.37;
	koff_tnchca = 0.032e-3;
	kon_tnchmg = 3.0e-3;
	koff_tnchmg = 3.33e-3;
	kon_cam = 34.0;
	Bmax_CaM = 24.0e-3;
	koff_cam = 238.0e-3;
	kon_myoca = 13.8;
	Bmax_myosin = 140.0e-3;
	koff_myoca = 0.46e-3;
	kon_myomg = 0.0157;
	koff_myomg = 0.057e-3;
	kon_sr = 100.0;
	Bmax_SR = 19.0 * 0.9e-3;
	koff_sr = 60.0e-3;
	fcaCaj = 0.0;
	fcaCaMSL = 0.0;
	Qpow = (Temp - 310.0) / 10.0;




	ecl = 1.0 / FoRT * log(Cli / Clo);
	gks_junc = 0.0035;
	gks_sl = 0.0035;

	/* epi ventricle parameters */
	GtoSlow = 0.0156; // 0.0156 (mS/microF)
	GtoFast = 0.165;   // 0.1144  (mS/microF) Haibo
	iksmultiply = 1.0;  // gks  0.0035  (mS/microF)
	ik1multiply = 1.0; //  gk1 0.35     (mS/microF)
	ina_c = 0.0; // ina ss
	/**************************************/

	kon_csqn = 100.0;
	Bmax_Csqn = 140.0e-3 * Vmyo / Vsr;
	koff_csqn = 65.0;
	MaxSR = 15.0;
	MinSR = 1.0;

	J_ca_juncsl = 8.2413e-13;
	J_ca_slmyo = 3.7243e-12;
	J_na_juncsl = 1.8313e-14;
	//  J_na_juncsl = 1.8313e-10; // CHANGE: Simulink and C
	J_na_slmyo = 1.6386e-12;
	//   J_na_slmyo = 1.6386e-8; // CHANGE: Simulink and C

	//
	//  Nernst Potentials
	// time: time (msec)

	ibarca_j = pCa * 4.0 * Vm * Frdy * FoRT * (0.341 * Y[1] * exp(2.0 * Vm * FoRT) - 0.341 * Cao) / (exp(2.0 * Vm * FoRT) - 1.0);
	I_Ca_junc = Fjunc_CaL * ibarca_j * Y[10] * Y[11] * (1.0 - Y[12] + fcaCaj) * pow(Q10CaL, Qpow) * 0.45 * 1.0;
	eca_junc = 1.0 / FoRT / 2.0 * log(Cao / Y[1]);
	I_cabk_junc = Fjunc * GCaB * (Vm - eca_junc);
	I_pca_junc = Fjunc * pow(Q10SLCaP, Qpow) * IbarSLCaP * pow(Y[1], 1.6) / (pow(KmPCa, 1.6) + pow(Y[1], 1.6));
	Ka_junc = 1.0 / (1.0 + pow(Kdact / Y[1], 2.0));
	s1_junc = exp(nu * Vm * FoRT) * pow(Y[31], 3.0) * Cao;
	s2_junc = exp((nu - 1.0) * Vm * FoRT) * pow(Nao, 3.0) * Y[1];
	s3_junc = KmCai * pow(Nao, 3.0) * (1.0 + pow(Y[31] / KmNai, 3.0)) + pow(KmNao, 3.0) * Y[1] * (1.0 + Y[1] / KmCai) + KmCao * pow(Y[31], 3.0) + pow(Y[31], 3.0) * Cao + pow(Nao, 3.0) * Y[1];
	I_ncx_junc = Fjunc * IbarNCX * pow(Q10NCX, Qpow) * Ka_junc * (s1_junc - s2_junc) / s3_junc / (1.0 + ksat * exp((nu - 1.0) * Vm * FoRT));
	I_Ca_tot_junc = I_Ca_junc + I_cabk_junc + I_pca_junc - 2.0 * I_ncx_junc;
	ibarca_sl = pCa * 4.0 * Vm * Frdy * FoRT * (0.341 * Y[2] * exp(2.0 * Vm * FoRT) - 0.341 * Cao) / (exp(2.0 * Vm * FoRT) - 1.0);
	I_Ca_sl = Fsl_CaL * ibarca_sl * Y[10] * Y[11] * (1.0 - Y[13] + fcaCaMSL) * pow(Q10CaL, Qpow) * 0.45 * 1.0;
	eca_sl = 1.0 / FoRT / 2.0 * log(Cao / Y[2]);
	I_cabk_sl = Fsl * GCaB * (Vm - eca_sl);
	I_pca_sl = Fsl * pow(Q10SLCaP, Qpow) * IbarSLCaP * pow(Y[2], 1.6) / (pow(KmPCa, 1.6) + pow(Y[2], 1.6));
	Ka_sl = 1.0 / (1.0 + pow(Kdact / Y[2], 2.0));
	s1_sl = exp(nu * Vm * FoRT) * pow(Y[32], 3.0) * Cao;
	s2_sl = exp((nu - 1.0) * Vm * FoRT) * pow(Nao, 3.0) * Y[2];
	s3_sl = KmCai * pow(Nao, 3.0) * (1.0 + pow(Y[32] / KmNai, 3.0)) + pow(KmNao, 3.0) * Y[2] * (1.0 + Y[2] / KmCai) + KmCao * pow(Y[32], 3.0) + pow(Y[32], 3.0) * Cao + pow(Nao, 3.0) * Y[2];
	I_ncx_sl = Fsl * IbarNCX * pow(Q10NCX, Qpow) * Ka_sl * (s1_sl - s2_sl) / s3_sl / (1.0 + ksat * exp((nu - 1.0) * Vm * FoRT));
	I_Ca_tot_sl = I_Ca_sl + I_cabk_sl + I_pca_sl - 2.0 * I_ncx_sl;
	J_CaB_junction = kon_sll * Y[1] * (Bmax_SLlowj - Y[25]) - koff_sll * Y[25] + (kon_slh * Y[1] * (Bmax_SLhighj - Y[23]) - koff_slh * Y[23]);
	J_SRCarel = ks * Y[36] / 1.0 * (Y[33] - Y[1]);

	J_SRleak = 5.348e-6 * (1.0 + 0.25 * AF) * (Y[33] - Y[1]);

	dY[1] = -I_Ca_tot_junc * Cmem / (Vjunc * 2.0 * Frdy) + J_ca_juncsl / Vjunc * (Y[2] - Y[1]) - J_CaB_junction + J_SRCarel * Vsr / Vjunc + J_SRleak * Vmyo / Vjunc;
	J_CaB_sl = kon_sll * Y[2] * (Bmax_SLlowsl - Y[26]) - koff_sll * Y[26] + (kon_slh * Y[2] * (Bmax_SLhighsl - Y[24]) - koff_slh * Y[24]);
	dY[2] = -I_Ca_tot_sl * Cmem / (Vsl * 2.0 * Frdy) + J_ca_juncsl / Vsl * (Y[1] - Y[2]) + J_ca_slmyo / Vsl * (Y[0] - Y[2]) - J_CaB_sl;
	J_serca = pow(Q10SRCaP, Qpow) * Vmax_SRCaP * (pow(Y[0] / Kmf, hillSRCaP) - pow(Y[33] / Kmr, hillSRCaP)) / (1.0 + pow(Y[0] / Kmf, hillSRCaP) + pow(Y[33] / Kmr, hillSRCaP));
	J_CaB_cytosol = kon_tncl * Y[0] * (Bmax_TnClow - Y[9]) - koff_tncl * Y[9] + kon_tnchca * Y[0] * (Bmax_TnChigh - Y[7] - Y[8]) - koff_tnchca * Y[7] + kon_tnchmg * Mgi * (Bmax_TnChigh - Y[7] - Y[8]) - koff_tnchmg * Y[8] + kon_cam * Y[0] * (Bmax_CaM - Y[3]) - koff_cam * Y[3] + kon_myoca * Y[0] * (Bmax_myosin - Y[4] - Y[5]) - koff_myoca * Y[4] + kon_myomg * Mgi * (Bmax_myosin - Y[4] - Y[5]) - koff_myomg * Y[5] + (kon_sr * Y[0] * (Bmax_SR - Y[6]) - koff_sr * Y[6]);
	dY[0] = -J_serca * Vsr / Vmyo - J_CaB_cytosol + J_ca_slmyo / Vmyo * (Y[2] - Y[0]);
	dY[9] = kon_tncl * Y[0] * (Bmax_TnClow - Y[9]) - koff_tncl * Y[9];
	dY[7] = kon_tnchca * Y[0] * (Bmax_TnChigh - Y[7] - Y[8]) - koff_tnchca * Y[7];
	dY[8] = kon_tnchmg * Mgi * (Bmax_TnChigh - Y[7] - Y[8]) - koff_tnchmg * Y[8];
	dY[3] = kon_cam * Y[0] * (Bmax_CaM - Y[3]) - koff_cam * Y[3];
	dY[4] = kon_myoca * Y[0] * (Bmax_myosin - Y[4] - Y[5]) - koff_myoca * Y[4];
	dY[5] = kon_myomg * Mgi * (Bmax_myosin - Y[4] - Y[5]) - koff_myomg * Y[5];
	dY[6] = kon_sr * Y[0] * (Bmax_SR - Y[6]) - koff_sr * Y[6];

	/* atrial ICaL */
	dss = 1.0 / (1.0 + exp(-(Vm + 3.0 * ISO + 9.0) / 6.0));
	taud = 1.0 * dss * (1.0 - exp(-(Vm + 3.0 * ISO + 9.0) / 6.0)) / (0.035 * (Vm + 3.0 * ISO + 9.0));
	fss = 1.0 / (1.0 + exp((Vm + 3.0 * ISO + 30.0) / 7.0)) + 0.2 / (1.0 + exp((50.0 - Vm - 3.0 * ISO) / 20.0));
	tauf = 1.0 / (0.0197 * exp(-pow(0.0337 * (Vm + 3.0 * ISO + 25.0), 2.0)) + 0.02);
	dY[10] = (dss - Y[10]) / taud;   // derivative of d
	dY[11] = (fss - Y[11]) / tauf;   // derivative of f
	dY[12] = 1.7 * Y[1] / 1.0 * (1.0 - Y[12]) - 11.9e-3 * Y[12]; // derivative of fca_bj
	dY[13] = 1.7 * Y[2] / 1.0 * (1.0 - Y[13]) - 11.9e-3 * Y[13]; // derivative of fca_bsl
	ibark = pK * Vm * Frdy * FoRT * (0.75 * Y[27] * exp(Vm * FoRT) - 0.75 * Ko) / (exp(Vm * FoRT) - 1.0);
	ibarna_j = pNa * Vm * Frdy * FoRT * (0.75 * Y[31] * exp(Vm * FoRT) - 0.75 * Nao) / (exp(Vm * FoRT) - 1.0);
	ibarna_sl = pNa * Vm * Frdy * FoRT * (0.75 * Y[32] * exp(Vm * FoRT) - 0.75 * Nao) / (exp(Vm * FoRT) - 1.0);
	I_Ca = I_Ca_junc + I_Ca_sl;
	/***********************/

	I_CaK = ibark * Y[10] * Y[11] * (Fjunc_CaL * (fcaCaj + (1.0 - Y[12])) + Fsl_CaL * (fcaCaMSL + (1.0 - Y[13]))) * pow(Q10CaL, Qpow) * 0.45 * 1.0;
	I_CaNa_junc = Fjunc_CaL * ibarna_j * Y[10] * Y[11] * (1.0 - Y[12] + fcaCaj) * pow(Q10CaL, Qpow) * 0.45 * 1.0;
	I_CaNa_sl = Fsl_CaL * ibarna_sl * Y[10] * Y[11] * (1.0 - Y[13] + fcaCaMSL) * pow(Q10CaL, Qpow) * 0.45 * 1.0;
	I_CaNa = I_CaNa_junc + I_CaNa_sl;
	I_Catot = I_Ca + I_CaK + I_CaNa;
	I_cabk = I_cabk_junc + I_cabk_sl;
	I_ClCa_junc = Fjunc * GClCa / (1.0 + KdClCa / Y[1]) * (Vm - ecl);
	I_ClCa_sl = Fsl * GClCa / (1.0 + KdClCa / Y[2]) * (Vm - ecl);
	I_ClCa = I_ClCa_junc + I_ClCa_sl;
	I_Clbk = GClB * (Vm - ecl);

	/* I_ki */
	ek = 1.0 / FoRT * log(Ko / Y[27]);
	aki = 1.02 / (1.0 + exp(0.2385 * (Vm - ek - 59.215)));
	bki = (0.49124 * exp(0.08032 * (Vm + 5.476 - ek)) + exp(0.06175 * (Vm - ek - 594.31))) / (1.0 + exp(-0.5143 * (Vm - ek + 4.753)));
	kiss = aki / (aki + bki);
	I_ki = (1.0 + 1.0 * AF) * 0.0525 * sqrt(Ko / 5.4) * kiss * (Vm - ek);



	kp_kp = 1.0 / (1.0 + exp(7.488 - Vm / 5.98));
	I_kp_junc = Fjunc * gkp * kp_kp * (Vm - ek);
	I_kp_sl = Fsl * gkp * kp_kp * (Vm - ek);
	I_kp = I_kp_junc + I_kp_sl;


	// I_kr
	gkr = 0.035 * sqrt(Ko / 5.4);
	xrss = 1.0 / (1.0 + exp(-(Vm + 10.0) / 5.0));
	tauxr = 550.0 / (1.0 + exp((-22.0 - Vm) / 9.0)) * 6.0 / (1.0 + exp((Vm - (-11.0)) / 9.0)) + 230.0 / (1.0 + exp((Vm - (-40.0)) / 20.0));
	dY[14] = (xrss - Y[14]) / tauxr;
	rkr = 1.0 / (1.0 + exp((Vm + 74.0) / 24.0));
	I_kr = gkr * Y[14] * rkr * (Vm - ek);

	/*I_ks HH model*/
	gks_junc = 1.0 * (1.0 + 1.0 * AF + 2.0 * ISO) * 0.0035;
	gks_sl = 1.0 * (1.0 + 1.0 * AF + 2.0 * ISO) * 0.0035; // %FRA
	eks = 1.0 / FoRT * log((Ko + pNaK * Nao) / (Y[27] + pNaK * Y[30]));
	xsss = 1.0 / (1.0 + exp(-(Vm + 40.0 * ISO + 3.8) / 14.25));
	tauxs = 990.1 / (1.0 + exp(-(Vm + 40.0 * ISO + 2.436) / 14.12));
	dY[15] = (xsss - Y[15]) / tauxs;
	I_ks_junc = Fjunc * gks_junc * pow(Y[15], 2.0) * (Vm - eks);
	I_ks_sl = Fsl * gks_sl * pow(Y[15], 2.0) * (Vm - eks);
	I_ks = iksmultiply * (I_ks_junc + I_ks_sl);
	I_ncx = I_ncx_junc + I_ncx_sl;


	/* I_Na */
	mss = 1.0 / pow(1.0 + exp(-(56.86 + Vm) / 9.03), 2.0);
	taum = 0.1292 * exp(-pow((Vm + 45.79) / 15.54, 2.0)) + 0.06487 * exp(-pow((Vm - 4.823) / 51.12, 2.0));

	/*if (Vm >= -40.0)
	    ah = 0.0;
	else
	    ah = 0.057 * exp(-(Vm + 80.0) / 6.8);*/

	ah = Vm >= -40.0 ? 0.0 : 0.057 * exp(-(Vm + 80.0) / 6.8);

	if (Vm >= -40.0)
		bh = 0.77 / (0.13 * (1.0 + exp(-(Vm + 10.66) / 11.1)));
	else
		bh = 2.7 * exp(0.079 * Vm) + 3.1e5 * exp(0.3485 * Vm);

	tauh = 1.0 / (ah + bh);
	//  hss = 1.0/pow(1.0+exp((Vm+71.55)/7.43), 2.0);
	hss = (1.0 - ina_c) / pow(1 + exp((Vm + 71.55) / 7.43), 2.0) + ina_c;

	if (Vm >= -40.0)
		aj = 0.0;
	else
		aj = (-2.5428e4 * exp(0.2444 * Vm) - 6.948e-6 * exp(-0.04391 * Vm)) * (Vm + 37.78) / (1.0 + exp(0.311 * (Vm + 79.23)));

	if (Vm >= -40.0)
		bj = 0.6 * exp(0.057 * Vm) / (1.0 + exp(-0.1 * (Vm + 32.0)));
	else
		bj = 0.02424 * exp(-0.01052 * Vm) / (1.0 + exp(-0.1378 * (Vm + 40.14)));

	tauj = 1.0 / (aj + bj);
	//   jss = 1.0/pow(1.0+exp((Vm+71.55)/7.43), 2.0);
	jss = (1.0 - ina_c) / pow(1 + exp((Vm + 71.55) / 7.43), 2.0) + ina_c;

	dY[18] = (mss - Y[18]) / taum;
	dY[16] = (hss - Y[16]) / tauh;
	dY[17] = (jss - Y[17]) / tauj;

	ena_junc = 1.0 / FoRT * log(Nao / Y[31]);
	I_Na_junc = Fjunc * GNa * pow(Y[18], 3.0) * Y[16] * Y[17] * (Vm - ena_junc);
	ena_sl = 1.0 / FoRT * log(Nao / Y[32]);
	I_Na_sl = Fsl * GNa * pow(Y[18], 3.0) * Y[16] * Y[17] * (Vm - ena_sl);
	I_Na = I_Na_junc + I_Na_sl;


	/*late I_Na, only in AF condition*/
	// Matlab code, change to C when needed.
	/*//% Late I_Na
	GNaL=0.0025*AF;
	aml = 0.32 * (Vm + 47.13) / (1.0 - exp(-0.1 * (Vm + 47.13)));
	bml = 0.08 * exp(-Vm) / 11.0);
	hlinf = 1.0 / (1.0 + exp((Vm + 91.0) / 6.1));
	tauhl = 600.0;
	ydot(60) = aml * (1.0 - y(60)) - bml * y(60);
	ydot(61) = (hlinf - y(61)) / tauhl;

	I_NaL_junc = Fjunc * GNaL * pow(y(60), 3.0) * y(61) * (y(39) - ena_junc);
	I_NaL_sl   = Fsl * GNaL * pow(y(60), 3.0) * y(61) * (y(39) - ena_sl);
	I_NaL      = I_NaL_junc + I_NaL_sl;

	if t<9050
	    ydot(62)=0;
	else
	    ydot(62)=I_NaL;
	end

	*/



	I_nabk_junc = Fjunc * GNaB * (Vm - ena_junc);
	I_nabk_sl = Fsl * GNaB * (Vm - ena_sl);
	I_nabk = I_nabk_junc + I_nabk_sl;



	// I_nak % I_nak: Na/K Pump Current
	sigma = (exp(Nao / 67.3) - 1.0) / 7.0;
	fnak = 1.0 / (1.0 + 0.1245 * exp(-0.1 * Vm * FoRT) + 0.0365 * sigma * exp(-Vm * FoRT));
	I_nak_junc = Fjunc * IbarNaK * fnak * Ko / (1.0 + pow(KmNaip / Y[31], 4.0)) / (Ko + KmKo);
	I_nak_sl = Fsl * IbarNaK * fnak * Ko / (1.0 + pow(KmNaip / Y[32], 4.0)) / (Ko + KmKo);
	I_nak = I_nak_junc + I_nak_sl;



	I_pca = I_pca_junc + I_pca_sl;



	xtoss = 1.0 / (1.0 + exp(-(Vm + 1.0) / 11.0)); // Haibo
	ytoss = 1.0 / (1.0 + exp((Vm + 40.5) / 11.5)); // Haibo
	tauxtos = 9.0 / (1.0 + exp((Vm + 3.0) / 15.0)) + 0.5;
	tauytos = 800.0 / (1.0 + exp((Vm + 60.0) / 10.0)) + 30.0;
	dY[20] = (xtoss - Y[20]) / tauxtos;
	dY[22] = (ytoss - Y[22]) / tauytos;
	I_tos = 0.0 * GtoSlow * Y[20] * Y[22] * (Vm - ek); // Haibo I_to slow has been removed


	// 11/12/09; changed Itof to that from maleckar/giles/2009; removed I_tos
	// atrium
	// equations for activation;
	// the Ito_f in atrium, time constant, by haibo
	// tauxtof = 3.5 / (1.0 + exp(- pow((Vm / 30.0), 2.0))) + 1.5;    // for atrium  // former version, not right, haibo
	tauxtof = 3.5 * ( exp(- pow((Vm / 30.0), 2.0))) + 1.5;    // for atrium
	// tauytof = 85.0 * exp(-pow(Vm + 52.45, 2.0) / 15.8) + 24.14; // check haibo wrong version
	tauytof = 25.635 * exp(-(pow((Vm + 52.45) / 15.8827, 2.0))) + 24.14; //%14.14  from matlab version
	GtoFast = (1.0 - 0.7 * AF) * 0.165; // to be careful
	dY[19] = (xtoss - Y[19]) / tauxtof; // activation of Itof
	dY[21] = (ytoss - Y[21]) / tauytof; // inactivation of Itof by haibo

	I_tof = GtoFast * Y[19] * Y[21] * (Vm - ek);
	I_to = I_tos + I_tof;
	dY[25] = kon_sll * Y[1] * (Bmax_SLlowj - Y[25]) - koff_sll * Y[25];
	dY[26] = kon_sll * Y[2] * (Bmax_SLlowsl - Y[26]) - koff_sll * Y[26];
	dY[23] = kon_slh * Y[1] * (Bmax_SLhighj - Y[23]) - koff_slh * Y[23];
	dY[24] = kon_slh * Y[2] * (Bmax_SLhighsl - Y[24]) - koff_slh * Y[24];



	// if (Ikur_model_type == 0)
	// {
	// 	//// I_kur: Ultra rapid delayed rectifier Outward K Current
	// 	/*  Equation for IKur; New CNZ model.
	// 	*   atrium equations for activation;
	// 	*   modified by haibo
	// 	*   */
	// 	// Gkur = 0.0073235294117647053;  // fit v = 40 mv max current
	// 	Gkur = 0.0086220472440944884;  // fit v = 20 mv max current
	// 	Gkur = 1.0 * (1.0 - 0.5 * AF) * (1.0 + 2.0 * ISO) * Gkur * (1.0 + 0.2 * RA);
	// 	double V = Vm;
	// 	double K_Q10 = 3.5308257834747638;
	// 	//CNZ_a
	// 	ac_Ikur_inf = ((IKur_ac1_mult * 1.0) / (1 + exp((V + 17.6684 + IKur_ac1_shift) / (-5.75418 * IKur_ac1_grad))) ) * ((IKur_ac2_mult * 1.0) / (1 + exp((V + 8.4153 + IKur_ac2_shift) / (-11.51037561 * IKur_ac2_grad)))) + IKur_ac_add;
	// 	tau_ac_Ikur = (45.6666746826 / (1 + exp((V + 11.2306497073) / 11.5254705962)) + 4.26753514993)
	// 	              * (0.262186042981 / (1 + exp((V + 35.8658312707) / (-3.87510627762))) + 0.291755017928); //
	// 	tau_ac_Ikur = tau_ac_Ikur / K_Q10;

	// 	// ac = inf + (ac - inf) * exp(-(dt) / tau);
	// 	dY[39] = (ac_Ikur_inf - Y[39]) / tau_ac_Ikur;

	// 	// CNZ_i
	// 	inac_Ikur_inf = (IKur_inac_mult * 0.52424) / (1.0 + exp((V + 15.1142 + IKur_inac_shift) / (7.567021 * IKur_inac_grad))) + 0.4580778 + IKur_inac_add;
	// 	tau_inac_Ikur = 2328 / (1 + exp(((V) - 9.435) / (3.5827))) + 1739.139;
	// 	tau_inac_Ikur = tau_inac_Ikur / K_Q10;
	// 	dY[40] = (inac_Ikur_inf - Y[40]) / tau_inac_Ikur;

	// 	Ikur = Gkur * Ikur_cnd_mult * (4.5128 + 1.899769 / (1.0 + exp((V - 20.5232) / (-8.26597)))) * Y[39] * Y[40] * (V - ek);

	// }
	// else if ( Ikur_model_type == 1)
	// {
	// 	//// I_kur: Ultra rapid delayed rectifier Outward K Current
	// 	/* Ikur formulation 2, based on the original grandi model Ikur(from mal..2009)
	// 	*   take mutation effects into consideration.
	// 	*   modified by haibo
	// 	*     */
	// 	// original comments
	// 	//Equation for IKur; from Maleckar et al. 2009 - EG
	// 	//atrium equations for activation;
	// 	Gkur = 1.0 * (1.0 - 0.5 * AF) * (1.0 + 2.0 * ISO) * 0.045 * (1.0 + 0.2 * RA); //*(1+0.2*RA); //nS/pF maleckar 0.045
	// 	// Gkur = 0.5 * Gkur;
	// 	tau_ac_Ikur = 9.0 / (1.0 + exp((Vm + 5.0) / 12.0)) + 0.5;
	// 	tau_inac_Ikur = 590.0 / (1.0 +  exp((Vm + 60.0) / 10.0)) + 3050.0;
	// 	// shift and change slope according to the effects of the mutations, by haibo
	// 	ac_Ikur_inf = 1.0 / (1.0 + exp(-(Vm + mal_IKur_ac_shift + 6.0) / (8.6 * mal_IKur_ac_grad)));
	// 	inac_Ikur_inf = 1.0 / (1.0 + exp((Vm + mal_IKur_inac_shift + 7.5) / (10.0 * mal_IKur_inac_grad)));
	// 	dY[39] = (ac_Ikur_inf - Y[39]) / tau_ac_Ikur;
	// 	dY[40] = (inac_Ikur_inf - Y[40]) / tau_inac_Ikur;
	// 	Ikur = mal_mult * Gkur * Y[39] * Y[40] * (Vm - ek);
	// 	// Ikur = Gkur * Y[39] * Y[40] * (Vm - ek);
	// }



	Gkur = 1.0 * (1.0 - 0.5 * AF) * (1.0 + 2.0 * ISO) * 0.045 * (1.0 + 0.2 * RA); //*(1+0.2*RA); //nS/pF maleckar 0.045
	// Gkur = 0.5 * Gkur;
	tau_ac_Ikur = 9.0 / (1.0 + exp((Vm + 5.0) / 12.0)) + 0.5;
	tau_inac_Ikur = 590.0 / (1.0 +  exp((Vm + 60.0) / 10.0)) + 3050.0;
	// shift and change slope according to the effects of the mutations, by haibo
	ac_Ikur_inf = 1.0 / (1.0 + exp(-(Vm + Data->Simple_Ikur_ac_shift + 6.0) / (8.6 * Data->Simple_Ikur_ac_grad)));
	inac_Ikur_inf = 1.0 / (1.0 + exp((Vm + Data->Simple_Ikur_inac_shift + 7.5) / (10.0 * Data->Simple_Ikur_inac_grad)));
	dY[39] = (ac_Ikur_inf - Y[39]) / tau_ac_Ikur;
	dY[40] = (inac_Ikur_inf - Y[40]) / tau_inac_Ikur;
	Ikur = Data->Simple_GKur * Gkur * Y[39] * Y[40] * (Vm - ek);


	I_K_tot = I_to + I_kr + I_ks + I_ki - 2.0 * I_nak + I_CaK + I_kp + Ikur;
	dY[27] = 0.0;
	dNa_Bj_dt = kon_na * Y[31] * (Bmax_Naj - Y[28]) - koff_na * Y[28];
	dY[28] = dNa_Bj_dt;
	dNa_Bsl_dt = kon_na * Y[32] * (Bmax_Nasl - Y[29]) - koff_na * Y[29];
	dY[29] = dNa_Bsl_dt;
	I_Na_tot_junc = I_Na_junc + I_nabk_junc + 3.0 * I_ncx_junc + 3.0 * I_nak_junc + I_CaNa_junc;
	I_Na_tot_sl = I_Na_sl + I_nabk_sl + 3.0 * I_ncx_sl + 3.0 * I_nak_sl + I_CaNa_sl;
	I_Na_tot_sl2 = 3.0 * I_ncx_sl + 3.0 * I_nak_sl + I_CaNa_sl;
	I_Na_tot_junc2 = 3.0 * I_ncx_junc + 3.0 * I_nak_junc + I_CaNa_junc;
	dY[31] = -I_Na_tot_junc * Cmem / (Vjunc * Frdy) + J_na_juncsl / Vjunc * (Y[32] - Y[31]) - dNa_Bj_dt;
	dY[32] = -I_Na_tot_sl * Cmem / (Vsl * Frdy) + J_na_juncsl / Vsl * (Y[31] - Y[32]) + J_na_slmyo / Vsl * (Y[30] - Y[32]) - dNa_Bsl_dt;
	dY[30] = J_na_slmyo / Vmyo * (Y[32] - Y[30]);
	dY[34] = kon_csqn * Y[33] * (Bmax_Csqn - Y[34]) - koff_csqn * Y[34];
	dY[33] = J_serca - (J_SRleak * Vmyo / Vsr + J_SRCarel) - (kon_csqn * Y[33] * (Bmax_Csqn - Y[34]) - koff_csqn * Y[34]);
	kCaSR = MaxSR - (MaxSR - MinSR) / (1.0 + pow(ec50SR / Y[33], 2.5));
	koSRCa = koCa / kCaSR;
	kiSRCa = kiCa * kCaSR;
	RI = 1.0 - Y[37] - Y[36] - Y[35];
	dY[37] = kim * RI - kiSRCa * Y[1] * Y[37] - (koSRCa * pow(Y[1], 2.0) * Y[37] - kom * Y[36]);
	dY[36] = koSRCa * pow(Y[1], 2.0) * Y[37] - kom * Y[36] - (kiSRCa * Y[1] * Y[36] - kim * Y[35]);
	dY[35] = kiSRCa * Y[1] * Y[36] - kim * Y[35] - (kom * Y[35] - koSRCa * pow(Y[1], 2.0) * RI);





	/* I_clcfTR */
	GClCFTR = 0; //4.9e-3 * ISO;
	I_ClCFTR = GClCFTR * (Vm - ecl);

	I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;
	I_Cl_tot = I_ClCa + I_Clbk + I_ClCFTR;
	I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
	I_tot = I_Na_tot + I_Cl_tot + I_Ca_tot + I_K_tot;
	Data->dV = -(I_tot - Data->Istim);

	// for (i = 0; i < NEQ; i++)
	//     Ith(ydot, i + 1) = dY[i];
	// return (0);
}