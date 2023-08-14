#ifndef BETAAR__HPP
#define BETAAR__HPP
#include  <iomanip>
#include "signalling_para.hpp"
#include <cmath>
class betaAR_Update
{
public:
	betaAR_Update() {

		ISO      = 0;
		LCCtot   = 0.0250000000000000;
		RyRtot   = 0.135000000000000;
		PLBtot   = 38;
		TnItot   = 70;
		IKstot   = 0.0250000000000000;
		ICFTRtot = 0.0250000000000000;
		PP1tot   = 0.890000000000000;
		PLMtot   = 48;
		Myotot   = 70;
		IKrtot   = 0.0250000000000000;
		IKurtot  = 0.0250000000000000;
		INatot   = 0.0250000000000000;
		IClCatot = 0.0250000000000000;
		Itotot   = 0.0250000000000000;
		IK1tot   = 0.0250000000000000;
		AF_index = 0;
	};
	~betaAR_Update() {
		y = nullptr;
		ydot = nullptr;
	};

	// double y[40], ydot[40];
	double *y, *ydot;
	double ISO;
	double LCCtot;   // pin(2);
	double RyRtot;   // pin(3);
	double PLBtot;   // pin(4);
	double TnItot;   // pin(5);
	double IKstot;   // pin(6);
	double ICFTRtot;   // pin(7);
	double PP1tot;   // pin(8);
	double PLMtot;   // pin(9);
	double Myotot;   // pin(10);
	double IKrtot;   // pin(11);
	double IKurtot;   // pin(12);
	double INatot;   // pin(13);
	double IClCatot;   // pin(14);
	double Itotot;   // pin(15);
	double IK1tot;   // pin(16);
	double AF_index;   // pin(17); // AF_index



	const int ODE_NUM = 41;


	double LCCtot_Scale        = 1.0;
	double RyRtot_Scale        = 1.0;
	double PLBtot_Scale        = 1.0;
	double TnItot_Scale        = 1.0;
	double PLMtot_Scale        = 1.0;
	double b1ARtot_Scale       = 1.0;
	double Gstot_Scale         = 1.0;
	double ACtot_Scale         = 1.0;
	double ATP_Scale           = 1.0;
	double PDE3tot_Scale       = 1.0;
	double PDE4tot_Scale       = 1.0;
	double PKItot_Scale        = 1.0;
	double PKAIItot_Scale      = 1.0;
	double I1tot_Scale         = 1.0;
	double PP1tot_Scale        = 1.0;
	double PKACII_LCCtot_Scale = 1.0;
	double PP1_LCC_Scale       = 1.0;
	double PP2A_LCC_Scale      = 1.0;
	double PKAIIryrtot_Scale   = 1.0;
	double PP1ryr_Scale        = 1.0;
	double PP2Aryr_Scale       = 1.0;
	// double PP1_ikstot_Scale    = 1.0;
	double PKAII_ikstot_Scale  = 1.0;
	double Km_AC_basal_Scale   = 1.0;
	double Kd_AC_Gsa_Scale     = 1.0;
	double Km_PDE3_cAMP_Scale  = 1.0;
	double Km_PDE4_cAMP_Scale  = 1.0;
	double Km_PKA_I1_Scale     = 1.0;
	double Km_PP2A_I1_Scale    = 1.0;
	double Km_PKA_PLB_Scale    = 1.0;
	double Km_PP1_PLB_Scale    = 1.0;
	double Km_PKA_LCC_Scale    = 1.0;
	double Km_PP1_LCC_Scale    = 1.0;
	double Km_PP2A_LCC_Scale   = 1.0;
	double Km_pka_ryr_Scale    = 1.0;
	double Km_pp1_ryr_Scale    = 1.0;
	double Km_pp2a_ryr_Scale   = 1.0;
	double Km_PKA_TnI_Scale    = 1.0;
	double Km_PP2A_TnI_Scale   = 1.0;
	double Km_pka_iks_Scale    = 1.0;
	double Km_pp1_iks_Scale    = 1.0;
	double K_rate_Scale = 1.0;

	double PP1_ik1tot_Scale = 1.0;
	double PP1_itotot_Scale = 1.0;
	double PP1_inatot_Scale = 1.0;
	double PP1_ClCatot_Scale = 1.0;
	double PP1_CFTRtot_Scale = 1.0;
	double PP1_KURtot_Scale = 1.0;
	double PP1_ikrtot_Scale = 1.0;
	double PP1_ikstot_Scale = 1.0;

	double PP2A_TnI_Scale = 1.0;

	double PP2A_I1_Vmax_Scale = 1.0;

	int update_single_time_step(betaAR_para & para) {

		//  This module describes the beta-adrenergic signaling pathway.

		// //  State variables

		ISO      = para.ISO      ;
		LCCtot   = LCCtot_Scale * para.LCCtot   ; // pin(2);
		RyRtot   = RyRtot_Scale * para.RyRtot   ; // pin(3);
		PLBtot   = PLBtot_Scale * para.PLBtot   ; // pin(4);
		TnItot   = TnItot_Scale * para.TnItot   ; // pin(5);
		IKstot   = para.IKstot   ;   // pin(6);
		ICFTRtot = para.ICFTRtot ;   // pin(7);
		PP1tot   = PP1tot_Scale * para.PP1tot   ; // pin(8);
		PLMtot   = PLMtot_Scale * para.PLMtot   ; // pin(9);
		Myotot   = para.myotot   ;   // pin(10);
		IKrtot   = para.IKrtot   ;   // pin(11);
		IKurtot  = para.IKurtot  ;   // pin(12);
		INatot   = para.INatot   ;   // pin(13);
		IClCatot = para.IClCatot ;   // pin(14);
		Itotot   = para.Itotot   ;   // pin(15);
		IK1tot   = para.IK1tot   ;   // pin(16);
		AF_index = para.AF_index ;   // pin(17); // AF_index

		// double LR            = y[0];
		// double LRG           = y[1];
		// double RG            = y[2];
		// double b1AR_S464     = y[3];
		// double b1AR_S301     = y[4];
		// double GsaGTPtot     = y[5];
		// double GsaGDP        = y[6];
		// double Gsby          = y[7];
		// double AC_GsaGTP     = y[8];
		// // PDEp              =y(10);9
		// double PDE3p         = y[9];
		// double PDE4p         = y[10];
		// double cAMPtot       = y[11];
		// double RC_I          = y[12];
		// double RCcAMP_I      = y[13];
		// double RCcAMPcAMP_I  = y[14];
		// double RcAMPcAMP_I   = y[15];
		// double PKACI         = y[16];
		// double PKACI_PKI     = y[17];
		// double RC_II         = y[18];
		// double RCcAMP_II     = y[19];
		// double RCcAMPcAMP_II = y[20];
		// double RcAMPcAMP_II  = y[21];
		// double PKACII        = y[22];
		// double PKACII_PKI    = y[23];
		// double I1p_PP1       = y[24]; //  output CaMKII
		// double I1ptot        = y[25];
		// double PLBp          = y[26]; //  output
		// double PLMp          = y[27]; //  output
		// double LCCap         = y[28]; //  output
		// double LCCbp         = y[29]; //  output
		// double RyRp          = y[30]; //  output
		// double TnIp          = y[31]; //  output
		// double Myop          = y[32]; //  output
		// double KSp           = y[33]; //  output
		// double KRp           = y[34]; //  output
		// double CFTRp         = y[35]; //  output
		// double ClCap         = y[36]; //  output

		// double NAp           = y[37]; // output
		// double TOp           = y[38]; // output
		// double K1p           = y[39]; // output
		// double KURp          = y[40]; // output
		double LR            = y[0]; // bound receptor
		double LRG           = y[1]; // G protein-associated, ligand bound receptor
		double RG            = y[2]; //
		double b1AR_S464     = y[3];
		double b1AR_S301     = y[4];
		double GsaGTPtot     = y[5];
		double GsaGDP        = y[6];
		double Gsby          = y[7];
		double AC_GsaGTP     = y[8];
		double PDE4p          = y[9];
		double cAMPtot       = y[10];
		double RC_I          = y[11];
		double RCcAMP_I      = y[12];
		double RCcAMPcAMP_I  = y[13];
		double RcAMPcAMP_I   = y[14];
		double PKACI         = y[15];  // pka_1
		double PKACI_PKI     = y[16];
		double RC_II         = y[17];
		double RCcAMP_II     = y[18];
		double RCcAMPcAMP_II = y[19];
		double RcAMPcAMP_II  = y[20];
		double PKACII        = y[21];  // pka_2
		double PKACII_PKI    = y[22];
		double I1p_PP1       = y[23]; // output CaMKII
		double I1ptot        = y[24];
		double PLBp          = y[25]; // output
		double PLMp          = y[26]; // output
		double LCCap         = y[27]; // output
		double LCCbp         = y[28]; // output
		double RyRp          = y[29]; // output
		double TnIp          = y[30]; // output
		double KSp           = y[31]; // output
		double KRp           = y[32]; // output
		double KURp          = y[33]; // output
		double CFTRp         = y[34]; // output
		double ClCap         = y[35]; // output
		double Myop          = y[36]; // output
		double NAp           = y[37]; // output
		double TOp           = y[38]; // output
		double K1p           = y[39]; // output
		// double PDE3p         = y[9];
		double PDE3p         = y[40];  // added PDE3p for rabbit and human mdoels // 16:59:21, Wed, 18-September-2019, By Haibo
		// //  Drug Concentrations

		// ISO = pin(1); //  (uM) isoproterenol concentration - Ltot
		double FSK = 0; //  (uM) forskolin concentration
		double IBMX = 0; //  (uM) IBMX concentration
		// //  b-AR module

		// b1ARtot = 0.00528;        //  (uM) total b1-AR protein //  MOUSE
		double b1ARtot = b1ARtot_Scale * 0.028; //  RABBIT

		double kf_LR           = K_rate_Scale * 1;              //  (1/[uM ms]) forward rate for ISO binding to b1AR
		double kr_LR           = K_rate_Scale * 0.285;          //  (1/ms) reverse rate for ISO binding to b1AR
		double kf_LRG          = K_rate_Scale * 1;              //  (1/[uM ms]) forward rate for ISO:b1AR association with Gs
		double kr_LRG          = K_rate_Scale * 0.062;          //  (1/ms) reverse rate for ISO:b1AR association with Gs
		double kf_RG           = K_rate_Scale * 1;              //  (1/[uM ms]) forward rate for b1AR association with Gs
		double kr_RG           = K_rate_Scale * 33;             //  (1/ms) reverse rate for b1AR association with Gs

		double Gstot           = Gstot_Scale * 3.83;         //  (uM) total Gs protein
		double k_G_act         = K_rate_Scale * 16e-3;        //  (1/ms) rate constant for Gs activation
		double k_G_hyd         = K_rate_Scale * 0.8e-3;       //  (1/ms) rate constant for G-protein hydrolysis
		double k_G_reassoc     = K_rate_Scale * 1.21;         //  (1/[uM ms]) rate constant for G-protein reassociation

		double kf_bARK         = K_rate_Scale * 1.1e-6;       //  (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
		double kr_bARK         = K_rate_Scale * 2.2e-6;       //  (1/ms) reverse rate for b1AR phosphorylation by b1ARK
		double kf_PKA          = K_rate_Scale * 3.6e-6;       //  (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
		double kr_PKA          = K_rate_Scale * 2.2e-6;       //  (1/ms) reverse rate for b1AR phosphorylation by PKA

		double b1ARact = b1ARtot - b1AR_S464 - b1AR_S301;
		double b1AR = b1ARact - LR - LRG - RG;
		double Gs = Gstot - LRG - RG - Gsby;

		double dLR = kf_LR * ISO * b1AR - kr_LR * LR + kr_LRG * LRG - kf_LRG * LR * Gs;
		double dLRG = kf_LRG * LR * Gs - kr_LRG * LRG - k_G_act * LRG;
		double dRG = kf_RG * b1AR * Gs - kr_RG * RG - k_G_act * RG;

		double bARK_desens = kf_bARK * (LR + LRG);
		double bARK_resens = kr_bARK * b1AR_S464;
		double PKA_desens = kf_PKA * PKACI * b1ARact;
		double PKA_resens = kr_PKA * b1AR_S301;
		double db1AR_S464 = bARK_desens - bARK_resens; //  ydot(5)
		double db1AR_S301 = PKA_desens - PKA_resens; //  ydot(6)

		double G_act = k_G_act * (RG + LRG);
		double G_hyd = k_G_hyd * GsaGTPtot;
		double G_reassoc = k_G_reassoc * GsaGDP * Gsby;
		double dGsaGTPtot = G_act - G_hyd; //  ydot(7)
		double dGsaGDP = G_hyd - G_reassoc; //  ydot(8)
		double dGsby = G_act - G_reassoc; //  ydot(9)
		//  end b-AR module
		// //  cAMP module

		// ACtot = 70.57e-3;       //  (uM) total adenylyl cyclase //  MOUSE
		double ACtot = ACtot_Scale * 47e-3; //  RABBIT

		double ATP             = ATP_Scale * 5e3;          //  (uM) total ATP
		double k_AC_basal      = K_rate_Scale * 0.2e-3;       //  (1/ms) basal cAMP generation rate by AC
		double Km_AC_basal     = Km_AC_basal_Scale * 1.03e3;       //  (uM) basal AC affinity for ATP

		double Kd_AC_Gsa       = 0.4;            //  (uM) Kd for AC association with Gsa
		double kf_AC_Gsa       = K_rate_Scale * 1;            //  (1/[uM ms]) forward rate for AC association with Gsa
		double kr_AC_Gsa       = K_rate_Scale * Kd_AC_Gsa;    //  (1/ms) reverse rate for AC association with Gsa

		double k_AC_Gsa        = K_rate_Scale * 8.5e-3;       //  (1/ms) basal cAMP generation rate by AC:Gsa
		double Km_AC_Gsa       = Kd_AC_Gsa_Scale * 315.0;        //  (uM) AC:Gsa affinity for ATP

		double Kd_AC_FSK       = 44.0;           //  (uM) Kd for FSK binding to AC
		double k_AC_FSK        = K_rate_Scale * 7.3e-3;       //  (1/ms) basal cAMP generation rate by AC:FSK
		double Km_AC_FSK       = 860.0;          //  (uM) AC:FSK affinity for ATP

		//  MOUSE
		// PDEtot          = 22.85e-3;       //  (uM) total phosphodiesterase
		// k_cAMP_PDE      = 5e-3;           //  (1/ms) cAMP hydrolysis rate by PDE
		// k_cAMP_PDEp     = 2*k_cAMP_PDE;   //  (1/ms) cAMP hydrolysis rate by phosphorylated PDE
		// Km_PDE_cAMP     = 1.3;            //  (uM) PDE affinity for cAMP

		//  RABBIT -> Human
		double PDE3tot = PDE3tot_Scale * 0.75 * 0.036;   //  (uM) total phosphodiesterase
		double PDE4tot = PDE4tot_Scale * 0.75 * 0.036;   //  (uM) total phosphodiesterase
		double k_cAMP_PDE3 = 3.5e-3;               //  k_pde3        [1/ms]
		double k_cAMP_PDE3p = 2 * k_cAMP_PDE3; //  (1/ms) cAMP hydrolysis rate by phosphorylated PDE
		double Km_PDE3_cAMP = Km_PDE3_cAMP_Scale * 0.15;           //  Km_pde3       [uM]
		double k_cAMP_PDE4 = K_rate_Scale * 5.0e-3;             //  k_pde4        [1/ms]
		double k_cAMP_PDE4p = K_rate_Scale * 2 * 5.0e-3/*k_cAMP_PDE4*/; //  (1/ms) cAMP hydrolysis rate by phosphorylated PDE
		double Km_PDE4_cAMP = Km_PDE4_cAMP_Scale * 1.3;            //  Km_pde4       [uM]

		double Kd_PDE_IBMX     = 30.0;           //  (uM) Kd_R2cAMP_C for IBMX binding to PDE
		double k_PKA_PDE       = K_rate_Scale * 7.5e-3;       //  (1/ms) rate constant for PDE phosphorylation by type 1 PKA
		double k_PP_PDE        = K_rate_Scale * 1.5e-3;       //  (1/ms) rate constant for PDE dephosphorylation by phosphatases

		double cAMP = cAMPtot - (RCcAMP_I + 2 * RCcAMPcAMP_I + 2 * RcAMPcAMP_I) - (RCcAMP_II + 2 * RCcAMPcAMP_II + 2 * RcAMPcAMP_II);
		double AC = ACtot - AC_GsaGTP;
		double GsaGTP = GsaGTPtot - AC_GsaGTP;
		double dAC_GsaGTP = kf_AC_Gsa * GsaGTP * AC - kr_AC_Gsa * AC_GsaGTP;

		double AC_FSK = FSK * AC / Kd_AC_FSK;
		double AC_ACT_BASAL = k_AC_basal * AC * ATP / (Km_AC_basal + ATP);
		double AC_ACT_GSA = k_AC_Gsa * AC_GsaGTP * ATP / (Km_AC_Gsa + ATP);
		double AC_ACT_FSK = k_AC_FSK * AC_FSK * ATP / (Km_AC_FSK + ATP);

		//  MOUSE
		//  // PDE_IBMX = PDEtot*IBMX/Kd_PDE_IBMX;
		//  PDE_IBMX = PDEtot*IBMX/(Kd_PDE_IBMX+IBMX);
		//  PDE = PDEtot - PDE_IBMX - PDEp;
		//  dPDEp = k_PKA_PDE*PKACII*PDE - k_PP_PDE*PDEp;
		//  PDE_ACT = k_cAMP_PDE*PDE*cAMP/(Km_PDE_cAMP+cAMP) + k_cAMP_PDEp*PDEp*cAMP/(Km_PDE_cAMP+cAMP);
		//
		//  dcAMPtot = AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE_ACT; //  ydot(11)

		//  RABBIT                                //  Add constrain on total IBMX?
		// PDE3_IBMX = PDE3tot*IBMX/Kd_PDE_IBMX;
		double PDE3_IBMX = PDE3tot * IBMX / (Kd_PDE_IBMX + IBMX);
		double PDE3 = PDE3tot - PDE3_IBMX - PDE3p;
		double dPDE3p = k_PKA_PDE * PKACII * PDE3 - k_PP_PDE * PDE3p; //  ydot(10)
		double PDE3_ACT = k_cAMP_PDE3 * PDE3 * cAMP / (Km_PDE3_cAMP + cAMP) + k_cAMP_PDE3p * PDE3p * cAMP / (Km_PDE3_cAMP + cAMP);

		// PDE4_IBMX = PDE4tot*IBMX/Kd_PDE_IBMX;
		double PDE4_IBMX = PDE4tot * IBMX / (Kd_PDE_IBMX + IBMX);
		double PDE4 = PDE4tot - PDE4_IBMX - PDE4p;
		double dPDE4p = k_PKA_PDE * PKACII * PDE4 - k_PP_PDE * PDE4p; //  ydot() - NEW STATE VARIABLE NEEDED
		double PDE4_ACT = k_cAMP_PDE4 * PDE4 * cAMP / (Km_PDE4_cAMP + cAMP) + k_cAMP_PDE4p * PDE4p * cAMP / (Km_PDE4_cAMP + cAMP);

		double dcAMPtot = AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE3_ACT - PDE4_ACT; //  ydot(11)
		//  end cAMP module
		// //  PKA module

		double PKItot          = PKItot_Scale * 0.18;         //  (uM) total PKI
		double kf_RC_cAMP      = K_rate_Scale * 1;            //  (1/[uM ms]) Kd for PKA RC binding to cAMP
		double kf_RCcAMP_cAMP  = K_rate_Scale * 1;            //  (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
		double kf_RcAMPcAMP_C  = K_rate_Scale * 4.375;        //  (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
		double kf_PKA_PKI      = K_rate_Scale * 1;            //  (1/[uM ms]) Ki for PKA inhibition by PKI
		double kr_RC_cAMP      = K_rate_Scale * 1.64;         //  (1/ms) Kd for PKA RC binding to cAMP
		double kr_RCcAMP_cAMP  = K_rate_Scale * 9.14;         //  (1/ms) Kd for PKA RC:cAMP binding to cAMP
		double kr_RcAMPcAMP_C  = K_rate_Scale * 1;            //  (1/ms) Kd for PKA R:cAMPcAMP binding to C
		double kr_PKA_PKI      = K_rate_Scale * 2e-4;         //  (1/ms) Ki for PKA inhibition by PKI

		double epsilon         = 10;             //  (-) AKAP-mediated scaling factor

		// PKAIItot = 0.059;          //  (uM) total type 2 PKA //  MOUSE
		double PKAIItot = PKAIItot_Scale * 0.084; //  RABBIT

		double PKI = PKItot - PKACI_PKI - PKACII_PKI;

		double dRC_I = - kf_RC_cAMP * RC_I * cAMP + kr_RC_cAMP * RCcAMP_I;
		double dRCcAMP_I = - kr_RC_cAMP * RCcAMP_I + kf_RC_cAMP * RC_I * cAMP - kf_RCcAMP_cAMP * RCcAMP_I * cAMP + kr_RCcAMP_cAMP * RCcAMPcAMP_I;
		double dRCcAMPcAMP_I = - kr_RCcAMP_cAMP * RCcAMPcAMP_I + kf_RCcAMP_cAMP * RCcAMP_I * cAMP - kf_RcAMPcAMP_C * RCcAMPcAMP_I + kr_RcAMPcAMP_C * RcAMPcAMP_I * PKACI;
		double dRcAMPcAMP_I = - kr_RcAMPcAMP_C * RcAMPcAMP_I * PKACI + kf_RcAMPcAMP_C * RCcAMPcAMP_I;
		double dPKACI = - kr_RcAMPcAMP_C * RcAMPcAMP_I * PKACI + kf_RcAMPcAMP_C * RCcAMPcAMP_I - kf_PKA_PKI * PKACI * PKI + kr_PKA_PKI * PKACI_PKI; //  ydot(17)
		double dPKACI_PKI = - kr_PKA_PKI * PKACI_PKI + kf_PKA_PKI * PKACI * PKI;

		double dRC_II = - kf_RC_cAMP * RC_II * cAMP + kr_RC_cAMP * RCcAMP_II;
		double dRCcAMP_II = - kr_RC_cAMP * RCcAMP_II + kf_RC_cAMP * RC_II * cAMP - kf_RCcAMP_cAMP * RCcAMP_II * cAMP + kr_RCcAMP_cAMP * RCcAMPcAMP_II;
		double dRCcAMPcAMP_II = - kr_RCcAMP_cAMP * RCcAMPcAMP_II + kf_RCcAMP_cAMP * RCcAMP_II * cAMP - kf_RcAMPcAMP_C * RCcAMPcAMP_II + kr_RcAMPcAMP_C * RcAMPcAMP_II * PKACII;
		double dRcAMPcAMP_II = - kr_RcAMPcAMP_C * RcAMPcAMP_II * PKACII + kf_RcAMPcAMP_C * RCcAMPcAMP_II;
		double dPKACII = - kr_RcAMPcAMP_C * RcAMPcAMP_II * PKACII + kf_RcAMPcAMP_C * RCcAMPcAMP_II - kf_PKA_PKI * PKACII * PKI + kr_PKA_PKI * PKACII_PKI; //  ydot(18)
		double dPKACII_PKI = - kr_PKA_PKI * PKACII_PKI + kf_PKA_PKI * PKACII * PKI;
		//  end PKA module
		// //  I-1/PP1 module

		// double PP1tot = pin(2); //  PP1tot = 0.89; //  (uM) total phosphatase 1
		double I1tot           = I1tot_Scale * 0.3;          //  (uM) total inhibitor 1
		double k_PKA_I1        = K_rate_Scale * 60e-3;        //  (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
		double Km_PKA_I1       = Km_PKA_I1_Scale * 1.0;          //  (uM) Km for I-1 phosphorylation by type 1 PKA
		double Vmax_PP2A_I1    = PP2A_I1_Vmax_Scale * 14.0e-3;      //  (uM/ms) Vmax for I-1 dephosphorylation by PP2A
		double Km_PP2A_I1      = Km_PP2A_I1_Scale * 1.0;          //  (uM) Km for I-1 dephosphorylation by PP2A

		double Ki_PP1_I1       = 1.0e-3;         //  (uM) Ki for PP1 inhibition by I-1
		double kf_PP1_I1       = K_rate_Scale * 1;            //  (uM) Ki for PP1 inhibition by I-1
		double kr_PP1_I1       = K_rate_Scale * Ki_PP1_I1;    //  (uM) Ki for PP1 inhibition by I-1

		double I1 = I1tot - I1ptot;
		double PP1 = PP1tot - I1p_PP1;
		double I1p = I1ptot - I1p_PP1;
		double I1_phosph = k_PKA_I1 * PKACI * I1 / (Km_PKA_I1 + I1);
		double I1_dephosph = Vmax_PP2A_I1 * I1ptot / (Km_PP2A_I1 + I1ptot);

		double dI1p_PP1 = kf_PP1_I1 * PP1 * I1p - kr_PP1_I1 * I1p_PP1;
		double dI1ptot = I1_phosph - I1_dephosph; //  ydot
		//  end I-1/PP1 module
		// //  PLB module

		// double PLBtot = pin(3); // p(41) = PLBtot; //  PLBtot    [uM]
		double k_PKA_PLB = K_rate_Scale * 54e-3;    // p(44) = 54;     //  k_pka_plb     [1/ms]
		double Km_PKA_PLB = Km_PKA_PLB_Scale * 21;   // p(45) = 21;     //  Km_pka_plb    [uM]
		double k_PP1_PLB = K_rate_Scale * 8.5e-3;   // p(46) = 8.5;    //  k_pp1_plb     [1/ms]
		double Km_PP1_PLB = Km_PP1_PLB_Scale * 7.0;  // p(47) = 7.0;    //  Km_pp1_plb    [uM]

		double PLB = PLBtot - PLBp;
		double PLB_phosph = k_PKA_PLB * PKACI * PLB / (Km_PKA_PLB + PLB);
		double PLB_dephosph = k_PP1_PLB * PP1 * PLBp / (Km_PP1_PLB + PLBp);
		double dPLBp = PLB_phosph - PLB_dephosph; //  ydot
		//  end PLB module
		// //  PLM module (from PLB, different total concentration)

		// double PLMtot = pin(4); //  p(102) = PLMtot; //  PLMtot    [uM]
		double k_PKA_PLM = K_rate_Scale * 54e-3; //  p(103) = 54;     //  k_pka_plb     [1/ms]
		double Km_PKA_PLM = 21; //  p(104) = 21;     //  Km_pka_plb    [uM]
		double k_PP1_PLM = K_rate_Scale * 8.5e-3; //  p(105) = 8.5;    //  k_pp1_plb     [1/ms]
		double Km_PP1_PLM = 7.0; //  p(106) = 7.0;    //  Km_pp1_plb    [uM]

		double PLM = PLMtot - PLMp;
		double PLM_phosph = k_PKA_PLM * PKACI * PLM / (Km_PKA_PLM + PLM);
		double PLM_dephosph = k_PP1_PLM * PP1 * PLMp / (Km_PP1_PLM + PLMp);
		double dPLMp = PLM_phosph - PLM_dephosph; //  ydot
		//  end PLM module
		// //  LCC module

		// double LCCtot = pin(5); // p(53) = LCCtot; //  LCCtot        [uM]
		double PKACII_LCCtot = PKACII_LCCtot_Scale * 0.025; // p(54) = 0.025;  //  PKAIIlcctot   [uM]
		double PP1_LCC = PP1_LCC_Scale * 0.025; // p(55) = 0.025;  //  PP1lcctot     [uM]
		double PP2A_LCC = PP2A_LCC_Scale * 0.025; // p(56) = 0.025;  //  PP2Alcctot    [uM]
		double k_PKA_LCC = K_rate_Scale * 54e-3;    // p(57) = 54;     //  k_pka_lcc     [1/ms]
		double Km_PKA_LCC = Km_PKA_LCC_Scale * 21;   // p(58) = 21;// *1.6;     //  Km_pka_lcc    [uM]
		double k_PP1_LCC = K_rate_Scale * 8.52e-3;  // p(59) = 8.52;   //  k_pp1_lcc     [1/ms] RABBIT, MOUSE
		// p(59) = 8.5;   //  k_pp1_lcc     [1/sec] RAT
		double Km_PP1_LCC = Km_PP1_LCC_Scale * 3;    // p(60) = 3;      //  Km_pp1_lcc    [uM]
		double k_PP2A_LCC = K_rate_Scale * 10.1e-3;  // p(61) = 10.1;   //  k_pp2a_lcc    [1/ms]
		double Km_PP2A_LCC = Km_PP2A_LCC_Scale * 3;    // p(62) = 3;      //  Km_pp2a_lcc   [uM]

		double PKACII_LCC = (PKACII_LCCtot / PKAIItot) * PKACII;
		double LCCa = LCCtot - LCCap;
		double LCCa_phosph = epsilon * k_PKA_LCC * PKACII_LCC * LCCa / (Km_PKA_LCC + epsilon * LCCa);
		double LCCa_dephosph = epsilon * k_PP2A_LCC * PP2A_LCC * LCCap / (Km_PP2A_LCC + epsilon * LCCap);
		double dLCCap = LCCa_phosph - LCCa_dephosph; //  ydot
		double LCCb = LCCtot - LCCbp;
		double LCCb_phosph = epsilon * k_PKA_LCC * PKACII_LCC * LCCb / (Km_PKA_LCC + epsilon * LCCb);
		double LCCb_dephosph = epsilon * k_PP1_LCC * PP1_LCC * LCCbp / (Km_PP1_LCC + epsilon * LCCbp);
		double dLCCbp = LCCb_phosph - LCCb_dephosph; //  ydot
		//  end LCC module
		// //  RyR module

		// double RyRtot = pin(6); // p(63) = RyRtot; //  RyRtot        [uM]
		double PKAIIryrtot = PKAIIryrtot_Scale * 0.034; // p(64) = 0.034;  //  PKAIIryrtot   [uM]
		double PP1ryr = PP1ryr_Scale * 0.034; // p(65) = 0.034;  //  PP1ryr        [uM]
		double PP2Aryr = PP2Aryr_Scale * 0.034; // p(66) = 0.034;  //  PP2Aryr       [uM]
		double kcat_pka_ryr = K_rate_Scale *  54e-3;   // p(67) = 54;     //  kcat_pka_ryr  [1/ms]
		double Km_pka_ryr = Km_pka_ryr_Scale * 21;   // p(68) = 21;     //  Km_pka_ryr    [uM]
		double kcat_pp1_ryr = K_rate_Scale * 8.52e-3;  // p(69) = 8.52;   //  kcat_pp1_ryr  [1/ms]
		double Km_pp1_ryr = Km_pp1_ryr_Scale * 7;    // p(70) = 7;      //  Km_pp1_ryr    [uM]
		double kcat_pp2a_ryr = K_rate_Scale * 10.1e-3;  // p(71) = 10.1;   //  kcat_pp2a_ryr [1/ms]
		double Km_pp2a_ryr = Km_pp2a_ryr_Scale * 4.1;  // p(72) = 4.1;    //  Km_pp2a_ryr   [uM]

		double PKACryr = (PKAIIryrtot / PKAIItot) * PKACII;
		double RyR = RyRtot - RyRp;
		double RyRPHOSPH = epsilon * kcat_pka_ryr * PKACryr * RyR / (Km_pka_ryr + epsilon * RyR);
		double RyRDEPHOSPH1 = epsilon * kcat_pp1_ryr * PP1ryr * RyRp / (Km_pp1_ryr + epsilon * RyRp);
		double RyRDEPHOSPH2A = epsilon * kcat_pp2a_ryr * PP2Aryr * RyRp / (Km_pp2a_ryr + epsilon * RyRp);
		double dRyRp = RyRPHOSPH - RyRDEPHOSPH1 - RyRDEPHOSPH2A; //  ydot
		//  end RyR module
		// //  TnI module

		// double TnItot = pin(7); // p(73) = TnItot; //  TnItot        [uM]
		double PP2A_TnI = PP2A_TnI_Scale * 0.67;   //  PP2Atni       [uM]
		double k_PKA_TnI = K_rate_Scale * 54e-3;    //  kcat_pka_tni  [1/ms]
		double Km_PKA_TnI = 21;     //  Km_pka_tni    [uM]
		double k_PP2A_TnI = K_rate_Scale * 10.1e-3;  //  kcat_pp2a_tni [1/ms]
		double Km_PP2A_TnI = 4.1;    //  Km_pp2a_tni   [uM]

		double TnIn = TnItot - TnIp;
		double TnI_phosph = k_PKA_TnI * PKACI * TnIn / (Km_PKA_TnI + TnIn);
		double TnI_dephosph = k_PP2A_TnI * PP2A_TnI * TnIp / (Km_PP2A_TnI + TnIp);
		double dTnIp = TnI_phosph - TnI_dephosph; //  ydot(26)
		//  end TnI module
		// //  Myofilament module (from TnI)



		// not used in atria.
		// double Myotot = pin(8);        //  Myotot_bar    [uM]
		double PP2A_myo = 0.67;             //  PP2Amyo       [uM]
		double kcat_pka_myo = K_rate_Scale * 54e-3;         //  kcat_pka_myo  [1/ms]
		double Km_pka_myo = 21;            //  Km_pka_myo    [uM]
		double kcat_pp2a_myo = K_rate_Scale * 10.1e-3;      //  kcat_pp2a_myo [1/ms]
		double Km_pp2a_myo = 4.1;          //  Km_pp2a_myo   [uM]

		double Myon = Myotot - Myop; //  Non-phos = tot - phos
		double MyoPHOSPH = kcat_pka_myo * PKACI * Myon / (Km_pka_myo + Myon);
		double MyoDEPHOSPH = kcat_pp2a_myo * PP2A_myo * Myop / (Km_pp2a_myo + Myop);
		double dMyop = MyoPHOSPH - MyoDEPHOSPH; //  ydot
		//  end Myo module

		////  IKs module

		// double IKstot = pin(9);        //  IKstot_bar    [uM]
		double Yotiao_tot = 0.025;         //  Yotiao_tot    [uM]
		double K_yotiao = K_rate_Scale * 0.1e-3;         //  K_yotiao      [uM] ** apply G589D mutation here: 1e2 **
		double PKAII_ikstot = PKAII_ikstot_Scale * 0.025;     //  PKAII_ikstot  [uM]
		double PP1_ikstot = PP1_ikstot_Scale * 0.025;       //  PP1_ikstot    [uM]
		double k_pka_iks = K_rate_Scale * 1.87e-3; // 54;      //  k_pka_iks     [1/ms] //  adjusted as in Xie et al 2013
		double Km_pka_iks = Km_pka_iks_Scale * 21;          //  Km_pka_iks    [uM]
		double k_pp1_iks = K_rate_Scale * 0.19e-3; // 8.52;    //  k_pp1_iks     [1/ms] //  adjusted as in Xie et al 2013
		double Km_pp1_iks = Km_pp1_iks_Scale * 7;           //  Km_pp1_iks    [uM]

		//  Effect of G589D mutation (IKs-Yotiao)
		double y1 = (-(K_yotiao + IKstot - Yotiao_tot) + sqrt((K_yotiao + IKstot - Yotiao_tot) * (K_yotiao + IKstot - Yotiao_tot) + 4 * K_yotiao * Yotiao_tot)) / 2;
		double x1 = IKstot / (1 + y1 / K_yotiao);
		double y2 = (-(K_yotiao + IKstot - Yotiao_tot) - sqrt((K_yotiao + IKstot - Yotiao_tot) * (K_yotiao + IKstot - Yotiao_tot) + 4 * K_yotiao * Yotiao_tot)) / 2;
		double x2 = IKstot / (1 + y2 / K_yotiao);

		double free_IKs;// = x1 * (y1 > 0) + x2 * (y1 <= 0);
		double free_Yotiao;// = y1 * (y1 > 0) + y2 * (y1 <= 0);
		free_IKs = y1 > 0 ? x1 : x2;
		free_Yotiao = y1 > 0 ? y1 : y2;


		double IksYot = free_IKs * free_Yotiao / K_yotiao; //  [uM] //  IKs-Yotiao

		double PKACiks = IksYot / IKstot * (PKAII_ikstot / PKAIItot) * PKACII;
		double PP1iks = IksYot / IKstot * PP1_ikstot;

		double KSn = IKstot - KSp; //  Non-phos = tot - phos
		double IKS_PHOSPH = epsilon * k_pka_iks * PKACiks * KSn / (Km_pka_iks + epsilon * KSn);
		double IKS_DEPHOSPH = epsilon * k_pp1_iks * PP1iks * KSp / (Km_pp1_iks + epsilon * KSp);
		double dKSp = IKS_PHOSPH - IKS_DEPHOSPH; //  ydot
		//  end Iks module
		// //  IKr module (from IKs, without mutation)

		// double Krtot = pin(10);        //  IKrtot_bar    [uM]
		double PKAII_ikrtot = 0.025;       //  PKAII_ikrtot  [uM]
		double PP1_ikrtot = PP1_ikrtot_Scale * 0.025;         //  PP1_ikrtot    [uM]
		double k_pka_ikr = K_rate_Scale * 1.87e-3; // 54;      //  k_pka_ikr     [1/ms] //  adjusted as in Xie et al 2013
		double Km_pka_ikr = 21;            //  Km_pka_ikr    [uM]
		double k_pp1_ikr = K_rate_Scale * 0.19e-3; // 8.52;    //  k_pp1_ikr     [1/ms] //  adjusted as in Xie et al 2013
		double Km_pp1_ikr = 7;             //  Km_pp1_ikr    [uM]

		double KRn = IKrtot - KRp; //  Non-phos = tot - phos
		double PKACikr = (PKAII_ikrtot / PKAIItot) * PKACII;
		double IKR_PHOSPH = epsilon * k_pka_ikr * PKACikr * KRn / (Km_pka_ikr + epsilon * KRn);
		double IKR_DEPHOSPH = epsilon * k_pp1_ikr * PP1_ikrtot * KRp / (Km_pp1_ikr + epsilon * KRp);
		double dKRp = IKR_PHOSPH - IKR_DEPHOSPH;
		//  end IKr module


		//// Ikur module (from Ikr)
		//IKurtot = pin(9);// p(95) = IKurtot;   // Ikur_tot      [uM]
		double PKAII_KURtot = 0.025; // p (96) = 0.025;  // PKAII_KURtot [uM]
		double PP1_KURtot = PP1_KURtot_Scale * 0.025; // p(97) = 0.025;  // PP1_KURtot   [uM]
		double k_pka_KUR = K_rate_Scale * 1.87e-3; // 54e-3; // k_pka_ikur [1/ms] // adjusted as in Xie et al 2013
		double Km_pka_KUR = 21; // p(99) = 21;    // Km_pka_KUR   [uM]
		double k_pp1_KUR = K_rate_Scale * 0.19e-3; // 8.52e-3; // k_pp1_ikur [1/ms] // adjusted as in Xie et al 2013
		double Km_pp1_KUR = 7; // p(101) = 7;     // Km_pp1_KUR   [uM]

		double KURn = IKurtot - KURp;  // Non-phos = tot - phos
		double PKACII_KUR = (PKAII_KURtot / PKAIItot) * PKACII; // (PKA_KURtot/PKAIItot)*PKAIIact
		double KURphos = epsilon * PKACII_KUR * k_pka_KUR * KURn / (Km_pka_KUR + epsilon * KURn);
		double KURdephos = epsilon * PP1_KURtot * k_pp1_KUR * KURp / (Km_pp1_KUR + epsilon * KURp);
		double dKURp = KURphos - KURdephos;
		// end Ikur module


		// //  ICFTR module


		// not used in atria.
		// double CFTRtot = pin(11); //  ICFTRtot_bar  [uM]
		double PKAII_CFTRtot = 0.025;      //  PKAII_CFTRtot [uM]
		double PP1_CFTRtot = PP1_CFTRtot_Scale * 0.025;        //  PP1_CFTRtot   [uM]
		double k_pka_CFTR = K_rate_Scale * 54e-3;           //  k_pka_CFTR    [1/ms]
		double Km_pka_CFTR = 8.5;          //  Km_pka_CFTR   [uM]
		double k_pp1_CFTR = K_rate_Scale * 8.52e-3;         //  k_pp1_CFTR    [1/ms]
		double Km_pp1_CFTR = 7;            //  Km_pp1_CFTR   [uM]


		// not used in atria
		double CFTRn = ICFTRtot - CFTRp;  //  Non-phos = tot - phos
		double PKAC_CFTR = (PKAII_CFTRtot / PKAIItot) * PKACII; //  (PKACFTRtot/PKAIItot)*PKAIIact
		double CFTRphos = epsilon * CFTRn * PKAC_CFTR * k_pka_CFTR / (Km_pka_CFTR + epsilon * CFTRn);
		double CFTRdephos = PP1_CFTRtot * k_pp1_CFTR * epsilon * CFTRp / (Km_pp1_CFTR + epsilon * CFTRp);
		double dCFTRp = CFTRphos - CFTRdephos;
		//  end ICFTR module
		// //  IClCa module (from CFTR)

		// double ClCatot = pin(12);	//  IClCatot_bar  [uM]
		double PKAII_ClCatot = 0.025;      //  PKAII_ClCatot [uM]
		double PP1_ClCatot = PP1_ClCatot_Scale * 0.025;        //  PP1_ClCatot   [uM]
		double k_pka_ClCa = K_rate_Scale *  54e-3;          //  k_pka_ClCa    [1/ms]
		double Km_pka_ClCa = 8.5;          //  Km_pka_ClCa   [uM]
		double k_pp1_ClCa = K_rate_Scale * 8.52e-3;         //  k_pp1_ClCa    [1/ms]
		double Km_pp1_ClCa = 7;            //  Km_pp1_ClCa   [uM]

		double ClCan = IClCatot - ClCap;  //  Non-phos = tot - phos
		double PKAC_ClCa = (PKAII_ClCatot / PKAIItot) * PKACII; //  (PKACFTRtot/PKAIItot)*PKAIIact
		double ClCaphos = epsilon * ClCan * PKAC_ClCa * k_pka_ClCa / (Km_pka_ClCa + epsilon * ClCan);
		double ClCadephos = PP1_ClCatot * k_pp1_ClCa * epsilon * ClCap / (Km_pp1_ClCa + epsilon * ClCap);
		double dClCap = ClCaphos - ClCadephos;
		//  end ICl(Ca) module
		// //  ydot


		//// Ina module (from Ikr)

		//INatot = pin(6); // p(79) [uM]
		double PKAIIinatot = 0.025; // PKAII_inatot [uM]
		double PP1_inatot = PP1_inatot_Scale * 0.025; // PP1_inatot [uM]
		double k_pka_ina = K_rate_Scale * 1.87e-3; // 54e-3; // k_pka_ina [1/ms] // adjusted as in Xie et al 2013
		double Km_pka_ina = 21; // Km_pka_ina [uM]
		double k_pp1_ina = K_rate_Scale * 0.19e-3; // 8.52e-3; // k_pp1_ina [1/ms] // adjusted as in Xie et al 2013
		double Km_pp1_ina = 7; // Km_pp1_ina [uM]

		double NAn = INatot - NAp;
		double PKACII_ina = (PKAIIinatot / PKAIItot) * PKACII;
		double INA_PHOSPH = epsilon * k_pka_ina * PKACII_ina * NAn / (Km_pka_ina + epsilon * NAn);
		double INA_DEPHOSPH = epsilon * k_pp1_ina * PP1_inatot * NAp / (Km_pp1_ina + epsilon * NAp);
		double dNAp = INA_PHOSPH - INA_DEPHOSPH;
		// end Ina module
		//// Ito module (from Ikr)

		//Itotot = pin(6); // p(79) [uM]
		double PKAIIitotot = 0.025; // PKAII_ikrtot [uM]
		double PP1_itotot = PP1_itotot_Scale * 0.025; // PP1_ikrtot [uM]
		double k_pka_ito = K_rate_Scale * 1.87e-3; // 54e-3; // k_pka_ikr [1/ms] // adjusted as in Xie et al 2013
		double Km_pka_ito = 21; // Km_pka_ikr [uM]
		double k_pp1_ito = K_rate_Scale * 0.19e-3; // 8.52e-3; // k_pp1_ikr [1/ms] // adjusted as in Xie et al 2013
		double Km_pp1_ito = 7; // Km_pp1_ikr [uM]

		double TOn = Itotot - TOp;
		double PKACII_ito = (PKAIIitotot / PKAIItot) * PKACII;
		double Ito_PHOSPH = epsilon * k_pka_ito * PKACII_ito * TOn / (Km_pka_ito + epsilon * TOn);
		double Ito_DEPHOSPH = epsilon * k_pp1_ito * PP1_itotot * TOp / (Km_pp1_ito + epsilon * TOp);
		double dTOp = Ito_PHOSPH - Ito_DEPHOSPH;
		// end Ito module
		//// Ik1 module (from Ikr)

		//IK1tot = pin(6); // p(79) [uM]
		double PKAIIik1tot = 0.025; // PKAII_ikrtot [uM]
		double PP1_ik1tot = PP1_ik1tot_Scale * 0.025; // PP1_ikrtot [uM]
		double k_pka_ik1 = K_rate_Scale * 1.87e-3; // 54e-3; // k_pka_ikr [1/ms] // adjusted as in Xie et al 2013
		double Km_pka_ik1 = 21; // Km_pka_ikr [uM]
		double k_pp1_ik1 = K_rate_Scale * 0.19e-3; // 8.52e-3; // k_pp1_ikr [1/ms] // adjusted as in Xie et al 2013
		double Km_pp1_ik1 = 7; // Km_pp1_ikr [uM]

		double K1n = IK1tot - K1p;
		double PKACII_ik1 = (PKAIIik1tot / PKAIItot) * PKACII;
		double IK1_PHOSPH = epsilon * k_pka_ik1 * PKACII_ik1 * K1n / (Km_pka_ik1 + epsilon * K1n);
		double IK1_DEPHOSPH = epsilon * k_pp1_ik1 * PP1_ik1tot * K1p / (Km_pp1_ik1 + epsilon * K1p);
		double dK1p = IK1_PHOSPH - IK1_DEPHOSPH;


		// ydot(10)=dPDEp;
		// ydot[0] = dLR;
		// ydot[1 ] = dLRG;
		// ydot[2 ] = dRG;
		// ydot[3 ] = db1AR_S464;
		// ydot[4 ] = db1AR_S301;
		// ydot[5 ] = dGsaGTPtot;
		// ydot[6 ] = dGsaGDP;
		// ydot[7 ] = dGsby;
		// ydot[8 ] = dAC_GsaGTP;
		// ydot[9 ] = dPDE3p;
		// ydot[10] = dPDE4p;
		// ydot[11] = dcAMPtot;
		// ydot[12] = dRC_I;
		// ydot[13] = dRCcAMP_I;
		// ydot[14] = dRCcAMPcAMP_I;
		// ydot[15] = dRcAMPcAMP_I;
		// ydot[16] = dPKACI;
		// ydot[17] = dPKACI_PKI;
		// ydot[18] = dRC_II;
		// ydot[19] = dRCcAMP_II;
		// ydot[20] = dRCcAMPcAMP_II;
		// ydot[21] = dRcAMPcAMP_II;
		// ydot[22] = dPKACII;
		// ydot[23] = dPKACII_PKI;
		// ydot[24] = dI1p_PP1; //  output CaMKII
		// ydot[25] = dI1ptot;
		// ydot[26] = dPLBp; //  output
		// ydot[27] = dPLMp; //  output
		// ydot[28] = dLCCap; //  output
		// ydot[29] = dLCCbp; //  output
		// ydot[30] = dRyRp; //  output
		// ydot[31] = dTnIp; //  output
		// ydot[32] = dMyop; //  output
		// ydot[33] = dKSp; //  output
		// ydot[34] = dKRp; //  output
		// ydot[35] = dCFTRp; //  output
		// ydot[36] = dClCap; //  output
		// ydot[37] = dNAp; // output
		// ydot[38] = dTOp; // output
		// ydot[39] = dK1p; // output
		// ydot[40] = dKURp; // output

		ydot[0]  = dLR;
		ydot[1]  = dLRG;
		ydot[2]  = dRG;
		ydot[3]  = db1AR_S464;
		ydot[4]  = db1AR_S301;
		ydot[5]  = dGsaGTPtot;
		ydot[6]  = dGsaGDP;
		ydot[7]  = dGsby;
		ydot[8]  = dAC_GsaGTP;
		ydot[9]  = dPDE4p;
		ydot[10] = dcAMPtot;
		ydot[11] = dRC_I;
		ydot[12] = dRCcAMP_I;
		ydot[13] = dRCcAMPcAMP_I;
		ydot[14] = dRcAMPcAMP_I;
		ydot[15] = dPKACI;
		ydot[16] = dPKACI_PKI;
		ydot[17] = dRC_II;
		ydot[18] = dRCcAMP_II;
		ydot[19] = dRCcAMPcAMP_II;
		ydot[20] = dRcAMPcAMP_II;
		ydot[21] = dPKACII;
		ydot[22] = dPKACII_PKI;
		ydot[23] = dI1p_PP1; // output CaMKII
		ydot[24] = dI1ptot;
		ydot[25] = dPLBp; // output
		ydot[26] = dPLMp; // output
		ydot[27] = dLCCap; // output
		ydot[28] = dLCCbp; // output
		ydot[29] = dRyRp; // output
		ydot[30] = dTnIp; // output
		ydot[31] = dKSp; // output
		ydot[32] = dKRp; // output
		ydot[33] = dKURp; // output
		ydot[34] = dCFTRp; // output
		ydot[35] = dClCap; // output
		ydot[36] = dMyop; // output
		ydot[37] = dNAp; // output
		ydot[38] = dTOp; // output
		ydot[39] = dK1p; // output
		ydot[40]  = dPDE3p;  // add PDE3p for human model





		// ydot[0]  = dLR;
		// ydot[1]  = dLRG;
		// ydot[2]  = dRG;
		// ydot[3]  = db1AR_S464;
		// ydot[4]  = db1AR_S301;
		// ydot[5]  = dGsaGTPtot;
		// ydot[6]  = dGsaGDP;
		// ydot[7]  = dGsby;
		// ydot[8]  = dAC_GsaGTP;
		// ydot[9]  = dPDE3p;
		// ydot[10] = dcAMPtot;
		// ydot[11] = dRC_I;
		// ydot[12] = dRCcAMP_I;
		// ydot[13] = dRCcAMPcAMP_I;
		// ydot[14] = dRcAMPcAMP_I;
		// ydot[15] = dPKACI;
		// ydot[16] = dPKACI_PKI;
		// ydot[17] = dRC_II;
		// ydot[18] = dRCcAMP_II;
		// ydot[19] = dRCcAMPcAMP_II;
		// ydot[20] = dRcAMPcAMP_II;
		// ydot[21] = dPKACII;
		// ydot[22] = dPKACII_PKI;
		// ydot[23] = dI1p_PP1; // output CaMKII
		// ydot[24] = dI1ptot;
		// ydot[25] = dPLBp; // output
		// ydot[26] = dPLMp; // output
		// ydot[27] = dLCCap; // output
		// ydot[28] = dLCCbp; // output
		// ydot[29] = dRyRp; // output
		// ydot[30] = dTnIp; // output
		// ydot[31] = dKSp; // output
		// ydot[32] = dKRp; // output
		// ydot[33] = dKURp; // output
		// ydot[34] = dCFTRp; // output
		// ydot[35] = dClCap; // output
		// ydot[36] = dmyop; // output
		// ydot[37] = dNAp; // output
		// ydot[38] = dTOp; // output
		// ydot[39] = dK1p; // output
		// ydot[40] = dPDE4p;

		return 0;

	}

	int update_single_time_step( double t, betaAR_para & para) {

		//  This module describes the beta-adrenergic signaling pathway.

		// //  State variables

		ISO      = para.ISO      ;
		LCCtot   = LCCtot_Scale * para.LCCtot   ; // pin(2);
		RyRtot   = RyRtot_Scale * para.RyRtot   ; // pin(3);
		PLBtot   = PLBtot_Scale * para.PLBtot   ; // pin(4);
		TnItot   = TnItot_Scale * para.TnItot   ; // pin(5);
		IKstot   = para.IKstot   ;   // pin(6);
		ICFTRtot = para.ICFTRtot ;   // pin(7);
		PP1tot   = PP1tot_Scale * para.PP1tot   ; // pin(8);
		PLMtot   = PLMtot_Scale * para.PLMtot   ; // pin(9);
		Myotot   = para.myotot   ;   // pin(10);
		IKrtot   = para.IKrtot   ;   // pin(11);
		IKurtot  = para.IKurtot  ;   // pin(12);
		INatot   = para.INatot   ;   // pin(13);
		IClCatot = para.IClCatot ;   // pin(14);
		Itotot   = para.Itotot   ;   // pin(15);
		IK1tot   = para.IK1tot   ;   // pin(16);
		AF_index = para.AF_index ;   // pin(17); // AF_index

		// double LR            = y[0];
		// double LRG           = y[1];
		// double RG            = y[2];
		// double b1AR_S464     = y[3];
		// double b1AR_S301     = y[4];
		// double GsaGTPtot     = y[5];
		// double GsaGDP        = y[6];
		// double Gsby          = y[7];
		// double AC_GsaGTP     = y[8];
		// // PDEp              =y(10);9
		// double PDE3p         = y[9];
		// double PDE4p         = y[10];
		// double cAMPtot       = y[11];
		// double RC_I          = y[12];
		// double RCcAMP_I      = y[13];
		// double RCcAMPcAMP_I  = y[14];
		// double RcAMPcAMP_I   = y[15];
		// double PKACI         = y[16];
		// double PKACI_PKI     = y[17];
		// double RC_II         = y[18];
		// double RCcAMP_II     = y[19];
		// double RCcAMPcAMP_II = y[20];
		// double RcAMPcAMP_II  = y[21];
		// double PKACII        = y[22];
		// double PKACII_PKI    = y[23];
		// double I1p_PP1       = y[24]; //  output CaMKII
		// double I1ptot        = y[25];
		// double PLBp          = y[26]; //  output
		// double PLMp          = y[27]; //  output
		// double LCCap         = y[28]; //  output
		// double LCCbp         = y[29]; //  output
		// double RyRp          = y[30]; //  output
		// double TnIp          = y[31]; //  output
		// double Myop          = y[32]; //  output
		// double KSp           = y[33]; //  output
		// double KRp           = y[34]; //  output
		// double CFTRp         = y[35]; //  output
		// double ClCap         = y[36]; //  output

		// double NAp           = y[37]; // output
		// double TOp           = y[38]; // output
		// double K1p           = y[39]; // output
		// double KURp          = y[40]; // output
		double LR            = y[0]; // bound receptor
		double LRG           = y[1]; // G protein-associated, ligand bound receptor
		double RG            = y[2]; //
		double b1AR_S464     = y[3];
		double b1AR_S301     = y[4];
		double GsaGTPtot     = y[5];
		double GsaGDP        = y[6];
		double Gsby          = y[7];
		double AC_GsaGTP     = y[8];
		double PDE4p          = y[9];
		double cAMPtot       = y[10];
		double RC_I          = y[11];
		double RCcAMP_I      = y[12];
		double RCcAMPcAMP_I  = y[13];
		double RcAMPcAMP_I   = y[14];
		double PKACI         = y[15];  // pka_1
		double PKACI_PKI     = y[16];
		double RC_II         = y[17];
		double RCcAMP_II     = y[18];
		double RCcAMPcAMP_II = y[19];
		double RcAMPcAMP_II  = y[20];
		double PKACII        = y[21];  // pka_2
		double PKACII_PKI    = y[22];
		double I1p_PP1       = y[23]; // output CaMKII
		double I1ptot        = y[24];
		double PLBp          = y[25]; // output
		double PLMp          = y[26]; // output
		double LCCap         = y[27]; // output
		double LCCbp         = y[28]; // output
		double RyRp          = y[29]; // output
		double TnIp          = y[30]; // output
		double KSp           = y[31]; // output
		double KRp           = y[32]; // output
		double KURp          = y[33]; // output
		double CFTRp         = y[34]; // output
		double ClCap         = y[35]; // output
		double Myop          = y[36]; // output
		double NAp           = y[37]; // output
		double TOp           = y[38]; // output
		double K1p           = y[39]; // output
		// double PDE3p         = y[9];
		double PDE3p         = y[40];  // added PDE3p for rabbit and human mdoels // 16:59:21, Wed, 18-September-2019, By Haibo
		// //  Drug Concentrations

		// ISO = pin(1); //  (uM) isoproterenol concentration - Ltot
		double FSK = 0; //  (uM) forskolin concentration
		double IBMX = 0; //  (uM) IBMX concentration
		// //  b-AR module

		// b1ARtot = 0.00528;        //  (uM) total b1-AR protein //  MOUSE
		double b1ARtot = b1ARtot_Scale * 0.028; //  RABBIT

		double kf_LR           = K_rate_Scale * 1;              //  (1/[uM ms]) forward rate for ISO binding to b1AR
		double kr_LR           = K_rate_Scale * 0.285;          //  (1/ms) reverse rate for ISO binding to b1AR
		double kf_LRG          = K_rate_Scale * 1;              //  (1/[uM ms]) forward rate for ISO:b1AR association with Gs
		double kr_LRG          = K_rate_Scale * 0.062;          //  (1/ms) reverse rate for ISO:b1AR association with Gs
		double kf_RG           = K_rate_Scale * 1;              //  (1/[uM ms]) forward rate for b1AR association with Gs
		double kr_RG           = K_rate_Scale * 33;             //  (1/ms) reverse rate for b1AR association with Gs

		double Gstot           = Gstot_Scale * 3.83;         //  (uM) total Gs protein
		double k_G_act         = K_rate_Scale * 16e-3;        //  (1/ms) rate constant for Gs activation
		double k_G_hyd         = K_rate_Scale * 0.8e-3;       //  (1/ms) rate constant for G-protein hydrolysis
		double k_G_reassoc     = K_rate_Scale * 1.21;         //  (1/[uM ms]) rate constant for G-protein reassociation

		double kf_bARK         = K_rate_Scale * 1.1e-6;       //  (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
		double kr_bARK         = K_rate_Scale * 2.2e-6;       //  (1/ms) reverse rate for b1AR phosphorylation by b1ARK
		double kf_PKA          = K_rate_Scale * 3.6e-6;       //  (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
		double kr_PKA          = K_rate_Scale * 2.2e-6;       //  (1/ms) reverse rate for b1AR phosphorylation by PKA

		double b1ARact = b1ARtot - b1AR_S464 - b1AR_S301;
		double b1AR = b1ARact - LR - LRG - RG;
		double Gs = Gstot - LRG - RG - Gsby;

		double dLR = kf_LR * ISO * b1AR - kr_LR * LR + kr_LRG * LRG - kf_LRG * LR * Gs;
		double dLRG = kf_LRG * LR * Gs - kr_LRG * LRG - k_G_act * LRG;
		double dRG = kf_RG * b1AR * Gs - kr_RG * RG - k_G_act * RG;

		double bARK_desens = kf_bARK * (LR + LRG);
		double bARK_resens = kr_bARK * b1AR_S464;
		double PKA_desens = kf_PKA * PKACI * b1ARact;
		double PKA_resens = kr_PKA * b1AR_S301;
		double db1AR_S464 = bARK_desens - bARK_resens; //  ydot(5)
		double db1AR_S301 = PKA_desens - PKA_resens; //  ydot(6)

		double G_act = k_G_act * (RG + LRG);
		double G_hyd = k_G_hyd * GsaGTPtot;
		double G_reassoc = k_G_reassoc * GsaGDP * Gsby;
		double dGsaGTPtot = G_act - G_hyd; //  ydot(7)
		double dGsaGDP = G_hyd - G_reassoc; //  ydot(8)
		double dGsby = G_act - G_reassoc; //  ydot(9)
		//  end b-AR module
		// //  cAMP module

		// ACtot = 70.57e-3;       //  (uM) total adenylyl cyclase //  MOUSE
		double ACtot = ACtot_Scale * 47e-3; //  RABBIT

		double ATP             = ATP_Scale * 5e3;          //  (uM) total ATP
		double k_AC_basal      = K_rate_Scale * 0.2e-3;       //  (1/ms) basal cAMP generation rate by AC
		double Km_AC_basal     = Km_AC_basal_Scale * 1.03e3;       //  (uM) basal AC affinity for ATP

		double Kd_AC_Gsa       = 0.4;            //  (uM) Kd for AC association with Gsa
		double kf_AC_Gsa       = K_rate_Scale * 1;            //  (1/[uM ms]) forward rate for AC association with Gsa
		double kr_AC_Gsa       = K_rate_Scale * Kd_AC_Gsa;    //  (1/ms) reverse rate for AC association with Gsa

		double k_AC_Gsa        = K_rate_Scale * 8.5e-3;       //  (1/ms) basal cAMP generation rate by AC:Gsa
		double Km_AC_Gsa       = Kd_AC_Gsa_Scale * 315.0;        //  (uM) AC:Gsa affinity for ATP

		double Kd_AC_FSK       = 44.0;           //  (uM) Kd for FSK binding to AC
		double k_AC_FSK        = K_rate_Scale * 7.3e-3;       //  (1/ms) basal cAMP generation rate by AC:FSK
		double Km_AC_FSK       = 860.0;          //  (uM) AC:FSK affinity for ATP

		//  MOUSE
		// PDEtot          = 22.85e-3;       //  (uM) total phosphodiesterase
		// k_cAMP_PDE      = 5e-3;           //  (1/ms) cAMP hydrolysis rate by PDE
		// k_cAMP_PDEp     = 2*k_cAMP_PDE;   //  (1/ms) cAMP hydrolysis rate by phosphorylated PDE
		// Km_PDE_cAMP     = 1.3;            //  (uM) PDE affinity for cAMP

		//  RABBIT -> Human
		double PDE3tot = PDE3tot_Scale * 0.75 * 0.036;   //  (uM) total phosphodiesterase
		double PDE4tot = PDE4tot_Scale * 0.75 * 0.036;   //  (uM) total phosphodiesterase
		double k_cAMP_PDE3 = 3.5e-3;               //  k_pde3        [1/ms]
		double k_cAMP_PDE3p = 2 * k_cAMP_PDE3; //  (1/ms) cAMP hydrolysis rate by phosphorylated PDE
		double Km_PDE3_cAMP = Km_PDE3_cAMP_Scale * 0.15;           //  Km_pde3       [uM]
		double k_cAMP_PDE4 = K_rate_Scale * 5.0e-3;             //  k_pde4        [1/ms]
		double k_cAMP_PDE4p = K_rate_Scale * 2 * 5.0e-3/*k_cAMP_PDE4*/; //  (1/ms) cAMP hydrolysis rate by phosphorylated PDE
		double Km_PDE4_cAMP = Km_PDE4_cAMP_Scale * 1.3;            //  Km_pde4       [uM]

		double Kd_PDE_IBMX     = 30.0;           //  (uM) Kd_R2cAMP_C for IBMX binding to PDE
		double k_PKA_PDE       = K_rate_Scale * 7.5e-3;       //  (1/ms) rate constant for PDE phosphorylation by type 1 PKA
		double k_PP_PDE        = K_rate_Scale * 1.5e-3;       //  (1/ms) rate constant for PDE dephosphorylation by phosphatases

		double cAMP = cAMPtot - (RCcAMP_I + 2 * RCcAMPcAMP_I + 2 * RcAMPcAMP_I) - (RCcAMP_II + 2 * RCcAMPcAMP_II + 2 * RcAMPcAMP_II);
		double AC = ACtot - AC_GsaGTP;
		double GsaGTP = GsaGTPtot - AC_GsaGTP;
		double dAC_GsaGTP = kf_AC_Gsa * GsaGTP * AC - kr_AC_Gsa * AC_GsaGTP;

		double AC_FSK = FSK * AC / Kd_AC_FSK;
		double AC_ACT_BASAL = k_AC_basal * AC * ATP / (Km_AC_basal + ATP);
		double AC_ACT_GSA = k_AC_Gsa * AC_GsaGTP * ATP / (Km_AC_Gsa + ATP);
		double AC_ACT_FSK = k_AC_FSK * AC_FSK * ATP / (Km_AC_FSK + ATP);

		//  MOUSE
		//  // PDE_IBMX = PDEtot*IBMX/Kd_PDE_IBMX;
		//  PDE_IBMX = PDEtot*IBMX/(Kd_PDE_IBMX+IBMX);
		//  PDE = PDEtot - PDE_IBMX - PDEp;
		//  dPDEp = k_PKA_PDE*PKACII*PDE - k_PP_PDE*PDEp;
		//  PDE_ACT = k_cAMP_PDE*PDE*cAMP/(Km_PDE_cAMP+cAMP) + k_cAMP_PDEp*PDEp*cAMP/(Km_PDE_cAMP+cAMP);
		//
		//  dcAMPtot = AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE_ACT; //  ydot(11)

		//  RABBIT                                //  Add constrain on total IBMX?
		// PDE3_IBMX = PDE3tot*IBMX/Kd_PDE_IBMX;
		double PDE3_IBMX = PDE3tot * IBMX / (Kd_PDE_IBMX + IBMX);
		double PDE3 = PDE3tot - PDE3_IBMX - PDE3p;
		double dPDE3p = k_PKA_PDE * PKACII * PDE3 - k_PP_PDE * PDE3p; //  ydot(10)
		double PDE3_ACT = k_cAMP_PDE3 * PDE3 * cAMP / (Km_PDE3_cAMP + cAMP) + k_cAMP_PDE3p * PDE3p * cAMP / (Km_PDE3_cAMP + cAMP);

		// PDE4_IBMX = PDE4tot*IBMX/Kd_PDE_IBMX;
		double PDE4_IBMX = PDE4tot * IBMX / (Kd_PDE_IBMX + IBMX);
		double PDE4 = PDE4tot - PDE4_IBMX - PDE4p;
		double dPDE4p = k_PKA_PDE * PKACII * PDE4 - k_PP_PDE * PDE4p; //  ydot() - NEW STATE VARIABLE NEEDED
		double PDE4_ACT = k_cAMP_PDE4 * PDE4 * cAMP / (Km_PDE4_cAMP + cAMP) + k_cAMP_PDE4p * PDE4p * cAMP / (Km_PDE4_cAMP + cAMP);

		double dcAMPtot = AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE3_ACT - PDE4_ACT; //  ydot(11)
		//  end cAMP module
		// //  PKA module

		double PKItot          = PKItot_Scale * 0.18;         //  (uM) total PKI
		double kf_RC_cAMP      = K_rate_Scale * 1;            //  (1/[uM ms]) Kd for PKA RC binding to cAMP
		double kf_RCcAMP_cAMP  = K_rate_Scale * 1;            //  (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
		double kf_RcAMPcAMP_C  = K_rate_Scale * 4.375;        //  (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
		double kf_PKA_PKI      = K_rate_Scale * 1;            //  (1/[uM ms]) Ki for PKA inhibition by PKI
		double kr_RC_cAMP      = K_rate_Scale * 1.64;         //  (1/ms) Kd for PKA RC binding to cAMP
		double kr_RCcAMP_cAMP  = K_rate_Scale * 9.14;         //  (1/ms) Kd for PKA RC:cAMP binding to cAMP
		double kr_RcAMPcAMP_C  = K_rate_Scale * 1;            //  (1/ms) Kd for PKA R:cAMPcAMP binding to C
		double kr_PKA_PKI      = K_rate_Scale * 2e-4;         //  (1/ms) Ki for PKA inhibition by PKI

		double epsilon         = 10;             //  (-) AKAP-mediated scaling factor

		// PKAIItot = 0.059;          //  (uM) total type 2 PKA //  MOUSE
		double PKAIItot = PKAIItot_Scale * 0.084; //  RABBIT

		double PKI = PKItot - PKACI_PKI - PKACII_PKI;

		double dRC_I = - kf_RC_cAMP * RC_I * cAMP + kr_RC_cAMP * RCcAMP_I;
		double dRCcAMP_I = - kr_RC_cAMP * RCcAMP_I + kf_RC_cAMP * RC_I * cAMP - kf_RCcAMP_cAMP * RCcAMP_I * cAMP + kr_RCcAMP_cAMP * RCcAMPcAMP_I;
		double dRCcAMPcAMP_I = - kr_RCcAMP_cAMP * RCcAMPcAMP_I + kf_RCcAMP_cAMP * RCcAMP_I * cAMP - kf_RcAMPcAMP_C * RCcAMPcAMP_I + kr_RcAMPcAMP_C * RcAMPcAMP_I * PKACI;
		double dRcAMPcAMP_I = - kr_RcAMPcAMP_C * RcAMPcAMP_I * PKACI + kf_RcAMPcAMP_C * RCcAMPcAMP_I;
		double dPKACI = - kr_RcAMPcAMP_C * RcAMPcAMP_I * PKACI + kf_RcAMPcAMP_C * RCcAMPcAMP_I - kf_PKA_PKI * PKACI * PKI + kr_PKA_PKI * PKACI_PKI; //  ydot(17)
		double dPKACI_PKI = - kr_PKA_PKI * PKACI_PKI + kf_PKA_PKI * PKACI * PKI;

		double dRC_II = - kf_RC_cAMP * RC_II * cAMP + kr_RC_cAMP * RCcAMP_II;
		double dRCcAMP_II = - kr_RC_cAMP * RCcAMP_II + kf_RC_cAMP * RC_II * cAMP - kf_RCcAMP_cAMP * RCcAMP_II * cAMP + kr_RCcAMP_cAMP * RCcAMPcAMP_II;
		double dRCcAMPcAMP_II = - kr_RCcAMP_cAMP * RCcAMPcAMP_II + kf_RCcAMP_cAMP * RCcAMP_II * cAMP - kf_RcAMPcAMP_C * RCcAMPcAMP_II + kr_RcAMPcAMP_C * RcAMPcAMP_II * PKACII;
		double dRcAMPcAMP_II = - kr_RcAMPcAMP_C * RcAMPcAMP_II * PKACII + kf_RcAMPcAMP_C * RCcAMPcAMP_II;
		double dPKACII = - kr_RcAMPcAMP_C * RcAMPcAMP_II * PKACII + kf_RcAMPcAMP_C * RCcAMPcAMP_II - kf_PKA_PKI * PKACII * PKI + kr_PKA_PKI * PKACII_PKI; //  ydot(18)
		double dPKACII_PKI = - kr_PKA_PKI * PKACII_PKI + kf_PKA_PKI * PKACII * PKI;
		//  end PKA module
		// //  I-1/PP1 module

		// double PP1tot = pin(2); //  PP1tot = 0.89; //  (uM) total phosphatase 1
		double I1tot           = I1tot_Scale * 0.3;          //  (uM) total inhibitor 1
		double k_PKA_I1        = K_rate_Scale * 60e-3;        //  (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
		double Km_PKA_I1       = Km_PKA_I1_Scale * 1.0;          //  (uM) Km for I-1 phosphorylation by type 1 PKA
		double Vmax_PP2A_I1    = PP2A_I1_Vmax_Scale * 14.0e-3;      //  (uM/ms) Vmax for I-1 dephosphorylation by PP2A
		double Km_PP2A_I1      = Km_PP2A_I1_Scale * 1.0;          //  (uM) Km for I-1 dephosphorylation by PP2A

		double Ki_PP1_I1       = 1.0e-3;         //  (uM) Ki for PP1 inhibition by I-1
		double kf_PP1_I1       = K_rate_Scale * 1;            //  (uM) Ki for PP1 inhibition by I-1
		double kr_PP1_I1       = K_rate_Scale * Ki_PP1_I1;    //  (uM) Ki for PP1 inhibition by I-1

		double I1 = I1tot - I1ptot;
		double PP1 = PP1tot - I1p_PP1;
		double I1p = I1ptot - I1p_PP1;
		double I1_phosph = k_PKA_I1 * PKACI * I1 / (Km_PKA_I1 + I1);
		double I1_dephosph = Vmax_PP2A_I1 * I1ptot / (Km_PP2A_I1 + I1ptot);

		double dI1p_PP1 = kf_PP1_I1 * PP1 * I1p - kr_PP1_I1 * I1p_PP1;
		double dI1ptot = I1_phosph - I1_dephosph; //  ydot
		//  end I-1/PP1 module
		// //  PLB module

		// double PLBtot = pin(3); // p(41) = PLBtot; //  PLBtot    [uM]
		double k_PKA_PLB = K_rate_Scale * 54e-3;    // p(44) = 54;     //  k_pka_plb     [1/ms]
		double Km_PKA_PLB = Km_PKA_PLB_Scale * 21;   // p(45) = 21;     //  Km_pka_plb    [uM]
		double k_PP1_PLB = K_rate_Scale * 8.5e-3;   // p(46) = 8.5;    //  k_pp1_plb     [1/ms]
		double Km_PP1_PLB = Km_PP1_PLB_Scale * 7.0;  // p(47) = 7.0;    //  Km_pp1_plb    [uM]

		double PLB = PLBtot - PLBp;
		double PLB_phosph = k_PKA_PLB * PKACI * PLB / (Km_PKA_PLB + PLB);
		double PLB_dephosph = k_PP1_PLB * PP1 * PLBp / (Km_PP1_PLB + PLBp);
		double dPLBp = PLB_phosph - PLB_dephosph; //  ydot
		//  end PLB module
		// //  PLM module (from PLB, different total concentration)

		// double PLMtot = pin(4); //  p(102) = PLMtot; //  PLMtot    [uM]
		double k_PKA_PLM = K_rate_Scale * 54e-3; //  p(103) = 54;     //  k_pka_plb     [1/ms]
		double Km_PKA_PLM = 21; //  p(104) = 21;     //  Km_pka_plb    [uM]
		double k_PP1_PLM = K_rate_Scale * 8.5e-3; //  p(105) = 8.5;    //  k_pp1_plb     [1/ms]
		double Km_PP1_PLM = 7.0; //  p(106) = 7.0;    //  Km_pp1_plb    [uM]

		double PLM = PLMtot - PLMp;
		double PLM_phosph = k_PKA_PLM * PKACI * PLM / (Km_PKA_PLM + PLM);
		double PLM_dephosph = k_PP1_PLM * PP1 * PLMp / (Km_PP1_PLM + PLMp);
		double dPLMp = PLM_phosph - PLM_dephosph; //  ydot
		//  end PLM module
		// //  LCC module

		// double LCCtot = pin(5); // p(53) = LCCtot; //  LCCtot        [uM]
		double PKACII_LCCtot = PKACII_LCCtot_Scale * 0.025; // p(54) = 0.025;  //  PKAIIlcctot   [uM]
		double PP1_LCC = PP1_LCC_Scale * 0.025; // p(55) = 0.025;  //  PP1lcctot     [uM]
		double PP2A_LCC = PP2A_LCC_Scale * 0.025; // p(56) = 0.025;  //  PP2Alcctot    [uM]
		double k_PKA_LCC = K_rate_Scale * 54e-3;    // p(57) = 54;     //  k_pka_lcc     [1/ms]
		double Km_PKA_LCC = Km_PKA_LCC_Scale * 21;   // p(58) = 21;// *1.6;     //  Km_pka_lcc    [uM]
		double k_PP1_LCC = K_rate_Scale * 8.52e-3;  // p(59) = 8.52;   //  k_pp1_lcc     [1/ms] RABBIT, MOUSE
		// p(59) = 8.5;   //  k_pp1_lcc     [1/sec] RAT
		double Km_PP1_LCC = Km_PP1_LCC_Scale * 3;    // p(60) = 3;      //  Km_pp1_lcc    [uM]
		double k_PP2A_LCC = K_rate_Scale * 10.1e-3;  // p(61) = 10.1;   //  k_pp2a_lcc    [1/ms]
		double Km_PP2A_LCC = Km_PP2A_LCC_Scale * 3;    // p(62) = 3;      //  Km_pp2a_lcc   [uM]

		double PKACII_LCC = (PKACII_LCCtot / PKAIItot) * PKACII;
		double LCCa = LCCtot - LCCap;
		double LCCa_phosph = epsilon * k_PKA_LCC * PKACII_LCC * LCCa / (Km_PKA_LCC + epsilon * LCCa);
		double LCCa_dephosph = epsilon * k_PP2A_LCC * PP2A_LCC * LCCap / (Km_PP2A_LCC + epsilon * LCCap);
		double dLCCap = LCCa_phosph - LCCa_dephosph; //  ydot
		double LCCb = LCCtot - LCCbp;
		double LCCb_phosph = epsilon * k_PKA_LCC * PKACII_LCC * LCCb / (Km_PKA_LCC + epsilon * LCCb);
		double LCCb_dephosph = epsilon * k_PP1_LCC * PP1_LCC * LCCbp / (Km_PP1_LCC + epsilon * LCCbp);
		double dLCCbp = LCCb_phosph - LCCb_dephosph; //  ydot
		//  end LCC module
		// //  RyR module

		// double RyRtot = pin(6); // p(63) = RyRtot; //  RyRtot        [uM]
		double PKAIIryrtot = PKAIIryrtot_Scale * 0.034; // p(64) = 0.034;  //  PKAIIryrtot   [uM]
		double PP1ryr = PP1ryr_Scale * 0.034; // p(65) = 0.034;  //  PP1ryr        [uM]
		double PP2Aryr = PP2Aryr_Scale * 0.034; // p(66) = 0.034;  //  PP2Aryr       [uM]
		double kcat_pka_ryr = K_rate_Scale *  54e-3;   // p(67) = 54;     //  kcat_pka_ryr  [1/ms]
		double Km_pka_ryr = Km_pka_ryr_Scale * 21;   // p(68) = 21;     //  Km_pka_ryr    [uM]
		double kcat_pp1_ryr = K_rate_Scale * 8.52e-3;  // p(69) = 8.52;   //  kcat_pp1_ryr  [1/ms]
		double Km_pp1_ryr = Km_pp1_ryr_Scale * 7;    // p(70) = 7;      //  Km_pp1_ryr    [uM]
		double kcat_pp2a_ryr = K_rate_Scale * 10.1e-3;  // p(71) = 10.1;   //  kcat_pp2a_ryr [1/ms]
		double Km_pp2a_ryr = Km_pp2a_ryr_Scale * 4.1;  // p(72) = 4.1;    //  Km_pp2a_ryr   [uM]

		double PKACryr = (PKAIIryrtot / PKAIItot) * PKACII;
		double RyR = RyRtot - RyRp;
		double RyRPHOSPH = epsilon * kcat_pka_ryr * PKACryr * RyR / (Km_pka_ryr + epsilon * RyR);
		double RyRDEPHOSPH1 = epsilon * kcat_pp1_ryr * PP1ryr * RyRp / (Km_pp1_ryr + epsilon * RyRp);
		double RyRDEPHOSPH2A = epsilon * kcat_pp2a_ryr * PP2Aryr * RyRp / (Km_pp2a_ryr + epsilon * RyRp);
		double dRyRp = RyRPHOSPH - RyRDEPHOSPH1 - RyRDEPHOSPH2A; //  ydot
		//  end RyR module
		// //  TnI module

		// double TnItot = pin(7); // p(73) = TnItot; //  TnItot        [uM]
		double PP2A_TnI = PP2A_TnI_Scale * 0.67;   //  PP2Atni       [uM]
		double k_PKA_TnI = K_rate_Scale * 54e-3;    //  kcat_pka_tni  [1/ms]
		double Km_PKA_TnI = 21;     //  Km_pka_tni    [uM]
		double k_PP2A_TnI = K_rate_Scale * 10.1e-3;  //  kcat_pp2a_tni [1/ms]
		double Km_PP2A_TnI = 4.1;    //  Km_pp2a_tni   [uM]

		double TnIn = TnItot - TnIp;
		double TnI_phosph = k_PKA_TnI * PKACI * TnIn / (Km_PKA_TnI + TnIn);
		double TnI_dephosph = k_PP2A_TnI * PP2A_TnI * TnIp / (Km_PP2A_TnI + TnIp);
		double dTnIp = TnI_phosph - TnI_dephosph; //  ydot(26)
		//  end TnI module
		// //  Myofilament module (from TnI)



		// not used in atria.
		// double Myotot = pin(8);        //  Myotot_bar    [uM]
		double PP2A_myo = 0.67;             //  PP2Amyo       [uM]
		double kcat_pka_myo = K_rate_Scale * 54e-3;         //  kcat_pka_myo  [1/ms]
		double Km_pka_myo = 21;            //  Km_pka_myo    [uM]
		double kcat_pp2a_myo = K_rate_Scale * 10.1e-3;      //  kcat_pp2a_myo [1/ms]
		double Km_pp2a_myo = 4.1;          //  Km_pp2a_myo   [uM]

		double Myon = Myotot - Myop; //  Non-phos = tot - phos
		double MyoPHOSPH = kcat_pka_myo * PKACI * Myon / (Km_pka_myo + Myon);
		double MyoDEPHOSPH = kcat_pp2a_myo * PP2A_myo * Myop / (Km_pp2a_myo + Myop);
		double dMyop = MyoPHOSPH - MyoDEPHOSPH; //  ydot
		//  end Myo module

		////  IKs module

		// double IKstot = pin(9);        //  IKstot_bar    [uM]
		double Yotiao_tot = 0.025;         //  Yotiao_tot    [uM]
		double K_yotiao = K_rate_Scale * 0.1e-3;         //  K_yotiao      [uM] ** apply G589D mutation here: 1e2 **
		double PKAII_ikstot = PKAII_ikstot_Scale * 0.025;     //  PKAII_ikstot  [uM]
		double PP1_ikstot = PP1_ikstot_Scale * 0.025;       //  PP1_ikstot    [uM]
		double k_pka_iks = K_rate_Scale * 1.87e-3; // 54;      //  k_pka_iks     [1/ms] //  adjusted as in Xie et al 2013
		double Km_pka_iks = Km_pka_iks_Scale * 21;          //  Km_pka_iks    [uM]
		double k_pp1_iks = K_rate_Scale * 0.19e-3; // 8.52;    //  k_pp1_iks     [1/ms] //  adjusted as in Xie et al 2013
		double Km_pp1_iks = Km_pp1_iks_Scale * 7;           //  Km_pp1_iks    [uM]

		//  Effect of G589D mutation (IKs-Yotiao)
		double y1 = (-(K_yotiao + IKstot - Yotiao_tot) + sqrt((K_yotiao + IKstot - Yotiao_tot) * (K_yotiao + IKstot - Yotiao_tot) + 4 * K_yotiao * Yotiao_tot)) / 2;
		double x1 = IKstot / (1 + y1 / K_yotiao);
		double y2 = (-(K_yotiao + IKstot - Yotiao_tot) - sqrt((K_yotiao + IKstot - Yotiao_tot) * (K_yotiao + IKstot - Yotiao_tot) + 4 * K_yotiao * Yotiao_tot)) / 2;
		double x2 = IKstot / (1 + y2 / K_yotiao);

		double free_IKs;// = x1 * (y1 > 0) + x2 * (y1 <= 0);
		double free_Yotiao;// = y1 * (y1 > 0) + y2 * (y1 <= 0);
		free_IKs = y1 > 0 ? x1 : x2;
		free_Yotiao = y1 > 0 ? y1 : y2;


		double IksYot = free_IKs * free_Yotiao / K_yotiao; //  [uM] //  IKs-Yotiao

		double PKACiks = IksYot / IKstot * (PKAII_ikstot / PKAIItot) * PKACII;
		double PP1iks = IksYot / IKstot * PP1_ikstot;

		double KSn = IKstot - KSp; //  Non-phos = tot - phos
		double IKS_PHOSPH = epsilon * k_pka_iks * PKACiks * KSn / (Km_pka_iks + epsilon * KSn);
		double IKS_DEPHOSPH = epsilon * k_pp1_iks * PP1iks * KSp / (Km_pp1_iks + epsilon * KSp);
		double dKSp = IKS_PHOSPH - IKS_DEPHOSPH; //  ydot
		//  end Iks module
		// //  IKr module (from IKs, without mutation)

		// double Krtot = pin(10);        //  IKrtot_bar    [uM]
		double PKAII_ikrtot = 0.025;       //  PKAII_ikrtot  [uM]
		double PP1_ikrtot = PP1_ikrtot_Scale * 0.025;         //  PP1_ikrtot    [uM]
		double k_pka_ikr = K_rate_Scale * 1.87e-3; // 54;      //  k_pka_ikr     [1/ms] //  adjusted as in Xie et al 2013
		double Km_pka_ikr = 21;            //  Km_pka_ikr    [uM]
		double k_pp1_ikr = K_rate_Scale * 0.19e-3; // 8.52;    //  k_pp1_ikr     [1/ms] //  adjusted as in Xie et al 2013
		double Km_pp1_ikr = 7;             //  Km_pp1_ikr    [uM]

		double KRn = IKrtot - KRp; //  Non-phos = tot - phos
		double PKACikr = (PKAII_ikrtot / PKAIItot) * PKACII;
		double IKR_PHOSPH = epsilon * k_pka_ikr * PKACikr * KRn / (Km_pka_ikr + epsilon * KRn);
		double IKR_DEPHOSPH = epsilon * k_pp1_ikr * PP1_ikrtot * KRp / (Km_pp1_ikr + epsilon * KRp);
		double dKRp = IKR_PHOSPH - IKR_DEPHOSPH;
		//  end IKr module


		//// Ikur module (from Ikr)
		//IKurtot = pin(9);// p(95) = IKurtot;   // Ikur_tot      [uM]
		double PKAII_KURtot = 0.025; // p (96) = 0.025;  // PKAII_KURtot [uM]
		double PP1_KURtot = PP1_KURtot_Scale * 0.025; // p(97) = 0.025;  // PP1_KURtot   [uM]
		double k_pka_KUR = K_rate_Scale * 1.87e-3; // 54e-3; // k_pka_ikur [1/ms] // adjusted as in Xie et al 2013
		double Km_pka_KUR = 21; // p(99) = 21;    // Km_pka_KUR   [uM]
		double k_pp1_KUR = K_rate_Scale * 0.19e-3; // 8.52e-3; // k_pp1_ikur [1/ms] // adjusted as in Xie et al 2013
		double Km_pp1_KUR = 7; // p(101) = 7;     // Km_pp1_KUR   [uM]

		double KURn = IKurtot - KURp;  // Non-phos = tot - phos
		double PKACII_KUR = (PKAII_KURtot / PKAIItot) * PKACII; // (PKA_KURtot/PKAIItot)*PKAIIact
		double KURphos = epsilon * PKACII_KUR * k_pka_KUR * KURn / (Km_pka_KUR + epsilon * KURn);
		double KURdephos = epsilon * PP1_KURtot * k_pp1_KUR * KURp / (Km_pp1_KUR + epsilon * KURp);
		double dKURp = KURphos - KURdephos;
		// end Ikur module


		// //  ICFTR module


		// not used in atria.
		// double CFTRtot = pin(11); //  ICFTRtot_bar  [uM]
		double PKAII_CFTRtot = 0.025;      //  PKAII_CFTRtot [uM]
		double PP1_CFTRtot = PP1_CFTRtot_Scale * 0.025;        //  PP1_CFTRtot   [uM]
		double k_pka_CFTR = K_rate_Scale * 54e-3;           //  k_pka_CFTR    [1/ms]
		double Km_pka_CFTR = 8.5;          //  Km_pka_CFTR   [uM]
		double k_pp1_CFTR = K_rate_Scale * 8.52e-3;         //  k_pp1_CFTR    [1/ms]
		double Km_pp1_CFTR = 7;            //  Km_pp1_CFTR   [uM]


		// not used in atria
		double CFTRn = ICFTRtot - CFTRp;  //  Non-phos = tot - phos
		double PKAC_CFTR = (PKAII_CFTRtot / PKAIItot) * PKACII; //  (PKACFTRtot/PKAIItot)*PKAIIact
		double CFTRphos = epsilon * CFTRn * PKAC_CFTR * k_pka_CFTR / (Km_pka_CFTR + epsilon * CFTRn);
		double CFTRdephos = PP1_CFTRtot * k_pp1_CFTR * epsilon * CFTRp / (Km_pp1_CFTR + epsilon * CFTRp);
		double dCFTRp = CFTRphos - CFTRdephos;
		//  end ICFTR module
		// //  IClCa module (from CFTR)

		// double ClCatot = pin(12);	//  IClCatot_bar  [uM]
		double PKAII_ClCatot = 0.025;      //  PKAII_ClCatot [uM]
		double PP1_ClCatot = PP1_ClCatot_Scale * 0.025;        //  PP1_ClCatot   [uM]
		double k_pka_ClCa = K_rate_Scale *  54e-3;          //  k_pka_ClCa    [1/ms]
		double Km_pka_ClCa = 8.5;          //  Km_pka_ClCa   [uM]
		double k_pp1_ClCa = K_rate_Scale * 8.52e-3;         //  k_pp1_ClCa    [1/ms]
		double Km_pp1_ClCa = 7;            //  Km_pp1_ClCa   [uM]

		double ClCan = IClCatot - ClCap;  //  Non-phos = tot - phos
		double PKAC_ClCa = (PKAII_ClCatot / PKAIItot) * PKACII; //  (PKACFTRtot/PKAIItot)*PKAIIact
		double ClCaphos = epsilon * ClCan * PKAC_ClCa * k_pka_ClCa / (Km_pka_ClCa + epsilon * ClCan);
		double ClCadephos = PP1_ClCatot * k_pp1_ClCa * epsilon * ClCap / (Km_pp1_ClCa + epsilon * ClCap);
		double dClCap = ClCaphos - ClCadephos;
		//  end ICl(Ca) module
		// //  ydot


		//// Ina module (from Ikr)

		//INatot = pin(6); // p(79) [uM]
		double PKAIIinatot = 0.025; // PKAII_inatot [uM]
		double PP1_inatot = PP1_inatot_Scale * 0.025; // PP1_inatot [uM]
		double k_pka_ina = K_rate_Scale * 1.87e-3; // 54e-3; // k_pka_ina [1/ms] // adjusted as in Xie et al 2013
		double Km_pka_ina = 21; // Km_pka_ina [uM]
		double k_pp1_ina = K_rate_Scale * 0.19e-3; // 8.52e-3; // k_pp1_ina [1/ms] // adjusted as in Xie et al 2013
		double Km_pp1_ina = 7; // Km_pp1_ina [uM]

		double NAn = INatot - NAp;
		double PKACII_ina = (PKAIIinatot / PKAIItot) * PKACII;
		double INA_PHOSPH = epsilon * k_pka_ina * PKACII_ina * NAn / (Km_pka_ina + epsilon * NAn);
		double INA_DEPHOSPH = epsilon * k_pp1_ina * PP1_inatot * NAp / (Km_pp1_ina + epsilon * NAp);
		double dNAp = INA_PHOSPH - INA_DEPHOSPH;
		// end Ina module
		//// Ito module (from Ikr)

		//Itotot = pin(6); // p(79) [uM]
		double PKAIIitotot = 0.025; // PKAII_ikrtot [uM]
		double PP1_itotot = PP1_itotot_Scale * 0.025; // PP1_ikrtot [uM]
		double k_pka_ito = K_rate_Scale * 1.87e-3; // 54e-3; // k_pka_ikr [1/ms] // adjusted as in Xie et al 2013
		double Km_pka_ito = 21; // Km_pka_ikr [uM]
		double k_pp1_ito = K_rate_Scale * 0.19e-3; // 8.52e-3; // k_pp1_ikr [1/ms] // adjusted as in Xie et al 2013
		double Km_pp1_ito = 7; // Km_pp1_ikr [uM]

		double TOn = Itotot - TOp;
		double PKACII_ito = (PKAIIitotot / PKAIItot) * PKACII;
		double Ito_PHOSPH = epsilon * k_pka_ito * PKACII_ito * TOn / (Km_pka_ito + epsilon * TOn);
		double Ito_DEPHOSPH = epsilon * k_pp1_ito * PP1_itotot * TOp / (Km_pp1_ito + epsilon * TOp);
		double dTOp = Ito_PHOSPH - Ito_DEPHOSPH;
		// end Ito module
		//// Ik1 module (from Ikr)

		//IK1tot = pin(6); // p(79) [uM]
		double PKAIIik1tot = 0.025; // PKAII_ikrtot [uM]
		double PP1_ik1tot = PP1_ik1tot_Scale * 0.025; // PP1_ikrtot [uM]
		double k_pka_ik1 = K_rate_Scale * 1.87e-3; // 54e-3; // k_pka_ikr [1/ms] // adjusted as in Xie et al 2013
		double Km_pka_ik1 = 21; // Km_pka_ikr [uM]
		double k_pp1_ik1 = K_rate_Scale * 0.19e-3; // 8.52e-3; // k_pp1_ikr [1/ms] // adjusted as in Xie et al 2013
		double Km_pp1_ik1 = 7; // Km_pp1_ikr [uM]

		double K1n = IK1tot - K1p;
		double PKACII_ik1 = (PKAIIik1tot / PKAIItot) * PKACII;
		double IK1_PHOSPH = epsilon * k_pka_ik1 * PKACII_ik1 * K1n / (Km_pka_ik1 + epsilon * K1n);
		double IK1_DEPHOSPH = epsilon * k_pp1_ik1 * PP1_ik1tot * K1p / (Km_pp1_ik1 + epsilon * K1p);
		double dK1p = IK1_PHOSPH - IK1_DEPHOSPH;


		// ydot(10)=dPDEp;
		// ydot[0] = dLR;
		// ydot[1 ] = dLRG;
		// ydot[2 ] = dRG;
		// ydot[3 ] = db1AR_S464;
		// ydot[4 ] = db1AR_S301;
		// ydot[5 ] = dGsaGTPtot;
		// ydot[6 ] = dGsaGDP;
		// ydot[7 ] = dGsby;
		// ydot[8 ] = dAC_GsaGTP;
		// ydot[9 ] = dPDE3p;
		// ydot[10] = dPDE4p;
		// ydot[11] = dcAMPtot;
		// ydot[12] = dRC_I;
		// ydot[13] = dRCcAMP_I;
		// ydot[14] = dRCcAMPcAMP_I;
		// ydot[15] = dRcAMPcAMP_I;
		// ydot[16] = dPKACI;
		// ydot[17] = dPKACI_PKI;
		// ydot[18] = dRC_II;
		// ydot[19] = dRCcAMP_II;
		// ydot[20] = dRCcAMPcAMP_II;
		// ydot[21] = dRcAMPcAMP_II;
		// ydot[22] = dPKACII;
		// ydot[23] = dPKACII_PKI;
		// ydot[24] = dI1p_PP1; //  output CaMKII
		// ydot[25] = dI1ptot;
		// ydot[26] = dPLBp; //  output
		// ydot[27] = dPLMp; //  output
		// ydot[28] = dLCCap; //  output
		// ydot[29] = dLCCbp; //  output
		// ydot[30] = dRyRp; //  output
		// ydot[31] = dTnIp; //  output
		// ydot[32] = dMyop; //  output
		// ydot[33] = dKSp; //  output
		// ydot[34] = dKRp; //  output
		// ydot[35] = dCFTRp; //  output
		// ydot[36] = dClCap; //  output
		// ydot[37] = dNAp; // output
		// ydot[38] = dTOp; // output
		// ydot[39] = dK1p; // output
		// ydot[40] = dKURp; // output

		ydot[0]  = dLR;
		ydot[1]  = dLRG;
		ydot[2]  = dRG;
		ydot[3]  = db1AR_S464;
		ydot[4]  = db1AR_S301;
		ydot[5]  = dGsaGTPtot;
		ydot[6]  = dGsaGDP;
		ydot[7]  = dGsby;
		ydot[8]  = dAC_GsaGTP;
		ydot[9]  = dPDE4p;
		ydot[10] = dcAMPtot;
		ydot[11] = dRC_I;
		ydot[12] = dRCcAMP_I;
		ydot[13] = dRCcAMPcAMP_I;
		ydot[14] = dRcAMPcAMP_I;
		ydot[15] = dPKACI;
		ydot[16] = dPKACI_PKI;
		ydot[17] = dRC_II;
		ydot[18] = dRCcAMP_II;
		ydot[19] = dRCcAMPcAMP_II;
		ydot[20] = dRcAMPcAMP_II;
		ydot[21] = dPKACII;
		ydot[22] = dPKACII_PKI;
		ydot[23] = dI1p_PP1; // output CaMKII
		ydot[24] = dI1ptot;
		ydot[25] = dPLBp; // output
		ydot[26] = dPLMp; // output
		ydot[27] = dLCCap; // output
		ydot[28] = dLCCbp; // output
		ydot[29] = dRyRp; // output
		ydot[30] = dTnIp; // output
		ydot[31] = dKSp; // output
		ydot[32] = dKRp; // output
		ydot[33] = dKURp; // output
		ydot[34] = dCFTRp; // output
		ydot[35] = dClCap; // output
		ydot[36] = dMyop; // output
		ydot[37] = dNAp; // output
		ydot[38] = dTOp; // output
		ydot[39] = dK1p; // output
		ydot[40]  = dPDE3p;  // add PDE3p for human model





		// ydot[0]  = dLR;
		// ydot[1]  = dLRG;
		// ydot[2]  = dRG;
		// ydot[3]  = db1AR_S464;
		// ydot[4]  = db1AR_S301;
		// ydot[5]  = dGsaGTPtot;
		// ydot[6]  = dGsaGDP;
		// ydot[7]  = dGsby;
		// ydot[8]  = dAC_GsaGTP;
		// ydot[9]  = dPDE3p;
		// ydot[10] = dcAMPtot;
		// ydot[11] = dRC_I;
		// ydot[12] = dRCcAMP_I;
		// ydot[13] = dRCcAMPcAMP_I;
		// ydot[14] = dRcAMPcAMP_I;
		// ydot[15] = dPKACI;
		// ydot[16] = dPKACI_PKI;
		// ydot[17] = dRC_II;
		// ydot[18] = dRCcAMP_II;
		// ydot[19] = dRCcAMPcAMP_II;
		// ydot[20] = dRcAMPcAMP_II;
		// ydot[21] = dPKACII;
		// ydot[22] = dPKACII_PKI;
		// ydot[23] = dI1p_PP1; // output CaMKII
		// ydot[24] = dI1ptot;
		// ydot[25] = dPLBp; // output
		// ydot[26] = dPLMp; // output
		// ydot[27] = dLCCap; // output
		// ydot[28] = dLCCbp; // output
		// ydot[29] = dRyRp; // output
		// ydot[30] = dTnIp; // output
		// ydot[31] = dKSp; // output
		// ydot[32] = dKRp; // output
		// ydot[33] = dKURp; // output
		// ydot[34] = dCFTRp; // output
		// ydot[35] = dClCap; // output
		// ydot[36] = dmyop; // output
		// ydot[37] = dNAp; // output
		// ydot[38] = dTOp; // output
		// ydot[39] = dK1p; // output
		// ydot[40] = dPDE4p;

		return 0;

	}

	void initialiser() {
		// y[0]  = 8.14039232720744e-35;0
		// y[1]  = 4.53346233860292e-33;1.50000000000000e-323
		// y[2]  = 0.000480248667287148;0.00282631217659216
		// y[3]  = 4.31682189932130e-35;2.80092307054198e-34
		// y[4]  = 0.000648227269322366;0.000422049845562602
		// y[5]  = 0.00960497334574349;0.0565262435318454
		// y[6]  = 0.000621006098743802;0.000653599606880587
		// y[7]  = 0.0102259794442235;0.0571798431376496
		// y[8]  = 0.00141579105756161;0.00533273395766946
		// y[9]  = 0.00222180557941371;0.000455252770299857
		// y[10] = 1.02250237329510;0.492824467402742
		// y[11] = 0.804483262433527;0.989246288005767
		// y[12] = 0.141686904193383;0.0418093107110889
		// y[13] = 0.00447754605803494;0.000317058381423304
		// y[14] = 0.229043453816613;0.148318509390643
		// y[15] = 0.0855264085372521;0.00935237567061519
		// y[16] = 0.143516934727477;0.138966023176643
		// y[17] = 0.0510346588626229;0.0742118800151811
		// y[18] = 0.00898830735015808;0.00313647631295501
		// y[19] = 0.000284045730078040;2.37852785957148e-05
		// y[20] = 0.0576887979515687;0.0406236682838489
		// y[21] = 0.0215414450156286;0.00256157551132883
		// y[22] = 0.0361474568044903;0.0380621966453375
		// y[23] = 0.0727115692379701;0.00910389864720390
		// y[24] = 0.0728005360781959;0.00911423346266423
		// y[25] = 4.80816797556955;0.316908017800809
		// y[26] = 5.60341884069710;0.343822128706968
		// y[27] = 0.00548941510824692;0.000563675197467848
		// y[28] = 0.00626643428486632;0.000665675123735624
		// y[29] = 0.0275772651433693;0.00266714957326578
		// y[30] = 4.38884564916186;0.249491847096173
		// y[31] = 0.0137120366693140;0.00225809535533418
		// y[32] = 0.0137120366693140;0.00225809535533323
		// y[33] = 0.0137120366693140;0.0137120366693140
		// y[34] = 0.0164709839871346;0.00337311507850651
		// y[35] = 0.0164709839871346;0.00337311507850651
		// y[36] = 4.38884564916186;0.249491847096173
		// y[37] = 0.0137120366693140;0.0137120366693140
		// y[38] = 0.0137120366693140;0.0137000000000000
		// y[39] = 0.0137120366693140;0.0137000000000000
		// 0.000455252770299857
		y[0]  = 0;
		y[1]  = 1.50000000000000e-323;
		y[2]  = 0.00282631217659216;
		y[3]  = 2.80092307054198e-34;
		y[4]  = 0.000422049845562602;
		y[5]  = 0.0565262435318454;
		y[6]  = 0.000653599606880587;
		y[7]  = 0.0571798431376496;
		y[8]  = 0.00533273395766946;
		y[9]  = 0.000455252770299857;
		y[10] = 0.492824467402742;
		y[11] = 0.989246288005767;
		y[12] = 0.0418093107110889;
		y[13] = 0.000317058381423304;
		y[14] = 0.148318509390643;
		y[15] = 0.00935237567061519;
		y[16] = 0.138966023176643;
		y[17] = 0.0742118800151811;
		y[18] = 0.00313647631295501;
		y[19] = 2.37852785957148e-05;
		y[20] = 0.0406236682838489;
		y[21] = 0.00256157551132883;
		y[22] = 0.0380621966453375;
		y[23] = 0.00910389864720390;
		y[24] = 0.00911423346266423;
		y[25] = 0.316908017800809;
		y[26] = 0.343822128706968;
		y[27] = 0.000563675197467848;
		y[28] = 0.000665675123735624;
		y[29] = 0.00266714957326578;
		y[30] = 0.249491847096173;
		y[31] = 0.00225809535533418;
		y[32] = 0.00225809535533323;
		y[33] = 0.0137120366693140;
		y[34] = 0.00337311507850651;
		y[35] = 0.00337311507850651;
		y[36] = 0.249491847096173;
		y[37] = 0.0137120366693140;
		y[38] = 0.0137000000000000;
		y[39] = 0.0137000000000000;
		y[40] = 0.000455252770299857;



		// new baseline data
		y[0]  = 0.0;
		y[1]  = 0.0;
		y[2]  = 0.002826312;
		y[3]  = 9.324460999999999e-35;
		y[4]  = 0.0004220498;
		y[5]  = 0.05652624;
		y[6]  = 0.0006535996;
		y[7]  = 0.05717984;
		y[8]  = 0.005332734;
		y[9]  = 0.0004552528;
		y[10] = 0.4928245;
		y[11] = 0.9892463;
		y[12] = 0.04180931;
		y[13] = 0.0003170584;
		y[14] = 0.1483185;
		y[15] = 0.009352375999999999;
		y[16] = 0.138966;
		y[17] = 0.07421188;
		y[18] = 0.0031364759999999996;
		y[19] = 2.378528e-05;
		y[20] = 0.04062367;
		y[21] = 0.002561576;
		y[22] = 0.0380622;
		y[23] = 0.009103899;
		y[24] = 0.009114233000000001;
		y[25] = 0.31690799999999997;
		y[26] = 0.34382209999999996;
		y[27] = 0.0005636752;
		y[28] = 0.0006656750999999999;
		y[29] = 0.00266715;
		y[30] = 0.24949179999999996;
		y[31] = 0.002258095;
		y[32] = 0.002258095;
		y[33] = 0.002544429;
		y[34] = 0.003373115;
		y[35] = 0.003373115;
		y[36] = 0.24949179999999996;
		y[37] = 0.002544429;
		y[38] = 0.002544124;
		y[39] = 0.002544124;
		y[40] = 0.0004552528;


		y[0]  = 0;
		y[1]  = 4.90000000000000e-324;
		y[2]  =  0.00278719440918818;
		y[3]  =  -3.56383965716342e-38;
		y[4]  = 0.000809062546121539;
		y[5]  = 0.0557438881837657;
		y[6]  =  0.000653495970006231;
		y[7]  = 0.0563973841516459;
		y[8]  = 0.00526649609876001;
		y[9]  = 0.000649341565200468;
		y[10] = 0.569723858593256;
		y[11] = 0.961031412390004;
		y[12] = 0.0594156694446770;
		y[13] = 0.000659116309999281;
		y[14] = 0.158584968334472;
		y[15] = 0.0181835257561421;
		y[16] = 0.140401332055390;
		y[17] = 0.0705998038080935;
		y[18] = 0.00436482569855733;
		y[19] = 4.84203550193421e-05;
		y[20] = 0.0429827600251848;
		y[21] = 0.00492846557748966;
		y[22] = 0.0380543983245343;
		y[23] = 0.0174471587631342;
		y[24] = 0.0174671542957124;
		y[25] = 0.648149274513690;
		y[26] = 0.706550421700249;
		y[27] = 0.00106431129149780;
		y[28] = 0.00125265011818498;
		y[29] = 0.00506920915849830;
		y[30] = 0.514141026957991;
		y[31] = 0.00402115953635574;
		y[32] = 0.00402115953635397;
		y[33] = 0.00402115953635397;
		y[34] = 0.0058;//0.003373115;
		y[35] = 0.00579815110518361;
		y[36] = 0.514141026957991;
		y[37] = 0.00402115953635397;
		y[38] = 0.00402115953635397;
		y[39] = 0.00402115953635397;
		y[40] = 0.000649341565200468;
	}

	void print_to_file(double t, std::ofstream & output_file) {


		// output_file <<  std::setprecision(9);
		output_file <<  std::setprecision(9) << t << std::setprecision(4) << " ";


		for (int i = 0; i < 41; ++i)
		{
			output_file << y[i] << " ";
		}
		output_file << std::endl;
	}


	void print_to_file_percent(double t, std::ofstream &output_file) {


		double PLB_PKAn   = (PLBtot - y[25]) / PLBtot; // non-phosphorylated PLB targets
		double PLM_PKAp   = y[26] / PLMtot;   //
		double LCCa_PKAp  = y[27] / LCCtot;   //
		double LCCb_PKAp  = y[28] / LCCtot;   //
		double RyR_PKAp   = y[29] / RyRtot;   //
		double TnI_PKAp   = y[30] / TnItot;   //
		double IKs_PKAp   = y[31] / IKstot;   //
		double IKr_PKAp   = y[32] / IKrtot;   //
		double IKur_PKAp  = y[33] / IKurtot;   //
		double ICFTR_PKAp = y[34] / ICFTRtot;   //
		double IClCa_PKAp = y[35] / IClCatot;   //
		double Myo_PKAp   = y[36] / Myotot;   //
		double INa_PKAp   = y[37] / INatot;   //
		double Ito_PKAp   = y[38] / Itotot;   //
		double IK1_PKAp   = y[39] / IK1tot;   //

		output_file <<  std::setprecision(9)
		            << t << " "
		            <<  std::setprecision(4)

		            << PLB_PKAn   << " "
		            << PLM_PKAp   << " "
		            << LCCa_PKAp  << " "
		            << LCCb_PKAp  << " "
		            << RyR_PKAp   << " "
		            << TnI_PKAp   << " "
		            << IKs_PKAp   << " "
		            << IKr_PKAp   << " "
		            << IKur_PKAp  << " "
		            << ICFTR_PKAp << " "
		            << IClCa_PKAp << " "
		            << Myo_PKAp   << " "
		            << INa_PKAp   << " "
		            << Ito_PKAp   << " "
		            << IK1_PKAp   << std::endl;

	}


	void set_scale_parameters_long(double scale[]) {
		LCCtot_Scale        = scale[0];
		RyRtot_Scale        = scale[1];
		PLBtot_Scale        = scale[2];
		TnItot_Scale        = scale[3];
		PLMtot_Scale        = scale[4];
		b1ARtot_Scale       = scale[5];
		Gstot_Scale         = scale[6];
		ACtot_Scale         = scale[7];
		ATP_Scale           = scale[8];
		PDE3tot_Scale       = scale[9];
		PDE4tot_Scale       = scale[10];
		PKItot_Scale        = scale[11];
		PKAIItot_Scale      = scale[12];
		I1tot_Scale         = scale[13];
		PP1tot_Scale        = scale[14];
		PKACII_LCCtot_Scale = scale[15];
		PP1_LCC_Scale       = scale[16];
		PP2A_LCC_Scale      = scale[17];
		PKAIIryrtot_Scale   = scale[18];
		PP1ryr_Scale        = scale[19];
		PP2Aryr_Scale       = scale[20];
		PP1_ikstot_Scale    = scale[21];
		PKAII_ikstot_Scale  = scale[22];
		Km_AC_basal_Scale   = scale[23];
		Kd_AC_Gsa_Scale     = scale[24];
		Km_PDE3_cAMP_Scale  = scale[25];
		Km_PDE4_cAMP_Scale  = scale[26];
		Km_PKA_I1_Scale     = scale[27];
		Km_PP2A_I1_Scale    = scale[28];
		Km_PKA_PLB_Scale    = scale[29];
		Km_PP1_PLB_Scale    = scale[30];
		Km_PKA_LCC_Scale    = scale[31];
		Km_PP1_LCC_Scale    = scale[32];
		Km_PP2A_LCC_Scale   = scale[33];
		Km_pka_ryr_Scale    = scale[34];
		Km_pp1_ryr_Scale    = scale[35];
		Km_pp2a_ryr_Scale   = scale[36];
		Km_PKA_TnI_Scale    = scale[37];
		Km_PP2A_TnI_Scale   = scale[38];
		Km_pka_iks_Scale    = scale[39];
		Km_pp1_iks_Scale    = scale[40];
		K_rate_Scale = scale[41];
	}


	void set_scale_parameters(double scale[]) {
		b1ARtot_Scale      = scale[0];
		Gstot_Scale        = scale[1];
		ACtot_Scale        = scale[2];
		ATP_Scale          = scale[3];
		PDE3tot_Scale      = scale[4];
		PDE4tot_Scale      = scale[5];
		PKItot_Scale       = scale[6];
		I1tot_Scale        = scale[7];
		PP1tot_Scale       = scale[8];
		PP1_LCC_Scale      = scale[9];
		PP2A_LCC_Scale     = scale[10];
		PP1ryr_Scale       = scale[11];
		PP2Aryr_Scale      = scale[12];
		PP1_ik1tot_Scale   = scale[13];
		PP1_itotot_Scale   = scale[14];
		PP1_inatot_Scale   = scale[15];
		PP1_ClCatot_Scale  = scale[16];
		PP1_KURtot_Scale   = scale[17];
		PP1_ikrtot_Scale   = scale[18];
		PP1_ikstot_Scale   = scale[19];
		PP2A_TnI_Scale     = scale[20];
		PP2A_I1_Vmax_Scale = scale[21];
	}


	void print_scale_parameters() {

		std::cerr << b1ARtot_Scale  << std::endl
		          << Gstot_Scale        << std::endl
		          << ACtot_Scale        << std::endl
		          << ATP_Scale          << std::endl
		          << PDE3tot_Scale      << std::endl
		          << PDE4tot_Scale      << std::endl
		          << PKItot_Scale       << std::endl
		          << I1tot_Scale        << std::endl
		          << PP1tot_Scale       << std::endl
		          << PP1_LCC_Scale      << std::endl
		          << PP2A_LCC_Scale     << std::endl
		          << PP1ryr_Scale       << std::endl
		          << PP2Aryr_Scale      << std::endl
		          << PP1_ik1tot_Scale   << std::endl
		          << PP1_itotot_Scale   << std::endl
		          << PP1_inatot_Scale   << std::endl
		          << PP1_ClCatot_Scale  << std::endl
		          << PP1_KURtot_Scale   << std::endl
		          << PP1_ikrtot_Scale   << std::endl
		          << PP1_ikstot_Scale   << std::endl
		          << PP2A_TnI_Scale     << std::endl
		          << PP2A_I1_Vmax_Scale << std::endl;
	}

};

#endif

