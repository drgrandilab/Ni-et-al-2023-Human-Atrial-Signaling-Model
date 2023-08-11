
#ifndef SIGNALLING_PARA__HPP
#define SIGNALLING_PARA__HPP
#include <fstream>
#include <cmath>
#include <cstdlib>
#include  <iomanip>
#include <iostream>


class betaAR_para
{
public:
	betaAR_para() {

		ISO      = 0;
		LCCtot   = 0.0250000000000000;
		RyRtot   = 0.135000000000000;
		PLBtot   = 38;
		TnItot   = 70;
		IKstot   = 0.0250000000000000;
		ICFTRtot = 0.0250000000000000;
		PP1tot   = 0.890000000000000;
		PLMtot   = 48;
		myotot   = 70;
		IKrtot   = 0.0250000000000000;
		IKurtot  = 0.0250000000000000;
		INatot   = 0.0250000000000000;
		IClCatot = 0.0250000000000000;
		Itotot   = 0.0250000000000000;
		IK1tot   = 0.0250000000000000;
		AF_index = 0;
	};
	~betaAR_para() {};
	double ISO;
	double LCCtot;   // pin(2);
	double RyRtot;   // pin(3);
	double PLBtot;   // pin(4);
	double TnItot;   // pin(5);
	double IKstot;   // pin(6);
	double ICFTRtot;   // pin(7);
	double PP1tot;   // pin(8);
	double PLMtot;   // pin(9);
	double myotot;   // pin(10);
	double IKrtot;   // pin(11);
	double IKurtot;   // pin(12);
	double INatot;   // pin(13);
	double IClCatot;   // pin(14);
	double Itotot;   // pin(15);
	double IK1tot;   // pin(16);
	double AF_index;   // pin(17); // AF_index
	double Ca;




};



class CaM_para
{
public:
	CaM_para() {
		K         =  135.0000;
		Mg        =  1.0000;
		CaMtot    =  418.0000;
		Btot      =  0;
		CaMKIItot =  120.0000;
		CaNtot    =  3.6;//3.61750874231279   ;
		PP1tot    =  96.5000;
		Ca        =  0.3494;
		// cyclelength =  333.3333    ;
		// compartment =  2.0000;
	};
	~CaM_para() {};
	double K         ; // =  135.0000    ;
	double Mg        ; // =  1.0000  ;
	double CaMtot    ; // =  418.0000         ;
	double Btot      ; // =  0  ;
	double CaMKIItot ; // =  120.0000    ;
	double CaNtot    ; // =  3.6;//3.61750874231279   ;
	double PP1tot    ; // =  96.5000    ;
	double Ca        ; // =  0.3494  ;
	// cyclelength   =  333.3333    ;
	// compartment   =  2.0000;

	double CaMtot_Scale    = 1;
	double CaMKIItot_Scale = 1;
	double CaNtot_Scale    = 1;
	double PP1tot_Scale    = 1;   // de-phosphorylate CaMKII

	void set_scale_parameters(double scale []) {

		CaMtot_Scale    = scale[0];
		CaMKIItot_Scale = scale[1];
		CaNtot_Scale    = scale[2];
		// PP1tot_Scale    = scale[3];   // de-phosphorylate CaMKII, not having much effect yet becaused PP1 tot is too large


		CaMtot    = CaMtot * CaMtot_Scale;
		CaMKIItot = CaMKIItot * CaMKIItot_Scale;
		CaNtot    = CaNtot * CaNtot_Scale;
		// PP1tot    = PP1tot * PP1tot_Scale;
	}

	void print_scale_parameters() {

		std::cerr << CaMtot_Scale << std::endl
		          << CaMKIItot_Scale << std::endl
		          << CaNtot_Scale    << std::endl;
		          // << PP1tot_Scale    << std::endl;
	}

};




class CaMKII_para
{
public:
	CaMKII_para() {
		CaMKIIactDyad = 5.39034430012935;
		LCCtotDyad     = 28.2600000000000;
		RyRtot         = 382.600000000000;
		PP1_dyad       = 95.7000000000000;
		PP2A_dyad      = 95.7600000000000;
		OA             = 0;
		PLBtot         = 38;
		NaVtot         = 30;
		CaMKIIactSL   = 3.34842094343477e-07;
		LCCtotSL       = 0.0846000000000000;
		PP1_SL         = 0.570000000000000;
		PP1_PLB_avail  = 0.999999607597786;
	};

	~CaMKII_para() {};
	// double LCC_PKAp, LCC_CKdyadp, RyR2809p, RyR2815p, PLBT17p, LCC_CKslp;
	double CaMKIIactDyad, LCCtotDyad, RyRtot, PP1_dyad, PP2A_dyad, OA, PLBtot, NaVtot,
	       CaMKIIactSL, LCCtotSL, PP1_SL, PP1_PLB_avail;
};





class param_cell
{
public:
	param_cell() {};
	~param_cell() {};
	double CaMKIIactDyad, LCCtotDyad, RyRtot, PP1_dyad, PP2A_dyad, OA, PLBtot, NaVtot,
	       CaMKIIactSL, LCCtotSL, PP1_SL, PP1_PLB_avail;
};


class param_yCell {
public:
	param_yCell() {};
	~param_yCell() {};

	double LCC_PKAp, LCC_CKdyadp, RyR2809p, RyR2815p, PLBT17p, LCC_CKslp;
};

class signalling_para
{
public:

	CaM_para cell_CaM_para_dyad, cell_CaM_para_sl, cell_CaM_para_cyto;
	CaMKII_para cell_CaMKII_para;
	betaAR_para cell_betaAR_para;

	int AF = 0;

	double LCC_CKdyadp;
	double RyR_CKp;
	double PLB_CKp;
	double NaV_CKp;
	double Ikur_CKp;
	double Itof_CKp;
	double Gapjunct_CKp;

	double LCCa_PKAp ;
	double LCCb_PKAp ;
	double PLB_PKAn  ;
	double PLM_PKAp  ;
	double RyR_PKAp  ;
	double TnI_PKAp  ;
	double IKs_PKAp  ;
	double IKr_PKAp  ;
	double IKur_PKAp ;
	double ICFTR_PKAp;  // not present in HAM
	double IClCa_PKAp;
	double Myo_PKAp  ;
	double INa_PKAp  ;
	double Ito_PKAp  ;
	double IK1_PKAp ;



	// scaling parameters for creating populations of models // 15:55:10, Mon, 09-September-2019, By Haibo
	// 1 - ECC modules
	double INa_Scale = 1.0;
	double INaL_Scale = 1.0;
	double INab_Scale = 1.0;
	double INaK_Scale = 1.0;
	double Itof_Scale = 1.0;
	double IKur_Scale = 1.0;
	double IK2p_Scale = 1.0;
	double IKr_Scale = 1.0;
	double IKs_Scale = 1.0;
	double IK1_Scale = 1.0;
	double IKp_Scale = 1.0;
	double IKach_Scale = 1.0;
	double ISK_Scale = 1.0;
	double ICaL_Scale = 1.0;
	double ICab_Scale = 1.0;
	double ICap_Scale = 1.0;
	double INCX_Scale = 1.0;
	double IClCa_Scale = 1.0;
	double IClb_Scale = 1.0;
	double Jrel_Scale = 1.0;
	double Jserca_Scale = 1.0;
	double Jleak_Scale = 1.0;
	double Cleft_Buffer_Scale = 1.0;
	double Cytosol_Buffer_Scale = 1.0;
	double SR_Buffer_Scale = 1.0;

	// 2 - BetaAR modules


	// 3 - CaMKII modules


	// LTCC
	double G_LTCC = 1.0;
	double K_PKA_LTCC = 0.0;
	double LTCC_junc_mode2 = 0;
	double LTCC_sl_mode2 = 0;
	double fckiim2_j = 0;
	double fpkam2 = 0;
	// RyR
	double fCKII_RyR = 0;
	double fPKA_RyR = 0;
	double RyR_koSRCa_Scale = 1;
	// PLB
	double fCKII_PLB = 1;
	double fPKA_PLB = 1;
	double PLB_kmf_Scale = 1;
	// INa
	double kPKA_Ina_increase = 0;
	double kCKII_Nav_change = 0;
	// PLM
	double KmNaip_PKA = 0;
	// TnI
	double fPKA_TnI = 1.0;
	// IKs
	double kPKA_Iks = 1;
	// IKr
	double kPKA_Ikr = 1;
	double Vkr_PKA = 0;
	double dGkr_PKA = 0;
	// IKur
	double kCKII_Ikur = 1;
	double kPKA_Ikur = 0;
	// IClCa
	double kPKA_IClCa = 0;
	// Myo
	double kPKA_Myo = 0;
	// Ito
	double kPKA_Itof = 0;
	double kCKII_Itof_tau = 1;
	double kCKII_Itof_G = 1;
	double kCKII_Itof_vshift = 0;
	// Ik1
	double kPKA_IK1 = 1;
	double kCKII_IK1_G = 1;
	double IK1_CKp = 0;

	// Gap junctions:
	double kCKII_Gapjunct_G = 1;


	double PLB_PKAn_Scale    = 1.0;
	double PLM_PKAp_Scale    = 1.0;
	double LCCa_PKAp_Scale   = 1.0;
	double LCCb_PKAp_Scale   = 1.0;
	double RyR_PKAp_Scale    = 1.0;
	double TnI_PKAp_Scale    = 1.0;
	double IKs_PKAp_Scale    = 1.0;
	double IKr_PKAp_Scale    = 1.0;
	double IKur_PKAp_Scale   = 1.0;
	double IClCa_PKAp_Scale  = 1.0;
	double Myo_PKAp_Scale    = 1.0;
	double INa_PKAp_Scale    = 1.0;
	double Ito_PKAp_Scale    = 1.0;
	double IK1_PKAp_Scale    = 1.0;
	double LCC_CKdyadp_Scale = 1.0;
	double NaV_CKp_Scale     = 1.0;
	double RyR_CKp_Scale     = 1.0;
	double PLB_CKp_Scale     = 1.0;
	double Ikur_CKp_Scale    = 1.0;
	double Itof_CKp_Scale    = 1.0;
	double IK1_CKp_Scale     = 1.0;
	double ICFTR_PKAp_Scale  = 1.0;


	bool CaMKII_Inhibition = false;

	bool CaMKII_Double = false;


	bool No_CaMKII_K1 = false;
	bool No_CaMKII_NaV = false;
	bool No_CaMKII_GapG = false;


	// simulation mode control:
	int Ca_clamp = 0;//p(12);
	int Na_clamp = 0; //p(13);

	void set_CaMKII_Inhibition();

	void update_signaling_modulators_for_new_bAR_Module();
	void print_to_file(double t, std::ofstream & output_file);
	void set_scaling_for_phos_levels(double *Scale);

};

#endif
