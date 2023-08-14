#include "SingleCellParameter.hpp"

SingleCellPara::SingleCellPara(AtriaCellType cell_type, IKurMuationType mutation_type, AFType AF_model)
	:    m_cell_type(cell_type),
	     m_mutation_type(mutation_type),
	     m_AF_model(AF_model)
{
	DefaultParameters();
	SetCellTypePara(cell_type);
	SetIKurMutationPara(mutation_type);
	SetAFTypePara(AF_model);
}

void SingleCellPara::DefaultParameters() {

	// Default scaling factors
	GNa = GCaL = GKur = GKr = GK1 = GKs = Gto =  GNaCa = GNaK = GCap = fIRel = RyR = GSR_leak = GCaT = GbCa = 1.0;
	IKur_ac_grad = IKur_inac_grad = IKs_grad = IKr_ac_grad = 1;
	BULK_CONST = 1.25;
	fRyR = IKur_ac_shift = IKur_inac_shift = IKs_shift =  IK1_v_shift = IKr_ac_shift  =  0.0;
	IKur_ac1_grad = IKur_ac2_grad = 1;
	IKur_ac1_mult = IKur_inac_mult = IKur_ac2_mult = 1;
	IKur_ac_add = IKur_inac_add = 0;
	IKur_ac1_shift = IKur_ac2_shift = 0;
	Ito_Af = 1.0;
	If_vshift = 0;
	If_grad = 1;
	Gf = 0.4;
	Istim = 0;
	Fcell = 0;
	Zindex = 0;
	phase = 0;
	SAN_Gna = 0;
	SAN_Gto = 0;
	SAN_GcaL = 0;
	SAN_Gks = 0;
	SAN_Gkr = 0;
	SAN_Gf = 0;
	SAN_Gbna = 0;
	SAN_Gbca = 0;
	SAN_periph = 0;
	IKur_type_CNZ = 1;
	Simple_Ikur_ac_shift = 0;
	Simple_Ikur_ac_grad = 1;
	Simple_Ikur_inac_shift = 0;
	Simple_Ikur_inac_grad = 1;
	Simple_GKur = 1;
	IKr_grad = 1;

	drug_INa_concen = 0.0; // uM
	drug_INa_Ka = 0.1 / 1000.0; // ms-1*uM-1
	drug_INa_Ki = 0.1 / 1000.0; // ms-1 * uM -1
	drug_INa_La = 0.1; // ms-1
	drug_INa_Li = 0.01; // ms-1

	drug_IKur_KO = 0;
	drug_IKur_ZKO = 0;
	drug_IKur_KC = 0;
	drug_IKur_ZKC = 0;

	K2p_Vh_shift = 0.0;

	drug_IKur_LO = 0;
	drug_IKur_ZLO = 0;
	drug_IKur_LC = 0;
	drug_IKur_ZLC = 0;

	drug_IKur_concen = 0.0;

	RyRTauScale = 1.0;
	INa_Act_Shift = 0;
	Ito_Act_Shift = 0.0;
	SERCA_Scale = 1.0;

	INa = IKur = Ito = 0.0;
	S_GKv14 = 0.0;
	IK2p_Shift = 0;
	GNaL = 1.0;

	GSK = 1.0;
	GK2P = 1.0;
	GbNa = 1.0;
}

void SingleCellPara::SetINaBlockPara(double drug_con, double Ka, double Ki, double La, double Li) {
	drug_INa_concen = drug_con; // uM
	drug_INa_Ka = Ka; // ms-1*uM-1
	drug_INa_Ki = Ki; // ms-1 * uM -1
	drug_INa_La = La; // ms-1
	drug_INa_Li = Li; // ms-1

}


void SingleCellPara::SetCellTypePara(AtriaCellType cell_type) {
	m_cell_type = cell_type;
	If_vshift = 0;
	If_grad = 1;
	Gf = 0.4;

	if (cell_type == RA) {
		BULK_CONST = 1.25;
	}
	else  if (cell_type == PM) { // PM
		GCaL = 0.94 * GCaL;
	}
	else if (cell_type == CT) { // CT
		GCaL = 1.68;//2;  // note 1.68 = Biophysically based (Feng et al), 2 better replicates APD differences (1.68 is still within expt range)
		Gto = 1.35 * Gto;
	}
	else if (cell_type == RAA || cell_type == AS) { // RAA / AS
		BULK_CONST = 1.5;
		Gto = 0.53;//0.4; // or 0.53
		if (cell_type == AS) { // AS
			GCaL = 0.4;//0.25;
			GKur = 0.667;
			Gto = 0.4 * 0.4;
			IK1_v_shift = -6;
			GNa = 1.3;
		}
		/*if (IKur_type_CNZ == 0){
		    IKur_ac_shift = -14;
		    IKur_inac_shift = 25;
		    IKur_ac_grad = 2.1317;
		    IKur_inac_grad = 1.769;
		}*/
	}
	else if (cell_type == AVR) { // AVR
		GCaL = 0.67;
		GKr = 1.63;
		Gto = 0.6;
	}
	else if (cell_type == BB) { // BB
		GCaL = 2.32;//2.2;//1.6;
		Gto = 1.17;//0.8;
		GK1 = 0.85;
		BULK_CONST = 1.6;//1.22;

	}
	else if (cell_type == LA) { // LA
		Gf = 1;
		GKr = 1.6;

	}
	else if (cell_type == LAA) { // LAA
		Gf = 1;
		GKr = 1.6;
		Gto = 0.53;

		if (IKur_type_CNZ == 0) {
			IKur_ac_shift = -14;
			IKur_inac_shift = 25;
			IKur_ac_grad = 2.1317;
			IKur_inac_grad = 1.769;
		}
		GKur = 0.8;
	}
	/*else if (cell_type == PV) { // PV

	    IKr_ac_shift = -4.2;
	    IKr_ac_grad = 0.938;
	    IKs_shift = -9.6;
	    IKs_grad = 1.047;
	    IK1_v_shift = -15;
	    GK1 = 0.65;

	}
	else if (cell_type == 10) { // PV new (rough ehrlich) -- note wrong on Gkr and Gks
	    GKr = 0.75 * 1.6;
	    Gto = 0.8 * 0.53;
	    GKs = 0.7 * 1.8;
	    IKs_shift = -13.43;
	    IKs_grad = 0.5;//0.75;
	    BULK_CONST = 1.3;
	    GK1 = 0.7;
	    GCaL = 0.65;//0.75;
	    IK1_v_shift = -10;
	}*/
	else if (cell_type == PV) { // PV new Ehrlich + Cha (Primarily Ehrlich) - USE THIS PV MODEL
		GKr = 1.5 * 1.6; //1.6*1.6;
		Gto = 0.75;//0.72*0.53;
		GKs = 1.5;//1.75*1.8;
		GCaL = 0.7;
		GK1 = 0.62;//0.69;
	}
	else if (cell_type == SAN_C) { // SAN C
		GNa = 0.06;
		GKur = 0.26;
		GK1 = 0.42;
		GCaT = 4.5 * 17; //5.2*17;
		Gf = 4.4;
		GNaCa = 0.5;
		Gto = 0.4;
		GCaL = 0.68;
		GKs = 0.69;
		GKr = 0.45;
		BULK_CONST = 3;
		If_grad = 1.1;
	}
	else if (cell_type == SAN_P) { // SAN P

		GNa = 0.3;
		GKur = 0.26;
		GK1 = 0.3;//0.42;
		GCaT = 4.5 * 17;
		Gf = 3.93 * 4.4;
		GNaCa = 0.5;
		Gto = 0.8;//0.92;
		GCaL = 2.28;
		GKs = 2;//4.19;
		GKr = 2.75;
		BULK_CONST = 3;
		If_grad = 1.075;
	}

} // end cell_type Params


void SingleCellPara::SetIKurMutationPara(IKurMuationType mutation_type) {
	m_mutation_type = mutation_type;

	if (mutation_type == WT) { /*WT*/
		IKur_ac1_mult   = 1.0;
		IKur_ac1_shift  = 0.0;
		IKur_ac1_grad   = 1.0;
		IKur_ac2_mult   = 1.0;
		IKur_ac2_shift  = 0.0;
		IKur_ac2_grad   = 1.0;
		IKur_ac_add     = 0.0;
		IKur_inac_mult  = 1.0;
		IKur_inac_shift = 0.0;
		IKur_inac_grad  = 1.0;
		IKur_inac_add   = 0.0;
		GKur *= 1.0;
	}
	else if (mutation_type == D322H) {  /*D322H*/
		IKur_ac1_mult   = 1.0;
		IKur_ac1_shift  = -9.3932825747999988;
		IKur_ac1_grad   = 1.3785501205210822;
		IKur_ac2_mult   = 1.0;
		IKur_ac2_shift  = 15.229362644299998;
		IKur_ac2_grad   = 0.25822502847411422;
		IKur_ac_add     = 0.0;
		// Inactivation
		IKur_inac_mult  = 0.87254304449870279;
		IKur_inac_shift = -9.6150807654200001;
		IKur_inac_grad  = 0.80075124328979652;
		IKur_inac_add   = 0.073931087206000057;

		GKur *= 1.720232;


		Simple_Ikur_ac_shift += -3.2641;
		Simple_Ikur_ac_grad *= 0.8;
		Simple_Ikur_inac_shift += 9.61;
		Simple_Ikur_inac_grad *= 0.8;
		Simple_GKur *= 1.7995;

	}
	else if (mutation_type == E48G) { /*E48G*/
		IKur_ac1_mult   = 1.0;
		IKur_ac1_shift  =  -8.914978286109998;
		IKur_ac1_grad   = 1.687979145508135 ;
		IKur_ac2_mult   = 1.0;
		IKur_ac2_shift  = 14.395569738899999;
		IKur_ac2_grad   = 0.312429308790384;
		IKur_ac_add     = 0.0;
		// // Inactivation
		IKur_inac_mult  = 0.935206196986113;
		IKur_inac_shift = -4.02829;
		IKur_inac_grad  = 0.940280253452448;
		IKur_inac_add   = 0.039197146604;

		GKur *= 1.286906;


		Simple_Ikur_ac_shift += -3.02;
		Simple_Ikur_ac_grad *= 0.95;
		Simple_Ikur_inac_shift += 4.02;
		Simple_Ikur_inac_grad *= 0.94;
		Simple_GKur *= 1.32;



	}
	else if (mutation_type == A305T) {  /*A305T*/
		IKur_ac1_mult   = 1.0;
		IKur_ac1_shift  = -5.135031219399998 ;
		IKur_ac1_grad   = 2.40452882662690 ;
		IKur_ac2_mult   = 1.0;
		IKur_ac2_shift  = 13.643107896399998;
		IKur_ac2_grad   = 0.387170393885174;
		IKur_ac_add     = 0.0;
		// Inactivation
		IKur_inac_mult  = 0.874602534570044;
		IKur_inac_shift = -6.032311251379999;
		IKur_inac_grad  = 0.748490912918043;
		IKur_inac_add   = 0.070656047487;

		GKur *= 1.388228;

		Simple_Ikur_ac_shift += -4.068;
		Simple_Ikur_ac_grad *= 1.102;
		Simple_Ikur_inac_shift += 6.03;
		Simple_Ikur_inac_grad *= 0.74;
		Simple_GKur *= 1.445;

	}
	else if (mutation_type == Y155C) {   /*Y155C*/
		IKur_ac1_mult   = 1.0;
		IKur_ac1_shift  = -4.430626802899999 ;
		IKur_ac1_grad   = 1.651529981648123 ;
		IKur_ac2_mult   = 1.0;
		IKur_ac2_shift  = 3.3770968185;
		IKur_ac2_grad   =  0.517356087786261;
		IKur_ac_add     = 0.0;
		// Inactivation
		IKur_inac_mult  = 0.936476668495346;
		IKur_inac_shift = -5.012779270500001;
		IKur_inac_grad  = 0.781775743434569;
		IKur_inac_add   = 0.034583585385;

		GKur *= 0.461283;
		Simple_Ikur_ac_shift += -.89;
		Simple_Ikur_ac_grad *= 0.75;
		Simple_Ikur_inac_shift += 5.01;
		Simple_Ikur_inac_grad *= 0.78;
		Simple_GKur *= 0.4745;


	}
	else if (mutation_type == D469E) {  /*D469E*/

		IKur_ac1_mult   = 1.0;
		IKur_ac1_shift  = 0.0 ;
		IKur_ac1_grad   = 1.0;
		IKur_ac2_mult   = 1.0;
		IKur_ac2_shift  = 0.0;
		IKur_ac2_grad   =  1.0;
		IKur_ac_add     = 0.0;
		// Inactivation
		IKur_inac_mult  =  1.174004519012284;
		IKur_inac_shift = -4.5454272743;
		IKur_inac_grad  = 0.842432589996777;
		IKur_inac_add   = -0.09291551303;

		GKur *= 0.563448;

		Simple_Ikur_inac_shift += 4.5455;
		Simple_Ikur_inac_grad *= 0.84;
		Simple_GKur *= 0.5457;

	}
	else if (mutation_type == P488S) {
		/*P488S*/   /*note, P488s and D469E used the same activation if as the WT, and
		                                                P488S used the inactivation of P488S/WT, due to the absense of
		                                                P488S inactvation data*/

		IKur_ac1_mult   = 1.0;
		IKur_ac1_shift  = 0.0 ;
		IKur_ac1_grad   = 1.0;
		IKur_ac2_mult   = 1.0;
		IKur_ac2_shift  = 0.0;
		IKur_ac2_grad   =  1.0;
		IKur_ac_add     = 0.0;
		// Inactivation
		IKur_inac_mult  =  1.212942827258507;
		IKur_inac_shift = 1.4151025432;
		IKur_inac_grad  = 0.687940393974062;
		IKur_inac_add   = -0.110508855073;

		GKur *= 0.038593;


		Simple_Ikur_ac_grad *= 0.83;
		Simple_Ikur_inac_shift += -1.4251;
		Simple_Ikur_inac_grad *= 0.68;
		Simple_GKur *= 0.0384;

	} else {
		std::cerr << "wrong IKur mutation type \n";
		std::exit(0);
	}
} // end mutation_type function


void SingleCellPara::SetAFTypePara(AFType AF_model) {
	m_AF_model = AF_model;
	if (AF_model == NONE) { // SR condition, do not change anything
		GCaL = GCaL;
		Gto = Gto;
		GK1 = GK1;

	} else if (AF_model == AF1) { // Bosch based model  // based on Henggui's paper. in cardiovascular research 2005.
		GCaL = 0.26 * GCaL;
		Gto = 0.3 * Gto;
		// GK1 = 3.35 * GK1;
		GK1 = 2.06 * GK1;
		Ito_Act_Shift = 16;
		INa_Act_Shift = 1.6;
	}
	else if (AF_model == AF2) { // Workman based model
		GCaL = 0.35 * GCaL;
		Gto = 0.35 * Gto;
		GK1 = 1.75 * GK1;
	}
	else if (AF_model == AF3) { // Wagoner based model
		GCaL = 0.37 * GCaL;
		Gto = 0.34 * Gto;
		GKur = 0.51 * GKur;
		GK1 = 2.06 * GK1;
	}
	else if (AF_model == AF4) { // Generic
		// GCaL = 0.3 * GCaL;
		GCaL = 0.35 * GCaL;  //
		Gto = 0.35 * Gto;
		GK1 = 2.1 * GK1;

		// GNaCa = 1.55 * GNaCa;
		GNaCa = 1.40 * GNaCa;
		BULK_CONST = 1.5; // Increase according to shanmugam et al, AF_SR_updatke
		// GSR_leak = 1.25;
		GSR_leak = 1.50; // Voigt et al. 2012
		// GSR_leak = 1.50; // Voigt et al. 2012
		fIRel = 3 * fIRel;
		GKs = 2 * GKs;
		RyR = 2.5;
		fRyR = 2.5; // sensitivity to Cai/Cass, // 14:49:55, Wed, 03-January-2018, By Haibo
		RyRTauScale = 2.7;//3;
		SERCA_Scale = 0.6;//0.5;  // voigt 2012 online figure III
		Ito_Act_Shift = 16;
		INa_Act_Shift = 1.6;

		IK2p_Shift = -15.0;  // Vhalf more negative, see Schmidt et al. 2014 model

		GNaL *= 1.5; // double INaL with cAF remodelling. // 15:19:07, Sat, 11-November-2017, By Haibo
		// GK2P *= 3.0;   // implemented in the main function now. // 17:35:46, Wed, 03-January-2018, By Haibo
		// GSK *= 2.0;
		// GKur = 0.5 * GKur;
	}
	else {
		std::cerr << "wrong AF model type!!! \n";
		std::exit(0);
	}
}



void SingleCellPara::SetAcacetinEffectPara(double Acacetin_Con) {

	if (Acacetin_Con <= 1.0e-14)
	{
		Acacetin_Con = 1.0e-14;   // just want to make sure that this data be positive;
	}

	double Aca_IC50_IKur = 3.2; // um
	double Aca_IC50_Ito = 9.3; // um
	double Aca_IC50_IKr = 32.4; // um
	double Aca_IC50_IKs = 81.4; // um
	double Aca_hill_IKur = 0.8; // um
	double Aca_hill_Ito = 0.9; // um
	double Aca_hill_IKr = 0.9; // um
	double Aca_hill_IKs = 0.8; // um
	double IKur_Acacetin_B =  1.0 - 1.0 / (1 + pow((Aca_IC50_IKur / Acacetin_Con), Aca_hill_IKur));
	double Ito_Acacetin_B =  1.0 - 1.0 / (1 + pow((Aca_IC50_Ito / Acacetin_Con), Aca_hill_Ito));
	double IKr_Acacetin_B =  1.0 - 1.0 / (1 + pow((Aca_IC50_IKr / Acacetin_Con), Aca_hill_IKr));
	double IKs_Acacetin_B =  1.0 - 1.0 / (1 + pow((Aca_IC50_IKs / Acacetin_Con), Aca_hill_IKs));
	// GKur = IKur_Acacetin_B * GKur;


	drug_IKur_KO  =  0.00019445959;
	drug_IKur_LO  =  0.00029099869;
	drug_IKur_ZKO = -0.25731953173;
	drug_IKur_ZLO = -0.01280158014;
	drug_IKur_KC  =  0.00247871659;
	drug_IKur_LC  =  0.00028636715;
	drug_IKur_ZKC = -0.32702357338;
	drug_IKur_ZLC =  0.83707302269;

	drug_IKur_concen = Acacetin_Con;

	Gto  = Ito_Acacetin_B * Gto;
	GKr  = IKr_Acacetin_B * GKr;
	GKs  = IKs_Acacetin_B * GKs;
}


void SingleCellPara::SetISOParameters(double ISO) {

	double    fICaL = 1 / ( 1 + exp((-4.362 * (log(ISO) + 3.27)) ) );
	double    fRyR = 1 / ( 1 + exp((-2.55 * (log(ISO) + 1.2)) ) );
	double    fIKs = 1 / ( 1 + exp((-4 * (log(ISO) + 1.4)) ) );
	double    fIKur = 1 / ( 1 + exp((-1.8 * (log(ISO) + 2.2)) ) );
	double    fINa = 1 / ( 1 + exp((-2.55 * (log(ISO) + 2.5)) ) );
	double    fPLB = 1 / ( 1 + exp((-3 * (log(ISO) + 2.2)) ) );

	GCaL = (1 + (fICaL * 2)) * GCaL;  // 3 fold increase
	GKur = (1 + (fIKur * 0.6)) * GKur; // 1.6 fold
	GKs = (1 + (fIKs * 1.5)) * GKs;     // 2.5 fold
	BULK_CONST = (BULK_CONST + (fPLB * 0.1)); // plus a bit 0.1 ,,,,
	GNa = (1 + (fINa * 0.25)) * GNa;        // 1.25 fold
}
