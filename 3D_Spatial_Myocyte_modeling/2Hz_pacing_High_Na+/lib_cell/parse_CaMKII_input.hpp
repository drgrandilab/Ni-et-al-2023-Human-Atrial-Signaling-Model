
#ifndef PARSE_CAMKII_INPUT__HPP
#define PARSE_CAMKII_INPUT__HPP
#include <iostream>
#include "input_output.h"
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cstdlib>

#include <fstream>
#include <cmath>
#include <cstdlib>
#include  <iomanip>
#include <iostream>

class CaMKII_para
{
public:
	CaMKII_para() {
		G_LTCC            = 1;//0.997;
		fpkam2            = 0.0;
		fPKA_RyR          = 1;//0.997;
		fPKA_PLB          = 1.0;
		kPKA_Ina_increase = 0;// -0.017760065783109222;
		KmNaip_PKA        = 0;// -0.003097165199445358;
		fPKA_TnI          = 1.0;
		kPKA_Iks          = 0;// -0.018192831575892424;
		kPKA_Ikr          = 0;// -0.017760065783109222;
		Vkr_PKA           = 0;// -0.1776006578310922;
		kPKA_Ikur         = 0;// -0.017760065783109222;
		kPKA_IClCa        = 0;// 0.14709854567734026;
		kPKA_Myo          = 0;// -0.0004982281126052046;
		kPKA_Itof         = 0;// -0.017760065783109222;
		kPKA_IK1          = 0;// -0.017760065783109222;
	};
	~CaMKII_para() {};

	double LCC_CKdyadp = 0.0001;
	double RyR_CKp = 0.0001;
	double PLB_CKp = 0.0001;
	double NaV_CKp = 0.0001;
	double Ikur_CKp = 0.0001;
	double Itof_CKp = 0.0001;
	double Gapjunct_CKp = 0.0001;

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
	double fCKII_RyR = 1;
	double fPKA_RyR = 1;
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
	double kPKA_Iks = 0;
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
	double kPKA_IK1 = 0;
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
	void print_to_file(double t, std::string file);
	void read_CaMKII_levels_from_txt(std::string filename);

	void set_CaMKII_signaling_effects();

	void set_ISO();
};

#endif