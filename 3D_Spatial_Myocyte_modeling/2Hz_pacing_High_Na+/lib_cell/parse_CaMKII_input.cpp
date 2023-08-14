


#include "parse_CaMKII_input.hpp"

void CaMKII_para::print_to_file(double t, std::string file) {

	std::ofstream output_file(file);
	output_file <<  std::setprecision(9)
	            << t << " "  << std::setprecision(3)
	            << LCC_CKdyadp << " "
	            << RyR_CKp << " "
	            << PLB_CKp << " "
	            << NaV_CKp << " "
	            << Ikur_CKp << " "
	            << Itof_CKp << " "
	            << fckiim2_j         << " "
	            << fCKII_RyR         << " "
	            << fCKII_PLB         << " "
	            << kCKII_Nav_change  << " "
	            << kCKII_Ikur        << " "
	            << kCKII_Itof_vshift << " "
	            << kCKII_Itof_tau    << " "
	            << kCKII_IK1_G       << " "
	            << kCKII_Gapjunct_G  << " "
	            << fPKA_RyR          << " "
	            << RyR_koSRCa_Scale  << " "
	            << PLB_kmf_Scale
	            << std::endl;
	output_file.close();
}


void CaMKII_para::read_CaMKII_levels_from_txt(std::string filename) {
	double * data_vec = new double [6] ;
	read_num_data_file(filename, data_vec, 6);
	LCC_CKdyadp = data_vec[0];
	RyR_CKp = data_vec[1];
	PLB_CKp = data_vec[2];
	NaV_CKp = data_vec[3];
	Ikur_CKp = data_vec[4];
	Itof_CKp = data_vec[5];
	delete [] data_vec;
}


void CaMKII_para::set_CaMKII_signaling_effects() {

	double Gapjunct_CKp = NaV_CKp; // dct['NaV_CKp']
	double IK1_CKp = NaV_CKp; // dct['NaV_CKp']
	double Nav_CKp_wt = 2.5654 / 30; // # 8.5# 1X, 1 Hz
	double Nav_CKp_oe = 23.9069 / 30; // # 80# 6X, 1 Hz
	double Ikur_CKp_ko = 0;
	double Ikur_CKp_wt = 2.5654 / 30;
	double Itof_CKp_wt = 0.11;// #//2.5654 / 30.0;



	fckiim2_j         = LCC_CKdyadp * 0.1;
	fCKII_RyR         = 1 + 3 * (RyR_CKp - 0.21);
	fCKII_PLB         = (1 - 0.5 * PLB_CKp);
	kCKII_Nav_change  = (NaV_CKp - Nav_CKp_wt) / (Nav_CKp_oe - Nav_CKp_wt); // # 0 WT, 1 OE
	kCKII_Ikur        = 0.5 + 1.0 / ( 1 + exp( (Ikur_CKp - Ikur_CKp_wt) / -0.08 ));
	kCKII_Itof_vshift = -5 + 10 / ( 1 + exp( (Itof_CKp - Itof_CKp_wt) / -0.1 ));
	kCKII_Itof_tau    = 0.5 + 1.0 / ( 1 + exp( (Itof_CKp - Itof_CKp_wt) / -0.1 )) ;
	kCKII_IK1_G       = 0.5 + 1.0 / ( 1 + exp( (IK1_CKp  - 0.1) / -0.15 ));
	kCKII_Gapjunct_G  = 1.35 - (0.9 / (1 + exp( -(Gapjunct_CKp - 0.16) / 0.12 ) ) );
	// fPKA_RyR          = 1.0;
	// RyR_koSRCa_Scale  = fCKII_RyR + fPKA_RyR - 1.0;
	PLB_kmf_Scale     = fCKII_PLB;
}


void CaMKII_para::set_CaMKII_Inhibition() {
}


void CaMKII_para::set_ISO() {

	G_LTCC = 1.41;
	fpkam2 = 0.041047270323433624;
	fPKA_RyR = 1.4;
	fPKA_PLB = 0.6984523233691271;

	kPKA_Ina_increase = 0.9259086130727806;

	KmNaip_PKA = 2.4192638096159436;
	fPKA_TnI = 1.3000000000000003;
	kPKA_Iks = 0.9293274644480991;
	kPKA_Ikr = 0.9259086130727806;
	Vkr_PKA = 9.259086130727807;
	dGkr_PKA = 0.27777258392183417;
	kPKA_Ikur = 0.9259086130727806;
	kPKA_IClCa = 1.0199999999999998;
	kPKA_Myo = 0.7963056012382703;
	kPKA_Itof = 0.9259086130727806;
	kPKA_IK1 = 0.9259086130727806;


}