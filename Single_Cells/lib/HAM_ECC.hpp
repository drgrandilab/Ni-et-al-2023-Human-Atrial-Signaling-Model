#ifndef HAM_ECC_HPP
#define HAM_ECC_HPP
#include  <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "HAM_Constants.h"
#include "LTCC_Markov.hpp"
#include "signalling_para.hpp"
#include "Herg.hpp"

class HAM_ECC
{
public:
	HAM_ECC() {
		allow_stimulus_app = true;
		LTCC_mode_1_out = 0.0;
		LTCC_mode_2_out = 0.0;
	};
	~HAM_ECC() {};

	double V, dV;
	double I_Na_junc, I_Na_sl, I_Na;

	double I_NaL_junc, I_NaL_sl, I_NaL;
	double I_nabk_junc, I_nabk_sl, I_nabk;
	double I_nak_junc, I_nak_sl, I_nak;
	double I_to, I_tof;
	double I_kr, I_kur;
	double I_k2p;
	double I_ks_sl, I_ks_junc, I_ks;
	double I_ki_sl, I_ki_j, I_ki;
	double I_kp;
	double I_kach_j, I_kach_sl, I_kach;

	double I_sk_junc, I_sk_sl, I_sk, /*I_ClCa_jun, I_ClCa_s,*/I_ClCa, I_Clbk, I_ClCFTR;
	double I_Ca_junc, I_Ca_sl, I_Ca, I_CaK, I_CaNa_junc, I_CaNa_sl, I_CaNa, I_Catot;
	double I_cabk_junc, I_cabk_sl, I_cabk;
	double I_pca_junc, I_pca_sl, I_pca;
	double I_ncx_junc, I_ncx_sl, I_ncx, J_SRCarel, J_serca, J_SRleak,  I_Na_tot_junc, J_CaB_junction, J_CaB_sl;
	double J_CaB_cytosol,  I_Na_tot_sl;
	double I_ClCa_junc, I_ClCa_sl, I_Cl_tot;
	double I_app, I_Ca_tot, I_Ca_tot_junc, I_Ca_tot_sl, I_tot, I_Na_tot, I_K_tot;
	double CaSR, Caj, Casl, Cai;
	double Nai, Naj, Nasl;

	bool allow_stimulus_app;

	LTCC_Markov LTCC_m2, LTCC;
	LTCC_Markov LTCC_sl_m2, LTCC_sl;

	double LTCC_mode_1_out, LTCC_mode_2_out;

	const int ODE_NUM = 85 + 7;
	Herg IKr_markov;
	double *y, *ydot;
	// double ydot[65];

	void HAM_ECC_update_ODE(double t, signalling_para & para);
	void print_to_file(double t, std::ofstream & output_file) ;
	void print_to_file_Vm_only(double t, std::ofstream & output_file);
	void print_to_file_Vm_Ca(double t, std::ofstream & output_file);

	void initialiser();

};






// int human_atrial_ECC_ODE(double t, double *y, double *ydot, Cell * cell_par);


#endif