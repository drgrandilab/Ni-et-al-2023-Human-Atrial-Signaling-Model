#ifndef LTCC_Markov_HPP
#define LTCC_Markov_HPP

// #include "HAM_Cell.hpp"
#include <fstream>

#include "HAM_Constants.h"



class LTCC_Markov
{
public:
	LTCC_Markov();
	~LTCC_Markov() {};
	double state[10];  // seven states for model 1
	double ydot[10];  // seven states for model 1
	// LTCC Current - Fixed Parameters
	double Fjunc_CaL = 1.0;
	double Fsl_CaL = 1 - Fjunc_CaL;
	double I_Ca_junc_m1;
	double I_Na_junc_m1;
	double I_K_junc_m1;
	double k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12;
	double k1p, k2p, k3p, k4p, k5p, k6p, k7p, k8p, k9p, k10p, k11p;
	double r1, r2, s1, s2, s1p, s2p;
	double k13, k14;
	double Pna = 1.05* 0.65 * (1 /*- 0.5 * cAF*/) * (0.70 - 0.20) * 100 * 0.675e-9; // [cm/s]
	double Pk = 1.05* 0.65 * (1 /*- 0.5 * cAF*/) * (0.70 - 0.20) * 100 * 12.15e-9; // [cm/s]
	double Pca = 1.05* 0.65 * (1 /*- 0.5 * cAF*/) * (0.70 - 0.20) * 100 * 24.3e-6; // [cm/s] 10*0.45*5.4e-6

	double update_states(double Caj, double V, double najLCC, double kjLCC);
	double update_states_v2(double Caj, double V, double najLCC, double kjLCC);
	double update_states_v3(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v4(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v4_plus(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v4_plus_mode_2(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v4_plus_new(double Caj, double Vm, double najLCC, double kjLCC);

	double update_states_v4_plus_new_mode_2(double Caj, double Vm, double najLCC, double kjLCC);

	double update_states_v4_plus2(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v4_pluscopy(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v3_plus(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v3_plus_mode_2(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v1(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v4_plus2_mode_2(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v5(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v5_mode_2(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v6(double Caj, double Vm, double najLCC, double kjLCC);

	double update_states_v6_mode_2(double Caj, double Vm, double najLCC, double kjLCC);
	double update_states_v7(double Caj, double Vm, double najLCC, double kjLCC,double ISO_PKA=0);

	double update_states_v7_mode_2(double Caj, double Vm, double najLCC, double kjLCC,double ISO_PKA=0);
	
	double update_states_v7_BPS_2020(double Caj, double Vm, double najLCC, double kjLCC,double ISO_PKA=0);
	double update_states_v7_BPS_2020_mode_2(double Caj, double Vm, double najLCC, double kjLCC,double ISO_PKA=0);

	double update_states_v8(double Caj, double Vm, double najLCC, double kjLCC,double ISO_PKA=0);

	double update_states_v8_mode_2(double Caj, double Vm, double najLCC, double kjLCC,double ISO_PKA=0);

	double update_states_v9(double Caj, double Vm, double najLCC, double kjLCC,double ISO_PKA=0);

	double update_states_v9_mode_2(double Caj, double Vm, double najLCC, double kjLCC,double ISO_PKA=0);

	void print_to_file(double t, std::ofstream & output_file);
};


#endif
