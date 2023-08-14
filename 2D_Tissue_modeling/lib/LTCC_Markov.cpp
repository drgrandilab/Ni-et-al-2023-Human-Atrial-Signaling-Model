// #include "HAM_Cell.hpp"
#include <cmath>

#include "LTCC_Markov.hpp"

LTCC_Markov::LTCC_Markov() {
	for (int i = 0; i < 10; ++i)
	{
		state[i] = 0;
		ydot[i] = 0;
	}

	k1 = k2 = k3 = k4 = k5 = k6 = k7 = k8 = k9 = k10 = k11 = k12 = k1p = k2p = k3p = k4p = k5p = k6p = k7p = k8p = k9p = k10p = k11p = r1 = r2 = s1 = s2 = s1p = s2p = 0;
}




void LTCC_Markov::print_to_file(double t, std::ofstream & output_file) {

	output_file << t << " "  // 1
	            << I_Ca_junc_m1 << " "
	            << I_Na_junc_m1 << " "
	            << I_K_junc_m1 << " "
	            << state[0] << " "  // 5
	            << state[1] << " "  // 6
	            << state[2] << " "  // 7
	            << state[3] << " "  // 8
	            << state[4] << " "  // 9
	            << state[5] << " "  // 10
	            << state[6] << " "  // 11
	            << state[7] << " "  // 12
	            << 1 - (state[0] + state[1] + state[2] + state[3] + state[4] + state[5] + state[6] + state[7] + state[8 ] + state[9]  ) << " " // 13
	            << k1 << " "
	            << k2 << " " // 15
	            << k3 << " "
	            << k4 << " "
	            << k5 << " "
	            << k6 << " "
	            << k7 << " "  // 20
	            << k8 << " "
	            << k9 << " "
	            << k10 << " "
	            << k11 << " "
	            << k12 << " "  //25
	            << k1p << " "
	            << k2p << " "
	            << k3p << " "
	            << k4p << " "
	            << k5p << " " // 30
	            << k6p << " "
	            << k7p << " "
	            << k8p << " "
	            << k9p << " "
	            << k10p << " " //35
	            << k11p << " "
	            << r1 << " "
	            << r2 << " "
	            << s1 << " "
	            << s2 << " " //40
	            << s1p << " "
	            << s2p << " "



	            << std::endl;
}