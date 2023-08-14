#ifndef RyR_HPP
#define RyR_HPP
#include "xor_rand.hpp"

#include <random>



class RyR
{
public:
	RyR(int RYR_Num = 41, int ID = 0, double Jmax_set = 0.0147 * 18 * 0.5); // 02/08/2021
	~RyR();

	int RyR_1, RyR_2, RyR_3;
	xor_rand random_num_RyR;

	double Po;
	int N_RyR;
	double Jmax;

    double p12;
    double p14;
    double p21;
    double p23;
    double p43;
    double p41;
    double p34;
    double p32;
 
    double NRyR_1;
    double NRyR_2;
    double NRyR_3;
    double NRyR_4;
    int tmp1;
    int tmp2;

	double Update_RyR_stochastic(double dt, double Caj, double CaSRj);
	int Update_RyR_rates_stochastic(double num, double p);
	int Update_RyR_rates_stochastic_BINV(int num, double p);
	int Update_RyR_rates_stochastic_cpp11(int num, double p);
	int Update_RyR_rates_stochastic_python(int n, double p);
	void Update_RyR_rates_stochastic_MN(int n, double p1, double p2);
	int random_binomial(int n, double p);

	void set_Jmax(double Jmax_set) {
		Jmax = Jmax_set;
	}

	double static const MaxSR = 15;
	double static const MinSR = 1;
	double static const ec50SR = 450;
	double static const hkosrca = 2.5;

	double static const Ku = 5.0;
	double static const Kb = 0.005;
	double static const tauu = 125.0;
	double static const taub = 0.5;
	double static const tauc1 = 2.0;
	double static const tauc2 = 0.3;

	double static const taup = 0.022;
	double static const tautr = 5.0;
	double static const BCSQN0 = 400;
        double static const Kcp = 10.0; 

 	double static const pedk12 = 0.000001;
	double static const pedk43 = 0.000001;
	
	double static const BCSQN = 400;
	double static const Kc = 650.0;

	double static const nM = 15;
	double static const nD = 35;
	double static const rhoinf = 5000;
	double static const hh = 23.0; 
	double static const KK = 900; 
	std::mt19937 gen;
};



#endif
