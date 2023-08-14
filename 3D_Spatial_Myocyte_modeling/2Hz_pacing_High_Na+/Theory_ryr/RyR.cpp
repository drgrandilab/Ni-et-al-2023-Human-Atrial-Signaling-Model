#include "RyR.hpp"
#include <iostream>

RyR::RyR(int RYR_Num, int ID, double Jmax_set)
	: random_num_RyR(0, ID), Jmax(Jmax_set), gen(ID)
{

	N_RyR = RYR_Num;
	//RyR_1 = 0 + int (5.0 * random_num_RyR.gen_rand_uint() / (double)(UINT_MAX));
	RyR_1 = 0 + int (5.0 * random_num_RyR.gen_rand());
	RyR_2 = 0;
	RyR_3 = 0;

	tmp1 = 0;
	tmp2 = 0;
}

RyR::~RyR() {}

double RyR::Update_RyR_stochastic(double dt, double Caj, double CaSRj) {
	double rho = rhoinf * 1.0 / (1.0 + pow(KK / CaSRj, hh));
	if (rho < 0.0000001)rho = 0.0000001;
	double MM = (sqrt(1 + 8 * rho * BCSQN) - 1) / (4 * rho * BCSQN);
	
	double cp2 = Caj * Caj;
	double sgmd = cp2 * Caj / (Kcp * Kcp * Kcp + cp2 * Caj);
	double k12 = Ku * sgmd + pedk12;
	double k43 = Kb * sgmd + pedk43;

	double k14 = MM / taub * BCSQN / BCSQN0;
	double k21 = 1 / tauc1;
	double k23 = MM / taub * BCSQN / BCSQN0;
	double k41 = 1 / tauu;
	double k34 = 1 / tauc2;
	double k32 = k41 * k12 * k23 * k34 / k43 / k21 / k14;

	int RyR_4 = N_RyR - RyR_1 - RyR_2 - RyR_3;

	Update_RyR_rates_stochastic_MN(RyR_1, k12 * dt, k14 * dt);
	int ryr12 = tmp1;
	int ryr14 = tmp2;
	Update_RyR_rates_stochastic_MN(RyR_2, k23 * dt, k21 * dt);	
	int ryr23 = tmp1;
	int ryr21 = tmp2;
	Update_RyR_rates_stochastic_MN(RyR_3, k34 * dt, k32 * dt);	
	int ryr34 = tmp1;
	int ryr32 = tmp2;	
	Update_RyR_rates_stochastic_MN(RyR_4, k41 * dt, k43 * dt);
	int ryr41 = tmp1;
	int ryr43 = tmp2;
	// int ryr12 = Update_RyR_rates_stochastic_MN(RyR_1, k12 * dt);
	// int ryr14 = Update_RyR_rates_stochastic_MN(RyR_1, k14 * dt);
	// int ryr21 = Update_RyR_rates_stochastic_MN(RyR_2, k21 * dt);
	// int ryr23 = Update_RyR_rates_stochastic_MN(RyR_2, k23 * dt);
	// int ryr43 = Update_RyR_rates_stochastic_MN(RyR_4, k43 * dt);
	// int ryr41 = Update_RyR_rates_stochastic_MN(RyR_4, k41 * dt);
	// int ryr34 = Update_RyR_rates_stochastic_MN(RyR_3, k34 * dt);
	// int ryr32 = Update_RyR_rates_stochastic_MN(RyR_3, k32 * dt);

	p12 = k12 * dt;
	p14 = k14 * dt;
	p21 = k21 * dt;
	p23 = k23 * dt;
	p43 = k43 * dt;
	p41 = k41 * dt;
	p34 = k34 * dt;
	p32 = k32 * dt;

	RyR_1 = RyR_1 - (ryr12 + ryr14) + (ryr21 + ryr41);
	RyR_2 = RyR_2 - (ryr21 + ryr23) + (ryr12 + ryr32);
	RyR_3 = RyR_3 - (ryr34 + ryr32) + (ryr43 + ryr23);
	RyR_4 = N_RyR - RyR_1 - RyR_2 - RyR_3;

	if (RyR_1 < 0 )
	{
		std::cout << RyR_1 << "\t" << RyR_2 << "\t" << RyR_3 << "\t" << RyR_4 << std::endl;
		std::cout << ryr12 << "\t" << ryr14 << "\t" << ryr21 << "\t" << ryr41 << std::endl;
		std::exit(0);
		/* code */
	}

	// if (RyR_1 < 0 || RyR_2 < 0 || RyR_3 < 0 || RyR_1 + RyR_2 + RyR_3 > N_RyR)
	// {
	// 	if (RyR_1 < 0)
	// 	{
	// 		if (random_num_RyR.gen_rand_uint() % 2)
	// 			RyR_2 += RyR_1;
	// 		RyR_1 = 0;
	// 	}
	// 	if (RyR_2 < 0)
	// 	{
	// 		if (random_num_RyR.gen_rand_uint() % 2)
	// 		{
	// 			RyR_1 += RyR_2;
	// 			if (RyR_1 < 0)RyR_1 = 0;
	// 		}
	// 		else
	// 			RyR_3 += RyR_2;
	// 		RyR_2 = 0;
	// 	}
	// 	if (RyR_3 < 0)
	// 	{
	// 		if (random_num_RyR.gen_rand_uint() % 2)
	// 		{
	// 			RyR_2 += RyR_3;
	// 			if (RyR_2 < 0)RyR_2 = 0;
	// 		}
	// 		RyR_3 = 0;
	// 	}
	// 	if (RyR_1 + RyR_2 + RyR_3 > N_RyR)
	// 	{
	// 		RyR_4 = N_RyR - (RyR_1 + RyR_2 + RyR_3);
	// 		if (random_num_RyR.gen_rand_uint() % 2)
	// 		{
	// 			RyR_3 += RyR_4;
	// 			if (RyR_3 < 0)
	// 			{
	// 				RyR_3 -= RyR_4;
	// 				RyR_1 += RyR_4;
	// 				if (RyR_1 < 0)
	// 				{
	// 					RyR_1 -= RyR_4;
	// 					RyR_2 += RyR_4;
	// 				}
	// 			}
	// 		}
	// 		else
	// 		{
	// 			RyR_1 += RyR_4;
	// 			if (RyR_1 < 0)
	// 			{
	// 				RyR_1 -= RyR_4;
	// 				RyR_3 += RyR_4;
	// 				if (RyR_3 < 0)
	// 				{
	// 					RyR_3 -= RyR_4;
	// 					RyR_2 += RyR_4;
	// 				}
	// 			}
	// 		}
	// 	}
	// }

	NRyR_1 = RyR_1;
  NRyR_2 = RyR_2;
  NRyR_3 = RyR_3;
  NRyR_4 = N_RyR - RyR_1 - RyR_2 - RyR_3;

	Po = (RyR_2 + RyR_3) / 100.0;
	return Po;
}


void RyR::Update_RyR_rates_stochastic_MN(int n, double p1, double p2)
{
	if (p1+p2 >= 1){
		std::cout << "Need smaller dt Update_RyR_rates_stochastic_MN()!!!! "<< std::endl;
		std::exit(0);
	}
	if (n < 0){
		std::cout << "Negative trials Update_RyR_rates_stochastic_MN()!!!! "<< std::endl;
		std::exit(0);
	}
	if (n == 0){
		tmp1 = 0;
		tmp2 = 0;
	}else{
	  double remaining_p = 1.0;
	  double p [2] = {p1, p2};
	  int j;
	  int dn = n;
	  int d  = 3;
	  int X [3] = {0, 0, 0};
	  for (j = 0; j < (d - 1); j++) {
	    X[j] = random_binomial(dn, p[j] / remaining_p);
	    dn = dn - X[j];
	    if (dn  <= 0) {
	      break;
	    }
	    remaining_p -= p[j];
	  }
	  if (dn > 0) {
	      X[d - 1] = dn;
	  }
	  tmp1 = X[0];
	  tmp2 = X[1];
	}
}

int RyR::random_binomial(int n, double p) {
  double q;

  if ((n == 0) || (p == 0.0))
    return 0;

  if (p <= 0.5) {
    return Update_RyR_rates_stochastic_python(n, p);
  }
  else {
    q = 1.0 - p;
    return n - Update_RyR_rates_stochastic_python(n, q);
  }  
}


// this function simulates the transition of RyR in the model/ stochastically .
int RyR::Update_RyR_rates_stochastic(double num, double p)
{
	double res;
	double lambda = num * p;

	if (lambda > 12)
	{
		//Gaussian
		double x1, x2, w;
		do
		{
			x1 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			x2 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);
		w = sqrt((-2.0 * log(w)) / w);
		double y1 = x1 * w;
		//double y2=x2*w;
		res = y1 * sqrt(num * p * (1 - p)) + num * p; // *** ave=num*p , rho^2=num*p*(1-p)
		res = int(res + 0.5); //round
	}
	else if (100 * p < 6.6 + 52 * pow(num, double(-0.5)))
	{
		//Poisson
		double L = exp(-lambda);
		double k = 0;
		double pp = 1;
		do
		{
			k++;
			double u = random_num_RyR.gen_rand();
			pp *= u;
		} while (pp >= L);
		res = k - 1;
	}
	else
	{
		//Gaussian
		double x1, x2, w;
		do
		{
			x1 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			x2 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);
		w = sqrt((-2.0 * log(w)) / w);
		double y1 = x1 * w;
		//double y2=x2*w;
		res = y1 * sqrt(num * p * (1 - p)) + num * p; // *** ave=num*p , rho^2=num*p*(1-p)
		res = int(res + 0.5); //round
	}

	if (res < 0)
		res = 0;

	return res;
}

// 2022/9/9 BINV algorithm to generate random Binomial distributed numbers
// KACHITVICHYANUKUL and SCHMEISER, Binomial Random Variate Generation. 1980
int RyR:: Update_RyR_rates_stochastic_BINV(int num, double p)
{
	//std::cout << num << std::endl;
	if (num <= 0) return 0;
	if (p >= 1)
	{
		p = 1;
		return num;
	}
	if (num == 1) {
		if (random_num_RyR.gen_rand() <= p)
			return 1;
		else
			return 0;
	}

	double q = 1 - p;
	double s = p / q;
	double a = (num + 1) * s;
	double r = pow(q,num);
	double u = random_num_RyR.gen_rand();
	int    x = 0;

	while (u > r) {
		u = u - r;
		x = x + 1;
		r = (a / x - s) * r;
	}

	return x;
}


// using the C++11 library of binomial_distribution, this is very slow!!!
int RyR::Update_RyR_rates_stochastic_cpp11(int num, double p) {

	// perform 4 trials, each succeeds 1 in 2 times
	std::binomial_distribution<> d(num, p);

	return d(gen);

}

// https://github.com/numpy/numpy/blob/main/numpy/random/src/distributions/distributions.c
int RyR::Update_RyR_rates_stochastic_python( int n, double p) {
	if (n <= 0) return 0;
	if (p >= 1) {
		//std::cout << p << std::endl;
		p = 1;
		return n;
	}
	if (n == 1) {
		if (random_num_RyR.gen_rand() <= p)
			return 1;
		else
			return 0;
	}


  double q, qn, np, px, U;
  int X, bound;


  q     = 1.0 - p;
  qn    = exp(n * log(q));
  np    = n * p;
  bound = (int)fmin(n, np + 10.0 * sqrt(np * q + 1));

  X     = 0;
  px    = qn;
  U     = random_num_RyR.gen_rand();

  while (U > px) {
    X++;
    if (X > bound) {
      X = 0;
      px = qn;
      U = random_num_RyR.gen_rand();
    } else {
      U -= px;
      px = ((n - X + 1) * p * px) / (X * q);
    }
  }
  return X;
}





