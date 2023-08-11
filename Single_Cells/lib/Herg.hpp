

#ifndef HERG_Markov_HPP
#define HERG_Markov_HPP

// #include "HAM_Cell.hpp"
#include <fstream>

#include "HAM_Constants.h"
#include <math.h>

#include "signalling_para.hpp"

class Herg
{
public:
	Herg() {};
	~Herg() {};

	double state[10];  // seven states for model 1
	double ydot[10];  // seven states for model 1

	double IKr = 0;
	double GKr = 0.035 * sqrt(Ko / 5.4);

	double update_states(double Vm, double Drug_D, double Ek, signalling_para & para) {

		// model from https://www.frontiersin.org/articles/10.3389/fphys.2018.01402/full#supplementary-material
		double C1     = state[0];
		double C2     = state[1];
		double C3     = state[2];
		double I      = state[3];
		double O      = state[4];
		double I_star = state[5];
		double O_star = state[6];

		double alpha = 1.8e-1 * exp( 4.6e-2 * (Vm - 2.6e1 + para.Vkr_PKA) );
		double beta = 1.9e-1 * exp(-4.6e-2 * (Vm + para.Vkr_PKA));
		double alpha1 = 2.7e0;
		double beta1 = 1.6e0;
		double alpha2 = 7.4e-2 * exp(1.8e-2 * (Vm - 7.6e1 + para.Vkr_PKA));
		double beta2 = 2.7e-3 * exp(-6.5e-3 * (Vm + para.Vkr_PKA));
		double alphai = 1.7e-1 * exp(-1.7e-2 * (Vm + 1.6e1) * (4.5 / Ko));
		double betai = 8.5e-1 * exp(1.8e-2 * Vm * pow(4.5 / Ko, 0.3));
		double mu = alphai * beta2 / betai;


		double Ka = 2.01e1 / 1e3; //(1/uM/ms)
		double lA = 2.87e-4 / 1e3; // //(1/uM/ms)
		double Ki = 7.14e-1 / 1e3; //(1/uM/ms)
		double lI = 8.99e-6 / 1e3; //(1/uM/ms)


		double dC1     = beta * C2 - alpha * C1;
		double dC2     = alpha * C1 + beta1 * C3 - (beta + alpha1) * C2;
		double dC3     = alpha1 * C2 + beta2 * O + mu * I - (beta1 + alpha2 + alpha2) * C3;

		double dI      = alpha2 * C3 + betai * O + lI * I_star - (mu + alphai + Ki * Drug_D) * I;
		double dO      = alpha2 * C3 + alphai * I + lA * O_star - (beta2 + betai + Ka * Drug_D) * O;
		double dI_star = Ki * Drug_D * I - lI * I_star;
		double dO_star = Ka * Drug_D * O - lA * O_star;

		ydot [0] = dC1    ;
		ydot [1] = dC2    ;
		ydot [2] = dC3    ;
		ydot [3] = dI     ;
		ydot [4] = dO     ;
		ydot [5] = dI_star;
		ydot [6] = dO_star;

		// double

		IKr =  (1 + para.dGkr_PKA) * GKr * O * (Vm - Ek);


	}


};

#endif