#ifndef HAM_SIGNALLING_HPP
#define HAM_SIGNALLING_HPP

#include  <iomanip>
#include <iostream>
#include <fstream>


// #include "betaAR.hpp"
#include "betaAR_update.hpp"  // using the updated Beta AR module
#include "CaM.hpp"
#include "CaMKII.hpp"
#include "HAM_ECC.hpp"



// int human_atrial_ECC_ODE(double t, double *y, double *ydot, Cell * cell_par);

// int human_atrial_ECC_initialiser(double *y);

class HAM_Signalling
{
public:
	double * y;
	double * ydot;

	int ODE_NUM;

	// modules
	HAM_ECC ECC_Module;

	betaAR_Update betaAR_Module;

	CaM CaM_Module_dyad, CaM_Module_sl, CaM_Module_cyto;

	CaMKII CaMKII_Module;

	double BCL, S2;

	// parameters
	signalling_para cell_para;
	bool Stimuli_in_ODE;
	HAM_Signalling(int AF = 0, double ISO = 0, bool CaMKII_double = false);
	~HAM_Signalling();

	void Master_ODE_update_CVODE(double t/*, double *ydot*/);
	void Master_ODE_update(double t, double *ydot);

	void assign_AF_ISO_CaMKII_para(int AF = 0, double ISO = 0, bool CaMKII_double = false);


	void read_initial_condition(const char* filename);
	void output_inital_condition(const char* filename);
	void assign_cell_pop_para(double *in_para);
	void report_cell_pop_para();
};














#endif