// #include "Lsoda_function_wrap.h"
// #include <cvode/cvode.h>                          /* prototypes for CVODE fcts., consts.          */
#include "cvode_solver.hpp"
#include "cvode_function.h"


int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

	HAM_Signalling *Data = (HAM_Signalling*) user_data;
	int NEQ = Data->ODE_NUM;



	// std::cout << NEQ<< std::endl;
	// std::cout << t<< std::endl;

	// realtype Y[NEQ];
	// double dY[NEQ]={0};
	// double
	// #pragma novector
	for (int i = 0; i < NEQ; i++)
		Data->y[i] = Ith(y, i + 1);

	// std::cout << Ith(y, 38+1)<< std::endl;

	Data->Master_ODE_update_CVODE(t);


	for (int i = 0; i < NEQ; i++)
		Ith(ydot, i + 1) = Data->ydot[i];
	// std::cout << Ith(ydot, 1)<< std::endl;

	// Y = Data->y;
	return 0;

}




int fnew(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

	HAM_Signalling *Data = (HAM_Signalling*) user_data;
	int NEQ_m = Data->ODE_NUM;

	for (int i = 0; i < NEQ_m; i++)
		Data->y[i] = Ith(y, i + 1);
	// Data->V  = Data->y[38];
	Data->Master_ODE_update_CVODE(t);


	for (int i = 0; i < NEQ_m; i++)
		Ith(ydot, i + 1) = Data->ydot[i];

	return 0;

}



int fnew_vm_as_para(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

	HAM_Signalling *Data = (HAM_Signalling*) user_data;
	int NEQ_m = Data->ODE_NUM;

	for (int i = 0; i < NEQ_m; i++)
		Data->y[i] = Ith(y, i + 1);
	Data->ECC_Module.y[38] = Data->ECC_Module.V;// Vm updated outside solver;
	Data->Master_ODE_update_CVODE(t);
	Data->ydot[38] = 0;  // Vm updated outside solver;

	for (int i = 0; i < NEQ_m; i++)
		Ith(ydot, i + 1) = Data->ydot[i];
	// std::cout << Ith(ydot, 1)<< std::endl;

	// y = Data->y;
	return 0;

}


