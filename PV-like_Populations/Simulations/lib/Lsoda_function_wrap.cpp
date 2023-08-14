#include "Lsoda_function_wrap.h"



void lsoda_generic_ODE(double t, double *Y, double * dY, void *usr_data) {


	Cell *Data = (Cell*) usr_data;

	// in the grandi model, Y[38] is the membrane potential
	// Y[38] = Data->V;
	// Data -> V = Y[38];
	human_atrial_ECC_ODE(t, Y, dY, Data );
	// Data->dV = dY[38];  // get the dV value save in the single cell data, and leave to solve using forward euler.
	// dY[38] = 0.0;  // remember to set the derivative of vm in the single cell to be zero.
}



void lsoda_generic_ODE_Vm_as_para(double t, double *Y, double * dY, void *usr_data) {

	Cell *Data = (Cell*) usr_data;
	double V;

	/*if (t < 1000) V = -80;
	else if (t < 1000 + 500) V = 20;
	else
	{
		V=-80;
	}*/
	// in the grandi model, Y[38] is the membrane potential
	Y[38] = Data->V;
	// Data -> V = Y[38];
	human_atrial_ECC_ODE(t, Y, dY, Data );
	// Data->dV = dY[38];  // get the dV value save in the single cell data, and leave to solve using forward euler.
	dY[38] = 0.0;  // remember to set the derivative of vm in the single cell to be zero.
}




void lsoda_generic_HAM_Signalling_ODE(double t, double *Y, double * dY, void *usr_data) {


	HAM_Signalling *Data = (HAM_Signalling*) usr_data;

	// Y = Data->y;
	Data->Master_ODE_update(t, dY);



	// in the grandi model, Y[38] is the membrane potential
	// Y[38] = Data->V;
	// Data -> V = Y[38];
	// human_atrial_ECC_ODE(t, Y, dY, Data );
	// Data->dV = dY[38];  // get the dV value save in the single cell data, and leave to solve using forward euler.
	// dY[38] = 0.0;  // remember to set the derivative of vm in the single cell to be zero.
}
void lsoda_generic_HAM_Signalling_ODE_Vm_as_para(double t, double *Y, double * dY, void *usr_data) {


	HAM_Signalling *Data = (HAM_Signalling*) usr_data;
	Y[38] = Data->ECC_Module.V;


	// Y = Data->y;
	Data->Master_ODE_update(t, dY);

	dY[38] = 0.0;  // remember to set the derivative of vm in the single cell to be zero.


	// in the grandi model, Y[38] is the membrane potential
	// Y[38] = Data->V;
	// Data -> V = Y[38];
	// human_atrial_ECC_ODE(t, Y, dY, Data );
	// Data->dV = dY[38];  // get the dV value save in the single cell data, and leave to solve using forward euler.
	// dY[38] = 0.0;  // remember to set the derivative of vm in the single cell to be zero.
}