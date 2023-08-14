/*generic function 

	for holding c type ode functions, decoupling voltage from the code
*/
#ifndef LSODA_FUNCTION_WRAP_H
#define LSODA_FUNCTION_WRAP_H

// #include "parameters.h"
// #include "cvode_struct.h"
// #include "Human_PF_Iyer_2014.h"

// #define NEQ_Human_PF_Iyer_2014 82
// int f_Human_PF_Iyer_2014(realtype t, N_Vector y, N_Vector ydot, void *user_data);

#include "SingleCellParameter.hpp"
#include "lsoda_C.h"
#include "HAM_Cell.hpp"
#include "HAM_Signalling.hpp"

void lsoda_generic_ODE(double t, double *Y, double * dY, void *Data );
void lsoda_generic_ODE_Vm_as_para(double t, double *Y, double * dY, void *Data );

void lsoda_generic_HAM_Signalling_ODE(double t, double *Y, double * dY, void *usr_data) ;
void lsoda_generic_HAM_Signalling_ODE_Vm_as_para(double t, double *Y, double * dY, void *usr_data);
#endif