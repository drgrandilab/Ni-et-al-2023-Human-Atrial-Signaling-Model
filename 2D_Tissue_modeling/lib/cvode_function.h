// #include "Lsoda_function_wrap.h"
// #include <cvode/cvode.h>                          /* prototypes for CVODE fcts., consts.          */
#include "cvode_solver.hpp"
#include "HAM_Signalling.hpp"

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

int fnew(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int fnew_vm_as_para(realtype t, N_Vector y, N_Vector ydot, void *user_data);
