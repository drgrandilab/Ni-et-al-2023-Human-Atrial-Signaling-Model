#ifndef CVODE_SOLVER_HPP
#define CVODE_SOLVER_HPP


#include <stdio.h>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

/* Header files with a description of contents used */
#include <cvode/cvode.h>                          /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_serial.h>               /* access to serial N_Vector                    */
#include <sunmatrix/sunmatrix_dense.h>            /* access to dense SUNMatrix                    */
#include <sunlinsol/sunlinsol_dense.h>            /* access to dense SUNLinearSolver              */
#include <sunmatrix/sunmatrix_band.h>             /* access to band SUNMatrix                     */
#include <sunlinsol/sunlinsol_band.h>             /* access to band SUNLinearSolver               */
#include <cvode/cvode_diag.h>                     /* access to CVDIAG linear solver               */
#include "sunnonlinsol/sunnonlinsol_newton.h"     /* access to the newton SUNNonlinearSolver      */
#include "sunnonlinsol/sunnonlinsol_fixedpoint.h" /* access to the fixed point SUNNonlinearSolver */
#include <sundials/sundials_types.h>              /* definition of realtype                       */
#include <sundials/sundials_math.h>               /* contains the macros ABS, SUNSQR, and EXP     */
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ    */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ    */




class cvode_solver
{
public:
	cvode_solver(int NEQ, double max_dt = 1) {
		t = 0;
		m_NEQ = NEQ;
		dt_max = max_dt;
		// std::cout << m_NEQ << std::endl;

		y = N_VNew_Serial(m_NEQ); // for solver

		// ydot = N_VNew_Serial(NEQ); // for solver
	}


	int set_IC(double *state) {
		if (!state) {
			std::cout << " error in state cvode!!" << std::endl;
			std::exit(0);
		}
		for (int i = 0; i < m_NEQ; i++)
			Ith(y, i + 1) = state[i];
		return 0;
	}

	int set_IC(N_Vector state) {
		// if (!state) {
		// 	std::cout << " error in state cvode!!" << std::endl;
		// 	std::exit(0);
		// }
		for (int i = 0; i < m_NEQ; i++)
			Ith(y, i + 1) = Ith(state, i + 1);
		return 0;
	}


	int initialise_mem(CVRhsFn f) {
		cvode_mem = CVodeCreate(CV_BDF);//, CV_NEWTON);
		CVodeInit(cvode_mem, f, t, y);
		CVodeSStolerances(cvode_mem, 10e-7, 5e-7);

		CVodeSetMaxStep(cvode_mem, dt_max);
		A = SUNDenseMatrix(m_NEQ, m_NEQ);
		LS = SUNLinSol_Dense(y, A);
		// NLS = SUNNonlinSol_Newton(y);
		CVodeSetLinearSolver(cvode_mem, LS, A);
		// CVodeSetNonlinearSolver(cvode_mem, NLS);
		return 0;
	}



	int set_user_data(void *usr_data) {
		CVodeSetUserData(cvode_mem, usr_data);


		return (0);
	}
	int solve_single_step(double tout) {
		CVode(cvode_mem, tout, y, &t, CV_NORMAL); // 1 time step solution.
		return 0;
	}

	void destroy_cvode_mem() {

		/* Free the u vector          */
		CVodeFree(&cvode_mem);  /* Free the integrator memory */
		SUNLinSolFree(LS);      /* Free linear solver memory  */
		// SUNNonlinSolFree(NLS);      /* Free linear solver memory  */
		SUNMatDestroy(A);       /* Free the matrix memory     */
		// free(data);

	}
	~cvode_solver() {
		N_VDestroy(y);
		destroy_cvode_mem();
	}





	realtype t;// tout;
	N_Vector y, ydot; // for solve
	SUNMatrix A;
	SUNLinearSolver LS;
	SUNNonlinearSolver NLS;
	int m_NEQ = 0;
	void *cvode_mem; // for solver
	double dt_max;
};


#endif