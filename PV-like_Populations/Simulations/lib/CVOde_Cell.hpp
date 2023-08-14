#ifndef CVODE_CELL__HPP
#define CVODE_CELL__HPP

#include "cvode_solver.hpp"
// #include "stimulus.h"
// #include "APInfo.hpp"

#include "HAM_Signalling.hpp"
typedef int  (*cvode_func) (realtype t, N_Vector y, N_Vector ydot, void *user_data);

class CVOde_Cell
{
public:
	HAM_Signalling cell;
	cvode_solver cvode;// (NEQ, 0.2);


	cvode_func member_cvode_func;
	std::ofstream output_file;//( "HAM_wrap_out.dat");                      // output filename

	int m_cell_ID;
	int m_cell_type_ID;


	CVOde_Cell(int NEQ, double t_step_max, cvode_func ode_function,
	           bool output_data = false, int cellID = 0, int cell_type_ID = 0)
		: cell(HAM_Signalling())
		, cvode(cvode_solver(NEQ, t_step_max))
		, member_cvode_func(ode_function)
		, m_cell_ID(cellID)
		, m_cell_type_ID(cell_type_ID)
	{
		// cvode_func member_cvode_func = ode_function;
		initialise();
		// cell.V  = cell.y[38];  // assign Vm here

		/*if (output_data) {
			output_file.open( "HAM_wrap_out.dat");                      // output filename
		}*/

	}
	~CVOde_Cell() {
		/*if (output_file.is_open()) {
			output_file.close();
		}*/
	}

	void initialise() {

		cvode.set_IC(cell.y);
		cell.ECC_Module.V  = cell.ECC_Module.y[38];  // assign Vm here
		cvode.initialise_mem(member_cvode_func);
		cvode.set_user_data(&cell);
	}


	void solve_single_time_step_vm_para(double tout, double dt) {
		// CVode(cvode_mem, tout, y, &t, CV_NORMAL); // 1 time step solution.

		// std::cout << cell.ECC_Module.dV << std::endl;
		cvode.solve_single_step(tout);
		cell.ECC_Module.V += dt * cell.ECC_Module.dV;
	}


	void solve_single_time_step(double tout, double dt) {
		// CVode(cvode_mem, tout, y, &t, CV_NORMAL); // 1 time step solution.
		cvode.solve_single_step(tout);
		// cell.ECC_Module.V += dt * cell.ECC_Module.dV;
	}

	/*void print_to_file() {
		cell.print_to_file(cvode.t, output_file);
	}
	*/
	void read_initial_condition(const char * filename) {
		cvode.destroy_cvode_mem();
		cell.read_initial_condition(filename);
		initialise();
	}



};



#endif