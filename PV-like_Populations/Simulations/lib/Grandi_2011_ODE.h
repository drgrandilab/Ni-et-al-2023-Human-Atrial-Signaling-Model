// Header file for the ICs and RHS functions for the Colman 2013 model

#ifndef GRANDI_ODE_2011_H
#define GRANDI_ODE_2011_H

// #include "atrial_parameter.h"
#include "SingleCellParameter.hpp"
// #include <cmath>



#define Y0 /*RCONST*/(8.597401e-5)    // 0: Ca_i (mM) (in Ca_Concentrations)      
#define Y1 /*RCONST*/(1.737475e-4)   // 1: Ca_j (mM) (in Ca_Concentrations)      
#define Y2 /*RCONST*/(1.031812e-4)   // 2: Ca_sl (mM) (in Ca_Concentrations)      
#define Y3 /*RCONST*/(2.911916e-4)   // 3: CaM (mM) (in Cytosolic_Ca_Buffers)      
#define Y4 /*RCONST*/(1.298754e-3)    // 4: Myo_c (mM) (in Cytosolic_Ca_Buffers)      
#define Y5 /*RCONST*/(1.381982e-1)      // 5: Myo_m (mM) (in Cytosolic_Ca_Buffers)      
#define Y6 /*RCONST*/(2.143165e-3)    // 6: SRB (mM) (in Cytosolic_Ca_Buffers)      
#define Y7 /*RCONST*/(1.078283e-1)      // 7: Tn_CHc (mM) (in Cytosolic_Ca_Buffers)      
#define Y8 /*RCONST*/(1.524002e-2)     // 8: Tn_CHm (mM) (in Cytosolic_Ca_Buffers)      
#define Y9 /*RCONST*/(8.773191e-3)    // 9: Tn_CL (mM) (in Cytosolic_Ca_Buffers)      
#define Y10 /*RCONST*/(7.175662e-6)   // 10: d (dimensionless) (in I_Ca)      
#define Y11 /*RCONST*/(1.000681)     // 11: f (dimensionless) (in I_Ca)      
#define Y12 /*RCONST*/(2.421991e-2)    // 12: f_Ca_Bj (dimensionless) (in I_Ca)      
#define Y13 /*RCONST*/(1.452605e-2)    // 13: f_Ca_Bsl (dimensionless) (in I_Ca)      
#define Y14 /*RCONST*/(8.641386e-3)    // 14: x_kr (dimensionless) (in I_Kr)      
#define Y15 /*RCONST*/(5.412034e-3)   // 15: x_ks (dimensionless) (in I_Ks)      
#define Y16 /*RCONST*/(9.867005e-1)     // 16: h (dimensionless) (in I_Na)      
#define Y17 /*RCONST*/(9.915620e-1)     // 17: j (dimensionless) (in I_Na)      
#define Y18 /*RCONST*/(1.405627e-3)           // 18: m (dimensionless) (in I_Na)      
#define Y19 /*RCONST*/(9.945511e-1)  // 19: x_to_f (dimensionless) (in I_to)          
#define Y20 /*RCONST*/(4.051574e-3)  // 20: x_to_s (dimensionless) (in I_to)          
#define Y21 /*RCONST*/(9.945511e-1)     // 21: y_to_f (dimensionless) (in I_to)      
#define Y22 /*RCONST*/(9.945511e-1)     // 22: y_to_s (dimensionless) (in I_to)      
#define Y23 /*RCONST*/(7.347888e-3)    // 23: SLH_j (mM) (in Junctional_and_SL_Ca_Buffers)      
#define Y24 /*RCONST*/(7.297378e-2)     // 24: SLH_sl (mM) (in Junctional_and_SL_Ca_Buffers)      
#define Y25 /*RCONST*/(9.566355e-3)   // 25: SLL_j (mM) (in Junctional_and_SL_Ca_Buffers)      
#define Y26 /*RCONST*/(1.110363e-1)   // 26: SLL_sl (mM) (in Junctional_and_SL_Ca_Buffers)      
#define Y27 /*RCONST*/(120.0)                 // 27: K_i (mM) (in K_Concentration)
#define Y28 /*RCONST*/(3.539892)       // 28: Na_Bj (mM) (in Na_Buffers)  
#define Y29 /*RCONST*/(7.720854e-1)     // 29: Na_Bsl (mM) (in Na_Buffers)      
#define Y30 /*RCONST*/(9.136)     // 30: Na_i (mM) (in Na_Concentrations)      
#define Y31 /*RCONST*/(9.136)      // 31: Na_j (mM) (in Na_Concentrations)      
#define Y32 /*RCONST*/(9.136)      // 32: Na_sl (mM) (in Na_Concentrations)      
#define Y33 /*RCONST*/(0.01)     // 33: Ca_sr (mM) (in SR_Ca_Concentrations)      
#define Y34 /*RCONST*/(1.242988)      // 34: Csqn_b (mM) (in SR_Ca_Concentrations)      
#define Y35 /*RCONST*/(1.024274e-7)   // 35: Ry_Ri (mM) (in SR_Fluxes)      
#define Y36 /*RCONST*/(8.156628e-7)   // 36: Ry_Ro (mM) (in SR_Fluxes)      
#define Y37 /*RCONST*/(8.884332e-1)     // 37: Ry_Rr (mM) (in SR_Fluxes)      
#define Y38 /*RCONST*/(-8.09763e+1)     // 38: V_m (mV) (in membrane_potential)      
#define Y39 /*RCONST*/(0.0)                   // 39: ikur gate1
#define Y40 /*RCONST*/(1.0)                   // 40: ikur gate2                      




void Grandi_2011_ODE_Initialise(double *state, int celltype);
void Grandi_2011_ODE_VH_Update_Initialise(double *state, int celltype);
// double Grandi_2011_ODE(const double dt, const double Istim,  SingleCellPara &data, double *state);
// double Grandi_2011_ODE_New_INa(const double dt, const double Istim,  SingleCellPara &data, double *state);


void Grandi_2011_ODE(double, double *, double *, void *);
void Grandi_2011_ODE_VH_Update(double t, double *Y, double * dY, void *usr_data);
void Grandi_2011_ODE(double t, double *Y, double * dY, SingleCellPara *Data);




#endif



