#include <math.h>
#include "HAM_Signalling.hpp"
#include "stimulus.h"

// first to add all functions here.
// 9:37:15, Mon, 04-March-2019, By Haibo

HAM_Signalling::HAM_Signalling(int AF, double ISO, bool CaMKII_double) {

	Stimuli_in_ODE = true;

	// set totoal numbe of ODEs, preparing for allocating memory
	ODE_NUM = ECC_Module. ODE_NUM + CaM_Module_dyad. ODE_NUM + CaM_Module_sl. ODE_NUM + CaM_Module_cyto. ODE_NUM + CaMKII_Module. ODE_NUM + betaAR_Module. ODE_NUM;

	y = new double [ODE_NUM];
	ydot = new double [ODE_NUM];


	for (int i = 0; i < ODE_NUM; ++i)
	{
		ydot[i] = 0;
	}

	// ydot will be assigned in Master ODEs


	ECC_Module.y      = &y[0];
	CaM_Module_dyad.y = &y[0 + ECC_Module.ODE_NUM ];
	CaM_Module_sl.y   = &y[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM];
	CaM_Module_cyto.y = &y[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM + CaM_Module_sl.ODE_NUM];
	CaMKII_Module.y   = &y[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM + CaM_Module_sl.ODE_NUM + CaM_Module_cyto.ODE_NUM];
	betaAR_Module.y   = &y[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM + CaM_Module_sl.ODE_NUM + CaM_Module_cyto.ODE_NUM + CaMKII_Module.ODE_NUM ];

	ECC_Module.initialiser();
	CaM_Module_dyad.initialiser(0);
	CaM_Module_sl.initialiser(1);
	CaM_Module_cyto.initialiser(2);
	CaMKII_Module.initialiser();
	betaAR_Module.initialiser();



	// concentrations for each compartment in the CaM modules
	cell_para.cell_CaM_para_dyad.K         =  135.0000    ;
	cell_para.cell_CaM_para_dyad.Mg        =  1.0000  ;
	cell_para.cell_CaM_para_dyad.CaMtot    =  418.0000         ;
	cell_para.cell_CaM_para_dyad.Btot      =  0  ;  // note that here it's zero
	cell_para.cell_CaM_para_dyad.CaMKIItot =  120.0000    ;
	cell_para.cell_CaM_para_dyad.CaNtot    =  3.6;//3.61750874231279   ;
	cell_para.cell_CaM_para_dyad.PP1tot    =  96.5000    ;
	cell_para.cell_CaM_para_dyad.Ca        =  0.3494  ;

	cell_para.cell_CaM_para_sl.K         =  135;
	cell_para.cell_CaM_para_sl.Mg        =  1;
	cell_para.cell_CaM_para_sl.CaMtot    =  5.65000000000000;
	cell_para.cell_CaM_para_sl.Btot      =  24.2000000000000;
	cell_para.cell_CaM_para_sl.CaMKIItot =  0.0995160000000000;
	cell_para.cell_CaM_para_sl.CaNtot    =  0.00300000000000000;
	cell_para.cell_CaM_para_sl.PP1tot    =  0.570000000000000;
	cell_para.cell_CaM_para_sl.Ca        =  0.243283499183151;

	cell_para.cell_CaM_para_cyto.K         =  135;
	cell_para.cell_CaM_para_cyto.Mg        =  1;
	cell_para.cell_CaM_para_cyto.CaMtot    =  5.65000000000000;
	cell_para.cell_CaM_para_cyto.Btot      =  24.2000000000000;
	cell_para.cell_CaM_para_cyto.CaMKIItot =  0.0995160000000000;
	cell_para.cell_CaM_para_cyto.CaNtot    =  0.00300000000000000;
	cell_para.cell_CaM_para_cyto.PP1tot    =  0.570000000000000;
	cell_para.cell_CaM_para_cyto.Ca        =  0.219701936526104;


	assign_AF_ISO_CaMKII_para(AF, ISO, CaMKII_double);
}

void HAM_Signalling::assign_AF_ISO_CaMKII_para(int AF, double ISO, bool CaMKII_double) {
	// change ISO concentration here.
	cell_para.cell_betaAR_para.ISO = ISO;
	cell_para.AF = AF;

	if (AF == 1) {
		cell_para.cell_CaM_para_dyad.CaMKIItot =  120.0000 * 2.0;

		cell_para.cell_CaM_para_sl.CaMKIItot =  0.0995160000000000 * 2.0;

		cell_para.cell_CaM_para_cyto.CaMKIItot =  0.0995160000000000 * 2.0;

	}

	if (CaMKII_double) {
		cell_para.cell_CaM_para_dyad.CaMKIItot =  120.0000 * 2.0    ;

		cell_para.cell_CaM_para_sl.CaMKIItot =  0.0995160000000000 * 2.0;

		cell_para.cell_CaM_para_cyto.CaMKIItot =  0.0995160000000000 * 2.0;

	}

}


HAM_Signalling::~HAM_Signalling() {
	// hosue keeping here
	if (y != nullptr) {
		delete [] y;
		y = nullptr;
	}
	if (ydot != nullptr) {
		delete [] ydot;
		ydot = nullptr;
	}
}




// here, ydot from HAM_Signalling class
void HAM_Signalling::Master_ODE_update_CVODE(double t/*, double *ydot*/) {



	// cell_para.Istim  = S1S2(0.0, -12.5, BCL, 10000, S2, t, 5.0);
	if (Stimuli_in_ODE) {
		ECC_Module.I_app = -S1S2(0.0, -12.5, BCL, 10000, BCL, t, 5.0);
		// ECC_Module.I_app = S1S2_num(0.0, 12.5, 300, 100, BCL, 10000,t, 5.0);  // (double stim_start_time, double stim_current, double BCL_s1, int S1_num, double BCL_s2, int S2_num, double current_time, double stim_duration)

		// if (t > 401e3 + 10)   // stop pacing at time = 401e3 ms
		if (ECC_Module.allow_stimulus_app == false)
			ECC_Module.I_app = 0.;
	}


	// assign ydot addresses
	ECC_Module.ydot      = &ydot[0];
	CaM_Module_dyad.ydot = &ydot[0 + ECC_Module.ODE_NUM ];
	CaM_Module_sl.ydot   = &ydot[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM];
	CaM_Module_cyto.ydot = &ydot[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM + CaM_Module_sl.ODE_NUM];
	CaMKII_Module.ydot   = &ydot[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM + CaM_Module_sl.ODE_NUM + CaM_Module_cyto.ODE_NUM];
	betaAR_Module.ydot   = &ydot[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM + CaM_Module_sl.ODE_NUM + CaM_Module_cyto.ODE_NUM + CaMKII_Module.ODE_NUM ];


	/*ECC_Module.y[35] =  1.9953e-06*1000;// y[35];  (pCa = 5.7)
	ECC_Module.y[36] =  1.9953e-06*1000;// y[36];
	ECC_Module.y[37] = 1.9953e-06*1000;//  y[37];
	*/


	/*ECC_Module.y[35] =  1e-7 * 1000;// y[35];  (pCa = 7)
	ECC_Module.y[36] =  1e-7 * 1000;// y[36];
	ECC_Module.y[37] = 1e-7 * 1000;//  y[37];*/

	// update Cai as inputs for CaM modules
	cell_para.cell_CaM_para_dyad.Ca =  ECC_Module.y[35] * 1e3;   // changing units
	cell_para.cell_CaM_para_sl.Ca   =  ECC_Module.y[36] * 1e3;
	cell_para.cell_CaM_para_cyto.Ca =  ECC_Module.y[37] * 1e3;
	// update CaMKII parameters here// 14:21:22, Mon, 04-March-2019, By Haibo


	// CaMKII phosphorylation module
	cell_para.cell_CaMKII_para.CaMKIIactDyad = cell_para.cell_CaM_para_dyad.CaMKIItot * (CaM_Module_dyad.y[7] + CaM_Module_dyad.y[8] + CaM_Module_dyad.y[9] + CaM_Module_dyad.y[10]); // Multiply total by fraction
	cell_para.cell_CaMKII_para.CaMKIIactSL = cell_para.cell_CaM_para_sl.CaMKIItot * (CaM_Module_sl.y[7] + CaM_Module_sl.y[8] + CaM_Module_sl.y[9] + CaM_Module_sl.y[10]); // Multiply total by fraction

	// CaMKIIact_SL = cell_para.cell_CaM_para_sl.CaMKIItot * (y(ny_ECC + ny_cam + 8) + y(ny_ECC + ny_cam + 9) + y(ny_ECC + ny_cam + 10) + y(ny_ECC + ny_cam + 11));
	//PP1_PLB_avail = y(ny_ECC+3*ny_cam+ny_CaMKII+22)./PP1_PLBtot + .0091;  // Active PP1 near PLB / total PP1 conc + basal value
	// double PP1_PLB_avail = 1 - y(ny_ECC + 3 * ny_cam + ny_CaMKII + 24) / PP1_PLBtot + 0.081698; // NEW ODEs - Mouse model
	cell_para.cell_CaMKII_para.PP1_PLB_avail = 1.0 - betaAR_Module.y[23] / cell_para.cell_betaAR_para.PP1tot + 0.081698; // NEW ODEs - Mouse model



	cell_para.LCC_CKdyadp = CaMKII_Module.y[1] / cell_para.cell_CaMKII_para.LCCtotDyad; // fractional CaMKII-dependent LCC dyad phosphorylation
	cell_para.NaV_CKp = CaMKII_Module.y[2] / cell_para.cell_CaMKII_para.NaVtot; // fractional CaMKII-dependent NaV phosphorylation
	cell_para.RyR_CKp = CaMKII_Module.y[3] / cell_para.cell_CaMKII_para.RyRtot; // fractional CaMKII-dependent RyR phosphorylation
	cell_para.PLB_CKp = CaMKII_Module.y[4] / cell_para.cell_CaMKII_para.PLBtot; // fractional CaMKII-dependent PLB phosphorylation




	// update PKA-phosphorylated factions here

	// adding scaling factor for the PKA phoslevels;
	cell_para.PLB_PKAn   =  (cell_para.cell_betaAR_para.PLBtot - cell_para.PLB_PKAn_Scale * betaAR_Module.y[25]) / cell_para.cell_betaAR_para.PLBtot; // non-phosphorylated PLB targets
	cell_para.PLM_PKAp   = cell_para.PLM_PKAp_Scale * betaAR_Module.y[26] / cell_para.cell_betaAR_para.PLMtot;   //
	cell_para.LCCa_PKAp  = cell_para.LCCa_PKAp_Scale * betaAR_Module.y[27] / cell_para.cell_betaAR_para.LCCtot;   //
	cell_para.LCCb_PKAp  = cell_para.LCCb_PKAp_Scale * betaAR_Module.y[28] / cell_para.cell_betaAR_para.LCCtot;   //
	cell_para.RyR_PKAp   = cell_para.RyR_PKAp_Scale * betaAR_Module.y[29] / cell_para.cell_betaAR_para.RyRtot;   //
	cell_para.TnI_PKAp   = cell_para.TnI_PKAp_Scale * betaAR_Module.y[30] / cell_para.cell_betaAR_para.TnItot;   //
	cell_para.IKs_PKAp   = cell_para.IKs_PKAp_Scale * betaAR_Module.y[31] / cell_para.cell_betaAR_para.IKstot;   //
	cell_para.IKr_PKAp   = cell_para.IKr_PKAp_Scale * betaAR_Module.y[32] / cell_para.cell_betaAR_para.IKrtot;   //
	cell_para.IKur_PKAp  = cell_para.IKur_PKAp_Scale * betaAR_Module.y[33] / cell_para.cell_betaAR_para.IKurtot;   //
	cell_para.ICFTR_PKAp = cell_para.ICFTR_PKAp_Scale * betaAR_Module.y[34] / cell_para.cell_betaAR_para.ICFTRtot;   //
	cell_para.IClCa_PKAp = cell_para.IClCa_PKAp_Scale * betaAR_Module.y[35] / cell_para.cell_betaAR_para.IClCatot;   //
	cell_para.Myo_PKAp   = cell_para.Myo_PKAp_Scale * betaAR_Module.y[36] / cell_para.cell_betaAR_para.myotot;   //
	cell_para.INa_PKAp   = cell_para.INa_PKAp_Scale * betaAR_Module.y[37] / cell_para.cell_betaAR_para.INatot;   //
	cell_para.Ito_PKAp   = cell_para.Ito_PKAp_Scale * betaAR_Module.y[38] / cell_para.cell_betaAR_para.Itotot;   //
	cell_para.IK1_PKAp   = cell_para.IK1_PKAp_Scale * betaAR_Module.y[39] / cell_para.cell_betaAR_para.IK1tot;   //
	// applying same CaMKII-dep phophorylation factors to IKur and Itof,
	cell_para.Ikur_CKp = cell_para.NaV_CKp;
	cell_para.Itof_CKp = cell_para.NaV_CKp;
	cell_para.IK1_CKp = cell_para.NaV_CKp;

	// set Gap junction connexins phosp by CaMKII, same as INa
	cell_para.Gapjunct_CKp = cell_para.NaV_CKp;


	// adding scaling factor for the CaMKII phoslevels;
	cell_para.LCC_CKdyadp *= cell_para.LCC_CKdyadp_Scale;
	cell_para.NaV_CKp     *= cell_para.NaV_CKp_Scale    ;
	cell_para.RyR_CKp     *= cell_para.RyR_CKp_Scale    ;
	cell_para.PLB_CKp     *= cell_para.PLB_CKp_Scale    ;
	cell_para.Ikur_CKp    *= cell_para.Ikur_CKp_Scale   ;
	cell_para.Itof_CKp    *= cell_para.Itof_CKp_Scale   ;
	cell_para.IK1_CKp     *= cell_para.IK1_CKp_Scale    ;



	if (cell_para.CaMKII_Inhibition) {
		cell_para.set_CaMKII_Inhibition();
	}




	cell_para.update_signaling_modulators_for_new_bAR_Module(); // this calculates signalling modulating effects, producing scaling/shifting factors to be included in the ECC module;

	CaM_Module_dyad.update_single_time_step(cell_para.cell_CaM_para_dyad);
	CaM_Module_sl.update_single_time_step(cell_para.cell_CaM_para_sl);
	CaM_Module_cyto.update_single_time_step(cell_para.cell_CaM_para_cyto);
	CaMKII_Module.update_single_time_step(cell_para.cell_CaMKII_para);
	ECC_Module.HAM_ECC_update_ODE(t, cell_para);
	// adding beta AR later
	betaAR_Module.update_single_time_step(cell_para.cell_betaAR_para);


	// incorporate Ca buffering from CaM, convert JCaCyt from uM / msec to mM / msec
	// should remove CaM buffers in the ECC module later.
	ECC_Module.ydot[36 - 1] = ECC_Module.ydot[36 - 1] + 1e-3 * CaM_Module_dyad.JCa;
	ECC_Module.ydot[37 - 1] = ECC_Module.ydot[37 - 1] + 1e-3 * CaM_Module_sl.JCa;
	ECC_Module.ydot[38 - 1] = ECC_Module.ydot[38 - 1] + 1e-3 * CaM_Module_cyto.JCa;



	// Cell geometry // from ECC module
	// cellLength = 100;     // cell length [um]
	// cellRadius = 10.25;   // cell radius [um]
	// Vcell = pi * cellRadius *cellRadius  * cellLength * 1e-15; // [L]
	// Vmyo = 0.65 * Vcell; Vsl = 0.02 * Vcell; Vjunc = 1 * 0.0539 * .01 * Vcell;
	// incorporate CaM diffusion between compartments//   CaM diffuses between compartments 
	double kDyadSL = 3.6363e-16;	// [L/msec]
	double kSLmyo = 8.587e-15;     // [L/msec]
	double k0Boff = 0.0014;        // [s^-1]
	double k0Bon = k0Boff / 0.2;   // [uM^-1 s^-1] kon = koff/Kd
	double k2Boff = k0Boff / 100;  // [s^-1]
	double k2Bon = k0Bon;          // [uM^-1 s^-1]
	double k4Boff = k2Boff;        // [s^-1]
	double k4Bon = k0Bon;          // [uM^-1 s^-1]

	double CaMtotDyad = cell_para.cell_CaM_para_dyad.CaMKIItot  * (CaM_Module_dyad.y[6] + CaM_Module_dyad.y[7] + CaM_Module_dyad.y[8] + CaM_Module_dyad.y[9])
	                    + CaM_Module_dyad.y[0] + CaM_Module_dyad.y[1] + CaM_Module_dyad.y[2] + CaM_Module_dyad.y[3] + CaM_Module_dyad.y[4] + CaM_Module_dyad.y[5]
	                    + CaM_Module_dyad.y[12] + CaM_Module_dyad.y[13] + CaM_Module_dyad.y[14];


	// sum(y_camDyad(1: 6)) + CaMKIItotDyad * sum(y_camDyad(7: 10)) + sum(y_camDyad(13: 15));
	// adding comments: diffusion of CaM between 3 compartments.

	// ** NOTE: Btotdyad being sent to the dyad camODEfile is set to zero, but is used below for transfer between SL and dyad
	double BtotDyad = 1.54 / 8.293e-4; // [uM] from the matlab code. this value is also huge (compared to others)

	double Bdyad = BtotDyad - CaMtotDyad; // [uM dyad]
	double J_cam_dyadSL = 1e-3 * (k0Boff * CaM_Module_dyad.y[0] - k0Bon * Bdyad * CaM_Module_sl.y[0]); // [uM/msec dyad]
	double J_ca2cam_dyadSL = 1e-3 * (k2Boff * CaM_Module_dyad.y[1] - k2Bon * Bdyad * CaM_Module_sl.y[1]); // [uM/msec dyad]
	double J_ca4cam_dyadSL = 1e-3 * (k2Boff * CaM_Module_dyad.y[2] - k4Bon * Bdyad * CaM_Module_sl.y[2]); // [uM/msec dyad]
	double J_cam_SLmyo = kSLmyo * (CaM_Module_sl.y[0] - CaM_Module_cyto.y[0]); // [umol/msec]
	double J_ca2cam_SLmyo = kSLmyo * (CaM_Module_sl.y[1] - CaM_Module_cyto.y[1]); // [umol/msec]
	double J_ca4cam_SLmyo = kSLmyo * (CaM_Module_sl.y[2] - CaM_Module_cyto.y[2]); // [umol/msec]

	CaM_Module_dyad.ydot[0] = CaM_Module_dyad.ydot[0] - J_cam_dyadSL;
	CaM_Module_dyad.ydot[1] = CaM_Module_dyad.ydot[1] - J_ca2cam_dyadSL;
	CaM_Module_dyad.ydot[2] = CaM_Module_dyad.ydot[2] - J_ca4cam_dyadSL;
	CaM_Module_sl.ydot[0] = CaM_Module_sl.ydot[0] + J_cam_dyadSL * Vjunc / Vsl - J_cam_SLmyo / Vsl;
	CaM_Module_sl.ydot[1] = CaM_Module_sl.ydot[1] + J_ca2cam_dyadSL * Vjunc / Vsl - J_ca2cam_SLmyo / Vsl;
	CaM_Module_sl.ydot[2] = CaM_Module_sl.ydot[2] + J_ca4cam_dyadSL * Vjunc / Vsl - J_ca4cam_SLmyo / Vsl;
	CaM_Module_cyto.ydot[0] = CaM_Module_cyto.ydot[0] + J_cam_SLmyo / Vmyo;
	CaM_Module_cyto.ydot[1] = CaM_Module_cyto.ydot[1] + J_ca2cam_SLmyo / Vmyo;
	CaM_Module_cyto.ydot[2] = CaM_Module_cyto.ydot[2] + J_ca4cam_SLmyo / Vmyo;
}



// here, ydot will be from LSODA class
void HAM_Signalling::Master_ODE_update(double t, double *ydot) {



	// cell_para.Istim  = S1S2(0.0, -12.5, BCL, 10000, S2, t, 5.0);
	ECC_Module.I_app = -S1S2(0.0, -12.5, BCL, 10000, 1000, t, 5.0); // leaving stimulus function in the ODE speeds the implementation up.
	if (t > 401e3 + 10)
		ECC_Module.I_app = 0.;
	ECC_Module.ydot      = &ydot[0];
	CaM_Module_dyad.ydot = &ydot[0 + ECC_Module.ODE_NUM ];
	CaM_Module_sl.ydot   = &ydot[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM];
	CaM_Module_cyto.ydot = &ydot[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM + CaM_Module_sl.ODE_NUM];
	CaMKII_Module.ydot   = &ydot[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM + CaM_Module_sl.ODE_NUM + CaM_Module_cyto.ODE_NUM];
	betaAR_Module.ydot   = &ydot[0 + ECC_Module.ODE_NUM  + CaM_Module_dyad.ODE_NUM + CaM_Module_sl.ODE_NUM + CaM_Module_cyto.ODE_NUM + CaMKII_Module.ODE_NUM ];


	// update Cai as inputs for CaM modules
	cell_para.cell_CaM_para_dyad.Ca =  y[35] * 1e3;
	cell_para.cell_CaM_para_sl.Ca   =  y[36] * 1e3;
	cell_para.cell_CaM_para_cyto.Ca =  y[37] * 1e3;
	// update CaMKII parameters here


	// CaMKII phosphorylation module
	cell_para.cell_CaMKII_para.CaMKIIactDyad = cell_para.cell_CaM_para_dyad.CaMKIItot * (CaM_Module_dyad.y[7] + CaM_Module_dyad.y[8] + CaM_Module_dyad.y[9] + CaM_Module_dyad.y[10]); // Multiply total by fraction
	cell_para.cell_CaMKII_para.CaMKIIactSL = cell_para.cell_CaM_para_sl.CaMKIItot * (CaM_Module_sl.y[7] + CaM_Module_sl.y[8] + CaM_Module_sl.y[9] + CaM_Module_sl.y[10]); // Multiply total by fraction

	// CaMKIIact_SL = cell_para.cell_CaM_para_sl.CaMKIItot * (y(ny_ECC + ny_cam + 8) + y(ny_ECC + ny_cam + 9) + y(ny_ECC + ny_cam + 10) + y(ny_ECC + ny_cam + 11));
	//PP1_PLB_avail = y(ny_ECC+3*ny_cam+ny_CaMKII+22)./PP1_PLBtot + .0091;  // Active PP1 near PLB / total PP1 conc + basal value
	// double PP1_PLB_avail = 1 - y(ny_ECC + 3 * ny_cam + ny_CaMKII + 24) / PP1_PLBtot + 0.081698; // NEW ODEs - Mouse model
	cell_para.cell_CaMKII_para.PP1_PLB_avail = 1.0 - betaAR_Module.y[23] / cell_para.cell_betaAR_para.PP1tot + 0.081698; // NEW ODEs - Mouse model



	cell_para.LCC_CKdyadp = CaMKII_Module.y[1] / cell_para.cell_CaMKII_para.LCCtotDyad; // fractional CaMKII-dependent LCC dyad phosphorylation
	cell_para.NaV_CKp = CaMKII_Module.y[2] / cell_para.cell_CaMKII_para.NaVtot; // fractional CaMKII-dependent NaV phosphorylation
	cell_para.RyR_CKp = CaMKII_Module.y[3] / cell_para.cell_CaMKII_para.RyRtot; // fractional CaMKII-dependent RyR phosphorylation
	cell_para.PLB_CKp = CaMKII_Module.y[4] / cell_para.cell_CaMKII_para.PLBtot; // fractional CaMKII-dependent PLB phosphorylation




	// update PKA-phosphorylated factions here

	// adding scaling factor for the PKA phoslevels;
	cell_para.PLB_PKAn   =  (cell_para.cell_betaAR_para.PLBtot - cell_para.PLB_PKAn_Scale * betaAR_Module.y[25]) / cell_para.cell_betaAR_para.PLBtot; // non-phosphorylated PLB targets
	cell_para.PLM_PKAp   = cell_para.PLM_PKAp_Scale * betaAR_Module.y[26] / cell_para.cell_betaAR_para.PLMtot;   //
	cell_para.LCCa_PKAp  = cell_para.LCCa_PKAp_Scale * betaAR_Module.y[27] / cell_para.cell_betaAR_para.LCCtot;   //
	cell_para.LCCb_PKAp  = cell_para.LCCb_PKAp_Scale * betaAR_Module.y[28] / cell_para.cell_betaAR_para.LCCtot;   //
	cell_para.RyR_PKAp   = cell_para.RyR_PKAp_Scale * betaAR_Module.y[29] / cell_para.cell_betaAR_para.RyRtot;   //
	cell_para.TnI_PKAp   = cell_para.TnI_PKAp_Scale * betaAR_Module.y[30] / cell_para.cell_betaAR_para.TnItot;   //
	cell_para.IKs_PKAp   = cell_para.IKs_PKAp_Scale * betaAR_Module.y[31] / cell_para.cell_betaAR_para.IKstot;   //
	cell_para.IKr_PKAp   = cell_para.IKr_PKAp_Scale * betaAR_Module.y[32] / cell_para.cell_betaAR_para.IKrtot;   //
	cell_para.IKur_PKAp  = cell_para.IKur_PKAp_Scale * betaAR_Module.y[33] / cell_para.cell_betaAR_para.IKurtot;   //
	cell_para.ICFTR_PKAp = cell_para.ICFTR_PKAp_Scale * betaAR_Module.y[34] / cell_para.cell_betaAR_para.ICFTRtot;   //
	cell_para.IClCa_PKAp = cell_para.IClCa_PKAp_Scale * betaAR_Module.y[35] / cell_para.cell_betaAR_para.IClCatot;   //
	cell_para.Myo_PKAp   = cell_para.Myo_PKAp_Scale * betaAR_Module.y[36] / cell_para.cell_betaAR_para.myotot;   //
	cell_para.INa_PKAp   = cell_para.INa_PKAp_Scale * betaAR_Module.y[37] / cell_para.cell_betaAR_para.INatot;   //
	cell_para.Ito_PKAp   = cell_para.Ito_PKAp_Scale * betaAR_Module.y[38] / cell_para.cell_betaAR_para.Itotot;   //
	cell_para.IK1_PKAp   = cell_para.IK1_PKAp_Scale * betaAR_Module.y[39] / cell_para.cell_betaAR_para.IK1tot;   //
	// applying same CaMKII-dep phophorylation factors to IKur and Itof,
	cell_para.Ikur_CKp = cell_para.NaV_CKp;
	cell_para.Itof_CKp = cell_para.NaV_CKp;
	cell_para.IK1_CKp = cell_para.NaV_CKp;


	// adding scaling factor for the CaMKII phoslevels;
	cell_para.LCC_CKdyadp *= cell_para.LCC_CKdyadp_Scale;
	cell_para.NaV_CKp     *= cell_para.NaV_CKp_Scale    ;
	cell_para.RyR_CKp     *= cell_para.RyR_CKp_Scale    ;
	cell_para.PLB_CKp     *= cell_para.PLB_CKp_Scale    ;
	cell_para.Ikur_CKp    *= cell_para.Ikur_CKp_Scale   ;
	cell_para.Itof_CKp    *= cell_para.Itof_CKp_Scale   ;
	cell_para.IK1_CKp     *= cell_para.IK1_CKp_Scale    ;

	if (cell_para.CaMKII_Inhibition) {
		cell_para.set_CaMKII_Inhibition();
	}

	cell_para.update_signaling_modulators_for_new_bAR_Module(); // this calculates signalling modulating effects, producing scaling/shifting factors to be included in the ECC module;

	CaM_Module_dyad.update_single_time_step(cell_para.cell_CaM_para_dyad);
	CaM_Module_sl.update_single_time_step(cell_para.cell_CaM_para_sl);
	CaM_Module_cyto.update_single_time_step(cell_para.cell_CaM_para_cyto);
	CaMKII_Module.update_single_time_step(cell_para.cell_CaMKII_para);
	ECC_Module.HAM_ECC_update_ODE(t, cell_para);
	// adding beta AR later
	betaAR_Module.update_single_time_step(cell_para.cell_betaAR_para);


	// incorporate Ca buffering from CaM, convert JCaCyt from uM / msec to mM / msec
	// should remove CaM buffers in the ECC module later.
	ECC_Module.ydot[36 - 1] = ECC_Module.ydot[36 - 1] + 1e-3 * CaM_Module_dyad.JCa;
	ECC_Module.ydot[37 - 1] = ECC_Module.ydot[37 - 1] + 1e-3 * CaM_Module_sl.JCa;
	ECC_Module.ydot[38 - 1] = ECC_Module.ydot[38 - 1] + 1e-3 * CaM_Module_cyto.JCa;


	// Cell geometry // from ECC module
	// cellLength = 100;     // cell length [um]
	// cellRadius = 10.25;   // cell radius [um]
	// Vcell = pi * cellRadius *cellRadius  * cellLength * 1e-15; // [L]
	// Vmyo = 0.65 * Vcell; Vsl = 0.02 * Vcell; Vjunc = 1 * 0.0539 * .01 * Vcell;
	// incorporate CaM diffusion between compartments//   CaM diffuses between compartments 
	double kDyadSL = 3.6363e-16;	// [L/msec]
	double kSLmyo = 8.587e-15;     // [L/msec]
	double k0Boff = 0.0014;        // [s^-1]
	double k0Bon = k0Boff / 0.2;   // [uM^-1 s^-1] kon = koff/Kd
	double k2Boff = k0Boff / 100;  // [s^-1]
	double k2Bon = k0Bon;          // [uM^-1 s^-1]
	double k4Boff = k2Boff;        // [s^-1]
	double k4Bon = k0Bon;          // [uM^-1 s^-1]

	double CaMtotDyad = cell_para.cell_CaM_para_dyad.CaMKIItot  * (CaM_Module_dyad.y[6] + CaM_Module_dyad.y[7] + CaM_Module_dyad.y[8] + CaM_Module_dyad.y[9])
	                    + CaM_Module_dyad.y[0] + CaM_Module_dyad.y[1] + CaM_Module_dyad.y[2] + CaM_Module_dyad.y[3] + CaM_Module_dyad.y[4] + CaM_Module_dyad.y[5]
	                    + CaM_Module_dyad.y[12] + CaM_Module_dyad.y[13] + CaM_Module_dyad.y[14];


	// sum(y_camDyad(1: 6)) + CaMKIItotDyad * sum(y_camDyad(7: 10)) + sum(y_camDyad(13: 15));
	// adding comments: diffusion of CaM between 3 compartments.

	// ** NOTE: Btotdyad being sent to the dyad camODEfile is set to zero, but is used below for transfer between SL and dyad
	double BtotDyad = 1.54 / 8.293e-4; // [uM] from the matlab code. this value is also huge (compared to others)

	double Bdyad = BtotDyad - CaMtotDyad; // [uM dyad]
	double J_cam_dyadSL = 1e-3 * (k0Boff * CaM_Module_dyad.y[0] - k0Bon * Bdyad * CaM_Module_sl.y[0]); // [uM/msec dyad]
	double J_ca2cam_dyadSL = 1e-3 * (k2Boff * CaM_Module_dyad.y[1] - k2Bon * Bdyad * CaM_Module_sl.y[1]); // [uM/msec dyad]
	double J_ca4cam_dyadSL = 1e-3 * (k2Boff * CaM_Module_dyad.y[2] - k4Bon * Bdyad * CaM_Module_sl.y[2]); // [uM/msec dyad]
	double J_cam_SLmyo = kSLmyo * (CaM_Module_sl.y[0] - CaM_Module_cyto.y[0]); // [umol/msec]
	double J_ca2cam_SLmyo = kSLmyo * (CaM_Module_sl.y[1] - CaM_Module_cyto.y[1]); // [umol/msec]
	double J_ca4cam_SLmyo = kSLmyo * (CaM_Module_sl.y[2] - CaM_Module_cyto.y[2]); // [umol/msec]

	CaM_Module_dyad.ydot[0] = CaM_Module_dyad.ydot[0] - J_cam_dyadSL;
	CaM_Module_dyad.ydot[1] = CaM_Module_dyad.ydot[1] - J_ca2cam_dyadSL;
	CaM_Module_dyad.ydot[2] = CaM_Module_dyad.ydot[2] - J_ca4cam_dyadSL;
	CaM_Module_sl.ydot[0] = CaM_Module_sl.ydot[0] + J_cam_dyadSL * Vjunc / Vsl - J_cam_SLmyo / Vsl;
	CaM_Module_sl.ydot[1] = CaM_Module_sl.ydot[1] + J_ca2cam_dyadSL * Vjunc / Vsl - J_ca2cam_SLmyo / Vsl;
	CaM_Module_sl.ydot[2] = CaM_Module_sl.ydot[2] + J_ca4cam_dyadSL * Vjunc / Vsl - J_ca4cam_SLmyo / Vsl;
	CaM_Module_cyto.ydot[0] = CaM_Module_cyto.ydot[0] + J_cam_SLmyo / Vmyo;
	CaM_Module_cyto.ydot[1] = CaM_Module_cyto.ydot[1] + J_ca2cam_SLmyo / Vmyo;
	CaM_Module_cyto.ydot[2] = CaM_Module_cyto.ydot[2] + J_ca4cam_SLmyo / Vmyo;
}

void HAM_Signalling::read_initial_condition(const char* filename) {

	FILE *file;
	// open the file
	file = fopen(filename, "rb");
	if (!file) {
		// perror(filename);

		std::cerr << filename << " not opened!!! " << std::endl;
		// print_error_info_file_open_failure(filename);
		std::exit(0);
	}
// fread(array, sizeof(double), num, in);
	// output
	long int rw = fread(y, sizeof(double), ODE_NUM, file);
	if (rw != ODE_NUM) {
		std::cerr << rw << " / " << ODE_NUM << " doubles read from " << filename << ", exiting..." << std::endl;
		// exit(EXIT_FAILURE);
		std::exit(0);
	}
	if (file)
		fclose(file);
}


void HAM_Signalling::output_inital_condition(const char* filename) {

	FILE *file;
	// open the file
	file = fopen(filename, "wb");
	if (!file) {
		// perror(filename);

		std::cerr << filename << " not opened!!! " << std::endl;
		// print_error_info_file_open_failure(filename);
		std::exit(0);
	}

	// output
	long int rw = fwrite(y, sizeof(double), ODE_NUM, file);
	if (rw != ODE_NUM) {
		std::cerr << rw << " / " << ODE_NUM << " doubles written to " << filename << ", exiting..." << std::endl;
		// exit(EXIT_FAILURE);
		std::exit(0);
	}
	if (file)
		fclose(file);
}

void HAM_Signalling::assign_cell_pop_para(double *in_para) {
	// to be updated
	/*if (in_para != nullptr) {
		for (int i = 0; i < 18; ++i) {
			para[i] *= in_para[i];
		}
	}*/


	cell_para.INa_Scale            =  in_para[-4 + 4]; //     *  1; // ;// 0.991102;// 1;
	cell_para.INaL_Scale           =  in_para[-4 + 5]; //     *  1; // 1.;//0.948192;// 1.2;
	cell_para.INab_Scale           =  in_para[-4 + 6]; //     *  1; // 1.;//0.950898;// 1;//1.077453;
	cell_para.INaK_Scale           =  in_para[-4 + 7]; //     *  1; // 1;//0.932910;// 1;//0.9*1.075329;
	cell_para.Itof_Scale           =  in_para[-4 + 8]; //     *  1; // 1;//1.007522;// 1.041409;
	cell_para.IKur_Scale           =  in_para[-4 + 9]; //     *  1; // 1.07;//1.074665;// 1.028833;
	cell_para.IK2p_Scale           =  in_para[-4 + 10]; //      * 1; // 1.1;//0.984595;//  1.315084;
	cell_para.IKr_Scale            =  in_para[-4 + 11]; //      * 1; // 0.95 ; //0.934642;//  0.899033;
	cell_para.IKs_Scale            =  in_para[-4 + 12]; //      * 1; // 0.95 ; //0.946997;//  1.177648;
	cell_para.IK1_Scale            =  in_para[-4 + 13]; //      * 1; // 0.95 ; //0.960132;//  0.902799;
	cell_para.IKp_Scale            =  in_para[-4 + 14]; //      * 1; // 1.05;// 1.045922;//  1.095730;
	cell_para.IKach_Scale          =  in_para[-4 + 15]; //      * 1; // 0.95 ; //0.891492;//  0.895409;
	cell_para.ISK_Scale            =  in_para[-4 + 16]; //      * 1; // 1.05;//1.082985;//  1.099031;
	cell_para.ICaL_Scale           =  in_para[-4 + 17]; //      * 1; // 1.037031;//  1.037185;
	cell_para.ICab_Scale           =  in_para[-4 + 18]; //      * 1; // 1.073250;//  0.997408;
	cell_para.ICap_Scale           =  in_para[-4 + 19]; //      * 1; // 1.022257;//  1.020868;
	cell_para.INCX_Scale           =  in_para[-4 + 20]; //      * 1; // 1.05;//1.15;//1.101249;//  0.972194;
	cell_para.IClCa_Scale          =  in_para[-4 + 21]; //      * 1; // 1.178030;//  0.912827;
	cell_para.IClb_Scale           =  in_para[-4 + 22]; //      * 1; // 1.1;//0.956994;//  1.1;
	cell_para.Jrel_Scale           =  in_para[-4 + 23]; //      * 1; // 1.273543;//  0.980100;
	cell_para.Jserca_Scale         =  in_para[-4 + 24]; //      * 1; // 1.036176;//  1.013384;
	cell_para.Jleak_Scale          =  in_para[-4 + 25]; //      * 1;// 1.134138;//  0.899613;
	cell_para.Cleft_Buffer_Scale   =  in_para[-4 + 26]; //      *  1;//1.111454;//  1.301314;
	cell_para.Cytosol_Buffer_Scale =  in_para[-4 + 27]; //      *  1; // 1.069467;//  0.940860;
	cell_para.SR_Buffer_Scale      =  in_para[-4 + 28]; //      *  1; // 1.064127;//  0.958842;
}



void HAM_Signalling::report_cell_pop_para() {
	// to be updated
	/*if (in_para != nullptr) {
		for (int i = 0; i < 18; ++i) {
			para[i] *= in_para[i];
		}
	}*/

	std::cout
	        << cell_para.INa_Scale            << std::endl  //in_para[-4 + 4]; //     *  1; // ;// 0.991102;// 1;
	        << cell_para.INaL_Scale           << std::endl  //in_para[-4 + 5]; //     *  1; // 1.;//0.948192;// 1.2;
	        << cell_para.INab_Scale           << std::endl  //in_para[-4 + 6]; //     *  1; // 1.;//0.950898;// 1;//1.077453;
	        << cell_para.INaK_Scale           << std::endl  //in_para[-4 + 7]; //     *  1; // 1;//0.932910;// 1;//0.9*1.075329;
	        << cell_para.Itof_Scale           << std::endl  //in_para[-4 + 8]; //     *  1; // 1;//1.007522;// 1.041409;
	        << cell_para.IKur_Scale           << std::endl  //in_para[-4 + 9]; //     *  1; // 1.07;//1.074665;// 1.028833;
	        << cell_para.IK2p_Scale           << std::endl  //in_para[-4 + 10]; //      * 1; // 1.1;//0.984595;//  1.315084;
	        << cell_para.IKr_Scale            << std::endl  //in_para[-4 + 11]; //      * 1; // 0.95 ; //0.934642;//  0.899033;
	        << cell_para.IKs_Scale            << std::endl  //in_para[-4 + 12]; //      * 1; // 0.95 ; //0.946997;//  1.177648;
	        << cell_para.IK1_Scale            << std::endl  //in_para[-4 + 13]; //      * 1; // 0.95 ; //0.960132;//  0.902799;
	        << cell_para.IKp_Scale            << std::endl  //in_para[-4 + 14]; //      * 1; // 1.05;// 1.045922;//  1.095730;
	        << cell_para.IKach_Scale          << std::endl  //in_para[-4 + 15]; //      * 1; // 0.95 ; //0.891492;//  0.895409;
	        << cell_para.ISK_Scale            << std::endl  //in_para[-4 + 16]; //      * 1; // 1.05;//1.082985;//  1.099031;
	        << cell_para.ICaL_Scale           << std::endl  //in_para[-4 + 17]; //      * 1; // 1.037031;//  1.037185;
	        << cell_para.ICab_Scale           << std::endl  //in_para[-4 + 18]; //      * 1; // 1.073250;//  0.997408;
	        << cell_para.ICap_Scale           << std::endl  //in_para[-4 + 19]; //      * 1; // 1.022257;//  1.020868;
	        << cell_para.INCX_Scale           << std::endl  //in_para[-4 + 20]; //      * 1; // 1.05;//1.15;//1.101249;//  0.972194;
	        << cell_para.IClCa_Scale          << std::endl  //in_para[-4 + 21]; //      * 1; // 1.178030;//  0.912827;
	        << cell_para.IClb_Scale           << std::endl  //in_para[-4 + 22]; //      * 1; // 1.1;//0.956994;//  1.1;
	        << cell_para.Jrel_Scale           << std::endl  //in_para[-4 + 23]; //      * 1; // 1.273543;//  0.980100;
	        << cell_para.Jserca_Scale         << std::endl  //in_para[-4 + 24]; //      * 1; // 1.036176;//  1.013384;
	        << cell_para.Jleak_Scale          << std::endl  //in_para[-4 + 25]; //      * 1;// 1.134138;//  0.899613;
	        << cell_para.Cleft_Buffer_Scale   << std::endl  //in_para[-4 + 26]; //      *  1;//1.111454;//  1.301314;
	        << cell_para.Cytosol_Buffer_Scale << std::endl  //in_para[-4 + 27]; //      *  1; // 1.069467;//  0.940860;
	        << cell_para.SR_Buffer_Scale      << std::endl;  //in_para[-4 + 28]; //      *  1; // 1.064127;//  0.958842;
}
