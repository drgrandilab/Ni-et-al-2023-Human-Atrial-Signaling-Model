/*


    // 13:28:09, Mon, 11-February-2019, By Haibo
    // 17:12:20, Wed, 06-February-2019, By Haibo
    Update: to include more biomarkers of AP and CaT

    Update: accurate measurement of RMP. (or MDP most distolic Potential)


    Haibo Ni
    Fri 26 Aug 2016 17:51:50 BST


    updated: calculate plateau potentials by Haibo Ni Wed 06 Jul 2016 15:45:53 BST

    simple function to measure the APDs

    Haibo Ni
    qiangzi.ni@gmail.com
    update Jan. 27. 2015
    Wed 25 May 2016 15:10:59 BST
*/


#ifndef AP_INFOR_H
#define AP_INFOR_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
// 0.36787944117144233
//

class APInfor {
public:
	double Vmax, Vmin, timeAPDstart, timeAPDend, dVdtmax, APD_out_end[3], Vmin_last;
	double APD50_out[3], APD30_out[3], APD75_out[3], APD20_out[3];
	double Vmax_out[3], Vmin_out[3], dVdtmax_out[3], INa_max_out[3];
	double plateau_potential_out[3];
	double CaD50_out[3], CaD80_out[3], CaD90_out[3], CaD_peak_T_out[3], CaD_Tau_out[3];
	double dCadt_time_out[3];

	double simulation_interval[3];

	int APD_switch, APD_count, Vmax_switch, APD_out_swhich;
	int CaD_switch;
	double  CaD50, CaD80, CaD90, CaD_Tau, CaD_peak_T;
	double CaT_50, CaT_80, CaT_Tau, CaT_90;
	double dCadt_time, dCadt_max, CaT_prev;

	double APD90, APD50, APD30, APD20, APD75;
	double Vm_prev;
	double AMP, AMP_last, AMP_over_last;
	double APD90_prev;
	double t_since_up_stroke;
	double v90, v75, v50, v30, v20;
	double ICaL_in, RyR_in;
	double Istim_prev;
	double Diastolic_value;
	double Systolic_value_H, Systolic_value_L;
	unsigned int N_stim;
	double Diastolic_value_2, Systolic_value_L_2, Systolic_value_H_2;
	double INa_max;
	double plateau_potential;
	double dVdt_Time, V_m60_Time;
	double dVdt_Time_RP[3];
	double V_m60_Time_RP[3];
	double CaT_max, CaT_min, CaT_max_out[3], CaT_min_out[3];
	double dVdt_pre;
	double time_vmax;
	std::ofstream out;
	bool fileout;
	std::vector<float> APD_Vec;
	std::vector<float> INa_max_Vec;
	std::vector<float> dVdtmax_Vec;

	bool APD_record, CaD_record;
	int Stim_counter;


	int current_stim_counter[3];




	APInfor(const char *filename  = "APDmeasure.dat", bool file_app_mode = false, bool _fileout = false);
	~APInfor();

	void MeasureAPD90(double t, double Istim, double BCL, double dt, double Vm);
	void MeasureAPD90_INa(double t, double Istim, double BCL, double dt, double Vm, double INa, double value);
	void MeasureAPD90Using_dVdtMax(double t, double Istim, double BCL, double dt, double Vm) ;
	// void MeasureAPD90_INa(double t, double Istim, double BCL, double dt, double Vm, double INa);
	void MeasureAPD90andDSValue(double t, double Istim, double BCL, double dt, double Vm, double Value);
	void MeasureAPD90andTwoDSValue(double t, double Istim, double BCL, double dt, double Vm, double Value_1, double Value_2);
	void MeasureAPD90andDSValuewith_INa(double INa_threshold, double INa, double t, double Istim, double BCL, double dt, double Vm, double Value);
	void MeasureAPD90andDSValuewithStrokeTime(double StrokeTime, double t, double Istim, double BCL, double dt, double Vm, double Value);
	void CalciumAccumulation(double t, double Istim, double BCL, double dt, double ICaL_con, double J_RyR);
	void ReportAPD();
	void ReportLastTwo();
	void ReportLastTwo(double);
	void ReportLast();
	void ReportLastThree();
};


APInfor::APInfor(const char *filename, bool file_app_mode, bool _fileout)  {
	Vmax = 0.0;
	Vmin = 0.0;
	timeAPDstart = 0.0;
	timeAPDend = 0.0;
	dVdtmax = 0.0;
	APD_out_end[0] = 0.0;
	APD_out_end[1] = 0.0;
	dVdt_pre = 0.0;
	APD_switch = 0;
	APD_count = 0;
	Vmax_switch = 0;
	APD_out_swhich = 0;
	APD90 = APD75 = APD50 = APD30 = APD20 = 0.0;
	Vm_prev = 0.0;
	v90 = 0.0;
	v75 = 0.0;
	v50 = 0.0;
	v30 = 0.0;
	v20 = 0.0;
	ICaL_in = 0.0;
	RyR_in = 0.0;
	Istim_prev = 0.0;
	AMP = 0.0;
	AMP_last = 0.0;;
	AMP_over_last = 0.0;
	APD90_prev = 0.0;
	t_since_up_stroke = 0.0;
	plateau_potential = 0.0;
	Diastolic_value = 0.0;
	Systolic_value_L = 0.0;
	Systolic_value_H = 0.0;
	Diastolic_value_2 = 0.0;
	Systolic_value_L_2 = 0.0;
	Systolic_value_H_2 = 0.0;
	N_stim = 0;
	INa_max = 0.0;
	Vmin_last = 0.0;
	dVdt_Time = 0;
	fileout = _fileout;
	dCadt_time = 0;
	CaT_prev = -1;
	dCadt_max = 0;
	CaT_max =  CaT_min =  CaD50 =  CaD80 =  CaD90 =  CaD_Tau =  CaD_peak_T = 0;
	CaT_50 =  CaT_80 = CaT_Tau = CaT_90 = 0;

	CaD_switch = 0;

	V_m60_Time = -1;


	APD_record = false;
	CaD_record = false;
	Stim_counter = 0;


	for (int i = 0; i < 3; ++i)
	{
		APD_out_end[i] = 0;
		APD50_out[i] = APD30_out[i] = APD75_out[i] = APD20_out[i] = 0;
		Vmax_out[i] = Vmin_out[i] = dVdtmax_out[i] = INa_max_out[i] = 0;
		plateau_potential_out[i] = 0;
		CaD50_out[i] = CaD80_out[i] = CaD90_out[i] = CaD_peak_T_out[i] = CaD_Tau_out[i] = 0;
		dCadt_time_out[i] = 0;
		simulation_interval[i] = 0;

		dVdt_Time_RP[i] = 0;
		V_m60_Time_RP[i] = 0;
		CaT_max_out[i] = CaT_min_out[i] = 0;
		current_stim_counter[i] = 0;

	}


	/* are you going to output in appending mode? */
	if (fileout)
		if (file_app_mode)
		{
			out.open(filename, std::ios::out | std::ios::app);
		} else {
			out.open(filename, std::ios::out);
		}
}

APInfor::~APInfor() {
	if (out.is_open()) {
		out.close();
	}
}


void APInfor::MeasureAPD90Using_dVdtMax(double t, double Istim, double BCL, double dt, double Vm)  {

	double dVdt = (Vm - Vm_prev) / dt;

	if ((dVdt < dVdt_pre) and (dVdt_pre > 10.0) and APD_switch == 0) { // threshold
		APD_switch = 1;
		Vmax = Vm;
		timeAPDstart = t;
		// Vmin = Vm;

		/* if (Vmin > Vm)
		 {
		     Vmin = Vm; // for the first Beat... added Haibo
		 }*/
		dVdtmax = dVdt_pre;
		N_stim ++;
		Vmax_switch = 0;
		INa_max = 0.0;
		plateau_potential = 0.0;

		// std::cerr << "ttttt" << std::endl;
	}




	if (APD_switch == 0)
	{
		if (Vmin >= Vm) {
			Vmin = Vm;
		}
	}


	if (APD_switch == 1) {
		if (Vm > -60.0) {
			if (Vm >= Vmax) {
				Vmax = Vm;
				time_vmax = t;
			}
			else {

				if ( t - time_vmax  > 8.0)
				{
					Vmax_switch = 1;
				}
			}
			if (Vmax_switch == 1) {
				APD_switch = 2;

				AMP_last = AMP;
				AMP = Vmax - Vmin;
				AMP_over_last = AMP / AMP_last;
				AMP_over_last = t;
				/*printf("%f %f\n", Vmax, Vmin);
				printf("%f %f\n", AMP, AMP_last);*/
				v90 = Vmax - 0.9 * (Vmax - Vmin);
				v75 = Vmax - 0.75 * (Vmax - Vmin);
				v50 = Vmax - 0.50 * (Vmax - Vmin);
				v20 = Vmax - 0.20 * (Vmax - Vmin);

				Vmin_last = Vmin;
				Vmin = 0;


				Vmax_switch = 0;

				APD_Vec.push_back(-100);

				// INa_max_out.push_back(0);
				dVdtmax_Vec.push_back(0);
				// printf("asfde\n");
			}
		}
	}

	/*

	    plateau_potential is defined as the average potential
	   between 10 and 50 ms following application of stimulation.
	   Haibo Ni.
	*/

	if (t > timeAPDstart + 10 and t < timeAPDstart + 50 )
	{
		plateau_potential = plateau_potential + dt * Vm;
	}

	if (APD_switch == 2) {

		if ((Vm_prev >= v20) && (Vm <= v20) ) {
			APD20 = t - timeAPDstart ;
		}
		else if ((Vm_prev >= v50) && (Vm <= v50 )) {
			APD50 = t - timeAPDstart ;
		}
		else if (Vm_prev >= v75 && Vm <= v75 ) {
			APD75 = t - timeAPDstart ;
		}
		else if (Vm_prev >= v90 && Vm <= v90) {
			APD_switch = 0;
			APD_count ++;

			APD90 = t - timeAPDstart;
			plateau_potential /= (50.0 - 10.0); // average
			// printf("dddddddddddd\n");
			if (fileout)
				out << t << " "
				    << APD90 << " "
				    << APD75 << " "
				    << APD50 << " "
				    << APD20 << " "
				    << Vmax  << " "
				    << Vmin_last  << " "
				    << dVdtmax  << " "
				    << plateau_potential  << " "
				    << AMP_over_last  << " "
				    << std::endl;

			if (APD_out_swhich == 0) {
				APD_out_end[0] = t - timeAPDstart;
				APD75_out[0] = APD75;
				APD50_out[0] = APD50;
				APD20_out[0] = APD20;
				Vmax_out[0] = Vmax;
				Vmin_out[0] = Vmin_last;
				dVdtmax_out[0] = dVdtmax;
				APD_out_swhich = 1;
				plateau_potential_out[0] = plateau_potential;
			}
			else if (APD_out_swhich == 1) {
				APD_out_end[1] = t - timeAPDstart;
				APD75_out[1]   = APD75;
				APD50_out[1]   = APD50;
				APD20_out[1]   = APD20;
				Vmax_out[1]    = Vmax;
				Vmin_out[1]    = Vmin_last;
				dVdtmax_out[1] = dVdtmax;
				plateau_potential_out[1] = plateau_potential;
				APD_out_swhich = 0;
			}
			APD_Vec.back() = APD90;
			dVdtmax_Vec.back() = dVdtmax;
		}
	}


	if (dVdt > dVdtmax) dVdtmax = dVdt;

	Vm_prev = Vm;
	Istim_prev = Istim;
	dVdt_pre = dVdt;
}

void APInfor::MeasureAPD90(double t, double Istim, double BCL, double dt, double Vm)  {


	if (((Istim >= 1e-3 || Istim <= -1e-3 ) && (fabs(Istim_prev) < 1e-10))) {
		APD_switch = 1;
		Vmax = Vm;
		timeAPDstart = t;
		// Vmin = Vm;

		if (Vmin > Vm)
		{
			Vmin = Vm; // for the first Beat... added Haibo
		}
		dVdtmax = 0.0;
		N_stim ++;
		Vmax_switch = 0;
		INa_max = 0.0;
		plateau_potential = 0.0;
	}


	if (APD_switch == 0)
	{
		if (Vmin >= Vm) {
			Vmin = Vm;
		}
	}


	if (APD_switch == 1) {
		if (Vm > -60.0) {
			if (Vm >= Vmax) {
				Vmax = Vm;
			}
			else Vmax_switch = 1;
			if (Vmax_switch == 1) {
				APD_switch = 2;
				AMP_last = AMP;
				AMP = Vmax - Vmin;
				AMP_over_last = AMP / AMP_last;
				/*printf("%f %f\n", Vmax, Vmin);
				printf("%f %f\n", AMP, AMP_last);*/
				v90 = Vmax - 0.9 * (Vmax - Vmin);
				v75 = Vmax - 0.75 * (Vmax - Vmin);
				v50 = Vmax - 0.50 * (Vmax - Vmin);
				v30 = Vmax - 0.30 * (Vmax - Vmin);

				Vmin_last = Vmin;
				Vmin = 0;
				// printf("asfde\n");
			}
		}
	}

	/*

	    plateau_potential is defined as the average potential
	   between 10 and 50 ms following application of stimulation.
	   Haibo Ni.
	*/

	if (t > timeAPDstart + 10 and t < timeAPDstart + 50 )
	{
		plateau_potential = plateau_potential + dt * Vm;
	}

	if (APD_switch == 2) {

		if ((Vm_prev >= v30) && (Vm <= v30) ) {
			APD30 = t - timeAPDstart ;
		}
		else if ((Vm_prev >= v50) && (Vm <= v50 )) {
			APD50 = t - timeAPDstart ;
		}
		else if (Vm_prev >= v75 && Vm <= v75 ) {
			APD75 = t - timeAPDstart ;
		}
		else if (Vm_prev >= v90 && Vm <= v90) {
			APD_switch = 0;
			APD_count ++;
			Vmax_switch = 0;
			APD90 = t - timeAPDstart;
			plateau_potential /= (50.0 - 10.0); // average
			// printf("dddddddddddd\n");
			if (fileout)
				out << APD_count *BCL << " "
				    << APD90 << " "
				    << APD75 << " "
				    << APD50 << " "
				    << APD30 << " "
				    << Vmax  << " "
				    << Vmin_last  << " "
				    << dVdtmax  << " "
				    << plateau_potential  << " "
				    << AMP_over_last  << " "
				    << std::endl;

			if (APD_out_swhich == 0) {
				APD_out_end[0] = t - timeAPDstart;
				APD75_out[0] = APD75;
				APD50_out[0] = APD50;
				APD30_out[0] = APD30;
				Vmax_out[0] = Vmax;
				Vmin_out[0] = Vmin_last;
				dVdtmax_out[0] = dVdtmax;
				APD_out_swhich = 1;
				plateau_potential_out[0] = plateau_potential;
			}
			else if (APD_out_swhich == 1) {
				APD_out_end[1] = t - timeAPDstart;
				APD75_out[1]   = APD75;
				APD50_out[1]   = APD50;
				APD30_out[1]   = APD30;
				Vmax_out[1]    = Vmax;
				Vmin_out[1]    = Vmin_last;
				dVdtmax_out[1] = dVdtmax;
				plateau_potential_out[1] = plateau_potential;
				APD_out_swhich = 0;
			}
			APD_Vec.push_back(APD90);
		}
	}

	double dVdt = (Vm - Vm_prev) / dt;
	if (dVdt > dVdtmax) dVdtmax = dVdt;

	Vm_prev = Vm;
	Istim_prev = Istim;
}

void APInfor::MeasureAPD90_INa(double t, double Istim, double BCL, double dt, double Vm, double INa, double CaT) {

	if (((Istim >= 1e-3 || Istim <= -1e-3 ) && (fabs(Istim_prev) < 1e-10))) {


		// check if APD values have been recorded or not.
		// if the previous AP did not reach APD90, then they are not.
		// hence APD values are recorded here.
		plateau_potential /= (50.0 - 10.0); // average


		APD_out_swhich =  Stim_counter % 3;
		APD_out_end[APD_out_swhich] = APD90;
		APD75_out[APD_out_swhich] = APD75;
		APD50_out[APD_out_swhich] = APD50;
		APD30_out[APD_out_swhich] = APD30;
		APD20_out[APD_out_swhich] = APD20;
		Vmax_out[APD_out_swhich] = Vmax;
		Vmin_out[APD_out_swhich] = Vmin_last;
		dVdtmax_out[APD_out_swhich] = dVdtmax;
		plateau_potential_out[APD_out_swhich] = plateau_potential;
		INa_max_out[APD_out_swhich] = INa_max;
		dVdt_Time_RP[APD_out_swhich] = dVdt_Time;
		V_m60_Time_RP[APD_out_swhich] = V_m60_Time;
		CaT_max_out[APD_out_swhich] = CaT_max;
		CaT_min_out[APD_out_swhich] = CaT_min;
		CaD50_out[APD_out_swhich] = CaD50;
		CaD80_out[APD_out_swhich] = CaD80;
		CaD90_out[APD_out_swhich] = CaD90;
		CaD_Tau_out[APD_out_swhich] = CaD_Tau;
		CaD_peak_T_out[APD_out_swhich] = CaD_peak_T;
		dCadt_time_out[APD_out_swhich] = dCadt_time;
		current_stim_counter[APD_out_swhich] = Stim_counter;

		if (out.is_open())
			out << std::setprecision(9) << current_stim_counter[APD_out_swhich] << " "
			    << simulation_interval[APD_out_swhich] << " "
			    << APD_out_end[APD_out_swhich] << " "
			    << APD75_out[APD_out_swhich] << " "
			    << APD50_out[APD_out_swhich] << " "
			    << APD30_out[APD_out_swhich] << " "
			    << APD20_out[APD_out_swhich] << " "
			    << Vmax_out[APD_out_swhich] << " "
			    << Vmin_out[APD_out_swhich] << " "
			    << dVdtmax_out[APD_out_swhich] << " "
			    << plateau_potential_out[APD_out_swhich] << " "
			    << INa_max_out[APD_out_swhich] << " "
			    << dVdt_Time_RP[APD_out_swhich] << " "
			    << V_m60_Time_RP[APD_out_swhich] << " "
			    << CaT_max_out[APD_out_swhich] << " "
			    << CaT_min_out[APD_out_swhich] << " "
			    << CaD50_out[APD_out_swhich]       << " "
			    << CaD80_out[APD_out_swhich]       << " "
			    << CaD90_out[APD_out_swhich]       << " "
			    << CaD_Tau_out[APD_out_swhich]     << " "
			    << CaD_peak_T_out[APD_out_swhich]  << " "
			    << dCadt_time_out[APD_out_swhich]  << " "
			    << std::endl;;



		// if (out.is_open())
		// 	out << std::setprecision(9) << APD_count *BCL << " "
		// 	    << APD90 << " "
		// 	    << APD75 << " "
		// 	    << APD50 << " "
		// 	    << APD30 << " "
		// 	    << APD20 << " "
		// 	    << Vmax  << " "
		// 	    << Vmin_last  << " "
		// 	    << dVdtmax  << " "
		// 	    << plateau_potential  << " "
		// 	    << INa_max << " "
		// 	    << dVdt_Time  << " "
		// 	    << V_m60_Time << " "
		// 	    << AMP_over_last  << " "
		// 	    << std::endl;

		APD_record  = false;
		CaD_record  = false;
		Stim_counter ++;
		simulation_interval[Stim_counter % 3] = t - timeAPDstart;




		APD_switch = 1;
		Vmax = Vm;
		timeAPDstart = t;
		dCadt_time = 0;
		// Vmin = Vm;

		if (Vmin > Vm)
		{
			Vmin = Vm; // for the first Beat... added Haibo
		}

		CaT_max = -100;
		CaT_min = 100;
		if (CaT_min > CaT)
			CaT_min = CaT;
		dVdtmax = 0.0;
		N_stim ++;
		Vmax_switch = 0;
		INa_max = 0.0;
		plateau_potential = 0.0;
		APD90 = -1.0;
		APD75 = -1.0;
		APD50 = -1.0;
		APD30 = -.0;
		APD20 = -1.0;
		Vmax = -100;
		Vmin_last = 0.0;
		dVdtmax = 0.0;
		INa_max = 0.0;
		dVdt_Time = -1;
		V_m60_Time = -1;
		dCadt_max = 0;

		CaD50 =  CaD80 =  CaD90 =  CaD_Tau =  CaD_peak_T = -100;
		CaT_50 =  CaT_80 = CaT_Tau = CaT_90 = -100.0;

		CaD_switch = 1;

	}


	if (APD_switch == 0)
	{
		if (Vmin >= Vm) {
			Vmin = Vm;
		}
	}


	if (CaD_switch == 0)
	{
		if (CaT_min >= CaT) {
			CaT_min = CaT;
		}
	}




	/*if (CaT_max < CaT)
	{
	    CaT_max = CaT;
	}
	if (CaT_min > CaT)
	{
	    CaT_min = CaT;
	}*/


	if (APD_switch == 2) {

		if (Vm >= Vmax) {
			Vmax = Vm;
			// APD_switch =1;
		}
	}

	if (CaD_switch == 2) {

		if (CaT_max <= CaT) {
			CaT_max = CaT;
		}
	}



	if (CaD_switch == 1) {

		double dCadt = (CaT - CaT_prev) / dt;
		if (dCadt > dCadt_max) {
			dCadt_max = dCadt;
			dCadt_time = t;
		}

		if (CaT > 0.00001) {

			if (CaT >= CaT_max) {
				CaT_max = CaT;
				CaD_peak_T = t - timeAPDstart;

			}
			else if (t > timeAPDstart + 50) {
				CaD_switch = 2;
				// AMP_last = AMP;
				// AMP = Vmax - Vmin;

				// AMP_over_last = AMP / AMP_last;



				/*printf("%f %f\n", Vmax, Vmin);
				printf("%f %f\n", AMP, AMP_last);*/
				CaT_90 = CaT_max - 0.9 * (CaT_max - CaT_min);
				CaT_Tau = CaT_max - (1 - 0.36787944117144233) * (CaT_max - CaT_min);
				CaT_50 = CaT_max - 0.50 * (CaT_max - CaT_min);
				CaT_80 = CaT_max - 0.80 * (CaT_max - CaT_min);
				// v20 = Vmax - 0.20 * (Vmax - Vmin);
				// Vmin_last = Vmin;
				// Vmin = 0;
				// printf("asfde\n");
			}
		}
	}





	if (APD_switch == 1) {

		double dVdt = (Vm - Vm_prev) / dt;
		if (dVdt > dVdtmax) {
			dVdtmax = dVdt;
			dVdt_Time = t;

		}


		if (INa_max > INa)
		{
			INa_max = INa;
			// std::cout << "yes" << std::endl;
		}

		if (Vm > -75.0) {

			if (Vm >= Vmax) {
				Vmax = Vm;
			}
			else if (t > timeAPDstart + 3) { //
				Vmax_switch = 1;
			}
			if (Vmax_switch == 1) {
				APD_switch = 2;
				AMP_last = AMP;
				AMP = Vmax - Vmin;
				AMP_over_last = AMP / AMP_last;
				/*printf("%f %f\n", Vmax, Vmin);
				printf("%f %f\n", AMP, AMP_last);*/
				v90 = Vmax - 0.9 * (Vmax - Vmin);
				v75 = Vmax - 0.75 * (Vmax - Vmin);
				v50 = Vmax - 0.50 * (Vmax - Vmin);
				v30 = Vmax - 0.30 * (Vmax - Vmin);
				v20 = Vmax - 0.20 * (Vmax - Vmin);
				// std::cout << v20 << std::endl;
				Vmin_last = Vmin;
				Vmin = 0;
				// printf("asfde\n");
			}
		}
	}

	/*

	   plateau_potential is defined as the average potential
	   between 10 and 50 ms following application of stimulation.
	   Haibo Ni.
	*/

	if (t > timeAPDstart + 10 and t < timeAPDstart + 50 )
	{
		plateau_potential = plateau_potential + dt * Vm;
	}



	if ((CaT_prev >= CaT_50) && (CaT <= CaT_50)  and  (CaD_switch == 2)) {
		CaD50 = t - timeAPDstart ;
		CaD_switch = 3;
	} else if ((CaD_switch == 3) and (CaT_prev >= CaT_Tau) && (CaT <= CaT_Tau) ) {
		CaD_Tau = t - (CaD_peak_T + timeAPDstart) ;
		CaD_switch = 4;
	} else if ((CaD_switch == 4) and (CaT_prev >= CaT_80) && (CaT <= CaT_80) ) {
		CaD80 = t - timeAPDstart ;
		CaD_switch = 5;
	}
	else if ((CaD_switch == 5) and (CaT_prev >= CaT_90) && (CaT <= CaT_90) ) {
		CaD90 = t - timeAPDstart ;
		CaD_switch = 0;
	}

	// std::cout << 'a' << std::endl;
	if ( APD_switch == 2 and (Vm_prev >= v20) and (Vm <= v20) ) {
		APD20 = t - dVdt_Time;
		APD_switch = 3;
	} else if ((APD_switch == 3) and (Vm_prev >= v30) and (Vm <= v30) ) {
		APD30 = t - dVdt_Time;
		APD_switch = 4;
	}
	else if ((APD_switch == 4) and (Vm_prev >= v50) and (Vm <= v50 )) {
		APD50 = t - dVdt_Time ;
		APD_count ++;
		Vmax_switch = 0;
		APD_switch = 5;
	}
	else if ((APD_switch == 5) and (Vm_prev >= v75) and (Vm <= v75) ) {
		APD75 = t - dVdt_Time ;
	}
	else if ((APD_switch == 5) and (Vm_prev >= v90 and Vm <= v90)) {
		APD_switch = 0;
		APD90 = t - dVdt_Time;

		// simulation_interval[APD_out_swhich] = t - timeAPDstart;

		// CaD50_out[3], CaD80_out[3], CaD90_out[3], CaD_peak_T_out[3], CaD_Tau_out[3]


	}

	if (fabs(Vm - (-60)) < 0.1) {
		V_m60_Time = t - dVdt_Time;
	}

	Vm_prev = Vm;
	CaT_prev = CaT;
	Istim_prev = Istim;
}






void APInfor::MeasureAPD90andDSValue(double t, double Istim, double BCL, double dt, double Vm, double Value)  {

	if (((Istim >= 1e-8 || Istim <= -1e-8 ) && (fabs(Istim_prev) < 1e-10) && (APD_switch == 0))) {
		if (APD_count > 0)
		{
			if (fileout)
				out << APD_count *BCL << " "
				    << BCL << " "
				    << APD90 << " "
				    << APD75 << " "
				    << APD50 << " "
				    << APD30 << " "
				    << Systolic_value_L << " "
				    << Systolic_value_H << " "

				    << Vmax << " "
				    << Vmin << " "
				    << dVdtmax << std::endl;

		}
		APD_switch = 1;
		Vmax = -60;
		timeAPDstart = t;
		Vmin = Vm;
		dVdtmax = 0.0;
		Diastolic_value = Value;
		Systolic_value_L = Value;
		Systolic_value_H = Value;
		N_stim ++;
		Vmax_switch = 0;

		// printf("ssssssssss\n");
	}
	if (Systolic_value_L > Value)
	{
		Systolic_value_L = Value;
	}
	if (Systolic_value_H < Value)
	{
		Systolic_value_H = Value;
	}
	if (APD_switch == 1) {
		if (Vm > -80.0) {
			if (Vm >= Vmax) {
				Vmax = Vm;
			}
			else Vmax_switch = 1;
			if (Vmax_switch == 1) {
				APD_switch = 2;
				AMP_last = AMP;
				AMP = Vmax - Vmin;
				AMP_over_last = AMP / AMP_last;
				/*printf("%f %f\n", Vmax, Vmin);

				printf("%f %f\n", AMP, AMP_last);*/
				v90 = Vmax - 0.9 * (Vmax - Vmin);
				v75 = Vmax - 0.75 * (Vmax - Vmin);
				v50 = Vmax - 0.50 * (Vmax - Vmin);
				v30 = Vmax - 0.30 * (Vmax - Vmin);
				// printf("asfde\n");
			}
		}
	}
	if (APD_switch == 2) {

		if ((Vm_prev >= v30) && (Vm <= v30) ) {
			APD30 = t - timeAPDstart ;
			// printf("aaaa\n");

		}
		else if ((Vm_prev >= v50) && (Vm <= v50 )) {
			APD50 = t - timeAPDstart ;
			// printf("bbbbbbbbb\n");

		}
		else if (Vm_prev >= v75 && Vm <= v75 ) {
			APD75 = t - timeAPDstart ;
			// printf("ccccccccccc\n");

		}
		else if (Vm_prev >= v90 && Vm <= v90) {
			APD_switch = 0;
			APD_count ++;
			Vmax_switch = 0;
			APD90 = t - timeAPDstart;
			APD_Vec.push_back(APD90);
			// printf("dddddddddddd\n");
			// printf("%f %f %f\n", Diastolic_value, Systolic_value_L, Systolic_value_H);

			// fprintf (out, "%.2f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", APD_count * BCL, APD90, APD75, APD50, APD30, Vmax, Vmin, dVdtmax);
			if (APD_out_swhich == 0) {
				APD_out_end[0] = t - timeAPDstart;
				APD_out_swhich = 1;
			}
			else if (APD_out_swhich == 1) {
				APD_out_end[1] = t - timeAPDstart;
				APD_out_swhich = 0;
			}
		}
	}

	double dVdt = (Vm - Vm_prev) / dt;
	if (dVdt > dVdtmax) dVdtmax = dVdt;
	Istim_prev = Istim;
	Vm_prev = Vm;
}
void APInfor::MeasureAPD90andTwoDSValue(double t, double Istim, double BCL, double dt, double Vm, double Value_1, double Value_2)  {

	if (((Istim >= 1e-8 || Istim <= -1e-8 ) && (fabs(Istim_prev) < 1e-10) && (APD_switch == 0))) {
		if (APD_count > 70)
		{
			if (fileout)
				out << APD_count *BCL << " "
				    << BCL << " "
				    << APD90 << " "
				    << APD75 << " "
				    << APD50 << " "
				    << APD30 << " "
				    << Systolic_value_L << " "
				    << Systolic_value_H << " "
				    << Systolic_value_L_2 << " "
				    << Systolic_value_H_2 << " "
				    << Vmax << " "
				    << Vmin << " "
				    << dVdtmax << std::endl;
		}

		APD_switch = 1;
		Vmax = -60;
		timeAPDstart = t;
		Vmin = Vm;
		dVdtmax = 0.0;
		Diastolic_value = Value_1;
		Systolic_value_L = Value_1;
		Systolic_value_H = Value_1;
		Diastolic_value_2 = Value_2;
		Systolic_value_L_2 = Value_2;
		Systolic_value_H_2 = Value_2;
		N_stim ++;
		Vmax_switch = 0;

		// printf("ssssssssss\n");
	}
	if (Systolic_value_L > Value_1)
	{
		Systolic_value_L = Value_1;
	}
	if (Systolic_value_H < Value_1)
	{
		Systolic_value_H = Value_1;
	}
	if (Systolic_value_L_2 > Value_2)
	{
		Systolic_value_L_2 = Value_2;
	}
	if (Systolic_value_H_2 < Value_2)
	{
		Systolic_value_H_2 = Value_2;
	}
	if (APD_switch == 1) {
		if (Vm > -30.0) {
			if (Vm >= Vmax) {
				Vmax = Vm;
			}
			else Vmax_switch = 1;
			if (Vmax_switch == 1) {
				APD_switch = 2;
				AMP_last = AMP;
				AMP = Vmax - Vmin;
				AMP_over_last = AMP / AMP_last;
				/*printf("%f %f\n", Vmax, Vmin);

				printf("%f %f\n", AMP, AMP_last);*/
				v90 = Vmax - 0.9 * (Vmax - Vmin);
				v75 = Vmax - 0.75 * (Vmax - Vmin);
				v50 = Vmax - 0.50 * (Vmax - Vmin);
				v30 = Vmax - 0.30 * (Vmax - Vmin);
				// printf("asfde\n");
			}
		}
	}
	if (APD_switch == 2) {

		if ((Vm_prev >= v30) && (Vm <= v30) ) {
			APD30 = t - timeAPDstart ;
			// printf("aaaa\n");

		}
		else if ((Vm_prev >= v50) && (Vm <= v50 )) {
			APD50 = t - timeAPDstart ;
			// printf("bbbbbbbbb\n");

		}
		else if (Vm_prev >= v75 && Vm <= v75 ) {
			APD75 = t - timeAPDstart ;
			// printf("ccccccccccc\n");

		}
		else if (Vm_prev >= v90 && Vm <= v90) {
			APD_switch = 0;
			APD_count ++;
			Vmax_switch = 0;
			APD90 = t - timeAPDstart;
			APD_Vec.push_back(APD90);
			// printf("dddddddddddd\n");
			// printf("%f %f %f\n", Diastolic_value, Systolic_value_L, Systolic_value_H);
			// fprintf (out, "%.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", APD_count * BCL, APD90, APD75, APD50, APD30, Vmax, Vmin, dVdtmax);
			if (APD_out_swhich == 0) {
				APD_out_end[0] = t - timeAPDstart;
				APD_out_swhich = 1;
			}
			else if (APD_out_swhich == 1) {
				APD_out_end[1] = t - timeAPDstart;
				APD_out_swhich = 0;
			}
		}
	}

	double dVdt = (Vm - Vm_prev) / dt;
	if (dVdt > dVdtmax) dVdtmax = dVdt;
	Istim_prev = Istim;
	Vm_prev = Vm;
}


/*INa is used to measure the upstroke */
void APInfor::MeasureAPD90andDSValuewith_INa(double INa_threshold, double INa, double t, double Istim, double BCL, double dt, double Vm, double Value)  {

	if (((Istim >= 1e-3 || Istim <= -1e-3 ) && (fabs(Istim_prev) < 1e-10))) {
		if (APD_count > 0)
		{
			// printf("dddddddddddd\n");
			// printf("%f %f %f\n", Diastolic_value, Systolic_value_L, Systolic_value_H);

			// fprintf (out, "%.2f %.3f %.3f %.3f %.3f %.3f %.3f  %.3f %.3f %.3f %.3f\n", APD_count * BCL, BCL, APD90, APD75, APD50, APD30, Systolic_value_L, Systolic_value_H, Vmax, Vmin, dVdtmax);
			if (fileout)
				out << APD_count *BCL << " "
				    << BCL << " "
				    << APD90 << " "
				    << APD75 << " "
				    << APD50 << " "
				    << APD30 << " "
				    << Systolic_value_L << " "
				    << Systolic_value_H << " "
				    << Vmax << " "
				    << Vmin << " "
				    << dVdtmax << std::endl;
		}
		APD_switch = 1;
		Vmax = -60;
		timeAPDstart = t;
		Vmin = Vm;
		dVdtmax = 0.0;
		Diastolic_value = Value;
		Systolic_value_L = Value;
		Systolic_value_H = Value;
		N_stim ++;
		Vmax_switch = 0;

		// printf("ssssssssss\n");
	}
	if (Systolic_value_L > Value)
	{
		Systolic_value_L = Value;
	}
	if (Systolic_value_H < Value)
	{
		Systolic_value_H = Value;
	}
	if (APD_switch == 1) {
		if (Vm > -30.0) {
			if ((Vm > Vmax) && (INa <= INa_threshold)) {
				Vmax = Vm;
			}
			else Vmax_switch = 1;
			if (Vmax_switch == 1) {
				APD_switch = 2;
				AMP_last = AMP;
				AMP = Vmax - Vmin;
				AMP_over_last = AMP / AMP_last;
				/*printf("%f %f\n", AMP, AMP_last);
				printf("%f %f\n", Vmax, Vmin);*/

				v90 = Vmax - 0.9 * (Vmax - Vmin);
				v75 = Vmax - 0.75 * (Vmax - Vmin);
				v50 = Vmax - 0.50 * (Vmax - Vmin);
				v30 = Vmax - 0.30 * (Vmax - Vmin);
				// printf("asfde\n");
			}
		}
	}
	if (APD_switch == 2) {

		if ((Vm_prev >= v30) && (Vm <= v30) ) {
			APD30 = t - timeAPDstart ;
			// printf("aaaa\n");

		}
		else if ((Vm_prev >= v50) && (Vm <= v50 )) {
			APD50 = t - timeAPDstart ;
			// printf("bbbbbbbbb\n");

		}
		else if (Vm_prev >= v75 && Vm <= v75 ) {
			APD75 = t - timeAPDstart ;
			// printf("ccccccccccc\n");

		}
		else if (Vm_prev >= v90 && Vm <= v90) {
			APD_switch = 0;
			APD_count ++;
			Vmax_switch = 0;
			APD90 = t - timeAPDstart;
			if (APD_out_swhich == 0) {
				APD_out_end[0] = t - timeAPDstart;
				APD_out_swhich = 1;
			}
			else if (APD_out_swhich == 1) {
				APD_out_end[1] = t - timeAPDstart;
				APD_out_swhich = 0;
			}
		}
	}

	double dVdt = (Vm - Vm_prev) / dt;
	if (dVdt > dVdtmax) dVdtmax = dVdt;
	Istim_prev = Istim;
	Vm_prev = Vm;
}

void APInfor::MeasureAPD90andDSValuewithStrokeTime(double StrokeTime, double t, double Istim, double BCL, double dt, double Vm, double Value)  {

	if (((Istim >= 1e-3 || Istim <= -1e-3 ) && (fabs(Istim_prev) < 1e-10))) {
		APD_switch = 1;
		Vmax = -60;
		timeAPDstart = t;
		Vmin = Vm;
		dVdtmax = 0.0;
		Diastolic_value = Value;
		Systolic_value_L = Value;
		Systolic_value_H = Value;
		N_stim ++;
		t_since_up_stroke = 0.0;
		Vmax_switch = 0;

		// printf("ssssssssss\n");
	}
	if (Systolic_value_L > Value)
	{
		Systolic_value_L = Value;
	}
	if (Systolic_value_H < Value)
	{
		Systolic_value_H = Value;
	}
	t_since_up_stroke += dt;

	if (APD_switch == 1) {

		if (Vm > -10.0) {
			if ((Vm > Vmax) && (t_since_up_stroke <= StrokeTime)) {
				Vmax = Vm;
				// printf("%f\n", t_since_up_stroke);
				// printf("%f\n", Vmax);

			}
			else Vmax_switch = 1;
			if (Vmax_switch == 1) {
				APD_switch = 2;
				AMP_last = AMP;
				AMP = Vmax - Vmin;
				AMP_over_last = AMP / AMP_last;
				// printf("%f %f\n", AMP, AMP_last);
				// printf("%f %f\n", Vmax, Vmin);

				v90 = Vmax - 0.9 * (Vmax - Vmin);
				v75 = Vmax - 0.75 * (Vmax - Vmin);
				v50 = Vmax - 0.50 * (Vmax - Vmin);
				v30 = Vmax - 0.30 * (Vmax - Vmin);
				// printf("asfde\n");
			}
		}
	}
	if (APD_switch == 2) {

		if ((Vm_prev >= v30) && (Vm <= v30) ) {
			APD30 = t - timeAPDstart ;
			// printf("aaaa\n");

		}
		else if ((Vm_prev >= v50) && (Vm <= v50 )) {
			APD50 = t - timeAPDstart ;
			// printf("bbbbbbbbb\n");

		}
		else if (Vm_prev >= v75 && Vm <= v75 ) {
			APD75 = t - timeAPDstart ;
			// printf("ccccccccccc\n");

		}
		else if (Vm_prev >= v90 && Vm <= v90) {
			APD_switch = 0;
			APD_count ++;
			Vmax_switch = 0;
			APD90 = t - timeAPDstart;
			APD_Vec.push_back(APD90);
			if (fileout)
				out << APD_count *BCL << " "
				    << APD90 << " "
				    << APD75 << " "
				    << APD50 << " "
				    << APD30 << " "
				    << Vmax  << " "
				    << Vmin  << " "
				    << dVdtmax << std::endl;

			if (APD_out_swhich == 0) {
				APD_out_end[0] = t - timeAPDstart;
				APD_out_swhich = 1;
			}
			else if (APD_out_swhich == 1) {
				APD_out_end[1] = t - timeAPDstart;
				APD_out_swhich = 0;
			}
		}
	}

	double dVdt = (Vm - Vm_prev) / dt;
	if (dVdt > dVdtmax) dVdtmax = dVdt;
	Istim_prev = Istim;
	Vm_prev = Vm;
}

/* measure the accumulations of calcium build up from the ICaL and RyR */
void APInfor::CalciumAccumulation(double t, double Istim, double BCL, double dt, double ICaL_con, double J_RyR) {
	if ((Istim >= 1e-3 || Istim <= -1e-3 ) && fabs(Istim_prev) < 1e-3) {
		if (fileout)
			out << t << " "
			    << ICaL_in << " "
			    << RyR_in  << std::endl;
		timeAPDstart = t;
		ICaL_in = 0.0;
		RyR_in = 0.0;
	}

	ICaL_in += ICaL_con * dt;
	RyR_in += J_RyR * dt;
	Istim_prev = Istim;
}

void APInfor::ReportAPD() {
	std::cout <<  "APD_out_end[0] = " <<  APD_out_end[0] << std::endl;
	std::cout <<  "APD_out_end[1] = " <<  APD_out_end[1] << std::endl;
}


void APInfor::ReportLastTwo() {
	std::cout /*<< BCL << "\t"*/
	        << std::setprecision(9) << APD_out_end[0] << "\t"
	        << APD75_out[0] << "\t"
	        << APD50_out[0] << "\t"
	        << APD30_out[0] << "\t"
	        << Vmax_out[0] << "\t"
	        << Vmin_out[0] << "\t"
	        << dVdtmax_out[0] << "\t"
	        << plateau_potential_out[0] << "\t"
	        << INa_max_out[0] << "\t"
	        << dVdt_Time_RP[0] << "\t"
	        /* std::cout*/ /*<< BCL << "\t"*/
	        << APD_out_end[1] << "\t"
	        << APD75_out[1] << "\t"
	        << APD50_out[1] << "\t"
	        << APD30_out[1] << "\t"
	        << Vmax_out[1] << "\t"
	        << Vmin_out[1] << "\t"
	        << dVdtmax_out[1] << "\t"
	        << plateau_potential_out[1] << "\t"
	        << INa_max_out[1] << "\t"
	        << dVdt_Time_RP[1] << "\t"
	        << V_m60_Time_RP[1]
	        << std::endl;
}

void APInfor::ReportLastTwo(double BCL) {
	std::cout << std::setprecision(9) << simulation_interval[0] << " "
	          << APD_out_end[0] << " "
	          << APD75_out[0] << " "
	          << APD50_out[0] << " "
	          << APD30_out[0] << " "
	          << APD20_out[0] << " "
	          << Vmax_out[0] << " "
	          << Vmin_out[0] << " "
	          << dVdtmax_out[0] << " "
	          << plateau_potential_out[0] << " "
	          << INa_max_out[0] << " "
	          << dVdt_Time_RP[0] << " "
	          << V_m60_Time_RP[0] << " "
	          << CaT_max_out[0] << " "
	          << CaT_min_out[0] << " "
	          << CaD50_out[0]       << " "
	          << CaD80_out[0]       << " "
	          << CaD90_out[0]       << " "
	          << CaD_Tau_out[0]     << " "
	          << CaD_peak_T_out[0]  << " "
	          << dCadt_time_out[0]  << " "

	          /* std::cout*/ << simulation_interval[1] << " "
	          << APD_out_end[1] << " "
	          << APD75_out[1] << " "
	          << APD50_out[1] << " "
	          << APD30_out[1] << " "
	          << APD20_out[1] << " "
	          << Vmax_out[1] << " "
	          << Vmin_out[1] << " "
	          << dVdtmax_out[1] << " "
	          << plateau_potential_out[1] << " "
	          << INa_max_out[1] << " "
	          << dVdt_Time_RP[1] << " "
	          << V_m60_Time_RP[1] << " "
	          << CaT_max_out[1] << " "
	          << CaT_min_out[1] << " "
	          << CaD50_out[1]       << " "
	          << CaD80_out[1]       << " "
	          << CaD90_out[1]       << " "
	          << CaD_Tau_out[1]     << " "
	          << CaD_peak_T_out[1]  << " "
	          << dCadt_time_out[1]  << " "
	          << std::endl;
}




void APInfor::ReportLastThree() {


	plateau_potential /= (50.0 - 10.0); // average

	APD_out_swhich =  Stim_counter % 3;
	APD_out_end[APD_out_swhich] = APD90;
	APD75_out[APD_out_swhich] = APD75;
	APD50_out[APD_out_swhich] = APD50;
	APD30_out[APD_out_swhich] = APD30;
	APD20_out[APD_out_swhich] = APD20;
	Vmax_out[APD_out_swhich] = Vmax;
	Vmin_out[APD_out_swhich] = Vmin_last;
	dVdtmax_out[APD_out_swhich] = dVdtmax;
	plateau_potential_out[APD_out_swhich] = plateau_potential;
	INa_max_out[APD_out_swhich] = INa_max;
	dVdt_Time_RP[APD_out_swhich] = dVdt_Time;
	V_m60_Time_RP[APD_out_swhich] = V_m60_Time;
	CaT_max_out[APD_out_swhich] = CaT_max;
	CaT_min_out[APD_out_swhich] = CaT_min;
	CaD50_out[APD_out_swhich] = CaD50;
	CaD80_out[APD_out_swhich] = CaD80;
	CaD90_out[APD_out_swhich] = CaD90;
	CaD_Tau_out[APD_out_swhich] = CaD_Tau;
	CaD_peak_T_out[APD_out_swhich] = CaD_peak_T;
	dCadt_time_out[APD_out_swhich] = dCadt_time;
	current_stim_counter[APD_out_swhich] = Stim_counter;
	// simulation_interval[APD_out_swhich] = t - timeAPDstart;



	for (int i = 0; i < 3; ++i)
	{
		std::cout << std::setprecision(9) << current_stim_counter[i] << " "
		          << simulation_interval[i] << " "
		          << APD_out_end[i] << " "
		          << APD75_out[i] << " "
		          << APD50_out[i] << " "
		          << APD30_out[i] << " "
		          << APD20_out[i] << " "
		          << Vmax_out[i] << " "
		          << Vmin_out[i] << " "
		          << dVdtmax_out[i] << " "
		          << plateau_potential_out[i] << " "
		          << INa_max_out[i] << " "
		          << dVdt_Time_RP[i] << " "
		          << V_m60_Time_RP[i] << " "
		          << CaT_max_out[i] << " "
		          << CaT_min_out[i] << " "
		          << CaD50_out[i]       << " "
		          << CaD80_out[i]       << " "
		          << CaD90_out[i]       << " "
		          << CaD_Tau_out[i]     << " "
		          << CaD_peak_T_out[i]  << " "
		          << dCadt_time_out[i]  << " ";
	}
	std::cout << Stim_counter << std::endl;
}



void APInfor::ReportLast() {
	std::cout << APD90 << " "
	          << APD75 << " "
	          << APD50 << " "
	          << APD30 << " "
	          << Vmax  << " "
	          << Vmin  << " "
	          << dVdtmax << " "
	          << plateau_potential << std::endl;
}

#endif