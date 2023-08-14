#ifndef CAM_HPP
#define CAM_HPP
#include  <iomanip>

#include "signalling_para.hpp"

class CaM
{
public:
	CaM() {
		K =  135.0000    ;
		Mg =  1.0000  ;
		CaMtot =  418.0000         ;
		Btot =  0  ;
		CaMKIItot =  120.0000    ;
		CaNtot =  3.6;//3.61750874231279   ;
		PP1tot =  96.5000    ;
		// Ca =  0.3494  ;
		// cyclelength =  333.3333    ;
		// compartment =  2.0000;
	};
	~CaM() {
		y = nullptr;
		ydot = nullptr;
	};

// 135	1	418	0	120	3.61750874231279	96.5000000000000	1	333.333333333333	2
	double *y, *ydot;

	double K, Mg, CaMtot, Btot, CaMKIItot, CaNtot, PP1tot;
	double Ca;
	// double cyclelength, compartment;
	double JCa;
	const int ODE_NUM = 15;

	int update_single_time_step(CaM_para & para) {

		K         = para.K         ;
		Mg        = para.Mg        ;
		CaMtot    = para.CaMtot;
		Btot      = para.Btot;
		CaMKIItot = para.CaMKIItot;
		CaNtot    = para.CaNtot;
		PP1tot    = para.PP1tot;
		Ca        = para.Ca;


		double CaM          = y[0];
		double Ca2CaM       = y[1];
		double Ca4CaM       = y[2];
		double CaMB         = y[3];
		double Ca2CaMB      = y[4];
		double Ca4CaMB      = y[5];
		double Pb2          = y[6];
		double Pb           = y[7];
		double Pt           = y[8];
		double Pt2          = y[9];
		double Pa           = y[10];
		double Ca4CaN       = y[11];
		double CaMCa4CaN    = y[12];
		double Ca2CaMCa4CaN = y[13];
		double Ca4CaMCa4CaN = y[14];

		double Kd02, Kd24;
		if (Mg <= 1) {
			Kd02 = 0.0025 * (1 + K / 0.94 - Mg / 0.012) * (1 + K / 8.1 + Mg / 0.022); // [uM^2]
			Kd24 = 0.128 * (1 + K / 0.64 + Mg / 0.0014) * (1 + K / 13.0 - Mg / 0.153); // [uM^2]
		} else {
			Kd02 = 0.0025 * (1 + K / 0.94 - 1 / 0.012 + (Mg - 1) / 0.060) * (1 + K / 8.1 + 1 / 0.022 + (Mg - 1) / 0.068); // [uM^2]
			Kd24 = 0.128 * (1 + K / 0.64 + 1 / 0.0014 + (Mg - 1) / 0.005) * (1 + K / 13.0 - 1 / 0.153 + (Mg - 1) / 0.150); // [uM^2]
		}
		double k20 = 10;               // [s^-1]
		double k02 = k20 / Kd02;       // [uM^-2 s^-1]
		double k42 = 500;              // [s^-1]
		double k24 = k42 / Kd24;       // [uM^-2 s^-1]

		// CaM buffering (B) parameters
		double k0Boff = 0.0014;        // [s^-1]
		double k0Bon = k0Boff / 0.2; // [uM^-1 s^-1] kon = koff/Kd
		double k2Boff = k0Boff / 100;  // [s^-1]
		double k2Bon = k0Bon;          // [uM^-1 s^-1]
		double k4Boff = k2Boff;        // [s^-1]
		double k4Bon = k0Bon;          // [uM^-1 s^-1]
		// using thermodynamic constraints
		double k20B = k20 / 100; // [s^-1] thermo constraint on loop 1
		double k02B = k02;     // [uM^-2 s^-1]
		double k42B = k42;     // [s^-1] thermo constraint on loop 2
		double k24B = k24;     // [uM^-2 s^-1]

		// CaMKII parameters
		// Wi Wa Wt Wp
		double kbi = 2.2;      // [s^-1] (Ca4CaM dissocation from Wb)
		double kib = kbi / 33.5e-3; // [uM^-1 s^-1]
		double kib2 = kib;
		double kb2i = kib2 * 5;
		double kb24 = k24;
		double kb42 = k42 * 33.5e-3 / 5.0;
		double kpp1 = 1.72;    // [s^-1] (PP1-dep dephosphorylation rates)
		double Kmpp1 = 11.5;   // [uM]
		double kta = kbi / 1000.0; // [s^-1] (Ca4CaM dissociation from Wt)
		double kat = kib;      // [uM^-1 s^-1] (Ca4CaM reassociation with Wa)
		double kt42 = k42 * 33.5e-6 / 5.0;
		double kt24 = k24;
		double kat2 = kib;
		double kt2a = kib * 5;

		// CaN parameters
		double kcanCaoff = 1;              // [s^-1]
		double kcanCaon = kcanCaoff / 0.5; // [uM^-1 s^-1]
		double kcanCaM4on = 46;            // [uM^-1 s^-1]
		double kcanCaM4off = 1.3e-3;       // [s^-1]
		double kcanCaM2on = kcanCaM4on;
		double kcanCaM2off = 2508 * kcanCaM4off;
		double kcanCaM0on = kcanCaM4on;
		double kcanCaM0off = 165 * kcanCaM2off;
		double k02can = k02;
		double k20can = k20 / 165;
		double k24can = k24;
		double k42can = k20 / 2508;


		// CaM Reaction fluxes
		double rcn02 = k02 * (Ca * Ca) * CaM - k20 * Ca2CaM;
		double rcn24 = k24 * (Ca * Ca) * Ca2CaM - k42 * Ca4CaM;
		// CaM buffer fluxes
		double B = Btot - CaMB - Ca2CaMB - Ca4CaMB;
		double rcn02B = k02B * (Ca * Ca) * CaMB - k20B * Ca2CaMB;
		double rcn24B = k24B * (Ca * Ca) * Ca2CaMB - k42B * Ca4CaMB;
		double rcn0B = k0Bon * CaM * B - k0Boff * CaMB;
		double rcn2B = k2Bon * Ca2CaM * B - k2Boff * Ca2CaMB;
		double rcn4B = k4Bon * Ca4CaM * B - k4Boff * Ca4CaMB;
		// CaN reaction fluxes
		double Ca2CaN = CaNtot - Ca4CaN - CaMCa4CaN - Ca2CaMCa4CaN - Ca4CaMCa4CaN;
		double rcnCa4CaN = kcanCaon * (Ca * Ca) * Ca2CaN - kcanCaoff * Ca4CaN;
		double rcn02CaN = k02can * (Ca * Ca) * CaMCa4CaN - k20can * Ca2CaMCa4CaN;
		double rcn24CaN = k24can * (Ca * Ca) * Ca2CaMCa4CaN - k42can * Ca4CaMCa4CaN;
		double rcn0CaN = kcanCaM0on * CaM * Ca4CaN - kcanCaM0off * CaMCa4CaN;
		double rcn2CaN = kcanCaM2on * Ca2CaM * Ca4CaN - kcanCaM2off * Ca2CaMCa4CaN;
		double rcn4CaN = kcanCaM4on * Ca4CaM * Ca4CaN - kcanCaM4off * Ca4CaMCa4CaN;
		// CaMKII reaction fluxes
		double Pi = 1 - Pb2 - Pb - Pt - Pt2 - Pa;  // comment by haibo: pi is thus the remaining proportion of CaMKII occupying Pi;

		double rcnCKib2 = kib2 * Ca2CaM * Pi - kb2i * Pb2;
		double rcnCKb2b = kb24 * (Ca * Ca) * Pb2 - kb42 * Pb;
		double rcnCKib = kib * Ca4CaM * Pi - kbi * Pb;
		double T_state = Pb + Pt + Pt2 + Pa;
		double kbt = 0.055 * T_state + 0.0074 * T_state * T_state + 0.015 * T_state * T_state * T_state; // was ist das?
		double rcnCKbt = kbt * Pb - kpp1 * PP1tot * Pt / (Kmpp1 + CaMKIItot * Pt);
		double rcnCKtt2 = kt42 * Pt - kt24 * (Ca * Ca) * Pt2;
		double rcnCKta = kta * Pt - kat * Ca4CaM * Pa;
		double rcnCKt2a = kt2a * Pt2 - kat2 * Ca2CaM * Pa;
		double rcnCKt2b2 = kpp1 * PP1tot * Pt2 / (Kmpp1 + CaMKIItot * Pt2);
		double rcnCKai = kpp1 * PP1tot * Pa / (Kmpp1 + CaMKIItot * Pa);

		// CaM equations
		double dCaM = 1e-3 * (-rcn02 - rcn0B - rcn0CaN);
		double dCa2CaM = 1e-3 * (rcn02 - rcn24 - rcn2B - rcn2CaN + CaMKIItot * (-rcnCKib2 + rcnCKt2a) );
		double dCa4CaM = 1e-3 * (rcn24 - rcn4B - rcn4CaN + CaMKIItot * (-rcnCKib + rcnCKta) );
		double dCaMB = 1e-3 * (rcn0B - rcn02B);
		double dCa2CaMB = 1e-3 * (rcn02B + rcn2B - rcn24B);
		double dCa4CaMB = 1e-3 * (rcn24B + rcn4B);

		// CaMKII equations
		double dPb2 = 1e-3 * (rcnCKib2 - rcnCKb2b + rcnCKt2b2); // Pb2
		double dPb = 1e-3 * (rcnCKib + rcnCKb2b - rcnCKbt);  // Pb
		double dPt = 1e-3 * (rcnCKbt - rcnCKta - rcnCKtt2);  // Pt
		double dPt2 = 1e-3 * (rcnCKtt2 - rcnCKt2a - rcnCKt2b2); // Pt2
		double dPa = 1e-3 * (rcnCKta + rcnCKt2a - rcnCKai); // Pa

		// CaN equations
		double dCa4CaN = 1e-3 * (rcnCa4CaN - rcn0CaN - rcn2CaN - rcn4CaN);                     // Ca4CaN


		// std::cout << rcnCa4CaN <<  rcn0CaN << rcn2CaN - rcn4CaN
		double dCaMCa4CaN = 1e-3 * (rcn0CaN - rcn02CaN);         // CaMCa4CaN
		double dCa2CaMCa4CaN = 1e-3 * (rcn2CaN + rcn02CaN - rcn24CaN); // Ca2CaMCa4CaN
		double dCa4CaMCa4CaN = 1e-3 * (rcn4CaN + rcn24CaN);         // Ca4CaMCa4CaN


		// update ydot of states.
		ydot[0] = dCaM;
		ydot[1] = dCa2CaM;
		ydot[2] = dCa4CaM;
		ydot[3] = dCaMB;
		ydot[4] = dCa2CaMB;
		ydot[5] = dCa4CaMB;
		ydot[6] = dPb2;
		ydot[7] = dPb;
		ydot[8] = dPt;
		ydot[9] = dPt2;
		ydot[10] = dPa;
		ydot[11] = dCa4CaN;
		ydot[12] = dCaMCa4CaN;
		ydot[13] = dCa2CaMCa4CaN;
		ydot[14] = dCa4CaMCa4CaN;

		// write to global variables for adjusting Ca buffering in EC coupling model
		// commmented by  // 15:20:42, Wed, 08-November-2017, By Haibo, makes sense to me
		JCa = 1e-3 * (2 * CaMKIItot * (rcnCKtt2 - rcnCKb2b) - 2 * (rcn02 + rcn24 + rcn02B + rcn24B + rcnCa4CaN + rcn02CaN + rcn24CaN)); // [uM/msec]

		return 0;

	}




	// CaM modules should have different intitial values for different compartments
	void initialiser(int CaM_compartment) {
		// CaM_compartment = 0 -> Dyad
		// CaM_compartment = 1 -> sl
		// CaM_compartment = 2 -> cyto
		if (CaM_compartment == 0) {
			y[0] = 294.545479542006;
			y[1] = 7.29988762926394;
			y[2] = 0.00355467213776979;
			y[3] = 0;
			y[4] = 0;
			y[5] = 0;
			y[6] = 0.567125177279042;
			y[7] = 0.0449077438977342;
			y[8] = 1.17847480402579e-05;
			y[9] = 4.26934065935237e-09;
			y[10] = 2.91929618330730e-09;
			y[11] = 0.000126745509662228;
			y[12] = 0.00318217933682596;
			y[13] = 0.0133769033811970;
			y[14] = 3.60078218589155;
		} else if (CaM_compartment == 1) {
			y[0]  = 0.0312058720142851;
			y[1]  = 0.000181385913425923;
			y[2]  = 2.17447735854594e-08;
			y[3]  = 2.17580578777376;
			y[4]  = 9.21011974618692;
			y[5]  = 0.000948921668443340;
			y[6]  = 3.63437020164358e-05;
			y[7]  = 3.35748048060805e-06;
			y[8]  = 6.82002650906883e-09;
			y[9]  = 8.42907708787532e-14;
			y[10] = 4.05529652818466e-10;
			y[11] = 0.000737399755237718;
			y[12] = 1.97072137084833e-06;
			y[13] = 5.56246916279762e-06;
			y[14] = 0.000785295757929995;
		} else if (CaM_compartment == 2) {
			y[0]  = 0.0310311973548212;
			y[1]  = 0.000162957843604266;
			y[2]  = 1.38465477578005e-08;
			y[3]  = 2.87218731579093;
			y[4]  = 2.61542720124795;
			y[5]  = 0.000219725429795659;
			y[6]  = 3.26263796848757e-05;
			y[7]  = 5.83609249218633e-07;
			y[8]  = 2.32362914511299e-12;
			y[9]  = 2.80386050960474e-17;
			y[10] = 1.32943351757868e-13;
			y[11] = 0.000404922509723759;
			y[12] = 1.07523749597845e-06;
			y[13] = 1.42126244821333e-06;
			y[14] = 2.48334685738420e-05;
		}
	}


	void print_content() {

		std::cerr << "K 		"<<  K 			<< std::endl
		          << "Mg       " << Mg         << std::endl
		          << "CaMtot   " << CaMtot     << std::endl
		          << "Btot     " << Btot       << std::endl
		          << "CaMKIItot" << CaMKIItot  << std::endl
		          << "CaNtot   " << CaNtot     << std::endl
		          << "PP1tot   " << PP1tot     << std::endl;
	}

	void print_to_file(double t, std::ofstream & output_file) {


		output_file <<  std::setprecision(9) << t << " ";

		for (int i = 0; i < 15; ++i)
		{
			output_file << std::setprecision(4) << y[i] << " ";
		}
		output_file << std::endl;
	}

};




#endif