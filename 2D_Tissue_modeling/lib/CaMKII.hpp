
#ifndef CAMKII__HPP
#define CAMKII__HPP

#include  <iomanip>

#include "signalling_para.hpp"
class CaMKII
{
public:
	CaMKII() {};
	~CaMKII() {
		y = nullptr;
		ydot = nullptr;
	};

	// This function computes the CaMKII-dependent phosphorylation profiles for
	// LCCs (dyadic and subsarcolemmal), RyRs, and PLB.
	//// Description of state variables
	// LCCp-PKA = y(1);        // [LCCp] by PKA (currently unused anywhere else)
	// LCCp-CaMKIIdyad = y(2); // Dyadic [LCCp] by dyadic CaMKII
	// RyR-Ser2809p = y(3);    // [RyR-Ser2809p] by PKA (currently unused anywhere else)
	// RyR-Ser2815p = y(4);    // [RyR-Ser2815p] by CaMKII
	// PLB-Thr17p = y(5);      // [PLB-Thr17p] by CaMKII
	// LCCp-CaMKIIsl = y(6);   // Subsarcolemmal [LCCp] by subsarcolemmal CaMKII
	//// PARAMETERS
	// Default PKA level (currently unused elsewhere)
	//PKAc = 95.6*.54;
	// Okadaic Acid inhibition params (based on Huke/Bers [2008])
	// Want to treat OA as non-competitive inhibitor of PP1 and PP2A
	// Okadaic Acid inhibition params (based on Huke/Bers [2008])
	// Want to treat OA as non-competitive inhibitor of PP1 and PP2A

	double Ki_OA_PP1 = 0.78;        // [uM] - Values from fit
	double Ki_OA_PP2A = 0.037;      // [uM] - Values from fit
	double CaMKIIactDyad, LCCtotDyad, RyRtot, PP1_dyad, PP2A_dyad, OA, PLBtot, NaVtot,
	       CaMKIIactSL, LCCtotSL, PP1_SL, PP1_PLB_avail;
	// double y[6], ydot[6];

	double *y, *ydot;
	const int ODE_NUM = 6;

	int update_single_time_step(const CaMKII_para & para) {

		CaMKIIactDyad = para.CaMKIIactDyad;
		LCCtotDyad    = para.LCCtotDyad;
		RyRtot        = para.RyRtot;
		PP1_dyad      = para.PP1_dyad;
		PP2A_dyad     = para.PP2A_dyad;
		OA            = para.OA;
		PLBtot        = para.PLBtot;
		NaVtot       = para.NaVtot;
		CaMKIIactSL   = para.CaMKIIactSL;
		LCCtotSL      = para.LCCtotSL;
		PP1_SL        = para.PP1_SL;
		PP1_PLB_avail = para.PP1_PLB_avail;


		double LCC_PKAp    = y[0];
		double LCC_CKdyadp = y[1];
		double RyR2809p    = y[2];
		double RyR2815p    = y[3];
		double PLBT17p     = y[4];
		double LCC_CKslp   = y[5];

		double NaVp = RyR2809p;
		//// OA inhibition term (non-competitive) for PP1 and PP2A
		double OA_PP1 = 1 / (1 + (OA / Ki_OA_PP1) * (OA / Ki_OA_PP1) * (OA / Ki_OA_PP1) );
		double OA_PP2A = 1 / (1 + (OA / Ki_OA_PP2A) * (OA / Ki_OA_PP2A) * (OA / Ki_OA_PP2A));
		//// ODE EQUATIONS
		//// LCC states (note: PP2A is acting on PKA site and PP1 on CKII site)
		// L-Type Ca Channel (LCC) parameters

		double k_ckLCC = 0.4;                  // [s^-1]
		double k_pp1LCC = 0.1103;              // [s^-1]
		double k_pkaLCC = 13.5;                // [s^-1]
		double k_pp2aLCC = 10.1;               // [s^-1]

		double KmCK_LCC = 12.0;                  // [uM]
		double KmPKA_LCC = 21.0;                 // [uM]
		double KmPP2A_LCC = 47.0;                // [uM]
		double KmPP1_LCC = 9.0;                  // [uM]

		// CaMKII phosphorylation of Dyadic LCCs
		double LCC_CKdyadn = LCCtotDyad - LCC_CKdyadp;
		double LCCDyad_PHOS = (k_ckLCC * CaMKIIactDyad * LCC_CKdyadn) / (KmCK_LCC + LCC_CKdyadn);
		double LCCDyad_DEPHOS = (k_pp1LCC * PP1_dyad * LCC_CKdyadp) / (KmPP1_LCC + LCC_CKdyadp) * OA_PP1;
		double dLCC_CKdyadp = LCCDyad_PHOS - LCCDyad_DEPHOS;
		// std::cout << LCCDyad_PHOS << "  " << LCCDyad_DEPHOS << std::endl;

		// CaMKII phosphorylation of Sub-sarcolemmal LCCs
		double LCC_CKsln = LCCtotSL - LCC_CKslp;
		double LCCSL_PHOS = (k_ckLCC * CaMKIIactSL * LCC_CKsln) / (KmCK_LCC + LCC_CKsln);
		double LCCSL_DEPHOS = (k_pp1LCC * PP1_SL * LCC_CKslp) / (KmPP1_LCC + LCC_CKslp) * OA_PP1;
		double dLCC_CKslp = LCCSL_PHOS - LCCSL_DEPHOS;

		// PKA phosphorylation (currently unused elsewhere)
		// LCC_PKAn = LCCtotDyad - LCC_PKAp;
		// dLCC_PKAp = (k_pkaLCC*PKAc*LCC_PKAn)/(KmPKA_LCC+LCC_PKAn) - ...
		//             (k_pp2aLCC*PP2A_dyad*LCC_PKAp)/(KmPP2A_LCC+LCC_PKAp)*OA_PP2A;
		double dLCC_PKAp = 0;
		//// RyR states
		// Ryanodine Receptor (RyR) parameters
		double k_ckRyR = 0.4;                  // [s^-1]
		double k_pkaRyR = 1.35;                // [s^-1]
		double k_pp1RyR = 1.07;                // [s^-1]
		double k_pp2aRyR = 0.481;              // [s^-1]

		double KmCK_RyR = 12;                  // [uM]
		double KmPKA_RyR = 21;                 // [uM]
		double KmPP1_RyR = 9;                  // [uM]
		double KmPP2A_RyR = 47;                // [uM]

		// Basal RyR phosphorylation (numbers based on param estimation)
		double kb_2809 = 0.51;                 // [uM/s] - PKA site
		double kb_2815 = 0.35;                 // [uM/s] - CaMKII site

		double RyR2815n = RyRtot - RyR2815p;
		double RyR_BASAL = kb_2815 * RyR2815n;
		double RyR_PHOS = (k_ckRyR * CaMKIIactDyad * RyR2815n) / (KmCK_RyR + RyR2815n);
		double RyR_PP1_DEPHOS = (k_pp1RyR * PP1_dyad * RyR2815p) / (KmPP1_RyR + RyR2815p) * OA_PP1;
		double RyR_PP2A_DEPHOS = (k_pp2aRyR * PP2A_dyad * RyR2815p) / (KmPP2A_RyR + RyR2815p) * OA_PP2A;
		double dRyR2815p = RyR_BASAL + RyR_PHOS - RyR_PP1_DEPHOS - RyR_PP2A_DEPHOS;

		RyR2815n = RyRtot - RyR2815p;
		RyR_BASAL = kb_2815 * RyR2815n;
		RyR_PHOS = (k_ckRyR * CaMKIIactDyad * RyR2815n) / (KmCK_RyR + RyR2815n);
		RyR_PP1_DEPHOS = (k_pp1RyR * PP1_dyad * RyR2815p) / (KmPP1_RyR + RyR2815p) * OA_PP1;
		RyR_PP2A_DEPHOS = (k_pp2aRyR * PP2A_dyad * RyR2815p) / (KmPP2A_RyR + RyR2815p) * OA_PP2A;
		dRyR2815p = RyR_BASAL + RyR_PHOS - RyR_PP1_DEPHOS - RyR_PP2A_DEPHOS;


		// PKA phosphorylation of Ser 2809 on RyR (currently unused elsewhere)
		// RyR2809n = RyRtot - RyR2809p;
		// dRyR2809p = kb_2809*RyR2809n + (k_pkaRyR*PKAc*RyR2809n)/(KmPKA_RyR+RyR2809n) - ...
		//             (k_pp1RyR*PP1_dyad*RyR2809p)/(KmPP1_RyR+RyR2809p)*OA_PP1;
		//// NaV states (from original RyR)
		// // Na channel (NaV) parameters
		// k_ckNaV = 0.4;                  // [s^-1]
		// k_pp1NaV = 1.07;                // [s^-1]
		//  k_pp2aNaV = 0.481;              // [s^-1]
		//
		// KmCK_NaV = 12;                  // [uM]
		// KmPP1_NaV = 9;                  // [uM]
		//     KmPP2A_NaV = 47;                // [uM]
		//
		// NaVn = NaVtot - NaVp;
		// NaV_BASAL = kb_2815*NaVn;
		// NaV_PHOS = (k_ckNaV*CaMKIIactDyad*NaVn)/(KmCK_NaV+NaVn);
		// NaV_PP1_DEPHOS = (k_pp1NaV*PP1_dyad*NaVp)/(KmPP1_NaV+NaVp)*OA_PP1;
		// NaV_PP2A_DEPHOS = (k_pp2aNaV*PP2A_dyad*NaVp)/(KmPP2A_NaV+NaVp)*OA_PP2A;
		// dNaVp = NaV_BASAL + NaV_PHOS - NaV_PP1_DEPHOS - NaV_PP2A_DEPHOS;

		// Na channel (NaV) parameters
		double k_ckNaV = 0.08;                  // [s^-1]
		double k_pp1NaV = 0.08;               // [s^-1]

		double KmCK_NaV = 12;                  // [uM]
		double KmPP1_NaV = 9;                  // [uM]

		double NaVn = NaVtot - NaVp;
		double NaV_PHOS = (k_ckNaV * CaMKIIactDyad * NaVn) / (KmCK_NaV + NaVn);
		double NaV_PP1_DEPHOS = (k_pp1NaV * PP1_dyad * NaVp) / (KmPP1_NaV + NaVp) * OA_PP1;
		double dNaVp = NaV_PHOS - NaV_PP1_DEPHOS;
		//// PLB states
		// Phospholamban (PLB) parameters
		double k_ckPLB = 8e-3;                 // [s^-1]
		double k_pp1PLB = .0428;               // [s^-1]

		double KmCK_PLB = 12;
		double KmPP1_PLB = 9;

		double PP1_PLB = PP1_dyad * PP1_PLB_avail;  // Inhibitor-1 regulation of PP1_dyad included here
		double PLBT17n = PLBtot - PLBT17p;
		double PLB_PHOS = (k_ckPLB * PLBT17n * CaMKIIactDyad) / (KmCK_PLB + PLBT17n);
		double PLB_DEPHOS = (k_pp1PLB * PP1_PLB * PLBT17p) / (KmPP1_PLB + PLBT17p) * OA_PP1;
		double dPLBT17p = PLB_PHOS - PLB_DEPHOS;
		//// Collect ODEs and convert to uM/ms
		// ydot = [dLCC_PKAp;
		//         dLCC_CKdyadp;
		//         dNaVp;//dRyR2809p;
		//         dRyR2815p;
		//         dPLBT17p;
		//         dLCC_CKslp].*10 ^ -3; // Convert to uM/ms

		ydot[0] = dLCC_PKAp * 1e-3;
		ydot[1] = dLCC_CKdyadp * 1e-3;
		ydot[2] = dNaVp * 1e-3;
		ydot[3] = dRyR2815p * 1e-3;
		ydot[4] = dPLBT17p * 1e-3;
		ydot[5] = dLCC_CKslp * 1e-3;

	}


	void initialiser() {


		y[0] = 16.4543925700766;
		y[1] = 15.2602358771970;
		y[2] = 2.65618652708035;
		y[3] = 74.4599383708654;
		y[4] = 0.448855577894634;
		y[5] = 2.63869941309990e-06;
	}




	void print_to_file(double t, std::ofstream & output_file) {

		// output_file <<  std::setprecision(9);
		output_file <<  std::setprecision(9) << t << " ";


		for (int i = 0; i < 6; ++i)
		{
			output_file << std::setprecision(4) << y[i] << " ";
		}

		output_file << CaMKIIactDyad << " " << RyRtot <<  std::endl;
	}

};





#endif