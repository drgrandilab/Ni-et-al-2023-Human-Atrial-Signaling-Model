/* SINGLECELLPARAMETER_HPP

* Update, to include simple IKur mutation parameters.
Date: Sat 24 Oct 2015 15:32:03 BST
* Update: Set IKr_grad = 1.0; at the very begeninning


*   All the parameters for single cell simulations
*   A class containing information of regional cell parameters,
    IKur mutations and AF remodelling.
    and Acacetin drug effects
*
*   Haibo NI qiangzi.ni@gmail.com
*   Date: 21 January 2015 (Wednesday)
    15:37
*/

#ifndef SINGLECELLPARAMETER_HPP
#define SINGLECELLPARAMETER_HPP
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "EnumSimulationCtrl.hpp"  // parameter input enum types;

class SingleCellPara {

public:
	SingleCellPara(AtriaCellType cell_type = RA, IKurMuationType mutation_type = WT, AFType AF_model = NONE);
	~SingleCellPara() {};

	void DefaultParameters();
	void SetCellTypePara(AtriaCellType cell_type);
	void SetIKurMutationPara(IKurMuationType mutation_type_type);
	void SetAFTypePara(AFType AF_model);
	void SetAcacetinEffectPara(double acacetin_con);
	void SetISOParameters(double ISO);
	void SetINaBlockPara(double drug_con, double Ka, double Ki, double La, double Li);

	double V;
	double dV;
	double Istim;

	double GNa, GNaL;
	double tau_ina_act;
	double tau_ina_inact;
	double drug_INa_Ka, drug_INa_La, drug_INa_Ki, drug_INa_Li;
	double drug_IKur_KO,     drug_IKur_ZKO,     drug_IKur_KC,     drug_IKur_ZKC,     drug_IKur_LO,     drug_IKur_ZLO,     drug_IKur_LC,     drug_IKur_ZLC;
	double drug_IKur_concen;
	double drug_INa_concen;

	double K2p_Vh_shift;


	double INa_Act_Shift;
	double Gto;
	double ITo_vshift;
	double Ito_Act_Shift;
	double Ito_Af;
	double S_GKv14;

	double GCaL;
	double ICaL_inac_shift;
	double ICaL_ac_shift;
	double ICaL_tau_rate;

	double GKs;
	double IKs_ac_grad;
	double IKs_shift;
	double IKs_grad;


	double GKr;
	double IKr_ac_grad;
	double IKr_ac_shift;
	double IKr_grad;


	double GK1;
	double IK1_v_shift;


	double Gf;
	double If_grad;
	double If_vshift;

	double GNaCa;
	double GNaK;
	double GCaT;
	double GCap;
	double GbCa;

	double BULK_CONST;
	double GSR_leak;
	double RyR;
	double fIRel;
	double fRyR;
	double RyRTauScale;
	double SERCA_Scale;

	double GKur;
	double IKur_cnd_add;
	double IKur_cnd_mult;
	double IKur_cnd_shift;
	double IKur_cnd_grad;
	double IKur_ac_shift;
	double IKur_ac_grad;
	double IKur_ac1_mult;
	double IKur_ac1_shift;
	double IKur_ac1_grad;
	double IKur_ac2_mult;
	double IKur_ac2_shift;
	double IKur_ac2_grad;
	double IKur_ac_add;
	double IKur_inac_mult;
	double IKur_inac_shift;
	double IKur_inac_grad;
	double IKur_inac_add;
	double Simple_Ikur_ac_shift;
	double Simple_Ikur_ac_grad;
	double Simple_Ikur_inac_shift;
	double Simple_Ikur_inac_grad;
	double Simple_GKur;

	int IKur_Model_type;
	int IKur_type_CNZ;
	int tau_type;
	float Fcell;
	float Zindex;
	int phase;
	float SAN_Gna;
	float SAN_Gto;
	float SAN_GcaL;
	float SAN_Gks;
	float SAN_Gkr;
	float SAN_Gf;
	float SAN_Gbna;
	float SAN_Gbca;
	int   SAN_periph;
	int state_len;
	int state_len2;
	AtriaCellType m_cell_type;
	IKurMuationType m_mutation_type;
	AFType m_AF_model;



	double IK2p_Shift;


	double GbNa;
	double GSK, GK2P;

	double INa, IKr, IKs, ICaL, IK1, Iab, IbK, IbCa, IK2p, ISK;
	double IbNa, ICap, INaCa, INaK, Ito, IKur, If, ICaT;

};



#endif
