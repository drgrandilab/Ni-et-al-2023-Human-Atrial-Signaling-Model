

#include "EnumSimulationCtrl.hpp"

AtriaCellType GetCellTypeFromLabel(int tissue) {
	if (tissue == 10) { // SAN
		return SAN_C;
	}
	else if (tissue == 11) { // CT
		return CT;

	}
	else if (tissue == 12) { // PM
		return PM;
	}
	else if (tissue == 202) { // RAA / APG
		return RAA;
	}
	else if (tissue == 201) { //LAA
		return LAA;
	}
	else if (tissue == 13 || tissue == 14 || tissue == 200) // RA (no change)11
	{
		return RA;
	}
	else if (tissue == 15) { // BB
		return BB;
	}
	else if (tissue == 16) { // LA
		return LA;
	}
	else if (tissue == 17) { // AS
		return AS;
	}
	else if (tissue == 18) { // AVR
		return AVR;
	}
	else if (tissue == 106) { // AVN - compact node
		return AVN;
	}
	else if (tissue == 104) { // AVN - inferior nodal extension
		return AVN;
	}

	else if (tissue == 105) { // AVN - Penetrating bundle
		return AVN;
	}
	else if (tissue == 102 || tissue == 103) { // AVN - transitional area
		return AVN;

	}
	else if (tissue == 101) { // PV
		return PV;
	}
	else {

		printf("wrong cell type!!!\n");
		std::exit(0);
	}
}