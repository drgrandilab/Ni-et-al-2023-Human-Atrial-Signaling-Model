#include "signalling_para.hpp"

void signalling_para::set_CaMKII_Inhibition() {
	LCC_CKdyadp = 0.01;
	RyR_CKp = 0.15;
	PLB_CKp = 0.0003;
	NaV_CKp = 0.0022;
	Ikur_CKp = 0.0022;
	Itof_CKp = 0.0022;
	Gapjunct_CKp = 0.0022;
}
// this is for the updated bAR module, as the baseline levels of phospholation has been changed in the new module here.
void signalling_para::update_signaling_modulators_for_new_bAR_Module() {
	// PKA PHOSPHOREGULATION OF LCC AVAILABILITY (beta subunit phosph)
	double fracLCCbpo =  0.05022;//0.2506574; // Derived quantity - (LCCbp(baseline) / LCCbtot)
	// G_LTCC = 1 + 0.5 * (LCCb_PKAp - fracLCCbpo) / (1 - fracLCCbpo); // 1.5 with max phosphorylation
	// G_LTCC = 1 + 0.6 * (LCCb_PKAp - fracLCCbpo) / (1 - fracLCCbpo); // 1.6 with max phosphorylation
	// G_LTCC = 1 + 1.2 * (LCCb_PKAp - fracLCCbpo) / (1 - fracLCCbpo); // 1.6 with max phosphorylation  // version 2
	// G_LTCC = 1 + 1.05 * (LCCb_PKAp - fracLCCbpo) / (1 - fracLCCbpo); // 1.6 with max phosphorylation  // version 2
	// G_LTCC = 1 + 1.5 * (LCCb_PKAp - fracLCCbpo) / (1 - fracLCCbpo); // 1.6 with max phosphorylation  // version 2
	G_LTCC = 1 + 0.5 * (LCCb_PKAp - fracLCCbpo) / (0.41 - fracLCCbpo); // 1.6 with max phosphorylation  // version 2
	// G_LTCC = 1 + 1.8 * (LCCb_PKAp - fracLCCbpo) / (1 - fracLCCbpo); // 1.6 with max phosphorylation  // version 2
	// G_LTCC = 1 + 1.85 * (LCCb_PKAp - fracLCCbpo) / (1 - fracLCCbpo); // 1.6 with max phosphorylation  // version 2
	// G_LTCC = 1 + 2.0 * (LCCb_PKAp - fracLCCbpo) / (1 - fracLCCbpo); // 1.6 with max phosphorylation  // version 2
	// G_LTCC = 1 + 1.0 * (LCCb_PKAp - fracLCCbpo) / (1 - fracLCCbpo); // 1.6 with max phosphorylation  // version 2
	// ICa_scale =  ICa_scale * favail;
	K_PKA_LTCC = (LCCb_PKAp - fracLCCbpo) / (0.41 - fracLCCbpo);
	// K_PKA_LTCC = (LCCb_PKAp - fracLCCbpo) / (1 - fracLCCbpo);

	if (K_PKA_LTCC > 0.999)
		K_PKA_LTCC = 0.999;
	// PKA - and CaMKII - DEPENDENT SHIFTING OF DYADIC LCCS TO MODE - 2
	double fracLCCapo = 0.04267; //0.2195766; // Derived quantity - (LCCap(baseline) / LCCatot)
	double fpkam2 = 0.15 * (LCCa_PKAp - fracLCCapo) / (1 - fracLCCapo); // Assumes max phosphorylation results in 15 % mode - 2 channels
	if (fpkam2 < 0)
		fpkam2 = 0;
	fckiim2_j = LCC_CKdyadp * 0.1; // Assumes max phosphorylation results in 10%// mode - 2 channels
	double fckiim2_sl = 0; // Set to zero since SL LCCp by CaMKII is negligible
	// Sum up total fraction of CKII and PKA - shifted mode - 2 channels
	LTCC_junc_mode2 = fckiim2_j + fpkam2;
	LTCC_sl_mode2 = fckiim2_sl + fpkam2;

	// CaMKII and PKA - dependent phosphoregulation of RyR Po (from Negroni et al - RABBIT)
	// fCKII_RyR = 1+ 2.5*(RyR_CKp-0.21); // (20 * RyR_CKp / 3.0 - 1.0 / 3.0);
	fCKII_RyR = 1 + 3 * (RyR_CKp - 0.21); // (20 * RyR_CKp / 3.0 - 1.0 / 3.0);
	// fCKII_RyR = 1;//(40 * RyR_CKp / 3.0 - 3 / 3.0);
	// fPKA_RyR = RyR_PKAp * 1.025 + 0.9750;
	double fracRyRpo = 0.03764;// 0.2042761; // Derived quantity - (RyRp(baseline) / RyRtot)
	// fPKA_RyR = 1 + (RyR_PKAp - fracRyRpo) / (1 - fracRyRpo); // 2 with max phosphorylation
	// fPKA_RyR = 1 + (RyR_PKAp - fracRyRpo) / (0.27423362962963 - fracRyRpo); // 2 with 100 nM ISO  // was 0.55
	// fPKA_RyR = 1 + (RyR_PKAp - fracRyRpo) / (1 - fracRyRpo); // 2 with 100 nM ISO  // was 0.55
	fPKA_RyR = 1 + 0.5 * (RyR_PKAp - fracRyRpo) / (0.356 - fracRyRpo); // 2 with 100 nM ISO  // was 0.55
	// fPKA_RyR = 1 + 0.25*(RyR_PKAp - fracRyRpo) / (0.356 - fracRyRpo); // 2 with 100 nM ISO  // was 0.55
	// fPKA_RyR = 1 + 1*(RyR_PKAp - fracRyRpo) / (0.356 - fracRyRpo); // 2 with 100 nM ISO  // was 0.55
	// fPKA_RyR = 1 + (RyR_PKAp - fracRyRpo) / (1 - fracRyRpo); // 2 with 100 nM ISO  // was 0.55
	RyR_koSRCa_Scale = fCKII_RyR + fPKA_RyR - 1.0;

	// CaMKII and PKA-dependent phosphoregulation of PLB (changes to SERCA flux)
	fCKII_PLB = (1 - 0.5 * PLB_CKp);
	double fracPKA_PLBo = 1 - 0.0171; // 0.1265308; // Derived quantity - ((PLBtot - PLBp(baseline))/PLBtot)
	//fPKA_PLB = (PLB_PKAn/fracPKA_PLBo)*3/4 + 1/4; // 0.25 with max PKA phosphorylation
	// fPKA_PLB = (PLB_PKAn / fracPKA_PLBo) * 0.4 + 0.6; // 0.50 with max PKA phosphorylation
	// fPKA_PLB = (PLB_PKAn / fracPKA_PLBo) * 0.55 + 0.45; // 0.50 with max PKA phosphorylation
	fPKA_PLB = (PLB_PKAn / fracPKA_PLBo) * 0.5 + 0.5; // 0.50 with max PKA phosphorylation
	PLB_kmf_Scale = fCKII_PLB < fPKA_PLB ? fCKII_PLB : fPKA_PLB;  // utilize the smaller factor

	// PKA-dependent INa phosphoregulation (from IKr)
	double fracPKA_Inao = 0.1613;//0.5484815; // Derived quantity (INa_PKAp(baseline)/INatot)
	double fracPKA_Inaiso = 0.6951; //0.7858287; // Derived quantity (INa_PKAp(ISO)/INatot)
	kPKA_Ina_increase = (INa_PKAp - fracPKA_Inao) / (fracPKA_Inaiso - fracPKA_Inao);
	// GNa = (1 + 0.25 * kPKA_Ina) * GNa; // 25% increase w/ 100 nM ISO
	// CaMKII-dependent INa phosphoregulation (from RyR)
	//Nav_CKp_wt = 0.1949; // Derived quantity (// phosph Nav (baseline))
	//Nav_CKp_oe = 0.6370; // Derived quantity (// phosph Nav (CaMKII-OE 2X))
	double Nav_CKp_wt = 2.5654 / 30; // 8.5// 1X, 1 Hz
	//Nav_CKp_oe = 7.6588/30; // 25.5// 2X, 1 Hz
	double Nav_CKp_oe = 23.9069 / 30; // 80// 6X, 1 Hz
	kCKII_Nav_change = (NaV_CKp - Nav_CKp_wt) / (Nav_CKp_oe - Nav_CKp_wt); // 0 WT, 1 OE
	// effects are here...
	// P2a5 = (7.8 - P2a5) * kCKII_Nav_change + P2a5;
	// P1b5 = (0.00467 - P1b5) * kCKII_Nav_change + P1b5;
	// P1b6 = (6.5449e-7 - P1b6) * kCKII_Nav_change + P1b6;
	// P1a7 = (0.3377e-3 - P1a7) * kCKII_Nav_change + P1a7;
	// P1b7 = (1.868e-4 - P1b7) * kCKII_Nav_change + P1b7;
	// P1a8 = (6.5e-6 - P1a8) * kCKII_Nav_change + P1a8;
	// P1b8 = (3.8e-3 - P1b8) * kCKII_Nav_change + P1b8;

	double KmNaip = 11;
	// PKA-dependent PLM phosphoregulation (from Negroni et al - RABBIT)
	double fracPKA_PLMo = 0.01476;//  0.1167379; // Derived quantity (PLM_PKAp(baseline)/PLMtot)
	double fracPKA_PLMiso = 0.8204; //0.8591454; // Derived quantity (PLM_PKAp(ISO)/PLMtot)
	double kPKA_PLM = KmNaip * (1.0 - 13.6 / 18.8) / (fracPKA_PLMiso / fracPKA_PLMo - 1); // PLM_PKAp ISO  this value = 0.4784
	//kPKA_PLM=KmNaip*(0.5)*(1-13.6/18.8)/(fracPKA_PLMiso/fracPKA_PLMo-1); // PLM_PKAp ISO (50%)
	KmNaip_PKA = -kPKA_PLM + kPKA_PLM * (PLM_PKAp / fracPKA_PLMo);
	// KmNaip_PKA *=;
	// KmNaip = KmNaip - KmNaip_PKA; // 27.66% reduction w/ ISO  // effects here

	// PKA-dependent phosphoregulation of TnI (increases Kd of TnC)
	double fracTnIpo = 0.007365; //0.0626978;  // Derived quantity (TnI_PKAp(baseline)/TnItot)
	fPKA_TnI = (1.45 - 0.45 * (1 - TnI_PKAp) / (1 - fracTnIpo)); // 1.45 with maximal phosphorylation
	// koff_tncl = koff_tncl * fPKA_TnI; // not used w/ myofilament  /effects

	// PKA-dependent IKs phosphoregulation (from Negroni et al - RABBIT)
	double fracPKA_Ikso = 0.1613;//0.5484815; // Derived quantity (IKs_PKAp(baseline)/Ikstot)
	double fracPKA_Iksiso = 0.6933; //0.7858287; // Derived quantity (IKs_PKAp(ISO)/Ikstot)
	kPKA_Iks = (IKs_PKAp - fracPKA_Ikso) / (fracPKA_Iksiso - fracPKA_Ikso); // 0-1 @ 100 nM ISO

	// PKA-dependent IKr phosphoregulation (from Negroni et al - RABBIT)
	double fracPKA_Ikro = 0.1613;//0.5484815; // Derived quantity (IKr_PKAp(baseline)/Ikrtot)
	double fracPKA_Ikriso = 0.6951; //0.7858287; // Derived quantity (IKr_PKAp(ISO)/Ikrtot)
	kPKA_Ikr = (IKr_PKAp - fracPKA_Ikro) / (fracPKA_Ikriso - fracPKA_Ikro);
	Vkr_PKA = 10 * kPKA_Ikr; // 10 mV shift w/ ISO     // effects
	dGkr_PKA = 0.3 * kPKA_Ikr; // 30% increase w/ ISO  // effects
	// total effect: SSA shifted and peak increased by 37%


	//// I_kur: Ultra rapid delayed rectifier Outward K Current
	// CaMKII effect (NEW, on G IKur), from NaV phosphorylation
	double Ikur_CKp_ko = 0;
	double Ikur_CKp_wt = 2.5654 / 30;
	//Itof_CKp_oe = 25.8/30;
	// kCKII_Ikur = ((1 - 2.8 / 4.9) * Ikur_CKp + 2.8 / 4.9 * Ikur_CKp_wt - Ikur_CKp_ko) / (Ikur_CKp_wt - Ikur_CKp_ko);  // this the original
	// kCKII_Ikur = 0.5 + 1.0 / ( 1 + exp( (Ikur_CKp - Ikur_CKp_wt) / -0.04 ));
	kCKII_Ikur = 0.5 + 1.0 / ( 1 + exp( (Ikur_CKp - Ikur_CKp_wt) / -0.08 ));
	// linear fitting between CaMKII-KO (about 0.57) and normal CaMKII (1)

	// PKA effect taken from IKr (PKA-IKur not present bAR module)

	// PKA-dependent IKur phosphoregulation (from IKr/IKs)
	double fracPKA_Ikuro = 0.1613;//0.5484815; // Derived quantity (IKr_PKAp(baseline)/Ikrtot)
	double fracPKA_Ikuriso = 0.6951; //0.7858287; // Derived quantity (IKr_PKAp(ISO)/Ikrtot)
	kPKA_Ikur = (IKur_PKAp - fracPKA_Ikuro) / (fracPKA_Ikuriso - fracPKA_Ikuro);
	//fracIKuravail = 1+(3-1)*kPKA_Ikur; // 3-fold increase w/ ISO (see old model)
	// fracIKuravail = 1 + (1.2 - 1) * kPKA_Ikur; // 20% increase w/ 100 nM ISO (Li et al 1996) // effects

	// PKA-dependent IClCa phosphoregulation (from Negroni et al - RABBIT)
	double fracPKA_IClCao = 0.1349246; //0.6588394; // Derived quantity (IClCa_PKAp(baseline)/Ikrtot)
	double fracPKA_IClCaiso = 0.7362976; //0.8634406; // Derived quantity (IClCa_PKAp(ISO)/Ikrtot)
	kPKA_IClCa = (IClCa_PKAp - fracPKA_IClCao) / (fracPKA_IClCaiso - fracPKA_IClCao);

	// PKA-dependent myofilament phosphoregulation (100 nM ISO)
	double fracPKA_Myoo = 0.007365;// 0.0626978; // Derived quantity (Myo_PKAp(baseline)/Myotot)
	double fracPKA_Myoiso = 0.8325; //0.8687614; // Derived quantity (Myo_PKAp(ISO)/Myotot)
	kPKA_Myo = (Myo_PKAp - fracPKA_Myoo) / (fracPKA_Myoiso - fracPKA_Myoo);

	// CaMKII effect (NEW, on Tau Inact Itof), from NaV phosphorylation
	double Itof_CKp_ko = 0;
	double Itof_CKp_wt = 0.11;//2.5654 / 30.0;
	//Itof_CKp_oe = 25.8/30;
	// kCKII_Itof_tau = ((1 - 43 / 66.9) * Itof_CKp + 43 / 66.9 * Itof_CKp_wt - Itof_CKp_ko) / (Itof_CKp_wt - Itof_CKp_ko);

	// https://www.ahajournals.org/doi/pdf/10.1161/CIRCEP.108.842799
	kCKII_Itof_tau = 0.5 + 1.0 / ( 1 + exp( (Itof_CKp - Itof_CKp_wt) / -0.1 )) ; // CaMKII accerates Ito recovery
	// linear fitting between CaMKII-KO (about 0.64) and normal CaMKII (1)
	kCKII_Itof_G = 1.0;//((1 - 7 / 6.5) * Itof_CKp + 7 / 6.5 * Itof_CKp_wt - Itof_CKp_ko) / (Itof_CKp_wt - Itof_CKp_ko);

	// shift .. https://www.ahajournals.org/doi/full/10.1161/01.RES.85.9.810
	kCKII_Itof_vshift = -5 + 10 / ( 1 + exp( (Itof_CKp - Itof_CKp_wt) / -0.1 ));    // smaller shift here  // 12:00:59, Thu, 21-March-2019, By Haibo
	// linear fitting between CaMKII-KO (about 0.64) and normal CaMKII (1)
	// modified for human myocytes
	//GtoFast=(1.0-0.7*cAF)*0.165; %nS/pF maleckar; %human atrium
	// GtoFast = kCKII_Itof_G*(1.0-0.45*cAF*RA-0.8*cAF*(1-RA))*0.165; %nS/pF maleckar; %human atrium

	// PKA-dependent Itof phosphoregulation (from IKr/IKs)
	double fracPKA_Itof_o = 0.1613;//0.5484815; // Derived quantity (IKr_PKAp(baseline)/Ikrtot)
	double fracPKA_Itof_iso = 0.6951; //0.7858287; // Derived quantity (IKr_PKAp(ISO)/Ikrtot)
	kPKA_Itof = (Ito_PKAp - fracPKA_Itof_o) / (fracPKA_Itof_iso - fracPKA_Itof_o);
	// fracItof_avail = 1 + (0.60 - 1) * kPKA_Itof; // 40% decrease w/ 100 nM ISO (Gonzales De La Fuente et al 2013) //effects

	// PKA-dependent IK1 phosphoregulation (from IKr/IKs)
	double fracPKA_IK1_o = 0.1613;//0.5484815; // Derived quantity (IKr_PKAp(baseline)/Ikrtot)
	double fracPKA_IK1_iso = 0.6951; //0.7858287; // Derived quantity (IKr_PKAp(ISO)/Ikrtot)
	kPKA_IK1 = (IK1_PKAp - fracPKA_IK1_o) / (fracPKA_IK1_iso - fracPKA_IK1_o);
	// fracIK1_avail = 1+(0.55-1)*kPKA_IK1; // 45% decrease w/ 100 nM ISO (Gonzales De La Fuente et al 2013) //effects

	// https://www.ahajournals.org/doi/pdf/10.1161/CIRCEP.108.842799
	// kCKII_IK1_G = 0.5 + 1.0 / ( 1 + exp( (IK1_CKp - 0.1) / -0.2 )); // values require updates// IK1_CKp=0 -> 0.57; IK1_CKp = 0.1 -> 1; IK1_CKp = 0.4 -> 1.5
	kCKII_IK1_G = 0.5 + 1.0 / ( 1 + exp( (IK1_CKp - 0.1) / -0.15 )); // values require updates IK1_CKp=0 -> 0.57; IK1_CKp = 0.1 -> 1; IK1_CKp = 0.4 -> 1.5

	if (No_CaMKII_K1) {
		kCKII_IK1_G = 1.0;
	}

	// gap junctions
	// print([gap_scaling(i) for i in [0, 0.1, 0.2, 0.3, 0.5 ]]) ===>
	// [1.1622523254065595, 1.0102133980816692, 0.8256868141839169, 0.6637122255095125, 0.4999943341418507]
	kCKII_Gapjunct_G = 1.35 - (0.9 / (1 + exp( -(Gapjunct_CKp - 0.16) / 0.12 ) ) );

	if (No_CaMKII_GapG) {
		kCKII_Gapjunct_G = 1.0;
	}


}


void signalling_para::print_to_file(double t, std::ofstream & output_file) {

	output_file <<  std::setprecision(9)
	            << t << " "  << std::setprecision(3)
	            << LCC_CKdyadp << " "
	            << RyR_CKp << " "
	            << PLB_CKp << " "
	            << NaV_CKp << " "
	            << Ikur_CKp << " "
	            << Itof_CKp << " "
	            << LCCa_PKAp  << " "
	            << LCCb_PKAp  << " "
	            << PLB_PKAn   << " " // 10
	            << PLM_PKAp   << " "
	            << RyR_PKAp   << " "
	            << TnI_PKAp   << " "
	            << IKs_PKAp   << " "
	            << IKr_PKAp   << " "  // 15
	            << IKur_PKAp  << " "
	            << ICFTR_PKAp << " "    // not present in HA << " "
	            << IClCa_PKAp << " "
	            << Myo_PKAp   << " "
	            << INa_PKAp   << " "   // 20
	            << Ito_PKAp   << " "
	            << IK1_PKAp   << " "
	            << kCKII_IK1_G << " "
	            << kPKA_IK1 << " "
	            << kPKA_Itof << " "     // 25
	            << kCKII_Itof_vshift << " "
	            << kCKII_Itof_tau << " "
	            << kPKA_Myo << " "
	            << kPKA_IClCa << " "
	            << kPKA_Ikur << " "  // 30
	            << kCKII_Ikur << " "
	            << kPKA_Ikr << " "
	            << kPKA_Iks << " "
	            << fPKA_TnI << " "
	            << KmNaip_PKA << " "  // 35
	            << kCKII_Nav_change << " "
	            << kPKA_Ina_increase << " "
	            << PLB_kmf_Scale << " "
	            << fPKA_PLB << " "
	            << fCKII_PLB << " "  // 40
	            << RyR_koSRCa_Scale << " "
	            << fPKA_RyR << " "
	            << fCKII_RyR << " "
	            << LTCC_sl_mode2 << " "
	            << LTCC_junc_mode2 << " "  // 45
	            << fckiim2_j << " "
	            << fpkam2 << " "
	            << G_LTCC << " " //48
	            << Vkr_PKA << " "
	            << kCKII_Gapjunct_G // 50
	            << std::endl;
}

// setting up scaling factor for the phosphorylation levels by PKA and CaMKII in the model
// 11:12:49, Fri, 10-January-2020, By Haibo
void signalling_para::set_scaling_for_phos_levels(double *Scale) {

	PLB_PKAn_Scale    *= Scale[0];
	PLM_PKAp_Scale    *= Scale[1];
	LCCa_PKAp_Scale   *= Scale[2];
	LCCb_PKAp_Scale   *= Scale[3];
	RyR_PKAp_Scale    *= Scale[4];
	TnI_PKAp_Scale    *= Scale[5];
	IKs_PKAp_Scale    *= Scale[6];
	IKr_PKAp_Scale    *= Scale[7];
	IKur_PKAp_Scale   *= Scale[8];
	IClCa_PKAp_Scale  *= Scale[9];
	Myo_PKAp_Scale    *= Scale[10];
	INa_PKAp_Scale    *= Scale[11];
	Ito_PKAp_Scale    *= Scale[12];
	IK1_PKAp_Scale    *= Scale[13];
	LCC_CKdyadp_Scale *= Scale[14];
	NaV_CKp_Scale     *= Scale[15];
	RyR_CKp_Scale     *= Scale[16];
	PLB_CKp_Scale     *= Scale[17];
	Ikur_CKp_Scale    *= Scale[18];
	Itof_CKp_Scale    *= Scale[19];
	IK1_CKp_Scale     *= Scale[20];
	// ICFTR_PKAp_Scale  *= Scale[22]; // not used in the atrial model

}