function ydot = NH_ECC(t, y, pECC, pPhos)

%MOD1 is for K1 (MOD1: GNa, GNaL, IClbkg) - anything but 3
%MOD2 is for K2
%MOD3 is for K3 - 3
MOD_ind = 3; %
parameters_min_opt17 = [0.999763767766896, 1.02704368132844, 0.904232756929961, 0.912089674296620, 1.06595749235083, 1.00817013225314, 0.561056138848061, 0.354046530214812, 1.35707583730889, 2.84563198875882, 1.52172787019912, 1, ...
    0.569127790357487, 1, 1.06023886956084, 1.28931577624649, 1.02096444680901, 0.956205050080743, 0.814370961024797, 0.854520304154750, 0.701949192397789, 0.833439913899504, 1.29893009653586, 1.42202828590400, 0.951213166214512];

% 130 state variables
ydot = zeros(size(y));

%% Assign passed-in parameters

% pECC = [cycleLength,AF_index,prot_index CaMKII_inhit, CaMKII_double, INa_Scale, INaL_Scale... % 1-7
%     INab_Scale, INaK_Scale, Itof_Scale, IKur_Scale, IK2p_Scale, IKr_Scale ... % 8 -13
%     IKs_Scale, IK1_Scale, IKp_Scale, IKach_Scale, ISK_Scale, ICaL_Scale, INCX_Scale... % 14 - 20
%     IClCa_Scale, IClb_Scale, Jrel_Scale, Jserca_Scale, Jleak_Scale, Cleft_Buffer_Scale... % 21-26
%     SR_Buffer_Scale ICab_Scale ICap_Scale ....% 27-31
%     Cytosol_Buffer_Scale Na_clamp Currents_record Output_current Voltage_step...% 32-36
%     gender_specific_array]; %-NH updated version % 37-45


%pECC
para.cycleLength            = pECC(1);
para.AF                     = pECC(2);
para.prot_index             = pECC(3);
para.CaMKII_inhit           = pECC(4);
para.CaMKII_double          = pECC(5);
para.INa_Scale              = pECC(6);
para.INaL_Scale             = pECC(7);
para.INab_Scale             = pECC(8);
para.INaK_Scale             = pECC(9);
para.Itof_Scale             = pECC(10);
para.IKur_Scale             = pECC(11);
para.IK2p_Scale             = pECC(12);
para.IKr_Scale              = pECC(13);
para.IKs_Scale              = pECC(14);
para.IK1_Scale              = pECC(15);
para.IKp_Scale              = pECC(16);
para.IKach_Scale            = pECC(17);
para.ISK_Scale              = pECC(18);
para.ICaL_Scale             = pECC(19);
para.INCX_Scale             = pECC(20);
para.IClCa_Scale            = pECC(21);
para.IClb_Scale             = pECC(22);
para.Jrel_Scale             = pECC(23);
para.Jserca_Scale           = pECC(24);
para.Jleak_Scale            = pECC(25);
para.Cleft_Buffer_Scale     = pECC(26);
para.SR_Buffer_Scale        = pECC(27);
para.ICab_Scale             = pECC(28);
para.ICap_Scale             = pECC(29);
para.Cytosol_Buffer_Scale   = pECC(30);
para.Na_clamp               = pECC(31);
para.Currents_record        = pECC(32);
para.Output_current         = pECC(33);

%Other Params
para.NaV_CKp            = pPhos(1); 
para.INa_PKAp           = pPhos(2);
para.LCCb_PKAp          = pPhos(3);
para.LCC_CKdyadp        = pPhos(4);
para.RyR_CKp            = pPhos(5);
para.RyR_PKAp           = pPhos(6);
para.PLB_CKp            = pPhos(7);
para.PLB_PKAn           = pPhos(8);
para.PLM_PKAp           = pPhos(9);
para.TnI_PKAp           = pPhos(10);
para.IKs_PKAp           = pPhos(11);
para.IKr_PKAp           = pPhos(12);
para.Ikur_CKp           = pPhos(13);
para.IKur_PKAp          = pPhos(14);
para.IClCa_PKAp         = pPhos(15);
para.Myo_PKAp           = pPhos(16);
para.Itof_CKp           = pPhos(17);
para.Ito_PKAp           = pPhos(18);
para.IK1_PKAp           = pPhos(19);
para.IK1_CKp            = pPhos(20);
para.Gapjunct_CKp       = pPhos(21);
para.LCCa_PKAp          = pPhos(22);

%Rename (can fix later) - NH
LCCa_PKAp = para.LCCa_PKAp;
LCCb_PKAp = para.LCCb_PKAp;

LCC_CKdyadp = para.LCC_CKdyadp; 
RyR_CKp = para.RyR_CKp ;
RyR_PKAp = para.RyR_PKAp; 
PLB_CKp = para.PLB_CKp;
PLB_PKAn = para.PLB_PKAn;
INa_PKAp = para.INa_PKAp;
PLM_PKAp = para.PLM_PKAp;
TnI_PKAp = para.TnI_PKAp;
IKs_PKAp = para.IKs_PKAp;
IKr_PKAp = para.IKr_PKAp;
Ikur_CKp = para.Ikur_CKp;
IKur_PKAp = para.IKur_PKAp;
IClCa_PKAp = para.IClCa_PKAp;
Itof_CKp = para.Itof_CKp;
Ito_PKAp  = para.Ito_PKAp;
IK1_PKAp = para.IK1_PKAp;
IK1_CKp = para.IK1_CKp;

AF = para.AF;
%% Signalging.cpp

fracLCCbpo = 0.05022; % Derived quantity - (LCCbp(baseline) / LCCbtot)

G_LTCC = 1 + 0.5 * (LCCb_PKAp - fracLCCbpo) / (0.41 - fracLCCbpo); % 1.6 with max phosphorylation  % version 2

K_PKA_LTCC = (LCCb_PKAp - fracLCCbpo) / (0.41 - fracLCCbpo);

if K_PKA_LTCC > 0.999
    K_PKA_LTCC = 0.999;
end

fracLCCapo = 0.04267; % Derived quantity - (LCCap(baseline) / LCCatot)
fpkam2 = 0.15 * (LCCa_PKAp - fracLCCapo) / (1 - fracLCCapo); % Assumes max phosphorylation results in 15 % mode - 2 channels
if fpkam2 < 0
    fpkam2 = 0;
end

fckiim2_j = LCC_CKdyadp * 0.1; % Assumes max phosphorylation results in 10%// mode - 2 channels
fckiim2_sl = 0; % Set to zero since SL LCCp by CaMKII is negligible

LTCC_junc_mode2 = fckiim2_j + fpkam2; 
LTCC_sl_mode2 = fckiim2_sl + fpkam2;
para.LTCC_sl_mode2  = LTCC_sl_mode2;

fCKII_RyR = 1 + 3 * (RyR_CKp - 0.21);

fracRyRpo = 0.03764; % 0.2042761; % Derived quantity - (RyRp(baseline) / RyRtot)

fPKA_RyR = 1 + 0.5 * (RyR_PKAp - fracRyRpo) / (0.356 - fracRyRpo);

RyR_koSRCa_Scale = fCKII_RyR + fPKA_RyR - 1.0;

fCKII_PLB = (1 - 0.5 * PLB_CKp);
fracPKA_PLBo = 1 - 0.0171; % 0.1265308; % Derived quantity - ((PLBtot - PLBp(baseline))/PLBtot)

fPKA_PLB = (PLB_PKAn / fracPKA_PLBo) * 0.5 + 0.5; % 0.50 with max PKA phosphorylation
PLB_kmf_Scale = min(fCKII_PLB, fPKA_PLB); % utilize the smaller factor

fracPKA_Inao = 0.1613; % 0.5484815; % Derived quantity (INa_PKAp(baseline)/INatot)
fracPKA_Inaiso = 0.6951; % 0.7858287; % Derived quantity (INa_PKAp(ISO)/INatot)
kPKA_Ina_increase = (INa_PKAp - fracPKA_Inao) / (fracPKA_Inaiso - fracPKA_Inao);

% Nav_CKp_wt = 2.5654 / 30; % 8.5// 1X, 1 Hz


% Nav_CKp_oe = 23.9069 / 30; % 80// 6X, 1 Hz
% kCKII_Nav_change = (NaV_CKp - Nav_CKp_wt) / (Nav_CKp_oe - Nav_CKp_wt); % 0 WT, 1 OE

KmNaip = 11;

fracPKA_PLMo = 0.01476; % 0.1167379; % Derived quantity (PLM_PKAp(baseline)/PLMtot)
fracPKA_PLMiso = 0.8204; % 0.8591454; % Derived quantity (PLM_PKAp(ISO)/PLMtot)
kPKA_PLM = KmNaip * (1.0 - 13.6 / 18.8) / (fracPKA_PLMiso / fracPKA_PLMo - 1); % PLM_PKAp ISO this value = 0.4784
KmNaip_PKA = -kPKA_PLM + kPKA_PLM * (PLM_PKAp / fracPKA_PLMo);

fracTnIpo = 0.007365; % 0.0626978;  % Derived quantity (TnI_PKAp(baseline)/TnItot)
fPKA_TnI = (1.45 - 0.45 * (1 - TnI_PKAp) / (1 - fracTnIpo)); % 1.45 with maximal phosphorylation


fracPKA_Ikso = 0.1613; % 0.5484815; % Derived quantity (IKs_PKAp(baseline)/Ikstot)
fracPKA_Iksiso = 0.6933; % 0.7858287; % Derived quantity (IKs_PKAp(ISO)/Ikstot)
kPKA_Iks = (IKs_PKAp - fracPKA_Ikso) / (fracPKA_Iksiso - fracPKA_Ikso); % 0-1 @ 100 nM ISO

fracPKA_Ikro = 0.1613; % 0.5484815; % Derived quantity (IKr_PKAp(baseline)/Ikrtot)
fracPKA_Ikriso = 0.6951; % 0.7858287; % Derived quantity (IKr_PKAp(ISO)/Ikrtot)
kPKA_Ikr = (IKr_PKAp - fracPKA_Ikro) / (fracPKA_Ikriso - fracPKA_Ikro);
Vkr_PKA = 10 * kPKA_Ikr; % 10 mV shift w/ ISO     % effects
dGkr_PKA = 0.3 * kPKA_Ikr; % 30% increase w/ ISO  % effects


% Ikur_CKp_ko = 0;
Ikur_CKp_wt = 2.5654 / 30;
kCKII_Ikur = 0.5 + 1.0 ./ (1 + exp((Ikur_CKp - Ikur_CKp_wt) ./ -0.08));

fracPKA_Ikuro = 0.1613; % 0.5484815; % Derived quantity (IKr_PKAp(baseline)/Ikrtot)
fracPKA_Ikuriso = 0.6951; % 0.7858287; % Derived quantity (IKr_PKAp(ISO)/Ikrtot)
kPKA_Ikur = (IKur_PKAp - fracPKA_Ikuro) / (fracPKA_Ikuriso - fracPKA_Ikuro);

fracPKA_IClCao = 0.1349246; % 0.6588394; % Derived quantity (IClCa_PKAp(baseline)/Ikrtot)
fracPKA_IClCaiso = 0.7362976; % 0.8634406; % Derived quantity (IClCa_PKAp(ISO)/Ikrtot)
kPKA_IClCa = (IClCa_PKAp - fracPKA_IClCao) / (fracPKA_IClCaiso - fracPKA_IClCao);

% fracPKA_Myoo = 0.007365; % 0.0626978; % Derived quantity (Myo_PKAp(baseline)/Myotot)
% fracPKA_Myoiso = 0.8325; % 0.8687614; % Derived quantity (Myo_PKAp(ISO)/Myotot)
% kPKA_Myo = (Myo_PKAp - fracPKA_Myoo) / (fracPKA_Myoiso - fracPKA_Myoo);

% Itof_CKp_ko = 0;
Itof_CKp_wt = 0.11; % 2.5654 / 30.0;

kCKII_Itof_tau = 0.5 + 1.0 ./ (1 + exp((Itof_CKp - Itof_CKp_wt) / -0.1)); % CaMKII accelerates Ito recovery

% Linear fitting between CaMKII-KO (about 0.64) and normal CaMKII (1)
% kCKII_Itof_G = 1.0; % ((1 - 7 / 6.5) * Itof_CKp + 7 / 6.5 * Itof_CKp_wt - Itof_CKp_ko) / (Itof_CKp_wt - Itof_CKp_ko);

% Shift - https://www.ahajournals.org/doi/full/10.1161/01.RES.85.9.810
kCKII_Itof_vshift = -5 + 10 ./ (1 + exp((Itof_CKp - Itof_CKp_wt) / -0.1)); % Smaller shift here
para.kCKII_Itof_vshift = kCKII_Itof_vshift;
% PKA-dependent Itof phosphoregulation (from IKr/IKs)
fracPKA_Itof_o = 0.1613; % Derived quantity (IKr_PKAp(baseline)/Ikrtot)
fracPKA_Itof_iso = 0.6951; % Derived quantity (IKr_PKAp(ISO)/Ikrtot)
kPKA_Itof = (Ito_PKAp - fracPKA_Itof_o) / (fracPKA_Itof_iso - fracPKA_Itof_o);


fracPKA_IK1_o = 0.1613; % Derived quantity (IKr_PKAp(baseline)/Ikrtot)
fracPKA_IK1_iso = 0.6951; % Derived quantity (IKr_PKAp(ISO)/Ikrtot)
kPKA_IK1 = (IK1_PKAp - fracPKA_IK1_o) / (fracPKA_IK1_iso - fracPKA_IK1_o);
para.kPKA_IK1 = kPKA_IK1;

kCKII_IK1_G = 0.5 + 1.0 ./ (1 + exp((IK1_CKp - 0.1) ./ -0.15)); % values require updates IK1_CKp=0 -> 0.57; IK1_CKp = 0.1 -> 1; IK1_CKp = 0.4 -> 1.5

No_CaMKII_GapG = 0;

if No_CaMKII_GapG == 1
    kCKII_IK1_G = 1.0;
end

% kCKII_Gapjunct_G = 1.35 - (0.9 ./ (1 + exp(-(Gapjunct_CKp - 0.16) ./ 0.12)));

% if No_CaMKII_GapG == 1
%     kCKII_Gapjunct_G = 1.0;
% end

%% Params Continued
para.KmNaip_PKA         = KmNaip_PKA;
para.kPKA_Itof          = kPKA_Itof;
para.kCKII_Itof_tau     = kCKII_Itof_tau;
para.dGkr_PKA           = dGkr_PKA;
para.Vkr_PKA            = Vkr_PKA;
para.kPKA_Iks           = kPKA_Iks;
para.kCKII_Ikur         = kCKII_Ikur;
para.kPKA_Ikur          = kPKA_Ikur;
para.kCKII_IK1_G        = kCKII_IK1_G;
% para.kPKA_IClCa         = kPKA_IClCa;
para.K_PKA_LTCC         = K_PKA_LTCC;
para.G_LTCC             = G_LTCC;
para.LTCC_junc_mode2    = LTCC_junc_mode2;
para.RyR_CKp            = RyR_CKp;
para.fPKA_TnI           = fPKA_TnI;
para.kPKA_Ina_increase  = kPKA_Ina_increase;
%% Constants
% Physical Constants
R = 8314;          % [J/kmol*K]
Frdy = 96485;      % [C/mol]
Temp = 310;        % [K]
FoRT = Frdy / R / Temp;

% Geometry
% Capacitance
Acell = 11e3;      % [um^2]
Cmem = Acell * 1e-14; % [F] 110 pF membrane capacitance

% Fractional currents in compartments
Fjunc = 0.11;
Fsl = 1 - Fjunc;
% Fjunc_CaL = 0.9;
% Fsl_CaL = 1 - Fjunc_CaL;

% Cell dimensions and volume
cellLength = 100;  % cell length [um]
cellRadius = 10.25; % cell radius [um]
Vcell = 3.14159265 * cellRadius^2 * cellLength * 1e-15; % [L]
Vmyo = 0.65 * Vcell;
Vsr = 0.035 * Vcell;
Vsl = 0.02 * Vcell;
Vjunc = 0.0539 * 0.01 * Vcell;

% Diffusion rates
diffusion_scale = 1.0;
J_ca_juncsl = diffusion_scale * 1 / 1.2134e12; % [L/msec] 8.2413e-13
J_ca_slmyo = diffusion_scale * 1 / 2.68510e11; % [L/msec] 3.2743e-12
J_na_juncsl = 1 / (1.6382e12 / 3 * 100);      % [L/msec] 1.8313e-14
J_na_slmyo = 1 / (1.8308e10 / 3 * 100);       % [L/msec] 1.6386e-12

% Protocol index
para.prot_index  = pECC(3);

if para.prot_index  == 1,     protocol = 'pace_cc';
elseif para.prot_index  == 2, protocol = 'pace_cc_erp';
% elseif para.prot_index  == 3, protocol = 'v_step';
elseif para.prot_index  == 4, protocol = 'pace_cc_dad';
end

pacing_rate = 1000*1/para.cycleLength; % (Hz) used only with 'pace_cc' and 'v_step'
RA = 0;         % Right ATRIUM
% ISO = 0;        % ISO (boolean)
CCh = 0;        % [uM]
ACh = 0;        % [uM]
K2P_cond = 1;   % K2P_cond
SK_cond = 1;    % SK_cond
SK_shift = 0;   % SK_shift

SA_par = parameters_min_opt17; % Copy values from parameters_min_opt17 to SA_par

%% Nernst Potentials
% I added nerst portions from my model- NH
Cli = 15;   % Intracellular Cl  [mM]
Clo = 150;  % Extracellular Cl  [mM]
Ko = 5.4;   % Extracellular K   [mM]
Nao = 140; %140; % Extracellular Na  [mM]
Cao = 1.8;  % Extracellular Ca  [mM]
Mgi = 1; % Intracellular Mg  [mM]



ena_junc = (1.0 / FoRT) * log(Nao ./ y(32));        % [mV]
ena_sl = (1.0 / FoRT) * log(Nao ./ y(33));          % [mV]
ek = (1.0 / FoRT) * log(Ko ./ y(35));               % [mV]
eca_junc = (1.0 / FoRT / 2) * log(Cao ./ y(36));    % [mV]
eca_sl = (1.0 / FoRT / 2) * log(Cao ./ y(37));      % [mV]
ecl = (1.0 / FoRT) * log(Cli ./ Clo);               % [mV]

%% INa
para.No_CaMKII_NaV = 0;

if para.No_CaMKII_NaV == 1
    Ina_h_shift = 0.3239; % Set Ina_h_shift value based on condition
else
    Ina_h_shift = 3.25 - 2 * 3.25 * 1 / (1 + exp(-(para.NaV_CKp - 0.12) / 0.1));
end

GNa = (0.25 * para.kPKA_Ina_increase + 1) * (1 - 0.1 * AF) * SA_par(1) * 9; % [mS/uF] Courtemanche (adjusted for dV/dt max) MOD1

am = (y(39) == -47.13) * ( 3.2 ) + (y(39) ~= -47.13) * ( 0.32 * (y(39) + 47.13) / (1.0 - exp( -0.1 * (y(39) + 47.13) ) ) );
bm = 0.08 * exp( -y(39) / 11.0 );


Vm = y(39) - Ina_h_shift;

ah = (Vm >= -40.0) * ( 0.0 ) + (y(39) < -40.0) * ( 0.135 * exp( -(Vm + 80.0 ) / 6.8 ) );
bh = (Vm >= -40.0) * ( 1.0 / ( 0.13 * ( 1.0 + exp( -(Vm + 10.66) / 11.1 ) ) ) ) + ...
    (Vm < -40.0) * ((3.56 * exp( 0.079 * Vm) + 3.1e5 * exp(0.35 * Vm)));

aj = (Vm >= -40.0) * (0.0) +(Vm < -40.0) * (( ( -127140 * exp(0.2444*Vm) - 3.474e-5 * exp(-0.04391 * Vm)) * (Vm + 37.78)) / (1.0 + exp( 0.311 * (Vm + 79.23) ) ));
bj = (Vm >= -40.0) * ((0.3 * exp(-2.535e-7*y(39))) / (1.0 + exp( -0.1 * (Vm + 32.0) ))) + ...
    (Vm < -40.0) * ((0.1212 * exp( -0.01052 * Vm )) / (1.0 + exp( -0.1378 * (Vm + 40.14) )));

ydot(1) = am * (1 - y(1)) - bm * y(1);
ydot(2) = ah * (1 - y(2)) - bh * y(2);
ydot(3) = aj * (1 - y(3)) - bj * y(3);


I_Na_junc = para.INa_Scale * Fjunc * GNa  * y(1)^3 * y(2) * y(3) * (y(39) - ena_junc);
I_Na_sl = para.INa_Scale * Fsl * GNa  * y(1)^3 * y(2) * y(3) * (y(39) - ena_sl);
I_Na = I_Na_junc + I_Na_sl;

%% INaL

GNaL = para.INaL_Scale * ((14.0 / 11.0) ./ (1 + exp(-(para.NaV_CKp - 0.12) / 0.1))) * SA_par(2) * (1.2 * 0.003) * SA_par(1) * 9 * (1 + AF * 0.25); % [mS/uF] (adjusted) MOD1 % AF changed by NH 11/21
am = 0.32 * (y(39) + 47.13) / (1 - exp(-0.1 * (y(39) + 47.13)));
bm = 0.08 * exp(-y(39) / 11);
inf = 1 ./ (1 + exp((y(39) + 91) / 6.1));
tau = 200;
ydot(40) = am * (1 - y(40)) - bm * y(40);
ydot(41) = (inf - y(41)) / tau;

I_NaL_junc = Fjunc * GNaL * y(40)^3 * y(41) * (y(39) - ena_junc);
I_NaL_sl = Fsl * GNaL * y(40)^3 * y(41) * (y(39) - ena_sl);
I_NaL = I_NaL_junc + I_NaL_sl; % No need for //#ok<NASGU> in MATLAB

%% I_NaB Na Background Current

GNaB = para.INab_Scale * 0.597e-3; % [mS/uF]

I_nabk_junc = Fjunc * GNaB * (y(39) - ena_junc);
I_nabk_sl = Fsl * GNaB * (y(39) - ena_sl);
I_nabk = I_nabk_junc + I_nabk_sl; % No need for //#ok<NASGU> in MATLAB

%% I_NaK Na/K Pump Current

IbarNaK = para.INaK_Scale * 1.0 * 1.26; % [uA/uF]
KmNaip = 11; % [mM]
KmNaip = KmNaip - para.KmNaip_PKA;
KmKo = 1.5; % [mM]

sigma = (exp(Nao / 67.3) - 1) / 7.0;
fnak = 1 / (1 + 0.1245 * exp(-0.1 * y(39) * FoRT) + 0.0365 * sigma * exp(-y(39) * FoRT));

I_nak_junc = Fjunc * IbarNaK * fnak * Ko ./ (1 + (KmNaip ./ y(32))^4) ./ (Ko + KmKo);
I_nak_sl = Fsl * IbarNaK * fnak * Ko ./ (1 + (KmNaip ./ y(33))^4) ./ (Ko + KmKo);
I_nak = I_nak_junc + I_nak_sl;

%% I_to: Transient Outward K Current

GtoFast = 1.2 * 0.75 * para.Itof_Scale * 1.2 * (1 + (0.6 - 1) * para.kPKA_Itof) * SA_par(5) * (1 - 0.7 * AF) * 0.165;
% GtoFast = 1.2 * 0.75 * para.Itof_Scale * 1.2 * (1 + (0.6 - 1) * para.kPKA_Itof) * SA_par(5) * (1 - 0.3 * AF) * 0.165;


Vm = y(39) - para.kCKII_Itof_vshift;

xtoss = 1 ./ (1 + exp(-(Vm + 1) / 11));
tauxtof = 3.5 * exp(-((Vm / 30) * (Vm / 30))) + 1.5;

Vm = y(39);
ytoss = 1 ./ (1 + exp((Vm + 40.5) / 11.5));
tauytof = (para.kCKII_Itof_tau - 1.0) * 20 ./ (1.0 + exp((Vm + 40) / -5.0)) + ...
    25.635 * exp(-(((Vm + 52.45) / 15.8827) * ((Vm + 52.45) / 15.8827))) + 24.14 - ...
    (para.kCKII_Itof_tau - 1.0) * 20 ./ (1.0 + exp((Vm + 50) / 5.0));

ydot(10) = (xtoss - y(10)) / tauxtof;
ydot(11) = (ytoss - y(11)) / tauytof;


I_tof = GtoFast * y(10) * y(11) * (y(39) - ek);
I_to = I_tof;

%% I_Kr: Rapidly Activatig K Current

% gkr = para.IKr_Scale * (1 + para.dGkr_PKA) * SA_par(6) * 0.95 * 0.035 * sqrt(Ko / 5.4);
% inf = 1 ./ (1 + exp(-(y(39) + 10 + para.Vkr_PKA) / 5.0));
% tau = 550 ./ (1 + exp((-22 - y(39)) / 9)) .* 6 ./ (1 + exp((y(39) - (-11.0)) / 9.0)) + 230 ./ (1 + exp((y(39) - (-40)) / 20));
ydot(12) =0; % (inf - y(12)) ./ tau;

% rkr = 1 ./ (1 + exp((y(39) + 74.0) / 24.0));
% I_kr = 0.9 * gkr * y(12) * rkr * (y(39) - ek);

GKr = 0.035 * sqrt(Ko / 5.4) * (1 - 0.6 * AF); % Added by NH
Ikr_tmp_state = zeros(7,1);
for i = 0:6
    Ikr_tmp_state(i+1) = y(86 + i);
end
% IKr_markov.update_states(y(39), 0.0, ek, para);
[Ikr_tmp_ydot, IKr_out] = IKr_ODE(Ikr_tmp_state, para.Vkr_PKA, para.dGkr_PKA, 0.0, y(39), GKr, ek, Ko);
for i = 0:6
    ydot(86 + i) = Ikr_tmp_ydot(i+1);
end

I_kr = 1.1 * 0.9 * para.IKr_Scale * 0.35 * IKr_out;
%% I_ks: Slowly Activating K Current

pNaK = 0.01833;
eks = (1 / FoRT) * log((Ko + pNaK * Nao) / (y(35) + pNaK * y(34)));

% OLTIKS is not defined, so assuming the values directly
% gks_junc = 1 * SA_par(7) * (1 + 1 * AF + 0.6 * para.kPKA_Iks) * 0.0035;
% gks_sl = 1 * SA_par(7) * (1 + 1 * AF + 0.6 * para.kPKA_Iks) * 0.0035;

% inf = 1 ./ (1 + exp(-(y(39) + 40 * ISO + 3.8) / 14.25));
% tau = 990.1 ./ (1 + exp(-(y(39) + 40 * ISO + 2.436) / 14.12));
% ydot(13) = (inf - y(13)) ./ tau;
ydot(13) = 0;

%IF DEFINED
% if OLTIKS
%     gks_junc = 1 * SA_par(6) * (1 + 1 * AF + 0.6 * para.kPKA_Iks) * 0.0035;
%     gks_sl = 1 * SA_par(6) * (1 + 1 * AF + 0.6 * para.kPKA_Iks) * 0.0035;
%
%     inf = 1 ./ (1 + exp(-(y(39) + 40 * ISO + 3.8) / 14.25));
%     tau = 990.1 ./ (1 + exp(-(y(39) + 40 * ISO + 2.436) / 14.12));
%     ydot(12) = (inf - y(12)) ./ tau;
%
%     I_ks_junc = para.IKs_Scale * Fjunc * gks_junc * y(12)^2 * (y(39) - eks);
%     I_ks_sl = para.IKs_Scale * Fsl * gks_sl * y(12)^2 * (y(39) - eks);
% else
%     gks_junc = 0;
%     gks_sl = 0;
%     I_ks_junc = 0;
%     I_ks_sl = 0;
% end

gks_factor = 1.5 * 8 * 0.0035;
P_g_0 = gks_factor * (0.2 + 0.2 * para.kPKA_Iks);
P_g_max = gks_factor * (0.8 + 0.5 * para.kPKA_Iks);
P_vh_0 = -1 - 10 * para.kPKA_Iks;
P_vh_max = -12 - 9 * para.kPKA_Iks;
P_tau_0 = 26 + 9 * para.kPKA_Iks;
P_tau_max = 40 + 4 * para.kPKA_Iks;
caks_junc = y(36);
caks_sl = y(37); % normal simulation
gks_junc = (P_g_0 + (P_g_max - P_g_0) ./ (1 + (150e-6 ./ caks_junc).^1.3)) * (1+1*AF); % Regulated by PKA
gks_sl = (P_g_0 + (P_g_max - P_g_0) ./ (1 + (150e-6 ./ caks_sl).^1.3)) * (1+1*AF); % Regulated by PKA
VsXs_Ca_junc = P_vh_0 + (P_vh_max - P_vh_0) ./ (1 + (350e-6 ./ caks_junc).^4.0); % Regulated by PKA
xsss_junc = 1.0 ./ (1 + exp(-(y(39) - VsXs_Ca_junc) ./ 25.0));

VsTs_Ca_junc = P_tau_0 + (P_tau_max - P_tau_0) ./ (1 + (150e-6 ./ caks_junc).^3); % Regulated by PKA
tauxs_junc = 2 * (50 + (50 + 350 * exp(-((y(39) + 30)^2) / 4000.0)) * 1.0 ./ (1 + exp(-(y(39) + VsTs_Ca_junc) / 10.0)));
VsXs_Ca_sl = P_vh_0 + (P_vh_max - P_vh_0) ./ (1 + (350e-6 ./ caks_sl).^4.0); % Regulated by PKA
xsss_sl = 1 ./ (1 + exp(-(y(39) - VsXs_Ca_sl) / 25.0));
VsTs_Ca_sl = P_tau_0 + (P_tau_max - P_tau_0) ./ (1 + (150e-6 ./ caks_sl).^3); % Regulated by PKA
tauxs_sl = 2 * (50 + (50 + 350 * exp(-((y(39) + 30)^2) / 4000.0)) * 1 ./ (1 + exp(-(y(39) + VsTs_Ca_sl) / 10.0)));

ydot(4) = (xsss_junc - y(4)) / tauxs_junc;
ydot(5) = (xsss_sl - y(5)) / tauxs_sl;
I_ks_junc = para.IKs_Scale * Fjunc * gks_junc * y(4)^2 * (y(39) - eks);
I_ks_sl = para.IKs_Scale * Fsl * gks_sl * y(5)^2 * (y(39) - eks);
I_ks = I_ks_junc + I_ks_sl;

%% I_kur: Ultra-Rapid Delayed Rectifier Outward K Current

Gkur = 1.2 * para.IKur_Scale * para.kCKII_Ikur * (1.0 + 2 * para.kPKA_Ikur) * 1.8 * SA_par(8) * (1 - 0.5 * AF) * 0.045 * (1 + 0.2 * RA);
% Gkur = 1.2 * para.IKur_Scale * para.kCKII_Ikur * (1.0 + 2 * para.kPKA_Ikur) * 1.8 * SA_par(8) * (1 - 0.2 * AF) * 0.045 * (1 + 0.2 * RA);


xkurss = 1 ./ (1 + exp((y(39) + 6) / -8.6));
tauxkur = 9.0 ./ (1 + exp((y(39) + 5) / 12.0)) + 0.5;
ykurss = 1 ./ (1 + exp((y(39) + 7.5) / 10));
tauykur = 590.0 ./ (1 + exp((y(39) + 60) / 10)) + 3050.0;

ydot(8) = (xkurss - y(8)) / tauxkur;
ydot(9) = (ykurss - y(9)) / tauykur;
I_kur = Gkur * y(8) * y(9) * (y(39) - ek);


% New IKur model, but not used in this model;

% CNZ_gkur = 1.05 * para.IKur_Scale * para.kCKII_Ikur * (1.0 + 1 * para.kPKA_Ikur) * (1 - 0.5 * AF) * 1.07 * 1.2 * 0.006398 * 0.85;
% K_Q10 = 3.0;
%
% inf = 1.0 ./ (1 + exp(-(y(39) + 6) / (8.6))); % change to simple version, haibo.  removed IKur mutations here.
% tau = (45.67 ./ (1 + exp((y(39) + 9.23) / 10.53)) + 4.27) * ( 0.55 ./ (1 + exp((y(39) + 60.87) / -22.88)) + 0.03); % incorporate deactivation Thu 20 Jul 2017 23:24:01 BST Haibo
% tau = tau / K_Q10;

% ydot(43) = (inf - y(43)) / tau;
ydot(43) = 0;
% inf = 1.0 ./ (1.0 + exp( (y(39) + 7.5) / 6.67)); % Thu 20 Jul 2017 23:24:21 BST
% tau = 1300.0 ./ (exp((y(39) - 133.0) / 97.0) * 5.04 + 0.2 * exp((y(39) + 10.9233071) / (-12.0))) + 370.0; % at 37deg  Thu 20 Jul 2017 23:24:35 BST haibo Feng et al. 1998
% ydot(44) = (inf - y(44)) / tau;
ydot(44) = 0;

% tau = 1400.0 ./ (exp((y(39) + 0.1) / 8.86) + 1.0) + 478.0 ./ (exp((y(39) - 41.0) / -3.0) + 1.0) + 6500.0; % // at 37deg  Thu 20 Jul 2017 23:24:35 BST haibo Feng et al. 1998
% ydot(45) = (inf - y(45)) / tau;
ydot(45) = 0;
% IKur = CNZ_gkur * (4.5128 + 1.899769 ./ (1.0 + exp((y(39) - 20.5232) ./ (-8.26597)))) * y(43) * (0.65 * y(44) + 0.35 * y(45)) * (y(39) - ek);


%% I_ki: Time-Independent K Current; IK1

% Gki = 0.75 * SA_par(9) * (1 + 0.4 * AF) * (2.1 * 0.0525) * sqrt(Ko / 5.4); % NH Change 11/29, the rest is added to IKach
Gki = 0.75 * SA_par(9) * (1 + 1 * AF) * (2.1 * 0.0525) * sqrt(Ko / 5.4); % NH Change 11/29, the rest is added to IKach

fracIK1_avail = 1 + (0.55 - 1) * para.kPKA_IK1; % 45 % decrease w / 100 nM ISO (Gonzales De La Fuente et al 2013)

Gki = kCKII_IK1_G * fracIK1_avail * Gki;

I_ki_j = para.IK1_Scale * 1 / 1.31 * 0.65 * (0.15 + 0.85 ./ (1 + (y(32) / 10.0) .* (y(32) / 10.0))) * Fjunc * Gki  * (y(39) - ek) / (1 + exp(0.07 * (y(39) - (ek + 6.94)))); % CRN IK1
I_ki_sl = para.IK1_Scale * 1 / 1.31 * 0.65 * (0.15 + 0.85 ./ (1 + (y(33) / 10.0) .* (y(33) / 10.0))) * Fsl * Gki  * (y(39) - ek) / (1 + exp(0.07 * (y(39) - (ek + 6.94)))); % CRN IK1

I_ki = I_ki_j + I_ki_sl;

%% I_kp: Plateau K current

gkp = 0.95 * SA_par(10) * 0.002;

kp_kp = 1 ./ (1 + exp(7.488 - y(39) / 5.98));
I_kp = para.IKp_Scale * gkp * kp_kp * (y(39) - ek);

%% I_K2P

gk2p = 0.9 * para.IK2p_Scale * 0.95 * SA_par(11) * K2P_cond * 0.0050 * (1 + AF * 2); % AF changed by NH 11/21

inf = 0.2 + 0.8 ./ (1 + exp(-(y(39) - 10 + 15 * AF) / 14));
tau = 2 + 40 ./ (1 + exp((y(39) + 25) .* (y(39) + 25) / 80));
ydot(42) = (inf - y(42)) / tau;

I_k2p = gk2p * y(42) * (y(39) - ek); % NEW K2P

%% I_k,ach: Muscarinic-Receptor-Activated K Current

Gkach_sl = 0;
Gkach_j = 0;
d_kach = 0;
r_kach = 0;
% fkach_sl = 0;
% fkach_j = 0;

% NEED to reassess once I consider ACH and CCH
if CCh == 0 && ACh == 0
    Gkach_sl = SA_par(11) * 0;
    Gkach_j = SA_par(11) * 0;
    d_kach = 0;
    r_kach = 0;
elseif CCh > 0 && ACh == 0
    % Na-dependence
    fkach_sl = 1.5 / (1 + (9 / y(33))^4); % Na-dependence sl
    fkach_j = 1.5 / (1 + (9 / y(32))^4); % Na-dependence j
    if AF == 0
        Gkach_sl = SA_par(11) * (RA * 1 + (1 - RA) * 0.3) * 0.10 * (1 + fkach_sl) / 0.9411765 * sqrt(Ko / 5.4);
        Gkach_j = SA_par(11) * (RA * 1 + (1 - RA) * 0.3) * 0.10 * (1 + fkach_j) / 0.9411765 * sqrt(Ko / 5.4);
    else
        Gkach_sl = SA_par(11) * (RA * 0.5 + (1 - RA) * 0.3) * 0.10 / 0.9411765 * sqrt(Ko / 5.4);
        Gkach_j = SA_par(11) * (RA * 0.5 + (1 - RA) * 0.3) * 0.10 / 0.9411765 * sqrt(Ko / 5.4);
    end
    % V-dependence
    r_kach = 0.055 + 0.4 / (1 + exp((y(39) - ek + 9.53) / 17.18));
    % Dose-dependence
    d_kach = CCh / (CCh + 0.125);
elseif CCh == 0 && ACh > 0
    % Na-dependence
    fkach_sl = 1.5 / (1 + (9 / y(33))^4); % Na-dependence sl
    fkach_j = 1.5 / (1 + (9 / y(32))^4); % Na-dependence j
    if AF == 0
        Gkach_sl = SA_par(11) * (RA * 1 + (1 - RA) * 0.3) * 5 * 0.10 * (1 + fkach_sl) / 0.9411765 * sqrt(Ko / 5.4);
        Gkach_j = SA_par(11) * (RA * 1 + (1 - RA) * 0.3) * 5 * 0.10 * (1 + fkach_j) / 0.9411765 * sqrt(Ko / 5.4);
    else
        Gkach_sl = SA_par(11) * (RA * 0.5 + (1 - RA) * 0.3) * 5 * 0.10 / 0.9411765 * sqrt(Ko / 5.4);
        Gkach_j = SA_par(11) * (RA * 0.5 + (1 - RA) * 0.3) * 5 * 0.10 / 0.9411765 * sqrt(Ko / 5.4);
    end
    % V-dependence
    r_kach = 0.08 + 0.4 / (1 + exp((y(39) + 91) / 12));
    % Dose-dependence
    d_kach = 1 / (1 + (0.03 / ACh)^2.1);
end

I_kach_sl = (para.IKach_Scale * Fsl * Gkach_sl  * d_kach * r_kach * (y(39) - ek));
I_kach_j = (para.IKach_Scale * Fjunc * Gkach_j  * d_kach * r_kach * (y(39) - ek));

% if AF == 1
%     % Na-dependence
%     fkach_sl = 1.5 / (1 + (9 / y(33 - 1))^4); % Na-dependence sl
%     fkach_j = 1.5 / (1 + (9 / y(32 - 1))^4); % Na-dependence j
% %     if AF == 0
% %         Gkach_sl = SA_par(11) * (RA * 1 + (1 - RA) * 0.3) * 5 * 0.10 * (1 + fkach_sl) / 0.9411765 * sqrt(Ko / 5.4);
% %         Gkach_j = SA_par(11) * (RA * 1 + (1 - RA) * 0.3) * 5 * 0.10 * (1 + fkach_j) / 0.9411765 * sqrt(Ko / 5.4);
% %     else
%         Gkach_sl = SA_par(11) * (RA * 0.5 + (1 - RA) * 0.3) * 5 * 0.10 / 0.9411765 * sqrt(Ko / 5.4);
%         Gkach_j = SA_par(11) * (RA * 0.5 + (1 - RA) * 0.3) * 5 * 0.10 / 0.9411765 * sqrt(Ko / 5.4);
% %     end
%     % V-dependence
%     r_kach = 0.08 + 0.4 / (1 + exp((y(39) + 91) / 12));
%     % Dose-dependence
%     d_kach = 1 / (1 + (0.03 / 0.06)^2.1); % Ach - 0.04
% 
%     I_kach_sl_constiutively = para.IKach_Scale * Fsl * Gkach_sl  * d_kach * r_kach * (y(39) - ek);
%     I_kach_j_constiutively = para.IKach_Scale * Fjunc * Gkach_j  * d_kach * r_kach * (y(39) - ek);
% 
% else
I_kach_sl_constiutively = 0;
I_kach_j_constiutively = 0;
% end

I_kach_sl = (I_kach_sl + I_kach_sl_constiutively) ;
I_kach_j = (I_kach_j + I_kach_j_constiutively);
I_kach = I_kach_sl + I_kach_j;


%% I_sk: Small-Conductance Ca-Activated K Current

par_sk = [0.0506381114404388, 0.273335569451572, 2.96381060498817, 0.199981221802789, 0.279328126521496, -86.9289059836381, 0.00636311816933264, 5.22915055145375];

gsk = (0.95 * para.ISK_Scale * 0.85 * SA_par(13) * SK_cond * par_sk(1) * (1 + AF * 1)); % Updated AF, NH 11/21

kdsk = SA_par(14) * (10^(SK_shift - 3.45)); % kdsk = (10^(ISK_shift-3.3)); % (mM)
gsk_ca_junc = 1.0 / (1 + exp((log10(kdsk) - log10(y(36))) / 0.3));
gsk_ca_sl = 1.0 / (1 + exp((log10(kdsk) - log10(y(37))) / 0.3));

gsk_vm = par_sk(2) / (1 + exp((y(39) - ek + par_sk(3)) * par_sk(4))) + par_sk(5) / (1 + exp((-(y(39) - ek + par_sk(6)) * par_sk(7))));

I_sk_junc = Fjunc * gsk * gsk_ca_junc * gsk_vm * (y(39) - ek);
I_sk_sl = Fsl * gsk * gsk_ca_sl * gsk_vm * (y(39) - ek);
I_sk = I_sk_junc + I_sk_sl;

%% I_ClCa and I_Clbk: Ca-activated and Background Cl Currents

GClCa = 1.178030 * SA_par(15) * 0.0548;   % [mS/uF]
KdClCa = (1 + (0.704 - 1) * kPKA_IClCa) * 100e-3;       % [mM]
% GClB = 9e-3;    % [mS/uF] MOD1

I_ClCa_junc = para.IClCa_Scale * Fjunc * GClCa / (1 + KdClCa / y(36)) * (y(39) - ecl);
I_ClCa_sl = para.IClCa_Scale * Fsl * GClCa / (1 + KdClCa / y(37)) * (y(39) - ecl);

I_ClCa = I_ClCa_junc + I_ClCa_sl;

I_Clbk = para.IClb_Scale * 1.05 * 0.19e-3 * (y(39) - ecl) / (1 - 0.94 * exp(2.5e-4 * (y(39) - ecl)));

GClCFTR = 0; % 4.9e-3*ISO;     % [mS/uF]
I_ClCFTR = GClCFTR * (y(39) - ecl);

%% I_Ca: L-type Ca Current

jLTCC_tmp_state = zeros(10,1);
jLTCC_tmp_state_m2 = zeros(10,1);


slLTCC_tmp_state = zeros(10,1);
slLTCC_tmp_state_m2 = zeros(10,1);


for i = 0:9
    jLTCC_tmp_state(i+1) = y(46 + i);
    slLTCC_tmp_state(i+1) = y(46 + 10 + i);
    jLTCC_tmp_state_m2(i+1) = y(46 + 20 + i);
    slLTCC_tmp_state_m2(i+1) =  y(46 + 10 + 20 + i);
end

[jLTCC_tmp_ydot, jLTCC_out] =  LTCC_mode_1(jLTCC_tmp_state, y(36), y(39),  y(32), y(35), para.K_PKA_LTCC, Frdy, FoRT,  Cao, Nao, Ko);
[jLTCC_tmp_ydot_m2, jLTCC_out_m2] =  LTCC_mode_2(jLTCC_tmp_state_m2, y(36), y(39),  y(32), y(35), para.K_PKA_LTCC, Frdy, FoRT,  Cao, Nao, Ko);


[slLTCC_tmp_ydot, slLTCC_out] =  LTCC_mode_1(slLTCC_tmp_state, y(37), y(39),  y(33), y(35), para.K_PKA_LTCC, Frdy, FoRT,  Cao, Nao, Ko);
[slLTCC_tmp_ydot_m2, slLTCC_out_m2] =  LTCC_mode_2(slLTCC_tmp_state_m2, y(37), y(39),  y(33), y(35), para.K_PKA_LTCC, Frdy, FoRT,  Cao, Nao, Ko);

for i = 0:9
    ydot(46 + i) = jLTCC_tmp_ydot(i+1);
    ydot(46 + 10 + i) = slLTCC_tmp_ydot(i+1);
    ydot(46 + 20 + i) = jLTCC_tmp_ydot_m2(i+1);
    ydot(46 + 20 + 10 + i) = slLTCC_tmp_ydot_m2(i+1);
end


scale = para.ICaL_Scale * 1.02 * 1.209 * 1.03 * 0.85 * 0.9 * 1.037031 * 1.05; %1.3 * 0.93;


I_Ca_junc = scale  * para.G_LTCC * (0.9 * ( (1 - para.LTCC_junc_mode2) * jLTCC_out.I_Ca_junc_m1  + para.LTCC_junc_mode2 * jLTCC_out_m2.I_Ca_junc_m1) );
I_CaNa_junc = scale  * para.G_LTCC * (0.9 * ( (1 - para.LTCC_junc_mode2) * jLTCC_out.I_Na_junc_m1  + para.LTCC_junc_mode2 * jLTCC_out_m2.I_Na_junc_m1) );
I_Ca_sl =  scale  * para.G_LTCC * (0.1 * ( (1 - para.LTCC_sl_mode2) * slLTCC_out.I_Ca_junc_m1  + para.LTCC_sl_mode2 * slLTCC_out_m2.I_Ca_junc_m1) );
I_CaNa_sl =  scale  * para.G_LTCC * (0.1 * ( (1 - para.LTCC_sl_mode2) * slLTCC_out.I_Na_junc_m1  + para.LTCC_sl_mode2 * slLTCC_out_m2.I_Na_junc_m1) );
I_CaK =  scale  * para.G_LTCC * (0.9 * ( (1 - para.LTCC_junc_mode2) * jLTCC_out.I_K_junc_m1  + para.LTCC_junc_mode2 * jLTCC_out_m2.I_K_junc_m1) ...
    + 0.1 * ( (1 - para.LTCC_sl_mode2) * slLTCC_out.I_K_junc_m1  + para.LTCC_sl_mode2 * slLTCC_out_m2.I_K_junc_m1) );



I_Ca = I_Ca_junc + I_Ca_sl;
I_CaNa = I_CaNa_junc + I_CaNa_sl;
I_Catot = I_Ca + I_CaK + I_CaNa;

%% I_cabk: Ca Background Current

GCaB = para.ICab_Scale * 6.0643e-4; % [uA/uF]

if MOD_ind == 3
%     GCaB = para.ICab_Scale * 6.0643e-4 * (1 + 0.5 * AF); % [uA/uF] MOD3
      GCaB = para.ICab_Scale * 6.0643e-4 ; % [uA/uF] MOD3 % Removed this AF Flag NH 11/21; Not sure why it is here

end

I_cabk_junc = Fjunc * GCaB * (y(39) - eca_junc);
I_cabk_sl = Fsl * GCaB * (y(39) - eca_sl);
I_cabk = I_cabk_junc + I_cabk_sl;

%% I_pca: Sarcolemmal Ca Pump Current

IbarSLCaP = para.ICap_Scale * SA_par(19) * 0.0471; % [uA/uF]

if MOD_ind == 3
    IbarSLCaP = para.ICap_Scale * 1.0 * 1.022257 * 0.9 * SA_par(19) * 0.0471 * 2; % [uA/uF] MOD3
end

KmPCa = 0.5e-3;    % [mM]
% Q10SLCaP = 2.35;   % [none]


I_pca_junc = Fjunc * IbarSLCaP * (y(36)^1.6) / ((KmPCa^1.6) + (y(36)^1.6));  % Q10SLCaP^Qpow* was removed
I_pca_sl = Fsl * IbarSLCaP * (y(37)^1.6) / ((KmPCa^1.6) + (y(37)^1.6));  % Q10SLCaP^Qpow* was removed
I_pca = I_pca_junc + I_pca_sl;

%% I_ncx: Na/Ca Exchanger flux

IbarNCX =  1 * SA_par(20) * (1 + 0.4 * AF) * 3.15;

if (MOD_ind == 3)
%     IbarNCX = 1.1 * 1.05 * 1.2 * SA_par(20) * (1 + 0.8 * AF) * 3.15; % [uA/uF] MOD3
    IbarNCX = 1.1 * 1.05 * 1.2 * SA_par(20) * (1 + 0.4 * AF) * 3.15; % [uA/uF] MOD3; %NH Changes to 40% same as SK paper

end

KmCai = 3.59e-3;    % [mM]
KmCao = 1.3;        % [mM]
KmNai = 12.29;      % [mM]
KmNao = 87.5;       % [mM]
ksat = 0.27;        % [none]
nu = 0.35;          % [none]
Kdact = 0.384e-3;   % [mM]
% Q10NCX = 1.57;      % [none]

Ka_junc = 1 / (1 + (Kdact / y(36)) * (Kdact / y(36)));
Ka_sl = 1 / (1 + (Kdact / y(37)) * (Kdact / y(37)));
s1_junc = exp(nu * y(39) * FoRT) * y(32)^3 * Cao;
s1_sl = exp(nu * y(39) * FoRT) * y(33)^3 * Cao;
Nao_3 = Nao^3;
Naj_3 = y(32)^3;
NaSL_3 = y(33)^3;
s2_junc = exp((nu - 1) * y(39) * FoRT) * Nao_3 * y(36);
s3_junc = KmCai * Nao_3 * (1 + (y(32) / KmNai)^3) + KmNao^3 * y(36) * (1 + y(36) / KmCai) + KmCao * Naj_3 + Naj_3 * Cao + Nao_3 * y(36);
s2_sl = exp((nu - 1) * y(39) * FoRT) * Nao_3 * y(37);
s3_sl = KmCai * Nao_3 * (1 + (y(33) / KmNai)^3) + KmNao^3 * y(37) * (1 + y(37) / KmCai) + KmCao * NaSL_3 + NaSL_3 * Cao + Nao_3 * y(37);

I_ncx_junc = para.INCX_Scale * Fjunc * IbarNCX * Ka_junc * (s1_junc - s2_junc) / s3_junc / (1 + ksat * exp((nu - 1) * y(39) * FoRT)); %Q10NCX^Qpow*removed
I_ncx_sl = para.INCX_Scale * Fsl * IbarNCX * Ka_sl * (s1_sl - s2_sl) / s3_sl / (1 + ksat * exp((nu - 1) * y(39) * FoRT)); %Q10NCX^Qpow*removed
I_ncx = I_ncx_junc + I_ncx_sl;

%% SR fluxes: Calcium Uptake, Release, and Leak

Vmax_SRCaP = 1.5 * SA_par(21) * 5.3114e-3 ; % [mM/msec] (286 umol/L cytosol/sec)

if MOD_ind == 3
%     Vmax_SRCaP = 1.1 * 1.5 * SA_par(21) * 5.3114e-3 * (1 - 0.25 * AF); % [mM/msec] MOD3
        Vmax_SRCaP = 1.1 * 1.5 * SA_par(21) * 5.3114e-3; % [mM/msec] MOD3; %NH removed AF Flag not sure what it is for

end

% Q10SRCaP = 2.6;          % [none]
Kmf = PLB_kmf_Scale * (2.5) * 0.246e-3 ; % [mM]
Kmr = 1.7;               % [mM]L cytosol
hillSRCaP = 1.787;       % [mM]
ks = 1.273543 * 1.1 * 25 ; % [1/ms] %flux

koCa = 10; % [mM^-2 1/ms]  % ISO will be modeled by PKA
kom = 0.06;              % [1/ms]
kiCa = 0.5;              % [1/mM/ms]
kim = 0.005;             % [1/ms]

ec50SR = 0.45 ;           % [mM]

MaxSR = 15;              % [mM]
MinSR = 1;               % [mM]

kleak = (1.0 / 3.0 + 10 * para.RyR_CKp / 3.0) * SA_par(23) * (1.0 + 0.25 * AF) * 5.348e-6;

kCaSR = MaxSR - (MaxSR - MinSR) / (1 + ((ec50SR / y(31))^2.5));
koSRCa = RyR_koSRCa_Scale * koCa / kCaSR;
kiSRCa = kiCa * kCaSR;

KoSR_Ca_2 = koSRCa * y(36)^2; % Original calculation

%RyR
RI = 1 - y(14) - y(15) - y(16);
ydot(14) = (kim * RI - kiSRCa * y(36) * y(14)) - (KoSR_Ca_2 * y(14) - kom * y(15)); % R
ydot(15) = (KoSR_Ca_2 * y(14) - kom * y(15)) - (kiSRCa * y(36) * y(15) - kim * y(16)); % O
ydot(16) = (kiSRCa * y(36) * y(15) - kim * y(16)) - (kom * y(16) - KoSR_Ca_2 * RI); % 0
J_SRCarel = para.Jrel_Scale * ks * y(15) * (y(31) - y(36)); % [mM/ms] RYR Release

%SERCA
J_serca = para.Jserca_Scale * Vmax_SRCaP * ((y(38) / Kmf)^hillSRCaP - (y(31) / Kmr)^hillSRCaP) / (1 + (y(38) / Kmf)^hillSRCaP + (y(31) / Kmr)^hillSRCaP); % [mM/ms] Q10SRCaP^Qpow* was removed

% J_serca = 1*Q10SRCaP^Qpow*Vmax_SRCaP*((y(38)/Kmf)^hillSRCaP-(y(31)/Kmr)^hillSRCaP)/(1+(y(38)/Kmf)^hillSRCaP+(y(31)/Kmr)^hillSRCaP);
J_SRleak = para.Jleak_Scale * kleak * (y(31) - y(36)); % [mM/ms]
%% Sodium and Calcium Buffering

BufferScaling_cytosol =  1 * para.Cytosol_Buffer_Scale;

Bmax_Naj = 7.561;       % [mM] // Na buffering
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3 * BufferScaling_cytosol;  % [mM]                      % TnC low affinity
koff_tncl = para.fPKA_TnI * 19.6e-3; % [1/ms]
kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3 * BufferScaling_cytosol; % [mM]                      % TnC high affinity
koff_tnchca = 0.032e-3; % [1/ms]
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms]
kon_tnchmg = 3e-3;      % [1/mM/ms]
% Bmax_CaM = 24e-3 * BufferScaling_cytosol;     % [mM] **? about setting to 0 in c-code**   % CaM buffering
% koff_cam = 238e-3;      % [1/ms]
% kon_cam = 34;           % [1/mM/ms]
Bmax_myosin = 140e-3 * BufferScaling_cytosol; % [mM]                      % Myosin buffering
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 19 * .9e-3 * BufferScaling_cytosol; % [mM] (Bers text says 47e-3) 19e-3
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]

BufferScaling = 1 * para.Cleft_Buffer_Scale; % scale both cleft and sub-membrane
Bmax_SLlowsl = BufferScaling * 37.4e-3 * Vmyo / Vsl;    % [mM]   	% SL buffering
Bmax_SLlowj = BufferScaling * 4.6e-3 * Vmyo / Vjunc * 0.1; % [mM]
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = BufferScaling * 13.4e-3 * Vmyo / Vsl;   % [mM]
Bmax_SLhighj = BufferScaling * 1.65e-3 * Vmyo / Vjunc * 0.1; % [mM]
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]

Bmax_Csqn = 140e-3 * Vmyo / Vsr * para.SR_Buffer_Scale;        % [mM]      % Csqn buffering
koff_csqn = 65;         % [1/ms]
kon_csqn = 100;         % [1/mM/ms]

%% Junctional and SL Na Buffers

ydot(17) = kon_na * y(32) * (Bmax_Naj - y(17)) - koff_na * y(17); % NaBj      [mM/ms]
ydot(18) = kon_na * y(33) * (Bmax_Nasl - y(18)) - koff_na * y(18); % NaBsl     [mM/ms]

% Cytosolic Ca Buffers
ydot(19) = kon_tncl * y(38) * (Bmax_TnClow - y(19)) - koff_tncl * y(19);  % TnCL      [mM/ms]
ydot(20) = kon_tnchca * y(38) * (Bmax_TnChigh - y(20) - y(21)) - koff_tnchca * y(20); % TnCHc     [mM/ms]
ydot(21) = kon_tnchmg * Mgi * (Bmax_TnChigh - y(20) - y(21)) - koff_tnchmg * y(21); % TnCHm     [mM/ms]
ydot(22) = 0.0; % kon_cam * y(37) * (Bmax_CaM - y(21)) - koff_cam * y(21);       % CaM       [mM/ms]  % With signaling, this part is taken care of by CaM modules
ydot(23) = kon_myoca * y(38) * (Bmax_myosin - y(23) - y(24)) - koff_myoca * y(23); % Myosin_ca [mM/ms]
ydot(24) = kon_myomg * Mgi * (Bmax_myosin - y(23) - y(24)) - koff_myomg * y(24); % Myosin_mg [mM/ms]
ydot(25) = kon_sr * y(38) * (Bmax_SR - y(25)) - koff_sr * y(25);          % SRB       [mM/ms]
J_CaB_cytosol = ydot(19) + ydot(20) + ydot(22) + ydot(23) + ydot(25);

% Junctional and SL Ca Buffers
ydot(26) = kon_sll * y(36) * (Bmax_SLlowj - y(26)) - koff_sll * y(26); % SLLj      [mM/ms]
ydot(27) = kon_sll * y(37) * (Bmax_SLlowsl - y(27)) - koff_sll * y(27); % SLLsl     [mM/ms]
ydot(28) = kon_slh * y(36) * (Bmax_SLhighj - y(28)) - koff_slh * y(28); % SLHj      [mM/ms]
ydot(29) = kon_slh * y(37) * (Bmax_SLhighsl - y(29)) - koff_slh * y(29); % SLHsl     [mM/ms]
J_CaB_junction = ydot(26) + ydot(28);
J_CaB_sl = ydot(27) + ydot(29);

%% Ion concentrations
% SR Ca Concentrations
ydot(30) = kon_csqn * y(31) * (Bmax_Csqn - y(30)) - koff_csqn * y(30); % Csqn      [mM/ms]
ydot(31) = J_serca - (J_SRleak * Vmyo / Vsr + J_SRCarel) - ydot(30); % Ca_sr     [mM/ms]

% Sodium Concentrations
I_Na_tot_junc = I_Na_junc + I_nabk_junc + 3 * I_ncx_junc + 3 * I_nak_junc + I_CaNa_junc + I_NaL_junc; % [uA/uF]
I_Na_tot_sl = I_Na_sl + I_nabk_sl + 3 * I_ncx_sl + 3 * I_nak_sl + I_CaNa_sl + I_NaL_sl; % [uA/uF]

ydot(32) = -I_Na_tot_junc * Cmem / (Vjunc * Frdy) + J_na_juncsl / Vjunc * (y(33) - y(32)) - ydot(17); %Naj
ydot(33) = -I_Na_tot_sl * Cmem / (Vsl * Frdy) + J_na_juncsl / Vsl * (y(32) - y(33)) + J_na_slmyo / Vsl * (y(34) - y(33)) - ydot(18); %Nasl
ydot(34) = (J_na_slmyo / Vmyo * (y(33) - y(34))) * (1 - para.Na_clamp); % [mM/msec] %Nai

% ydot(32) = 0; %Nasl
% ydot(33) = 0; %Naj
% ydot(34) = 0; % [mM/msec] %Nai
%
% Potassium Concentration
I_K_tot = I_to + I_kr + I_ks + I_ki - 2 * I_nak + I_CaK + I_kp + I_kur + I_kach + I_k2p + I_sk; % [uA/uF] //SVP: added IKur
ydot(35) = 0; % -I_K_tot*Cmem/(Vmyo*Frdy);         % [mM/msec]

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc + I_cabk_junc + I_pca_junc - 2 * I_ncx_junc;           % [uA/uF]
I_Ca_tot_sl = I_Ca_sl + I_cabk_sl + I_pca_sl - 2 * I_ncx_sl;    % [uA/uF]
ydot(36) = -I_Ca_tot_junc * Cmem / (Vjunc * 2 * Frdy) + J_ca_juncsl / Vjunc * (y(37) - y(36)) - J_CaB_junction + (J_SRCarel) * Vsr / Vjunc + J_SRleak * Vmyo / Vjunc; % Ca_j
ydot(37) = -I_Ca_tot_sl * Cmem / (Vsl * 2 * Frdy) + J_ca_juncsl / Vsl * (y(36) - y(37)) + J_ca_slmyo / Vsl * (y(38) - y(37)) - J_CaB_sl; % Ca_sl
ydot(38) = -J_serca * Vsr / Vmyo - J_CaB_cytosol + J_ca_slmyo / Vmyo * (y(37) - y(38)); % [mM/msec]

para.Ca_clamp = 0;

if para.Ca_clamp == 1 || para.Ca_clamp == 2
    ydot(35) = 0;
    ydot(36) = 0;
    ydot(37) = 0;
elseif para.Ca_clamp == 3
    ydot(37) = 0;
end

switch lower(protocol)

    case 'pace_cc' % pace w/ current injection at rate 'rate'
        rate = pacing_rate*1e-3;
        if mod(t,1/rate) <= 5
            I_app = 12.5;
        else
            I_app = 0.0;
        end


    case 'pace_cc_erp' % pace ERP protocol
        if t <= 5
            I_app = 12.5;
        elseif t > 5 && t <= rec_interval
            I_app = 0.0;
        elseif t > rec_interval && t <= rec_interval+5
            %I_app = 12.5;
            %I_app = 12.5*0.18; % DTE (1 Hz, nSR)
            I_app = (12.5*0.18)*2; % 2*DTE
        else
            I_app = 0.0;
        end

    case 'pace_cc_dad', % pace w/ current injection at rate 'rate'
%         if t < 10e3,
        if t < 470e3,

            rate = pacing_rate*1e-3;
        else
            rate = 0;
        end
        if mod(t,1/rate) <= 5, %this means there is 5ms of stim
            I_app = 12.5;
        else
            I_app = 0.0;
        end

end

%% Membrane Potentials
I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;
I_Cl_tot = I_ClCa + I_Clbk + I_ClCFTR;
I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
I_tot = I_Na_tot + I_Cl_tot + I_Ca_tot + I_K_tot;

ydot(39) = -(I_tot - I_app);


%% Currents
if para.Currents_record == 1 || para.Output_current == 1
    global tStep tArray ICa_store Ito_store INa_store IK1_store
    global Jserca_store IKs_store Jleak_store Incx_store
    global INaK_store Ikr_store INabk_store
    global Lmyo_store Fmyo_store Vmax_store
    global IKp_store Iclca_store Iclbk_store Ipca_store Icabk_store Cai_store
    global I_J_SRCarelstore
    global I_app_store I_NaLstore I_kurstore I_k2pstore I_kistore I_kachstore
    global I_skstore I_ClCFTRstore I_pcastore


    % Time
    if t > tArray(tStep)    % Roughly eliminates data from rejected time steps
        tStep = tStep + 1;
    end
    tArray(tStep) = t;


    % Currents
    INa_store(tStep) = I_Na;
    INabk_store(tStep) = I_nabk;
    INaK_store(tStep) = I_nak;
    Ikr_store(tStep) = I_kr;
    IKs_store(tStep) = I_ks;
    IKp_store(tStep) = I_kp;
    Ito_store(tStep) = I_to;
    IK1_store(tStep) = I_ki;
    Iclca_store(tStep) = I_ClCa;     % Total I_ClCa
    Iclbk_store(tStep) = I_Clbk;     % I_Clbk
    ICa_store(tStep) = I_Catot;
    Incx_store(tStep) = I_ncx;
    Ipca_store(tStep) = I_pca;
    Icabk_store(tStep) = I_cabk;
    % Ca
    Jleak_store(tStep,1) = J_SRCarel*Vsr/Vmyo + J_SRleak;   % Total Jleak [mmol/L cyt/ms]
    Jleak_store(tStep,2) = J_SRleak;                        % Passive SR leak only [mmol/L cyt/ms]
    Jserca_store(tStep) = J_serca;
    Cai_store(tStep) = y(38);

    % % Vmax
    Vmax_store(tStep) = ydot(39);

    % I Added these
    I_app_store(tStep) = I_app;
    I_NaLstore(tStep) = I_NaL;
    I_kurstore(tStep) = I_kur;
    I_k2pstore(tStep) = I_k2p;
    I_kistore(tStep) = I_ki;
    I_kachstore(tStep) = I_kach;
    I_skstore(tStep) = I_sk;
    I_ClCFTRstore(tStep) = I_ClCFTR;
    I_pcastore(tStep) = 0;  % not sure what this is
    I_J_SRCarelstore(tStep) = J_SRCarel;
end




end
