function ydot = NH_BetaAR(~,y,pin)


%Nate params match
%int conditions match
%equations now should match ; made changes

%% Assign passed in params

ISO = pin(1); % Ltot
LCCtot = pin(2);
RyRtot = pin(3);
PLBtot = pin(4);
TnItot = pin(5);
IKstot = pin(6);
ICFTRtot = pin(7);
PP1tot = pin(8);
PLMtot = pin(9);
myotot = pin(10);
IKrtot = pin(11);
IKurtot = pin(12);
INatot = pin(13);
IClCatot = pin(14);
Itotot = pin(15);
IK1tot = pin(16);
AF_index = pin(17); % AF_index

%% State variables
LR=y(1); % bound receptor   
LRG=y(2); % G protein-associated, ligand bound receptor
RG=y(3); % 
b1AR_S464=y(4);
b1AR_S301=y(5);
GsaGTPtot=y(6);
GsaGDP=y(7);
Gsby=y(8);
AC_GsaGTP=y(9);
PDE4p=y(10);
cAMPtot=y(11);
RC_I=y(12);
RCcAMP_I=y(13);
RCcAMPcAMP_I=y(14);
RcAMPcAMP_I=y(15);
PKACI=y(16);
PKACI_PKI=y(17);
RC_II=y(18);
RCcAMP_II=y(19);
RCcAMPcAMP_II=y(20);
RcAMPcAMP_II=y(21);
PKACII=y(22);
PKACII_PKI=y(23);
I1p_PP1=y(24); % output CaMKII
I1ptot=y(25);
PLBp=y(26); % output
PLMp=y(27); % output
LCCap=y(28); % output
LCCbp=y(29); % output
RyRp=y(30); % output
TnIp=y(31); % output
KSp=y(32); % output
KRp=y(33); % output
KURp=y(34); % output
CFTRp=y(35); % output
ClCap=y(36); % output
myop=y(37); % output
NAp=y(38); % output
TOp=y(39); % output
K1p=y(40); % output

% Haibo has PDE3 but thats 41 params?
PDE3p = y(41);


%% Scale Values

% LCCtot_Scale        = 1.0;
% RyRtot_Scale        = 1.0;
% PLBtot_Scale        = 1.0;
% TnItot_Scale        = 1.0;
% PLMtot_Scale        = 1.0;
% b1ARtot_Scale       = 1.0;
% Gstot_Scale         = 1.0;
% ACtot_Scale         = 1.0;
% ATP_Scale           = 1.0;
% PDE3tot_Scale       = 1.0;
% PDE4tot_Scale       = 1.0;
% PKItot_Scale        = 1.0;
% PKAIItot_Scale      = 1.0;
% I1tot_Scale         = 1.0;
% PP1tot_Scale        = 1.0;
% PKACII_LCCtot_Scale = 1.0;
% PP1_LCC_Scale       = 1.0;
% PP2A_LCC_Scale      = 1.0;
% PKAIIryrtot_Scale   = 1.0;
% PP1ryr_Scale        = 1.0;
% PP2Aryr_Scale       = 1.0;
PKAII_ikstot_Scale  = 1.0;
% Km_AC_basal_Scale   = 1.0;
% Kd_AC_Gsa_Scale     = 1.0;
% Km_PDE3_cAMP_Scale  = 1.0;
% Km_PDE4_cAMP_Scale  = 1.0;
% Km_PKA_I1_Scale     = 1.0;
% Km_PP2A_I1_Scale    = 1.0;
% Km_PKA_PLB_Scale    = 1.0;
% Km_PP1_PLB_Scale    = 1.0;
% Km_PKA_LCC_Scale    = 1.0;
% Km_PP1_LCC_Scale    = 1.0;
% Km_PP2A_LCC_Scale   = 1.0;
% Km_pka_ryr_Scale    = 1.0;
% Km_pp1_ryr_Scale    = 1.0;
% Km_pp2a_ryr_Scale   = 1.0;
% Km_PKA_TnI_Scale    = 1.0;
% Km_PP2A_TnI_Scale   = 1.0;
Km_pka_iks_Scale    = 1.0;
Km_pp1_iks_Scale    = 1.0;
K_rate_Scale = 1.0;
% PP1_ik1tot_Scale = 1.0;
% PP1_itotot_Scale = 1.0;
% PP1_inatot_Scale = 1.0;
% PP1_ClCatot_Scale = 1.0;
% PP1_CFTRtot_Scale = 1.0;
% PP1_KURtot_Scale = 1.0;
% PP1_ikrtot_Scale = 1.0;
PP1_ikstot_Scale = 1.0;
% PP2A_TnI_Scale = 1.0;
% PP2A_I1_Vmax_Scale = 1.0;
%% Drug Concentrations

%ISO = Ltot; % (uM) isoproterenol concentration - Ltot
FSK = 0; % (uM) forskolin concentration
IBMX = 0; % (uM) IBMX concentration     % potential issue with the IBMX in this model
										% // 16:51:18, Wed, 08-November-2017, By Haibo
%% -------- SIGNALING MODEL -----------

ydot = zeros(size(y));
%% b-AR module

b1ARtot         = 0.028; % (uM) total b1-AR protein % RABBIT (not used in the original code)
% b1ARtot         = 0.00528; % MOUSE   

kf_LR           = 1;              % (1/[uM ms]) forward rate for ISO binding to b1AR
kr_LR           = 0.285;          % (1/ms) reverse rate for ISO binding to b1AR
kf_LRG          = 1;              % (1/[uM ms]) forward rate for ISO:b1AR association with Gs
kr_LRG          = 0.062;          % (1/ms) reverse rate for ISO:b1AR association with Gs
kf_RG           = 1;              % (1/[uM ms]) forward rate for b1AR association with Gs
kr_RG           = 33;             % (1/ms) reverse rate for b1AR association with Gs

Gstot           = 3.83;           % (uM) total Gs protein
k_G_act         = 16e-3;          % (1/ms) rate constant for Gs activation
k_G_hyd         = 0.8e-3;         % (1/ms) rate constant for G-protein hydrolysis: GTP->GDP
k_G_reassoc     = 1.21;           % (1/[uM ms]) rate constant for G-protein reassociation

kf_bARK         = 1.1e-6;         % (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
kr_bARK         = 2.2e-6;         % (1/ms) reverse rate for b1AR phosphorylation by b1ARK
kf_PKA          = 3.6e-6;         % (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
kr_PKA          = 2.2e-6;         % (1/ms) reverse rate for b1AR phosphorylation by PKA

%{
% update comments // 15:29:24, Mon, 13-November-2017, By Haibo
$ The dynamics of these LR, LRG, RG is given in detail in:
%Yang, Jason H., and Jeffrey J. Saucerman. 
“Phospholemman Is a Negative Feed-Forward Regulator of Ca2+ 
in β-Adrenergic Signaling, Accelerating β-Adrenergic Inotropy.” 
Journal of Molecular and Cellular Cardiology 52, no. 5 (May 1, 2012): 1048–55. 
https://doi.org/10.1016/j.yjmcc.2011.12.015.
%}



b1ARact = b1ARtot - b1AR_S464 - b1AR_S301;
b1AR = b1ARact - LR - LRG - RG;
Gs = Gstot - LRG - RG - Gsby;

dLR = kf_LR*ISO*b1AR - kr_LR*LR + kr_LRG*LRG - kf_LRG*LR*Gs;
dLRG = kf_LRG*LR*Gs - kr_LRG*LRG - k_G_act*LRG;
dRG = kf_RG*b1AR*Gs - kr_RG*RG - k_G_act*RG;

bARK_desens = kf_bARK*(LR+LRG);
bARK_resens = kr_bARK*b1AR_S464;
PKA_desens = kf_PKA*PKACI*b1ARact;
PKA_resens = kr_PKA*b1AR_S301;
db1AR_S464 = bARK_desens - bARK_resens; % ydot(5)
db1AR_S301 = PKA_desens - PKA_resens; % ydot(6)

G_act = k_G_act*(RG+LRG);
G_hyd = k_G_hyd*GsaGTPtot;   %   GTP->GDP // 14:55:26, Mon, 13-November-2017, By Haibo
G_reassoc = k_G_reassoc*GsaGDP*Gsby;
dGsaGTPtot = G_act - G_hyd; % ydot(7)
dGsaGDP = G_hyd - G_reassoc; % ydot(8)
dGsby = G_act - G_reassoc; % ydot(9)
% end b-AR module
%% cAMP module

ACtot           = 47e-3;          % (uM) total adenylyl cyclase % RABBIT
% ACtot           = 70.57e-3; % MOUSE

ATP             = 5e3;            % (uM) total ATP
k_AC_basal      = 0.2e-3;         % (1/ms) basal cAMP generation rate by AC
Km_AC_basal     = 1.03e3;         % (uM) basal AC affinity for ATP

Kd_AC_Gsa       = 0.4;            % (uM) Kd for AC association with Gsa
kf_AC_Gsa       = 1;              % (1/[uM ms]) forward rate for AC association with Gsa
kr_AC_Gsa       = Kd_AC_Gsa;      % (1/ms) reverse rate for AC association with Gsa

k_AC_Gsa        = 8.5e-3;         % (1/ms) basal cAMP generation rate by AC:Gsa
Km_AC_Gsa       = 315.0;          % (uM) AC:Gsa affinity for ATP

Kd_AC_FSK       = 44.0;           % (uM) Kd for FSK binding to AC
k_AC_FSK        = 7.3e-3;         % (1/ms) basal cAMP generation rate by AC:FSK
Km_AC_FSK       = 860.0;          % (uM) AC:FSK affinity for ATP

%PDEtot          = 2*36e-3;       % (uM) total phosphodiesterase (RABBIT - PDE3 and PDE4?)
PDEtot          = 22.85e-3;       % (uM) total phosphodiesterase (MOUSE only PDE4)
k_cAMP_PDE      = 5e-3;           % (1/ms) cAMP hydrolysis rate by PDE
k_cAMP_PDEp     = 2*k_cAMP_PDE;   % (1/ms) cAMP hydrolysis rate by phosphorylated PDE
Km_PDE_cAMP     = 1.3;            % (uM) PDE affinity for cAMP

Kd_PDE_IBMX     = 30.0;           % (uM) Kd_R2cAMP_C for IBMX binding to PDE
k_PKA_PDE       = 7.5e-3;         % (1/ms) rate constant for PDE phosphorylation by type 1 PKA
k_PP_PDE        = 1.5e-3;         % (1/ms) rate constant for PDE dephosphorylation by phosphatases

cAMP = cAMPtot - (RCcAMP_I+2*RCcAMPcAMP_I+2*RcAMPcAMP_I) - (RCcAMP_II+2*RCcAMPcAMP_II+2*RcAMPcAMP_II);
AC = ACtot-AC_GsaGTP;
GsaGTP = GsaGTPtot - AC_GsaGTP;
dAC_GsaGTP = kf_AC_Gsa*GsaGTP*AC - kr_AC_Gsa*AC_GsaGTP;

AC_FSK = FSK*AC/Kd_AC_FSK;
AC_ACT_BASAL = k_AC_basal*AC*ATP/(Km_AC_basal+ATP);
AC_ACT_GSA = k_AC_Gsa*AC_GsaGTP*ATP/(Km_AC_Gsa+ATP);
AC_ACT_FSK = k_AC_FSK*AC_FSK*ATP/(Km_AC_FSK+ATP);

%Mouse
% PDE_IBMX = PDEtot*IBMX/Kd_PDE_IBMX;
% PDE = PDEtot - PDE_IBMX - PDEp;
% dPDEp = k_PKA_PDE*PKACII*PDE - k_PP_PDE*PDEp;
% PDE_ACT = k_cAMP_PDE*PDE*cAMP/(Km_PDE_cAMP+cAMP) + k_cAMP_PDEp*PDEp*cAMP/(Km_PDE_cAMP+cAMP);

%% ADDED THIS-NATE
PDE3tot =  0.75 * 0.036;   % (uM) total phosphodiesterase
PDE4tot =  0.75 * 0.036;   % (uM) total phosphodiesterase
k_cAMP_PDE3 = 3.5e-3;               % k_pde3        [1/ms]
k_cAMP_PDE3p = 2 * k_cAMP_PDE3; % (1/ms) cAMP hydrolysis rate by phosphorylated PDE
Km_PDE3_cAMP = 0.15;           % Km_pde3       [uM]
k_cAMP_PDE4 = 5.0e-3;             % k_pde4        [1/ms]
k_cAMP_PDE4p = 2 * 5.0e-3; % (1/ms) cAMP hydrolysis rate by phosphorylated PDE
Km_PDE4_cAMP = 1.3;            % Km_pde4       [uM]

Kd_PDE_IBMX = 30.0;           % (uM) Kd_R2cAMP_C for IBMX binding to PDE
k_PKA_PDE = 7.5e-3;       % (1/ms) rate constant for PDE phosphorylation by type 1 PKA
k_PP_PDE = 1.5e-3;       % (1/ms) rate constant for PDE dephosphorylation by phosphatases

cAMP = cAMPtot - (RCcAMP_I + 2 * RCcAMPcAMP_I + 2 * RcAMPcAMP_I) - (RCcAMP_II + 2 * RCcAMPcAMP_II + 2 * RcAMPcAMP_II);
AC = ACtot - AC_GsaGTP;
GsaGTP = GsaGTPtot - AC_GsaGTP;
dAC_GsaGTP = kf_AC_Gsa * GsaGTP * AC - kr_AC_Gsa * AC_GsaGTP;

AC_FSK = FSK * AC / Kd_AC_FSK;
AC_ACT_BASAL = k_AC_basal * AC * ATP / (Km_AC_basal + ATP);
AC_ACT_GSA = k_AC_Gsa * AC_GsaGTP * ATP / (Km_AC_Gsa + ATP);
AC_ACT_FSK = k_AC_FSK * AC_FSK * ATP / (Km_AC_FSK + ATP);

% RABBIT                            %added to match stefano's and Haibos
PDE3_IBMX = PDE3tot * IBMX / (Kd_PDE_IBMX + IBMX);
PDE3 = PDE3tot - PDE3_IBMX - PDE3p;
dPDE3p = k_PKA_PDE * PKACII * PDE3 - k_PP_PDE * PDE3p; % ydot(10)
PDE3_ACT = k_cAMP_PDE3 * PDE3 * cAMP / (Km_PDE3_cAMP + cAMP) + k_cAMP_PDE3p * PDE3p * cAMP / (Km_PDE3_cAMP + cAMP);

% PDE4_IBMX = PDE4tot * IBMX / Kd_PDE_IBMX;
PDE4_IBMX = PDE4tot * IBMX / (Kd_PDE_IBMX + IBMX);
PDE4 = PDE4tot - PDE4_IBMX - PDE4p;
dPDE4p = k_PKA_PDE * PKACII * PDE4 - k_PP_PDE * PDE4p; % ydot() - NEW STATE VARIABLE NEEDED
PDE4_ACT = k_cAMP_PDE4 * PDE4 * cAMP / (Km_PDE4_cAMP + cAMP) + k_cAMP_PDE4p * PDE4p * cAMP / (Km_PDE4_cAMP + cAMP);

dcAMPtot = AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE3_ACT - PDE4_ACT; % ydot(11)


% dcAMPtot = AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE_ACT; % ydot(15)

% end cAMP module
%% PKA module

PKAIItot        = 0.084;          % (uM) total type 2 PKA % RABBIT
% PKAIItot        = 0.059; % MOUSE

PKItot          = 0.18;           % (uM) total PKI
kf_RC_cAMP      = 1;              % (1/[uM ms]) Kd for PKA RC binding to cAMP
kf_RCcAMP_cAMP  = 1;              % (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
kf_RcAMPcAMP_C  = 4.375;          % (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
kf_PKA_PKI      = 1;              % (1/[uM ms]) Ki for PKA inhibition by PKI
kr_RC_cAMP      = 1.64;           % (1/ms) Kd for PKA RC binding to cAMP
kr_RCcAMP_cAMP  = 9.14;           % (1/ms) Kd for PKA RC:cAMP binding to cAMP
kr_RcAMPcAMP_C  = 1;              % (1/ms) Kd for PKA R:cAMPcAMP binding to C
kr_PKA_PKI      = 2e-4;           % (1/ms) Ki for PKA inhibition by PKI
epsilon         = 10;             % (-) AKAP-mediated scaling factor

PKI = PKItot - PKACI_PKI - PKACII_PKI;

dRC_I = - kf_RC_cAMP*RC_I*cAMP + kr_RC_cAMP*RCcAMP_I;
dRCcAMP_I = - kr_RC_cAMP*RCcAMP_I + kf_RC_cAMP*RC_I*cAMP - kf_RCcAMP_cAMP*RCcAMP_I*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_I;
dRCcAMPcAMP_I = - kr_RCcAMP_cAMP*RCcAMPcAMP_I + kf_RCcAMP_cAMP*RCcAMP_I*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_I + kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI;
dRcAMPcAMP_I = - kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I;
dPKACI = - kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I - kf_PKA_PKI*PKACI*PKI + kr_PKA_PKI*PKACI_PKI; % ydot(17)
dPKACI_PKI = - kr_PKA_PKI*PKACI_PKI + kf_PKA_PKI*PKACI*PKI;

dRC_II = - kf_RC_cAMP*RC_II*cAMP + kr_RC_cAMP*RCcAMP_II;
dRCcAMP_II = - kr_RC_cAMP*RCcAMP_II + kf_RC_cAMP*RC_II*cAMP - kf_RCcAMP_cAMP*RCcAMP_II*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_II;
dRCcAMPcAMP_II = - kr_RCcAMP_cAMP*RCcAMPcAMP_II + kf_RCcAMP_cAMP*RCcAMP_II*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_II + kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII;
dRcAMPcAMP_II = - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II;
dPKACII = - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II - kf_PKA_PKI*PKACII*PKI + kr_PKA_PKI*PKACII_PKI; % ydot(18)
dPKACII_PKI = - kr_PKA_PKI*PKACII_PKI + kf_PKA_PKI*PKACII*PKI;
% end PKA module
%% I-1/PP1 module

%PP1tot = pin(8); % PP1tot = 0.89; % (uM) total phosphatase 1
I1tot           = 0.3;            % (uM) total inhibitor 1
k_PKA_I1        = 60e-3;          % (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
Km_PKA_I1       = 1.0;            % (uM) Km for I-1 phosphorylation by type 1 PKA
Vmax_PP2A_I1    = 14.0e-3;        % (uM/ms) Vmax for I-1 dephosphorylation by PP2A
Km_PP2A_I1      = 1.0;            % (uM) Km for I-1 dephosphorylation by PP2A

Ki_PP1_I1       = 1.0e-3;         % (uM) Ki for PP1 inhibition by I-1
kf_PP1_I1       = 1;              % (uM) Ki for PP1 inhibition by I-1
kr_PP1_I1       = Ki_PP1_I1;      % (uM) Ki for PP1 inhibition by I-1

I1 = I1tot - I1ptot;
PP1 = PP1tot - I1p_PP1;
I1p = I1ptot - I1p_PP1;
dI1p_PP1 = kf_PP1_I1*PP1*I1p - kr_PP1_I1*I1p_PP1;

I1_phosph = k_PKA_I1*PKACI*I1/(Km_PKA_I1+I1);
I1_dephosph = Vmax_PP2A_I1*I1ptot/(Km_PP2A_I1+I1ptot);
dI1p_PP1 = kf_PP1_I1 * PP1 * I1p - kr_PP1_I1 * I1p_PP1;
dI1ptot = I1_phosph - I1_dephosph; % ydot(20)
% end I-1/PP1 module
%% PLB module

%PLBtot = pin(4); %p(41) = PLBtot; % PLBtot    [uM]
k_PKA_PLB = 54e-3; %p(44) = 54;     % k_pka_plb     [1/ms]
Km_PKA_PLB = 21; %p(45) = 21;     % Km_pka_plb    [uM]
k_PP1_PLB = 8.5e-3; %p(46) = 8.5;    % k_pp1_plb     [1/ms]
Km_PP1_PLB = 7.0; %p(47) = 7.0;    % Km_pp1_plb    [uM]

PLB = PLBtot - PLBp;
PLB_phosph = k_PKA_PLB*PKACI*PLB/(Km_PKA_PLB+PLB);
PLB_dephosph = k_PP1_PLB*PP1*PLBp/(Km_PP1_PLB+PLBp);
dPLBp = PLB_phosph - PLB_dephosph; % ydot(19)
% end PLB module
%% PLM module (included 09/18/12) MOUSE

%PLMtot = pin(10); % p(102) = PLMtot; % PLMtot    [uM]
k_PKA_PLM = 54e-3; % p(103) = 54;     % k_pka_plb     [1/ms]
Km_PKA_PLM = 21; % p(104) = 21;     % Km_pka_plb    [uM]
k_PP1_PLM = 8.5e-3; % p(105) = 8.5;    % k_pp1_plb     [1/ms]
Km_PP1_PLM = 7.0; % p(106) = 7.0;    % Km_pp1_plb    [uM]

PLM = PLMtot - PLMp;
PLM_phosph = k_PKA_PLM*PKACI*PLM/(Km_PKA_PLM+PLM);
PLM_dephosph = k_PP1_PLM*PP1*PLMp/(Km_PP1_PLM+PLMp);
dPLMp = PLM_phosph - PLM_dephosph; % ydot(32)
% end PLM module
%% LTCC module

%LCCtot = pin(2); %p(53) = LCCtot; % LCCtot        [uM]
PKACII_LCCtot = 0.025; %p(54) = 0.025;  % PKAIIlcctot   [uM]
PP1_LCC = 0.025; %p(55) = 0.025;  % PP1lcctot     [uM]
PP2A_LCC = 0.025; %p(56) = 0.025;  % PP2Alcctot    [uM]
k_PKA_LCC = 54e-3; %p(57) = 54;     % k_pka_lcc     [1/ms]
Km_PKA_LCC = 21; %p(58) = 21;%*1.6;     % Km_pka_lcc    [uM]
k_PP1_LCC = 8.52e-3; %p(59) = 8.52;   % k_pp1_lcc     [1/ms] RABBIT, MOUSE
                %p(59) = 8.5;   % k_pp1_lcc     [1/sec] RAT
Km_PP1_LCC = 3; %p(60) = 3;      % Km_pp1_lcc    [uM]
k_PP2A_LCC = 10.1e-3; %p(61) = 10.1;   % k_pp2a_lcc    [1/ms]
Km_PP2A_LCC = 3; %p(62) = 3;      % Km_pp2a_lcc   [uM]

PKACII_LCC = (PKACII_LCCtot/PKAIItot)*PKACII;
LCCa = LCCtot - LCCap;
LCCa_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCa/(Km_PKA_LCC+epsilon*LCCa);
LCCa_dephosph = epsilon*k_PP2A_LCC*PP2A_LCC*LCCap/(Km_PP2A_LCC+epsilon*LCCap);
dLCCap = LCCa_phosph - LCCa_dephosph; % ydot(23)
LCCb = LCCtot - LCCbp;
LCCb_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCb/(Km_PKA_LCC+epsilon*LCCb);
LCCb_dephosph = epsilon*k_PP1_LCC*PP1_LCC*LCCbp/(Km_PP1_LCC+epsilon*LCCbp);
dLCCbp = LCCb_phosph - LCCb_dephosph; % ydot(24)
% end LCC module
%% RyR module (not included in Yang-Saucerman)

%RyRtot = pin(3); %p(63) = RyRtot; % RyRtot        [uM]
PKACII_ryrtot = 0.034;  %p(64) = 0.034;  % PKAIIryrtot   [uM]
PP1ryr = 0.034;  %p(65) = 0.034;  % PP1ryr        [uM]
PP2Aryr = 0.034;  %p(66) = 0.034;  % PP2Aryr       [uM]
kcat_pka_ryr = 54e-3;     %p(67) = 54;     % kcat_pka_ryr  [1/ms]
Km_pka_ryr = 21;     %p(68) = 21;     % Km_pka_ryr    [uM]
kcat_pp1_ryr = 8.52e-3;   %p(69) = 8.52;   % kcat_pp1_ryr  [1/ms]
Km_pp1_ryr = 7;      %p(70) = 7;      % Km_pp1_ryr    [uM]
kcat_pp2a_ryr = 10.1e-3;   %p(71) = 10.1;   % kcat_pp2a_ryr [1/ms]
Km_pp2a_ryr = 4.1;    %p(72) = 4.1;    % Km_pp2a_ryr   [uM]

PKACII_ryr = (PKACII_ryrtot/PKAIItot)*PKACII;
RyR = RyRtot-RyRp;
RyRPHOSPH = epsilon*kcat_pka_ryr*PKACII_ryr*RyR/(Km_pka_ryr+epsilon*RyR);
RyRDEPHOSPH1 = epsilon*kcat_pp1_ryr*PP1ryr*RyRp/(Km_pp1_ryr+epsilon*RyRp);
RyRDEPHOSPH2A = epsilon*kcat_pp2a_ryr*PP2Aryr*RyRp/(Km_pp2a_ryr+epsilon*RyRp);
dRyRp = RyRPHOSPH-RyRDEPHOSPH1-RyRDEPHOSPH2A; % ydot(25)
% end RyR module
%% TnI module

%TnItot = pin(5); %p(73) = TnItot; % TnItot        [uM]
PP2A_TnI = 0.67;   % PP2Atni       [uM]
k_PKA_TnI = 54e-3;     % kcat_pka_tni  [1/ms]
Km_PKA_TnI = 21;     % Km_pka_tni    [uM]
k_PP2A_TnI = 10.1e-3;   % kcat_pp2a_tni [1/ms]
Km_PP2A_TnI = 4.1;    % Km_pp2a_tni   [uM]

TnI = TnItot - TnIp;
TnI_phosph = k_PKA_TnI*PKACI*TnI/(Km_PKA_TnI+TnI);
TnI_dephosph = k_PP2A_TnI*PP2A_TnI*TnIp/(Km_PP2A_TnI+TnIp);
dTnIp = TnI_phosph - TnI_dephosph; % ydot(26)
% end TnI module
%% Iks module (now without equations for mutation)

%not used in atria (this is on line502)
% PP2A_myo = 0.67;                 % PP2Amyo       [uM]
% kcat_pka_myo = K_rate_Scale * 54e-3;     % kcat_pka_myo  [1/ms]
% Km_pka_myo = 21;                % Km_pka_myo    [uM]
% kcat_pp2a_myo = K_rate_Scale * 10.1e-3;  % kcat_pp2a_myo [1/ms]
% Km_pp2a_myo = 4.1;               % Km_pp2a_myo   [uM]
% 
% Myon = Myotot - Myop; % Non-phos = tot - phos
% MyoPHOSPH = kcat_pka_myo * PKACI * Myon / (Km_pka_myo + Myon);
% MyoDEPHOSPH = kcat_pp2a_myo * PP2A_myo * Myop / (Km_pp2a_myo + Myop);
% dMyop = MyoPHOSPH - MyoDEPHOSPH; % ydot

Yotiao_tot = 0.025; % Yotiao_tot    [uM]
K_yotiao = K_rate_Scale * 0.1e-3;  % K_yotiao      [uM] ** apply G589D mutation here: 1e2 **
PKAII_ikstot = PKAII_ikstot_Scale * 0.025;  %  PKAII_ikstot  [uM]
PP1_ikstot = PP1_ikstot_Scale * 0.025;   % PP1_ikstot    [uM]
k_pka_iks = K_rate_Scale * 1.87e-3;% 54; % k_pka_iks     [1/ms] //  adjusted as in Xie et al 2013
Km_pka_iks = Km_pka_iks_Scale * 21;  %Km_pka_iks    [uM]
k_pp1_iks = K_rate_Scale * 0.19e-3; % 8.52;  % k_pp1_iks     [1/ms] //  adjusted as in Xie et al 2013
Km_pp1_iks = Km_pp1_iks_Scale * 7;  % Km_pp1_iks    [uM]
% 
y1 = (-(K_yotiao+IKstot-Yotiao_tot)+sqrt((K_yotiao+IKstot-Yotiao_tot)^2+4*K_yotiao*Yotiao_tot))/2;
x1 = IKstot/(1+y1/K_yotiao);
y2 = (-(K_yotiao+IKstot-Yotiao_tot)-sqrt((K_yotiao+IKstot-Yotiao_tot)^2+4*K_yotiao*Yotiao_tot))/2;
x2 = IKstot/(1+y2/K_yotiao);

free_IKs = x1*(y1 > 0) + x2*(y1 <= 0);
free_Yotiao = y1*(y1 > 0) + y2*(y1 <= 0);
IksYot = free_IKs*free_Yotiao/K_yotiao; % [uM] % IKs-Yotiao

PKACiks = IksYot / IKstot * (PKAII_ikstot / PKAIItot) * PKACII;
PP1iks = IksYot / IKstot * PP1_ikstot;

KSn = IKstot - KSp; % Non-phos = tot - phos
IKS_PHOSPH = epsilon * k_pka_iks * PKACiks * KSn / (Km_pka_iks + epsilon * KSn);
IKS_DEPHOSPH = epsilon * k_pp1_iks * PP1iks * KSp / (Km_pp1_iks + epsilon * KSp);
dKSp = IKS_PHOSPH - IKS_DEPHOSPH; % ydot
% end Iks module
%% Ikr module (from Iks)

%IKrtot = pin(6); % p(79) [uM]
PKAIIikrtot = 0.025; % PKAII_ikrtot [uM]
PP1_ikrtot = 0.025; % PP1_ikrtot [uM]
k_pka_ikr = 1.87e-3; % 54e-3; % k_pka_ikr [1/ms] % adjusted as in Xie et al 2013
Km_pka_ikr = 21; % Km_pka_ikr [uM]
k_pp1_ikr = 0.19e-3; % 8.52e-3; % k_pp1_ikr [1/ms] % adjusted as in Xie et al 2013
Km_pp1_ikr = 7; % Km_pp1_ikr [uM]

KRn = IKrtot-KRp;
PKACII_ikr = (PKAIIikrtot/PKAIItot)*PKACII;
IKR_PHOSPH = epsilon*k_pka_ikr*PKACII_ikr*KRn/(Km_pka_ikr+epsilon*KRn);
IKR_DEPHOSPH = epsilon*k_pp1_ikr*PP1_ikrtot*KRp/(Km_pp1_ikr+epsilon*KRp);
dKRp = IKR_PHOSPH-IKR_DEPHOSPH;
% end Ikr module
%% Ikur module (from Ikr)

%IKurtot = pin(9);% p(95) = IKurtot;   % Ikur_tot      [uM]
PKAII_KURtot = 0.025; % p (96) = 0.025;  % PKAII_KURtot [uM]
PP1_KURtot = 0.025; % p(97) = 0.025;  % PP1_KURtot   [uM]
k_pka_KUR = 1.87e-3; % 54e-3; % k_pka_ikur [1/ms] % adjusted as in Xie et al 2013
Km_pka_KUR = 21; % p(99) = 21;    % Km_pka_KUR   [uM]
k_pp1_KUR = 0.19e-3; % 8.52e-3; % k_pp1_ikur [1/ms] % adjusted as in Xie et al 2013
Km_pp1_KUR = 7; % p(101) = 7;     % Km_pp1_KUR   [uM]

KURn = IKurtot - KURp;  % Non-phos = tot - phos
PKACII_KUR = (PKAII_KURtot/PKAIItot)*PKACII; % (PKA_KURtot/PKAIItot)*PKAIIact
KURphos = epsilon*PKACII_KUR*k_pka_KUR*KURn/(Km_pka_KUR+epsilon*KURn);
KURdephos = epsilon*PP1_KURtot*k_pp1_KUR*KURp/(Km_pp1_KUR+epsilon*KURp);
dKURp = KURphos - KURdephos;
% end Ikur module
%% CFTR module (included 04/30/10)

%ICFTRtot = pin(7); % p(88) = ICFTRtot;  % ICFTR_tot      [uM]
PKAII_CFTRtot = 0.025;  %p(89) = 0.025;  % PKAII_CFTRtot [uM]
PP1_CFTRtot = 0.025;  %p(90) = 0.025;  % PP1_CFTRtot   [uM]
k_pka_CFTR = 54e-3;     %p(91) = 54;     % k_pka_CFTR    [1/ms]
Km_pka_CFTR = 8.5;    %p(92) = 8.5;    % Km_pka_CFTR   [uM]
k_pp1_CFTR = 8.52e-3;   %p(93) = 8.52;   % k_pp1_CFTR    [1/ms]
Km_pp1_CFTR = 7;      %p(94) = 7;      % Km_pp1_CFTR   [uM]

CFTRn = ICFTRtot - CFTRp;  % Non-phos = tot - phos
PKACII_CFTR = (PKAII_CFTRtot/PKAIItot)*PKACII; % (PKACFTRtot/PKAIItot)*PKAIIact
CFTRphos = epsilon*PKACII_CFTR*k_pka_CFTR*CFTRn/(Km_pka_CFTR+epsilon*CFTRn);
CFTRdephos = epsilon*PP1_CFTRtot*k_pp1_CFTR*CFTRp/(Km_pp1_CFTR+epsilon*CFTRp);
dCFTRp = CFTRphos - CFTRdephos; %dCFTRp = 0;
% end CFTR module
%% ClCa module (from CFTR)

%IClCatot = pin(7); % p(88) = IClCatot;  % IClCa_tot      [uM]
PKAII_ClCatot = 0.025;  %p(89) = 0.025;  % PKAII_ClCatot [uM]
PP1_ClCatot = 0.025;  %p(90) = 0.025;  % PP1_ClCatot   [uM]
k_pka_ClCa = 54e-3;     %p(91) = 54;     % k_pka_ClCa    [1/ms]
Km_pka_ClCa = 8.5;    %p(92) = 8.5;    % Km_pka_ClCa   [uM]
k_pp1_ClCa = 8.52e-3;   %p(93) = 8.52;   % k_pp1_ClCa    [1/ms]
Km_pp1_ClCa = 7;      %p(94) = 7;      % Km_pp1_ClCa   [uM]

ClCan = IClCatot - ClCap;  % Non-phos = tot - phos
PKACII_ClCa = (PKAII_ClCatot/PKAIItot)*PKACII; % (PKAClCatot/PKAIItot)*PKAIIact
ClCaphos = epsilon*PKACII_ClCa*k_pka_ClCa*ClCan/(Km_pka_ClCa+epsilon*ClCan);
ClCadephos = epsilon*PP1_ClCatot*k_pp1_ClCa*ClCap/(Km_pp1_ClCa+epsilon*ClCap);
dClCap = ClCaphos - ClCadephos;
% end ClCa module
%% Myofilament module (from TnI)
%Haibos
%not used in atria (this is on line502)
% PP2A_myo = 0.67;                 % PP2Amyo       [uM]
% kcat_pka_myo = K_rate_Scale * 54e-3;     % kcat_pka_myo  [1/ms]
% Km_pka_myo = 21;                % Km_pka_myo    [uM]
% kcat_pp2a_myo = K_rate_Scale * 10.1e-3;  % kcat_pp2a_myo [1/ms]
% Km_pp2a_myo = 4.1;               % Km_pp2a_myo   [uM]
% 
% Myon = Myotot - Myop; % Non-phos = tot - phos
% MyoPHOSPH = kcat_pka_myo * PKACI * Myon / (Km_pka_myo + Myon);
% MyoDEPHOSPH = kcat_pp2a_myo * PP2A_myo * Myop / (Km_pp2a_myo + Myop);
% dMyop = MyoPHOSPH - MyoDEPHOSPH; % ydot


%myotot = pin(5); %p(73) = myotot; % myotot        [uM]
PP2A_myo = 0.67;   % PP2Amyo       [uM]
k_PKA_myo = 54e-3;     % kcat_pka_myo  [1/ms]
Km_PKA_myo = 21;     % Km_pka_myo    [uM]
k_PP2A_myo = 10.1e-3;   % kcat_pp2a_myo [1/ms]
Km_PP2A_myo = 4.1;    % Km_pp2a_myo   [uM]

myon = myotot - myop;
myo_phosph = k_PKA_myo*PKACI*myon/(Km_PKA_myo+myon);
myo_dephosph = k_PP2A_myo*PP2A_myo*myop/(Km_PP2A_myo+myop);
dmyop = myo_phosph - myo_dephosph; % ydot(26)
% end myofilament module
%% Ina module (from Ikr)

%INatot = pin(6); % p(79) [uM]
PKAIIinatot = 0.025; % PKAII_inatot [uM]
PP1_inatot = 0.025; % PP1_inatot [uM]
k_pka_ina = 1.87e-3; % 54e-3; % k_pka_ina [1/ms] % adjusted as in Xie et al 2013
Km_pka_ina = 21; % Km_pka_ina [uM]
k_pp1_ina = 0.19e-3; % 8.52e-3; % k_pp1_ina [1/ms] % adjusted as in Xie et al 2013
Km_pp1_ina = 7; % Km_pp1_ina [uM]

NAn = INatot-NAp;
PKACII_ina = (PKAIIinatot/PKAIItot)*PKACII;
INA_PHOSPH = epsilon*k_pka_ina*PKACII_ina*NAn/(Km_pka_ina+epsilon*NAn);
INA_DEPHOSPH = epsilon*k_pp1_ina*PP1_inatot*NAp/(Km_pp1_ina+epsilon*NAp);
dNAp = INA_PHOSPH-INA_DEPHOSPH;
% end Ina module
%% Ito module (from Ikr)

%Itotot = pin(6); % p(79) [uM]
PKAIIitotot = 0.025; % PKAII_ikrtot [uM]
PP1_itotot = 0.025; % PP1_ikrtot [uM]
k_pka_ito = 1.87e-3; % 54e-3; % k_pka_ikr [1/ms] % adjusted as in Xie et al 2013
Km_pka_ito = 21; % Km_pka_ikr [uM]
k_pp1_ito = 0.19e-3; % 8.52e-3; % k_pp1_ikr [1/ms] % adjusted as in Xie et al 2013
Km_pp1_ito = 7; % Km_pp1_ikr [uM]

TOn = Itotot-TOp;
PKACII_ito = (PKAIIitotot/PKAIItot)*PKACII;
Ito_PHOSPH = epsilon*k_pka_ito*PKACII_ito*TOn/(Km_pka_ito+epsilon*TOn);
Ito_DEPHOSPH = epsilon*k_pp1_ito*PP1_itotot*TOp/(Km_pp1_ito+epsilon*TOp);
dTOp = Ito_PHOSPH-Ito_DEPHOSPH;
% end Ito module
%% Ik1 module (from Ikr)

%IK1tot = pin(6); % p(79) [uM]
PKAIIik1tot = 0.025; % PKAII_ikrtot [uM]
PP1_ik1tot = 0.025; % PP1_ikrtot [uM]
k_pka_ik1 = 1.87e-3; % 54e-3; % k_pka_ikr [1/ms] % adjusted as in Xie et al 2013
Km_pka_ik1 = 21; % Km_pka_ikr [uM]
k_pp1_ik1 = 0.19e-3; % 8.52e-3; % k_pp1_ikr [1/ms] % adjusted as in Xie et al 2013
Km_pp1_ik1 = 7; % Km_pp1_ikr [uM]

K1n = IK1tot-K1p;
PKACII_ik1 = (PKAIIik1tot/PKAIItot)*PKACII;
IK1_PHOSPH = epsilon*k_pka_ik1*PKACII_ik1*K1n/(Km_pka_ik1+epsilon*K1n);
IK1_DEPHOSPH = epsilon*k_pp1_ik1*PP1_ik1tot*K1p/(Km_pp1_ik1+epsilon*K1p);
dK1p = IK1_PHOSPH-IK1_DEPHOSPH;
% end Ik1 module
%% ydot

ydot(1)=dLR;
ydot(2)=dLRG;
ydot(3)=dRG;
ydot(4)=db1AR_S464;
ydot(5)=db1AR_S301;
ydot(6)=dGsaGTPtot;
ydot(7)=dGsaGDP;
ydot(8)=dGsby;
ydot(9)=dAC_GsaGTP;
ydot(10)=dPDE4p;
ydot(11)=dcAMPtot;
ydot(12)=dRC_I;
ydot(13)=dRCcAMP_I;
ydot(14)=dRCcAMPcAMP_I;
ydot(15)=dRcAMPcAMP_I;
ydot(16)=dPKACI;
ydot(17)=dPKACI_PKI;
ydot(18)=dRC_II;
ydot(19)=dRCcAMP_II;
ydot(20)=dRCcAMPcAMP_II;
ydot(21)=dRcAMPcAMP_II;
ydot(22)=dPKACII;
ydot(23)=dPKACII_PKI;
ydot(24)=dI1p_PP1; % output CaMKII
ydot(25)=dI1ptot;
ydot(26)=dPLBp; % output
ydot(27)=dPLMp; % output
ydot(28)=dLCCap; % output
ydot(29)=dLCCbp; % output
ydot(30)=dRyRp; % output
ydot(31)=dTnIp; % output
ydot(32)=dKSp; % output
ydot(33)=dKRp; % output
ydot(34)=dKURp; % output
ydot(35)=dCFTRp; % output
ydot(36)=dClCap; % output
ydot(37)=dmyop; % output
ydot(38)=dNAp; % output
ydot(39)=dTOp; % output
ydot(40)=dK1p; % output
ydot(41)=dPDE3p; % output