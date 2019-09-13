#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <cmath>
#include <map>
#include <string>
#include <algorithm>
#include <vector>
extern std::map<std::string,float> initializeParamMap(std::map<std::string,float> Plist,int,char**, bool);
extern std::map<std::string,float> initializeICSMap(std::map<std::string,float> Vlist,int ,char** , bool);
extern char names[25][25];
typedef std::vector <double> state_type;
class NutSig{

//Parameters
float w_rtg,
w_sak,
w_pka_camp,
w_mig_snf,
k_acc_glu,
w_sch9_torc,
k_mRNA_degr,
w_torc_egoin,
k_acc_pro,
sigma_gcn2,
w_sch9,
w_eif_gcn2,
w_gln1_gln3,
w_snf_sak,
w_tps_pka,
w_pka,
gamma_gcn2,
sigma_gis1,
w_dot,
w_tre,
w_torc,
w_pka_sch9,
gammatre,
w_gis,
sigma_gln1,
sigma_tps,
sigma_dot,
sigma_eif,
w_gln_snf,
w_torc_glut,
w_ego,
sigma_ras,
sigma_tor,
w_mig,
gammacyr,
gammagap,
sigma_rtg,
gammasnf,
Cyr1_T,
k_camp_deg,
sigma_mig1,
gammasch9,
w_snf_glc,
w_gcn4_gcn2_trna,
Gln1_T,
sigma_ego,
w_mig_pka,
Gcn2_T,
gammatps,
k_pr,
tRNA_total,
w_tps,
w_tre_pka,
gammator,
w_gcn4,
gammapka,
w_rtg_torc,
sigma_gln,
Sak_T,
w_ego_basal,
w_gln3,
k_acc_nh4,
PKA_T,
k_camp_cyr,
gammagln3,
w_torc_snf,
gammaeif,
Trehalase_T,
ATP,
gammapde,
k_degr,
Rtg13_T,
Ras_T,
tRNA_sensitivity,
w_gln_sit,
w_cyr_glu,
w_pde_pka,
Carbon,
Gcn4_T,
k_camp_pde,
gammaras,
sigma_snf,
w_snf,
sigma_gap,
w_gap_torc,
Gln3_T,
eIF_T,
w_cyr_snf,
w_dot_sch_pka,
k_transcription,
Dot6_T,
Glutamine_ext,
PDE_T,
w_sak_pka,
gammaego,
sigma_cyr,
w_gcn_torc,
gamma_mig,
Mig1_T,
sigma_sch9,
w_ras_glu,
Tps1_T,
Proline,
TORC1_T,
w_ras,
w_cyr,
w_ras_pka,
sigma_trehalase,
sigma_gcn4,
w_gis_pka,
w_gcn,
w_gap_N,
Gis1_T,
sigma_pde,
w_ego_gap,
NH4,
gammagln1,
w_gln1,
Snf1_T,
w_pde,
EGO_T,
w_eif,
sigma_pka,
sigma_sak,
Sch9_T,
EGOGAP_T,
w_torc_ego,
w_gis_sch;

public:
NutSig(std::map<std::string, float> Plist){
w_rtg = Plist["w_rtg"];
w_sak = Plist["w_sak"];
w_pka_camp = Plist["w_pka_camp"];
w_mig_snf = Plist["w_mig_snf"];
k_acc_glu = Plist["k_acc_glu"];
w_sch9_torc = Plist["w_sch9_torc"];
k_mRNA_degr = Plist["k_mRNA_degr"];
w_torc_egoin = Plist["w_torc_egoin"];
k_acc_pro = Plist["k_acc_pro"];
sigma_gcn2 = Plist["sigma_gcn2"];
w_sch9 = Plist["w_sch9"];
w_eif_gcn2 = Plist["w_eif_gcn2"];
w_gln1_gln3 = Plist["w_gln1_gln3"];
w_snf_sak = Plist["w_snf_sak"];
w_tps_pka = Plist["w_tps_pka"];
w_pka = Plist["w_pka"];
gamma_gcn2 = Plist["gamma_gcn2"];
sigma_gis1 = Plist["sigma_gis1"];
w_dot = Plist["w_dot"];
w_tre = Plist["w_tre"];
w_torc = Plist["w_torc"];
w_pka_sch9 = Plist["w_pka_sch9"];
gammatre = Plist["gammatre"];
w_gis = Plist["w_gis"];
sigma_gln1 = Plist["sigma_gln1"];
sigma_tps = Plist["sigma_tps"];
sigma_dot = Plist["sigma_dot"];
sigma_eif = Plist["sigma_eif"];
w_gln_snf = Plist["w_gln_snf"];
w_torc_glut = Plist["w_torc_glut"];
w_ego = Plist["w_ego"];
sigma_ras = Plist["sigma_ras"];
sigma_tor = Plist["sigma_tor"];
w_mig = Plist["w_mig"];
gammacyr = Plist["gammacyr"];
gammagap = Plist["gammagap"];
sigma_rtg = Plist["sigma_rtg"];
gammasnf = Plist["gammasnf"];
Cyr1_T = Plist["Cyr1_T"];
k_camp_deg = Plist["k_camp_deg"];
sigma_mig1 = Plist["sigma_mig1"];
gammasch9 = Plist["gammasch9"];
w_snf_glc = Plist["w_snf_glc"];
w_gcn4_gcn2_trna = Plist["w_gcn4_gcn2_trna"];
Gln1_T = Plist["Gln1_T"];
sigma_ego = Plist["sigma_ego"];
w_mig_pka = Plist["w_mig_pka"];
Gcn2_T = Plist["Gcn2_T"];
gammatps = Plist["gammatps"];
k_pr = Plist["k_pr"];
tRNA_total = Plist["tRNA_total"];
w_tps = Plist["w_tps"];
w_tre_pka = Plist["w_tre_pka"];
gammator = Plist["gammator"];
w_gcn4 = Plist["w_gcn4"];
gammapka = Plist["gammapka"];
w_rtg_torc = Plist["w_rtg_torc"];
sigma_gln = Plist["sigma_gln"];
Sak_T = Plist["Sak_T"];
w_ego_basal = Plist["w_ego_basal"];
w_gln3 = Plist["w_gln3"];
k_acc_nh4 = Plist["k_acc_nh4"];
PKA_T = Plist["PKA_T"];
k_camp_cyr = Plist["k_camp_cyr"];
gammagln3 = Plist["gammagln3"];
w_torc_snf = Plist["w_torc_snf"];
gammaeif = Plist["gammaeif"];
Trehalase_T = Plist["Trehalase_T"];
ATP = Plist["ATP"];
gammapde = Plist["gammapde"];
k_degr = Plist["k_degr"];
Rtg13_T = Plist["Rtg13_T"];
Ras_T = Plist["Ras_T"];
tRNA_sensitivity = Plist["tRNA_sensitivity"];
w_gln_sit = Plist["w_gln_sit"];
w_cyr_glu = Plist["w_cyr_glu"];
w_pde_pka = Plist["w_pde_pka"];
Carbon = Plist["Carbon"];
Gcn4_T = Plist["Gcn4_T"];
k_camp_pde = Plist["k_camp_pde"];
gammaras = Plist["gammaras"];
sigma_snf = Plist["sigma_snf"];
w_snf = Plist["w_snf"];
sigma_gap = Plist["sigma_gap"];
w_gap_torc = Plist["w_gap_torc"];
Gln3_T = Plist["Gln3_T"];
eIF_T = Plist["eIF_T"];
w_cyr_snf = Plist["w_cyr_snf"];
w_dot_sch_pka = Plist["w_dot_sch_pka"];
k_transcription = Plist["k_transcription"];
Dot6_T = Plist["Dot6_T"];
Glutamine_ext = Plist["Glutamine_ext"];
PDE_T = Plist["PDE_T"];
w_sak_pka = Plist["w_sak_pka"];
gammaego = Plist["gammaego"];
sigma_cyr = Plist["sigma_cyr"];
w_gcn_torc = Plist["w_gcn_torc"];
gamma_mig = Plist["gamma_mig"];
Mig1_T = Plist["Mig1_T"];
sigma_sch9 = Plist["sigma_sch9"];
w_ras_glu = Plist["w_ras_glu"];
Tps1_T = Plist["Tps1_T"];
Proline = Plist["Proline"];
TORC1_T = Plist["TORC1_T"];
w_ras = Plist["w_ras"];
w_cyr = Plist["w_cyr"];
w_ras_pka = Plist["w_ras_pka"];
sigma_trehalase = Plist["sigma_trehalase"];
sigma_gcn4 = Plist["sigma_gcn4"];
w_gis_pka = Plist["w_gis_pka"];
w_gcn = Plist["w_gcn"];
w_gap_N = Plist["w_gap_N"];
Gis1_T = Plist["Gis1_T"];
sigma_pde = Plist["sigma_pde"];
w_ego_gap = Plist["w_ego_gap"];
NH4 = Plist["NH4"];
gammagln1 = Plist["gammagln1"];
w_gln1 = Plist["w_gln1"];
Snf1_T = Plist["Snf1_T"];
w_pde = Plist["w_pde"];
EGO_T = Plist["EGO_T"];
w_eif = Plist["w_eif"];
sigma_pka = Plist["sigma_pka"];
sigma_sak = Plist["sigma_sak"];
Sch9_T = Plist["Sch9_T"];
EGOGAP_T = Plist["EGOGAP_T"];
w_torc_ego = Plist["w_torc_ego"];
w_gis_sch = Plist["w_gis_sch"];
}
float shs(float sigma, float omega){
return 1/(1+exp(-sigma*omega));}

float tRNA(float tRNA_tot,float AmAc){
return std::min(tRNA_tot, AmAc);}

float pRib(float rib_comp, float init_factor){
return std::min(rib_comp,init_factor);}

void operator() (const state_type &x, state_type &xdot, const double ){

//Variables
float EGOGAP = x[0],
cAMP = x[1],
Sch9 = x[2],
Gcn2 = x[3],
Gcn4 = x[4],
Glutamine = x[5],
Snf1 = x[6],
PKA = x[7],
eIF = x[8],
Gln3 = x[9],
Gln1 = x[10],
Rtg13 = x[11],
Gis1 = x[12],
Mig1 = x[13],
Cyr1 = x[14],
PDE = x[15],
Sak = x[16],
Dot6 = x[17],
TORC1 = x[18],
Tps1 = x[19],
Trehalase = x[20],
Ras = x[21],
EGO = x[22],
Protein = x[23],
Rib = x[24];

//Equations
xdot[0] = gammagap*(EGOGAP_T*shs(sigma_gap,w_gap_N*(1-Glutamine)-w_gap_torc*TORC1)-EGOGAP);
xdot[1] = k_camp_cyr*Cyr1*ATP-k_camp_pde*PDE*cAMP-k_camp_deg*cAMP;
xdot[2] = gammasch9*(Sch9_T*shs(sigma_sch9,w_sch9_torc*TORC1-w_sch9)-Sch9);
xdot[3] = gamma_gcn2*(Gcn2_T*shs(sigma_gcn2,w_gcn-w_gcn_torc*Sch9)-Gcn2) ;
xdot[4] = Gcn4_T*shs(sigma_gcn4,w_gcn4_gcn2_trna*std::min(Gcn2,tRNA_sensitivity*(tRNA_total-tRNA(tRNA_total,Glutamine)))-w_gcn4)-Gcn4;
xdot[5] = (k_acc_glu*Glutamine_ext + k_acc_pro*Proline + k_acc_nh4*NH4*Gln1*Carbon) - k_degr*Glutamine ;
xdot[6] = gammasnf*(Snf1_T*shs(sigma_snf, -w_snf_glc*Carbon + w_snf_sak*Sak - w_snf)-Snf1) ;
xdot[7] = gammapka*(PKA_T*shs(sigma_pka,w_pka_camp*cAMP-w_pka-w_pka_sch9*Sch9)-PKA)        ;
xdot[8] = gammaeif*(eIF_T*shs(sigma_eif,w_eif-w_eif_gcn2*Gcn2)-eIF);
xdot[9] = gammagln3*(Gln3_T*shs(sigma_gln,-w_gln3+w_gln_snf*Snf1+w_gln_sit*(1-TORC1))-Gln3);
xdot[10] = gammagln1*(Gln1_T*shs(sigma_gln1,w_gln1_gln3*Gln3-w_gln1)- Gln1);
xdot[11] = Rtg13_T*shs(sigma_rtg,-w_rtg_torc*TORC1+w_rtg)-Rtg13;
xdot[12] = Gis1_T*shs(sigma_gis1,-w_gis_pka*PKA-w_gis_sch*Sch9+w_gis)-Gis1;
xdot[13] = gamma_mig*(Mig1_T*shs(sigma_mig1, w_mig_pka*PKA - w_mig_snf*Snf1 + w_mig) - Mig1);
xdot[14] = gammacyr*(Cyr1_T*shs(sigma_cyr,w_cyr_glu*Carbon*Ras-w_cyr-w_cyr_snf*Snf1)-Cyr1);
xdot[15] = gammapde*(PDE_T*shs(sigma_pde,w_pde_pka*PKA-w_pde)-PDE);
xdot[16] = Sak_T*shs(sigma_sak,w_sak-w_sak_pka*PKA)-Sak;
xdot[17] = Dot6_T*shs(sigma_dot,-w_dot_sch_pka*Sch9*PKA+w_dot)-Dot6;
xdot[18] = gammator*(TORC1_T*shs(sigma_tor,w_torc_glut*Glutamine + w_torc_ego*EGO-w_torc_egoin*(1-EGO)-w_torc - w_torc_snf*Snf1)-TORC1);
xdot[19] = gammatps*(Tps1_T*shs(sigma_tps,w_tps_pka*(PKA_T-PKA)-w_tps)-Tps1);
xdot[20] = gammatre*(Trehalase_T*shs(sigma_trehalase,w_tre_pka*PKA-w_tre)-Trehalase);
xdot[21] = gammaras*(Ras_T*shs(sigma_ras,-w_ras_pka*PKA+w_ras_glu*Carbon+w_ras)-Ras);
xdot[22] = gammaego*(EGO_T*shs(sigma_ego, w_ego_gap*EGOGAP*(Glutamine_ext+0.5*NH4+0.01*Proline) - w_ego*(1-Glutamine) - w_ego_basal )-EGO);
xdot[23] = k_pr*ATP*std::min(pRib(Rib,eIF),tRNA(tRNA_total,Glutamine))*Protein;
xdot[24] = k_transcription*(1-Dot6) - k_mRNA_degr*Rib ;
}};
#endif