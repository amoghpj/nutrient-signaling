#include "model.h"
char names[25][25] = {"EGOGAP","cAMP","Sch9","Gcn2","Gcn4","Glutamine","Snf1","PKA","eIF","Gln3","Gln1","Rtg13","Gis1","Mig1","Cyr1","PDE","Sak","Dot6","TORC1","Tps1","Trehalase","Ras","EGO","Protein","Rib"};

std::map<std::string,float> initializeParamMap(std::map<std::string,float> Plist, int argc,char** argv, bool verb){
bool paramflag = false;
bool icsflag = false;
char ListOfPars[200][25] = {
"w_rtg",
"w_sak",
"w_pka_camp",
"w_mig_snf",
"k_acc_glu",
"w_sch9_torc",
"k_mRNA_degr",
"w_torc_egoin",
"k_acc_pro",
"sigma_gcn2",
"w_sch9",
"w_eif_gcn2",
"w_gln1_gln3",
"w_snf_sak",
"w_tps_pka",
"w_pka",
"gamma_gcn2",
"sigma_gis1",
"w_dot",
"w_tre",
"w_torc",
"w_pka_sch9",
"gammatre",
"w_gis",
"sigma_gln1",
"sigma_tps",
"sigma_dot",
"sigma_eif",
"w_gln_snf",
"w_torc_glut",
"w_ego",
"sigma_ras",
"sigma_tor",
"w_mig",
"gammacyr",
"gammagap",
"sigma_rtg",
"gammasnf",
"Cyr1_T",
"k_camp_deg",
"sigma_mig1",
"gammasch9",
"w_snf_glc",
"w_gcn4_gcn2_trna",
"Gln1_T",
"sigma_ego",
"w_mig_pka",
"Gcn2_T",
"gammatps",
"k_pr",
"tRNA_total",
"w_tps",
"w_tre_pka",
"gammator",
"w_gcn4",
"gammapka",
"w_rtg_torc",
"sigma_gln",
"Sak_T",
"w_ego_basal",
"w_gln3",
"k_acc_nh4",
"PKA_T",
"k_camp_cyr",
"gammagln3",
"w_torc_snf",
"gammaeif",
"Trehalase_T",
"ATP",
"gammapde",
"k_degr",
"Rtg13_T",
"Ras_T",
"tRNA_sensitivity",
"w_gln_sit",
"w_cyr_glu",
"w_pde_pka",
"Carbon",
"Gcn4_T",
"k_camp_pde",
"gammaras",
"sigma_snf",
"w_snf",
"sigma_gap",
"w_gap_torc",
"Gln3_T",
"eIF_T",
"w_cyr_snf",
"w_dot_sch_pka",
"k_transcription",
"Dot6_T",
"Glutamine_ext",
"PDE_T",
"w_sak_pka",
"gammaego",
"sigma_cyr",
"w_gcn_torc",
"gamma_mig",
"Mig1_T",
"sigma_sch9",
"w_ras_glu",
"Tps1_T",
"Proline",
"TORC1_T",
"w_ras",
"w_cyr",
"w_ras_pka",
"sigma_trehalase",
"sigma_gcn4",
"w_gis_pka",
"w_gcn",
"w_gap_N",
"Gis1_T",
"sigma_pde",
"w_ego_gap",
"NH4",
"gammagln1",
"w_gln1",
"Snf1_T",
"w_pde",
"EGO_T",
"w_eif",
"sigma_pka",
"sigma_sak",
"Sch9_T",
"EGOGAP_T",
"w_torc_ego",
"w_gis_sch",
};

Plist["w_rtg"] = 0.186412268402;
Plist["w_sak"] = 0.205140223044;
Plist["w_pka_camp"] = 102.109826697;
Plist["w_mig_snf"] = 1.21479427525;
Plist["k_acc_glu"] = 0.0492311292391;
Plist["w_sch9_torc"] = 1.96219039127;
Plist["k_mRNA_degr"] = 0.0730004044196;
Plist["w_torc_egoin"] = 0.30260360841;
Plist["k_acc_pro"] = 0.000214641490256;
Plist["sigma_gcn2"] = 20.0;
Plist["w_sch9"] = 0.565348248672;
Plist["w_eif_gcn2"] = 0.275514108409;
Plist["w_gln1_gln3"] = 0.519802074308;
Plist["w_snf_sak"] = 1.51831949656;
Plist["w_tps_pka"] = 0.573915723095;
Plist["w_pka"] = 0.0580860701396;
Plist["gamma_gcn2"] = 4.71028660913;
Plist["sigma_gis1"] = 10.0;
Plist["w_dot"] = 0.292969276177;
Plist["w_tre"] = 1.06748143104;
Plist["w_torc"] = 0.538903351393;
Plist["w_pka_sch9"] = 17.4971243642;
Plist["gammatre"] = 0.341746174318;
Plist["w_gis"] = 1.30388948487;
Plist["sigma_gln1"] = 1.0;
Plist["sigma_tps"] = 5.0;
Plist["sigma_dot"] = 20.0;
Plist["sigma_eif"] = 1.0;
Plist["w_gln_snf"] = 3.90445735339;
Plist["w_torc_glut"] = 0.863231820077;
Plist["w_ego"] = 0.284475372776;
Plist["sigma_ras"] = 1.0;
Plist["sigma_tor"] = 5.0;
Plist["w_mig"] = 10.6446018992;
Plist["gammacyr"] = 8.95550770947;
Plist["gammagap"] = 0.562145750461;
Plist["sigma_rtg"] = 10.0;
Plist["gammasnf"] = 0.819616489712;
Plist["Cyr1_T"] = 1.0;
Plist["k_camp_deg"] = 0.0838063704863;
Plist["sigma_mig1"] = 0.27;
Plist["gammasch9"] = 4.63314532353;
Plist["w_snf_glc"] = 1.15001639964;
Plist["w_gcn4_gcn2_trna"] = 1.53157062813;
Plist["Gln1_T"] = 1.0;
Plist["sigma_ego"] = 5.0;
Plist["w_mig_pka"] = 2.30555200295;
Plist["Gcn2_T"] = 1.0;
Plist["gammatps"] = 0.471347096514;
Plist["k_pr"] = 0.0201693698466;
Plist["tRNA_total"] = 2.47451146942;
Plist["w_tps"] = 0.0529949376897;
Plist["w_tre_pka"] = 3.07232519451;
Plist["gammator"] = 7.54843313097;
Plist["w_gcn4"] = 0.742819875259;
Plist["gammapka"] = 2.67588273674;
Plist["w_rtg_torc"] = 0.876816147213;
Plist["sigma_gln"] = 10.0;
Plist["Sak_T"] = 1.0;
Plist["w_ego_basal"] = 0.0109926702877;
Plist["w_gln3"] = 0.639194021909;
Plist["k_acc_nh4"] = 0.00146932708863;
Plist["PKA_T"] = 1.0;
Plist["k_camp_cyr"] = 10.8655884228;
Plist["gammagln3"] = 0.0809218108437;
Plist["w_torc_snf"] = 0.437340183291;
Plist["gammaeif"] = 0.471028772049;
Plist["Trehalase_T"] = 1.0;
Plist["ATP"] = 0.5;
Plist["gammapde"] = 0.281811406645;
Plist["k_degr"] = 0.0898253686393;
Plist["Rtg13_T"] = 1.0;
Plist["Ras_T"] = 1.0;
Plist["tRNA_sensitivity"] = 74.507463086;
Plist["w_gln_sit"] = 0.861343608332;
Plist["w_cyr_glu"] = 5.12925553067;
Plist["w_pde_pka"] = 2.89230085203;
Plist["Carbon"] = 0.5;
Plist["Gcn4_T"] = 1.0;
Plist["k_camp_pde"] = 14.1162421901;
Plist["gammaras"] = 1.82228382029;
Plist["sigma_snf"] = 3.0;
Plist["w_snf"] = 0.537671836477;
Plist["sigma_gap"] = 1.0;
Plist["w_gap_torc"] = 88.3263828822;
Plist["Gln3_T"] = 1.0;
Plist["eIF_T"] = 1.0;
Plist["w_cyr_snf"] = 0.119130742303;
Plist["w_dot_sch_pka"] = 0.162969752116;
Plist["k_transcription"] = 0.235950807701;
Plist["Dot6_T"] = 1.0;
Plist["Glutamine_ext"] = 0.0;
Plist["PDE_T"] = 1.0;
Plist["w_sak_pka"] = 0.375154791579;
Plist["gammaego"] = 50.6551762575;
Plist["sigma_cyr"] = 3.5;
Plist["w_gcn_torc"] = 1.29443701855;
Plist["gamma_mig"] = 0.655857797953;
Plist["Mig1_T"] = 1.0;
Plist["sigma_sch9"] = 8.0;
Plist["w_ras_glu"] = 0.206713370722;
Plist["Tps1_T"] = 1.0;
Plist["Proline"] = 0.0;
Plist["TORC1_T"] = 1.0;
Plist["w_ras"] = 0.0207813759533;
Plist["w_cyr"] = 1.35141380018;
Plist["w_ras_pka"] = 1.87230225794;
Plist["sigma_trehalase"] = 10.0;
Plist["sigma_gcn4"] = 5.0;
Plist["w_gis_pka"] = 3.30131718395;
Plist["w_gcn"] = 0.11547502832;
Plist["w_gap_N"] = 7.75731921965;
Plist["Gis1_T"] = 1.0;
Plist["sigma_pde"] = 1.9;
Plist["w_ego_gap"] = 2.20798971462;
Plist["NH4"] = 0.0;
Plist["gammagln1"] = 0.0635115135427;
Plist["w_gln1"] = 0.222322780374;
Plist["Snf1_T"] = 1.0;
Plist["w_pde"] = 0.382693037748;
Plist["EGO_T"] = 1.0;
Plist["w_eif"] = 3.72843982973;
Plist["sigma_pka"] = 1.0;
Plist["sigma_sak"] = 20.0;
Plist["Sch9_T"] = 1.0;
Plist["EGOGAP_T"] = 1.0;
Plist["w_torc_ego"] = 0.877473800695;
Plist["w_gis_sch"] = 0.842195608342;
 if (verb == true){
   std::cout<<"\nIn initializeParamMap()";   std::cout<<"\nWe have "<<argc<<" arguments\n";}
 
 for (int i=1;i<argc;i++){
    std::string arg = argv[i];
    if (arg == "--pars"){
      paramflag = true;
      icsflag = false;}
    if (i+1 < argc){
      if ((paramflag==true) && (icsflag==false)){
        arg = argv[i];
        for (int j = 0; j<200;j++){
          if (arg == ListOfPars[j]){
            if (verb ==true){
            std::cout<<arg<<" found in ListofPars!\n";
            std::cout<<"It will be assigned the value="<<atof(argv[i+1])<<"\n";}
            Plist[arg] = atof(argv[i+1]);}
        }}}}
return Plist;}

std::map<std::string,float> initializeICSMap(std::map<std::string,float> Vlist, int argc,char** argv, bool verb){
  bool paramflag =false;
  bool icsflag =false;
  char ListOfVars[200][25] = {
"EGOGAP",
"cAMP",
"Sch9",
"Gcn2",
"Gcn4",
"Glutamine",
"Snf1",
"PKA",
"eIF",
"Gln3",
"Gln1",
"Rtg13",
"Gis1",
"Mig1",
"Cyr1",
"PDE",
"Sak",
"Dot6",
"TORC1",
"Tps1",
"Trehalase",
"Ras",
"EGO",
"Protein",
"Rib",
};

Vlist["EGOGAP"] = 0.0;
Vlist["cAMP"] = 0.0;
Vlist["Sch9"] = 0.0;
Vlist["Gcn2"] = 1.0;
Vlist["Gcn4"] = 1.0;
Vlist["Glutamine"] = 0.001;
Vlist["Snf1"] = 1.0;
Vlist["PKA"] = 0.0;
Vlist["eIF"] = 0.0;
Vlist["Gln3"] = 1.0;
Vlist["Gln1"] = 0.0;
Vlist["Rtg13"] = 1.0;
Vlist["Gis1"] = 1.0;
Vlist["Mig1"] = 1.0;
Vlist["Cyr1"] = 0.0;
Vlist["PDE"] = 0.0;
Vlist["Sak"] = 1.0;
Vlist["Dot6"] = 1.0;
Vlist["TORC1"] = 0.0;
Vlist["Tps1"] = 1.0;
Vlist["Trehalase"] = 0.0;
Vlist["Ras"] = 0.0;
Vlist["EGO"] = 0.0;
Vlist["Protein"] = 0.1;
Vlist["Rib"] = 0.0;
 if (verb ==true){
   std::cout<<"In initializeICSMap()";   std::cout<<"\nWe have "<<argc<<" arguments";}

 
 for (int i=1;i<argc;i++){
    std::string arg = argv[i];
    if (arg == "--ics"){
      paramflag = false;
      icsflag = true;}
    if (i+1 < argc){
      if ((paramflag==false) && (icsflag==true)){
        arg = argv[i];

        for (int j = 0; j<200;j++){
          if (arg == ListOfVars[j]){

            if (verb == true){
            std::cout<<arg<<" found in ListofVars!\n";
            std::cout<<"It will be assigned the value="<<atof(argv[i+1])<<"\n";}

            Vlist[arg] = atof(argv[i+1]);

          }
          
        }}}}
    return Vlist;
}
