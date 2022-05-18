#include "../FitterMichal.h"
#include "../ReadFile.h"


void traditional_fit_data_michal_2dfit_step1_hel(int choose = 1){


  FitterMichal ang;
  ReadFile cm;
  ReadFile cm_mc_rec;
  cm.ReadFileTreeData("../../../david_data/download.root",0);

  cm_mc_rec.ReadFileTreeData("../../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1);
//  cm_mc_gen.ReadFileTreeGen("../../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root");

  //cout << "mass before setting" << cm.DiMuMass.size() << endl;


  vector <float> DiMuMass;
  vector <float> costheta;
  vector <float> phi;
  vector <float> DiMuMassRec;
  vector <float> costheta_rec;
  vector <float> phi_rec;
  for(int i=0 ; i< cm.daughter1.size();i++){
    AngleCalc * calc = new AngleCalc();
    calc->CosThetaHelicityFrame(cm.daughter_pos[i],cm.daughter_neg[i]);
    calc->PhiHelicityFrame(cm.daughter_pos[i],cm.daughter_neg[i]);
    DiMuMass.push_back(cm.parent[i].M());
    costheta.push_back(calc->ct_he);
    phi.push_back(calc->phi_he);
  }

  for(int i=0 ; i< cm_mc_rec.daughter1.size();i++){
    AngleCalc * calc = new AngleCalc();
    calc->CosThetaHelicityFrame(cm_mc_rec.daughter_pos[i],cm_mc_rec.daughter_neg[i]);
    calc->PhiHelicityFrame(cm_mc_rec.daughter_pos[i],cm_mc_rec.daughter_neg[i]);


    DiMuMassRec.push_back(cm_mc_rec.parent[i].M());
    costheta_rec.push_back(calc->ct_he);
    phi_rec.push_back(calc->phi_he);

  }
  ang.setyieldfname("Yield_helicity.root");

  ang.setmassthetaphilist(DiMuMass,costheta,phi);
  ang.setmcmassthetaphilist(DiMuMassRec,costheta_rec,phi_rec);
  //ang.setmassthetaphilist(cm.DiMuMass,cm.costheta,cm.phi);


  //ang.setmcmassthetaphilist(cm_mc_rec.DiMuMass,cm_mc_rec.costheta,cm_mc_rec.phi);
  ang.FitJPsiMass(6, 1, 0.2, 0, 0.5*TMath::Pi(), -0.6, -0.3, 116, 2.2, 4.5);
  ang.FitJPsiMass(7, 1, 0.2, 0.5*TMath::Pi(), TMath::Pi(), -0.6, -0.3, 116, 2.2, 4.5);
  ang.FitJPsiMass(8, 1, 0.2, TMath::Pi(), 1.5*TMath::Pi(), -0.6, -0.3, 116, 2.2, 4.5);
  ang.FitJPsiMass(9, 1, 0.2, 1.5*TMath::Pi(),2*TMath::Pi(), -0.6, -0.3, 116, 2.2, 4.5);

  ang.FitJPsiMass(12, 1, 0.2, 0, 0.5*TMath::Pi(), -0.3, 0.3, 116, 2.2, 4.5);
  ang.FitJPsiMass(13, 1, 0.2, 0.5*TMath::Pi(), TMath::Pi(), -0.3, 0.3, 116, 2.2, 4.5);
  ang.FitJPsiMass(14, 1, 0.2, TMath::Pi(), 1.5*TMath::Pi(), -0.3, 0.3, 116, 2.2, 4.5);
  ang.FitJPsiMass(15, 1, 0.2, 1.5*TMath::Pi(),2*TMath::Pi(), -0.3, 0.3, 116, 2.2, 4.5);

  ang.FitJPsiMass(18, 1, 0.2, 0., 0.5*TMath::Pi(), 0.3, 0.6, 116, 2.2, 4.5);
  ang.FitJPsiMass(19, 1, 0.2, 0.5*TMath::Pi(), TMath::Pi(), 0.3, 0.6, 116, 2.2, 4.5);
  ang.FitJPsiMass(20, 1, 0.2, TMath::Pi(), 1.5*TMath::Pi(), 0.3, 0.6, 116, 2.2, 4.5);
  ang.FitJPsiMass(21, 1, 0.2, 1.5*TMath::Pi(),2*TMath::Pi(), 0.3, 0.6, 116, 2.2, 4.5);







  }
