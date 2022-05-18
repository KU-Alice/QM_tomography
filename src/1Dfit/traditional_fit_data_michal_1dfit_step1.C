#include "../Fitter.h"
#include "../ReadFile.h"
//#include "generator.h"

void traditional_fit_data_michal_1dfit_step1(){
  Fitter ang;
  ReadFile cm;
  ReadFile cm_mc_rec;
  //ReadFile cm_mc_gen;

  cm.ReadFileTreeData("../../../david_data/download.root",0);

  cm_mc_rec.ReadFileTreeData("../../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1);
  //cm_mc_gen.ReadFileTreeGen("../../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root");

  //cout << "mass before setting" << cm.DiMuMass.size() << endl;
  vector <float> DiMuMass;
  vector <float> costheta;
  vector <float> phi;
  vector <float> DiMuMassRec;
  vector <float> costheta_rec;
  vector <float> phi_rec;
  for(int i=0 ; i< cm.daughter1.size();i++){
    AngleCalc * calc = new AngleCalc();
    calc->AngleCalculator_DT(cm.daughter1[i],cm.daughter2[i]);
    DiMuMass.push_back(cm.parent[i].M());
    costheta.push_back(calc->ct_dt);
    phi.push_back(calc->phi_dt);
  }

  for(int i=0 ; i< cm_mc_rec.daughter1.size();i++){
    AngleCalc * calc = new AngleCalc();
    calc->AngleCalculator_DT(cm_mc_rec.daughter1[i],cm_mc_rec.daughter2[i]);
    DiMuMassRec.push_back(cm_mc_rec.parent[i].M());
    costheta_rec.push_back(calc->ct_dt);
    phi_rec.push_back(calc->phi_dt);

  }

  ang.setmassthetaphilist(DiMuMass,costheta,phi);
  ang.setmcmassthetaphilist(DiMuMassRec,costheta_rec,phi_rec);
  // phi bins
  ang.FitJPsiMass(0, 1, 0.2, 0, 0.4*TMath::Pi(), -1, 1, 116, 2.2, 4.5,0);
  ang.FitJPsiMass(1, 1, 0.2, 0.4*TMath::Pi(), 0.8*TMath::Pi(), -1, 1, 116, 2.2, 4.5,0);
  ang.FitJPsiMass(2, 1, 0.2, 0.8*TMath::Pi(), 1.2*TMath::Pi(), -1, 1, 116, 2.2, 4.5,0);
  ang.FitJPsiMass(3, 1, 0.2, 1.2*TMath::Pi(),1.6*TMath::Pi(), -1, 1, 116, 2.2, 4.5,0);
  ang.FitJPsiMass(4, 1, 0.2, 1.6*TMath::Pi(),2*TMath::Pi(), -1, 1, 116, 2.2, 4.5,0);

  //costheta bins

  ang.FitJPsiMass(0, 1, 0.2, 0.0, 2*TMath::Pi(), -1, -0.5, 116, 2.2, 4.5,1);
  ang.FitJPsiMass(1, 1, 0.2, 0.0, 2*TMath::Pi(), -0.5, 0., 116, 2.2, 4.5,1);
  ang.FitJPsiMass(2, 1, 0.2, 0.0, 2*TMath::Pi(), 0., 0.5, 116, 2.2, 4.5,1);
  ang.FitJPsiMass(3, 1, 0.2, 0.0, 2*TMath::Pi(), 0.5, 1, 116, 2.2, 4.5,1);









  }
