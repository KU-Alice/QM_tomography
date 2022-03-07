#include "Fitter.h"
#include "Correlation.h"
#include "generator.h"

void traditional_fit_data_michal_1dfit_step1(){


  Fitter ang;
  Correlation cm;
  Correlation cm_mc_rec;
  Correlation cm_mc_gen;
  cm.ReadTree("../../david_data/download.root",0,0);

  cm_mc_rec.ReadTree("../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1,0);
  cm_mc_gen.ReadTree("../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1,1);

  //cout << "mass before setting" << cm.DiMuMass.size() << endl;



  ang.setmassthetaphilist(cm.DiMuMass,cm.costheta,cm.phi);
  ang.setmcmassthetaphilist(cm_mc_rec.DiMuMass,cm_mc_rec.costheta,cm_mc_rec.phi);
  // phi bins
  ang.FitJPsiMass(0, 1, 0.2, 0, 0.5*TMath::Pi(), -1, 1, 116, 2.2, 4.5);
  ang.FitJPsiMass(1, 1, 0.2, 0.5*TMath::Pi(), TMath::Pi(), -1, 1, 116, 2.2, 4.5);
  ang.FitJPsiMass(2, 1, 0.2, TMath::Pi(), 1.5*TMath::Pi(), -1, 1, 116, 2.2, 4.5);
  ang.FitJPsiMass(3, 1, 0.2, 1.5*TMath::Pi(),2*TMath::Pi(), -1, 1, 116, 2.2, 4.5);

  //costheta bins
  ang.FitJPsiMass(0, 1, 0.2, 0.0, 2*TMath::Pi(), -1, -0.5, 116, 2.2, 4.5);
  ang.FitJPsiMass(1, 1, 0.2, 0.0, 2*TMath::Pi(), -0.5, 0.5, 116, 2.2, 4.5);
  ang.FitJPsiMass(2, 1, 0.2, 0.0, 2*TMath::Pi(), 0.5, 1, 116, 2.2, 4.5);









  }
