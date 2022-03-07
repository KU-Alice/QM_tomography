#include "Fitter.h"
#include "Correlation.h"
#include "generator.h"

//#include "gROOt.h"
void traditional_fit_data_michal_1dfit_step2(int choose = 1){

  gStyle->SetOptStat(0);
  Fitter ang;
  Correlation cm;
  Correlation cm_mc_rec;
  Correlation cm_mc_gen;
  //cm.ReadTree("../../david_data/download.root",0,0);

  cm_mc_rec.ReadTree("../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1,0);
  cm_mc_gen.ReadTree("../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1,1);

  TFile* YieldFile = TFile::Open("Yield_JPsi_1DFit.root");
  YieldFile->ls();
  TH1F* yieldPhiStore_read = (TH1F*) YieldFile->Get("yieldPhiStore");
  TH1F* yieldCsTStore_read = (TH1F*) YieldFile->Get("yieldCsTStore");
  //yieldCsTStore_read->Draw();
  yieldPhiStore_read->Sumw2();
  yieldCsTStore_read->Sumw2();



  //TH1F* yieldPhiData;
  //TH1F* yieldCostData;
  TH1F* yieldPhiRec;
  TH1F* yieldPhiGen;

  yieldPhiRec = (TH1F*) yieldPhiStore_read->Clone("yieldPhiRec");
  yieldPhiGen = (TH1F*) yieldPhiStore_read->Clone("yieldPhiGen");

  yieldPhiGen->Reset();
  yieldPhiRec->Reset();


  TH1F* yieldCsTRec;
  TH1F* yieldCsTGen;
  yieldCsTRec = (TH1F*) yieldCsTStore_read->Clone("yieldCsTRec");
  yieldCsTGen = (TH1F*) yieldCsTStore_read->Clone("yieldCsTGen");

  yieldCsTGen->Reset();
  yieldCsTRec->Reset();




  TH1F* yieldPhiAXE;
  TH1F* yieldCsTAXE;
  yieldCsTAXE = (TH1F*) yieldCsTStore_read->Clone("yieldCsTAXE");
  yieldPhiAXE = (TH1F*) yieldPhiStore_read->Clone("yieldPhiAXE");
  yieldPhiAXE->Reset();
  yieldCsTAXE->Reset();


  TH1F* yieldPhiCorrected;
  TH1F* yieldCsTCorrected;
  yieldCsTCorrected = (TH1F*) yieldCsTStore_read->Clone("yieldCstCsCorrected");
  yieldPhiCorrected = (TH1F*) yieldPhiStore_read->Clone("yieldPhiCorrected");
  yieldCsTCorrected->Reset();
  yieldPhiCorrected->Reset();

  for(int i = 0;i<cm_mc_gen.costheta.size();i++){
    yieldCsTGen->Fill(cm_mc_gen.costheta[i]);
    yieldPhiGen->Fill(cm_mc_gen.phi[i]);


  }

  for(int i = 0;i<cm_mc_rec.costheta.size();i++){
    yieldCsTRec->Fill(cm_mc_gen.costheta[i]);
    yieldPhiRec->Fill(cm_mc_gen.phi[i]);


  }

  yieldPhiGen->Sumw2();
  yieldPhiRec->Sumw2();

  yieldCsTGen->Sumw2();
  yieldCsTRec->Sumw2();

  //yieldCsTGen->Draw();

  yieldPhiAXE->Divide(yieldPhiRec,yieldPhiGen);
  yieldCsTAXE->Divide(yieldCsTRec,yieldCsTGen);

  yieldPhiAXE->Sumw2();
  yieldCsTAXE->Sumw2();

  yieldCsTAXE->Draw();

  yieldPhiCorrected->Divide(yieldPhiStore_read,yieldPhiAXE);
  yieldCsTCorrected->Divide(yieldCsTStore_read,yieldCsTAXE);


  yieldPhiAXE->SetTitle("1D Phi A X E (Reco/Gen)");
  yieldCsTAXE->SetTitle("1D Costheta A X E (Reco/Gen)");
  yieldCsTCorrected->SetTitle("1D Costheta Corrected");
  yieldPhiCorrected->SetTitle("1D Phi Corrected");

  yieldCsTStore_read->SetTitle("Costheta Before Correction");
  yieldPhiStore_read->SetTitle("Phi Before Correction");

  yieldCsTGen->SetTitle("Generated Costheta");
  yieldPhiGen->SetTitle("Generated Phi");

  yieldCsTRec->SetTitle("Reconstructed Costheta");
  yieldPhiRec->SetTitle("Reconstructed Phi");





  TCanvas* canvas = new TCanvas("canvas","", 900,600);
  canvas->Divide(3,2);

  canvas->cd(1);
  yieldCsTStore_read->Draw();

  canvas->cd(2);

  yieldCsTAXE->Draw();


  canvas->cd(3);

  yieldCsTCorrected->Draw();


  canvas->cd(4);
  yieldPhiStore_read->Draw();

  canvas->cd(5);
  yieldPhiAXE->Draw();

  canvas->cd(6);
  yieldPhiCorrected->Draw();


  TCanvas* canvas2 = new TCanvas("canvas2","", 900,600);
  canvas2->Divide(3,2);

  canvas2->cd(1);
  yieldCsTStore_read->Draw();

  canvas2->cd(2);
  yieldCsTGen->Draw();
  canvas2->cd(3);
  yieldCsTRec->Draw();


  canvas2->cd(4);

  yieldPhiStore_read->Draw();
  canvas2->cd(5);

  yieldPhiGen->Draw();

  canvas2->cd(6);
  yieldPhiRec->Draw();
























  //yieldPhiCsTStore_read->Draw("COLZtext");


  }
