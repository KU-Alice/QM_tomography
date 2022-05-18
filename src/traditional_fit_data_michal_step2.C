#include "Fitter.h"
#include "Correlation.h"
#include "generator.h"

//#include "gROOt.h"
void traditional_fit_data_michal_step2(int choose = 1){

  gStyle->SetOptStat(0);
  Fitter ang;
  Correlation cm;
  Correlation cm_mc_rec;
  Correlation cm_mc_gen;
  //cm.ReadTree("../../david_data/download.root",0,0);

  cm_mc_rec.ReadTree("../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1,0);
  cm_mc_gen.ReadTree("../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1,1);

  TFile* YieldFile = TFile::Open("Yield_JPsi_2DFit.root");
  YieldFile->ls();
  TH2F* yieldPhiCsTStore_read = (TH2F*) YieldFile->Get("yieldPhiCsTStore");
  //yieldPhiCsTStore_read->Draw();

  TH2F* yieldPhiCsRec;
  TH2F* yieldPhiCsGen;
  yieldPhiCsRec = (TH2F*) yieldPhiCsTStore_read->Clone("yieldPhiCsRec");
  yieldPhiCsGen = (TH2F*) yieldPhiCsTStore_read->Clone("yieldPhiCsGen");

  TH2F* yieldPhiCsaxe;
  TH2F* yieldPhiCsCorrected;
  yieldPhiCsaxe = (TH2F*) yieldPhiCsTStore_read->Clone("yieldPhiCsaxe");
  yieldPhiCsCorrected = (TH2F*) yieldPhiCsTStore_read->Clone("yieldPhiCsCorrected");
  //yieldPhiCsGen = (TH2F*)
  //yieldPhiCsTStore_read->Copy(*yieldPhiCsGen);
  yieldPhiCsGen->Reset();
  yieldPhiCsRec->Reset();
  yieldPhiCsaxe->Reset();
  yieldPhiCsCorrected->Reset();
  //yieldPhiCsGen->Draw("");

  for(int i = 0;i<cm_mc_rec.costheta.size();i++){

    yieldPhiCsRec->Fill(cm_mc_rec.phi[i],cm_mc_rec.costheta[i]);




  }

  for(int i = 0;i<cm_mc_gen.costheta.size();i++){
    yieldPhiCsGen->Fill(cm_mc_gen.phi[i],cm_mc_gen.costheta[i]);


  }

  //yieldPhiCsRec->Draw("COLZtext");
  yieldPhiCsRec->Sumw2();
  yieldPhiCsGen->Sumw2();

  yieldPhiCsTStore_read->Sumw2();
  yieldPhiCsaxe->Divide(yieldPhiCsRec,yieldPhiCsGen);
  yieldPhiCsaxe->SetTitle("2D A X E (Reco/Gen)");

  yieldPhiCsTStore_read->SetTitle("Yield Before Correction");
  //yieldPhiCsTStore_read->Draw("COLZtext");

  yieldPhiCsaxe->Sumw2();
  yieldPhiCsCorrected->SetTitle("Corrected Yield");
  yieldPhiCsCorrected->Divide(yieldPhiCsTStore_read,yieldPhiCsaxe);

  //yieldPhiCsCorrected->Draw("COLZtext");

  TH1D* costheta;
  costheta = (TH1D*) yieldPhiCsCorrected->ProjectionY();
  costheta->SetTitle("2D Corrected CosTheta");
  //costheta->Draw("e");

  TH1D* phi;
  phi = (TH1D*) yieldPhiCsCorrected->ProjectionX();
  phi->SetTitle("2D Corrected Phi");
  //phi->Draw("e");


  TCanvas*TwoDCorrection = new TCanvas("TwoDCorrection", "2 Dimensional Correction", 1200,600 );
  TwoDCorrection->Divide(4,2);
  TwoDCorrection->cd(1);
  yieldPhiCsTStore_read->Draw("COLZtext");
  TwoDCorrection->cd(2);
  yieldPhiCsaxe->Draw("COLZtext");
  TwoDCorrection->cd(3);
  yieldPhiCsCorrected->Draw("COLZtext");

  TCanvas *TwoDCorrectionprojection = new TCanvas("TwoDCorrectionprojection", "2 Dimensional Correction projection into 1D", 600,300 );
  //TwoDCorrectionprojection->Divide(2,1);
  //TwoDCorrectionprojection->cd(1);
  TwoDCorrection->cd(4);
  costheta->Draw("e");
  //TwoDCorrectionprojection->cd(2);
  TwoDCorrection->cd(7);
  phi->Draw("e");


 // Comparison of Reocnstructed And Generated



  TCanvas* CosthetaPhi = new TCanvas("CosthetaPhi","Comparison of Reco and Data",2400,600);
  CosthetaPhi->Divide(4,2);


  CosthetaPhi->cd(1);
  TH1F* RecoCostheta = (TH1F*) yieldPhiCsRec->ProjectionY("RecoCostehta");
  RecoCostheta->SetTitle("ReConstructed CosTheta of Inco. J/#Psi");
  RecoCostheta->Draw("e");


  CosthetaPhi->cd(2);

  yieldPhiCsRec->SetTitle("2D Plot of ReConstructed CosThetaPhi of Inco. J/#Psi");
  yieldPhiCsRec->Draw("COLZtext");
  CosthetaPhi->cd(3);
  yieldPhiCsTStore_read->SetTitle("2DPlot of Costhetha Phi of Inco. J/#Psi Data");
  yieldPhiCsTStore_read->Draw("COLZtext");


  CosthetaPhi->cd(4);
  TH1F* DataCostheta = (TH1F*) yieldPhiCsTStore_read->ProjectionY("DataCostheta");
  DataCostheta->SetTitle("CosTheta of J/#Psi Data");
  DataCostheta->Draw("e");

  CosthetaPhi->cd(6);
  TH1F* ReCoPhi = (TH1F*) yieldPhiCsRec->ProjectionX("ReCoPhi");
  ReCoPhi->SetTitle("ReConstructed Phi of Inco. J/#Psi");
  ReCoPhi->Draw("e");

  CosthetaPhi->cd(7);
  TH1F* DataPhi = (TH1F*) yieldPhiCsTStore_read->ProjectionX("DataPhi");
  DataPhi->SetTitle(" Phi of Inco. J/#Psi Data");
  DataPhi->Draw("e");
















}
