#include "../ReadFile.h"

//#include "gROOt.h"
void traditional_fit_data_michal_2dfit_step2(int choose = 1){

  gStyle->SetOptStat(0);
  ReadFile cm;
  ReadFile cm_mc_rec;
  ReadFile cm_mc_gen;
  cm.ReadFileTreeData("../../../david_data/download.root",0);

  cm_mc_rec.ReadFileTreeData("../../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1);
  cm_mc_gen.ReadFileTreeGen("../../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root");

  //TFile* YieldFile = TFile::Open("Yield_helicity.root");
  TFile* YieldFile = TFile::Open("Yield_JPsi.root");
  //TFile* YieldFile = TFile::Open("Yield_JPsi.root");
  //YieldFile->ls();
  TH2F* yieldPhiCsTStore_read = (TH2F*) YieldFile->Get("yieldPhiCsTStore");
  yieldPhiCsTStore_read->Draw("COLZtext");

  TH2F* yieldPhiCsRec;
  TH2F* yieldPhiCsGen;
  yieldPhiCsRec = (TH2F*) yieldPhiCsTStore_read->Clone("yieldPhiCsRec");
  yieldPhiCsGen = (TH2F*) yieldPhiCsTStore_read->Clone("yieldPhiCsGen");

  //yieldPhiCsRec = new TH2F("yieldPhiCsRec","yieldPhiCsRec",1000,-1,1); //yieldPhiCsTStore_read->Clone("yieldPhiCsRec");
  //yieldPhiCsGen = (TH2F*) yieldPhiCsTStore_read->Clone("yieldPhiCsGen");

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


  vector <float> phi_rec;
  vector <float> costheta_rec;
  vector <float> phi_gen;
  vector <float> costheta_gen;

  for(int i = 0;i<cm_mc_gen.daughter1.size();i++){
    AngleCalc * calc = new AngleCalc();
    calc->AngleCalculator_DT(cm_mc_gen.daughter1[i],cm_mc_gen.daughter2[i]);
    yieldPhiCsGen->Fill(calc->phi_dt,calc->ct_dt);

    costheta_gen.push_back(calc->ct_dt);
    phi_gen.push_back(calc->phi_dt);
  }

  for(int i = 0;i<cm_mc_rec.daughter1.size();i++){
    AngleCalc * calc = new AngleCalc();
    calc->AngleCalculator_DT(cm_mc_rec.daughter1[i],cm_mc_rec.daughter2[i]);
    yieldPhiCsRec->Fill(calc->phi_dt,calc->ct_dt);
    costheta_rec.push_back(calc->ct_dt);
    phi_rec.push_back(calc->phi_dt);
  }


  //yieldPhiCsRec->Draw("COLZtext");
  //yieldPhiCsRec->Sumw2();
  //yieldPhiCsGen->Sumw2();

  yieldPhiCsTStore_read->Sumw2();
  yieldPhiCsaxe->Divide(yieldPhiCsRec,yieldPhiCsGen,1,1,"B");
  yieldPhiCsaxe->SetTitle("2D A X E (Reco/Gen)");

  yieldPhiCsTStore_read->SetTitle("Yield Before Correction");
  //yieldPhiCsTStore_read->Draw("COLZtext");

  yieldPhiCsaxe->Sumw2();
  yieldPhiCsCorrected->SetTitle("Corrected Yield");
  //yieldPhiCsCorrected->Divide(yieldPhiCsRec,yieldPhiCsaxe);
  yieldPhiCsCorrected->Divide(yieldPhiCsTStore_read,yieldPhiCsaxe);

  //yieldPhiCsCorrected->Draw("COLZtext");
  TH1D* costhetaaxe = (TH1D*) yieldPhiCsaxe->ProjectionY();
  TH1D* phiaxe = (TH1D*) yieldPhiCsaxe->ProjectionX();
  //costhetaaxe->Draw("text0");
  //phiaxe->Draw("text0");


  TH1D* costheta;
  costheta = (TH1D*) yieldPhiCsCorrected->ProjectionY();
  costheta->SetTitle("2D Corrected CosTheta");
  //costheta->Draw("e");

  TH1D* phi;
  phi = (TH1D*) yieldPhiCsCorrected->ProjectionX();
  phi->SetTitle("2D Corrected Phi");
  //phi->Draw("e");
  //yieldPhiCsaxe->Draw("COLZtext");
  //yieldPhiCsGen->Draw("COLZtext");
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
  costheta->Draw("text0e");
  //TwoDCorrectionprojection->cd(2);
  TwoDCorrection->cd(7);
  phi->Draw("text0e");


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

  Double_t Ngen1=  yieldPhiCsGen->GetBinContent(8);
  Double_t Nrec1=  yieldPhiCsRec->GetBinContent(8);

  Double_t eff = Nrec1/Ngen1;
  cout << "read value is: " << Ngen1 << ", " << Nrec1 << " , "<<eff << endl;
  Double_t Calc_Error = TMath::Sqrt((eff*(1-eff))/Ngen1);
  cout << "calculated error is : " << Calc_Error<< endl;
  cout << "error from plot is : " << yieldPhiCsaxe->GetBinError(8)<< endl;

  TFile * file = new TFile("2D_Corrected.root","recreate");
  //TFile * file = new TFile("2D_Corrected_helicity.root","recreate");
  costheta->Write();
  phi->Write();
  yieldPhiCsCorrected->Write();
  costhetaaxe->Write();
  phiaxe->Write();
  yieldPhiCsaxe->Write();
  yieldPhiCsGen->Write();
  yieldPhiCsRec->Write();
  yieldPhiCsTStore_read->Write();

  Plotter plot;
  plot.SetBinX(4);
  plot.SetBinY(3);
  plot.Plot2D(phi_rec,costheta_rec,"Generated;#Phi;Cos#Theta" , "COLZtext0","myanvas" );











}
