
#include "../ReadFile.h"

//#include "gROOt.h"
void traditional_fit_data_michal_1dfit_step2(int choose = 1){

  gStyle->SetOptStat(0);
  ReadFile cm;
  ReadFile cm_mc_rec;
  ReadFile cm_mc_gen;
  //cm.ReadTree("../../david_data/download.root",0,0);

  cm_mc_rec.ReadFileTreeData("../../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1);
  cm_mc_gen.ReadFileTreeGen("../../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root");

  TFile* YieldFile = TFile::Open("Yield_JPsi_1D_Step1.root");
//  YieldFile->ls();
  TH1F* yieldPhiStore_read = (TH1F*) YieldFile->Get("yieldPhiStore");
  TH1F* yieldCsTStore_read = (TH1F*) YieldFile->Get("yieldCsTStore");

  //yieldCsTStore_read->Draw("text0");
  //yieldPhiStore_read->Draw("text0");
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
  yieldCsTCorrected = (TH1F*) yieldCsTStore_read->Clone("yieldCstCorrected");
  yieldPhiCorrected = (TH1F*) yieldPhiStore_read->Clone("yieldPhiCorrected");
  yieldCsTCorrected->Reset();
  yieldPhiCorrected->Reset();

  for(int i = 0;i<cm_mc_gen.daughter1.size();i++){
    AngleCalc * calc = new AngleCalc();
    calc->AngleCalculator_DT(cm_mc_gen.daughter1[i],cm_mc_gen.daughter2[i]);
    yieldCsTGen->Fill(calc->ct_dt);
    yieldPhiGen->Fill(calc->phi_dt);


  }

  for(int i = 0;i<cm_mc_rec.daughter1.size();i++){
    AngleCalc * calc = new AngleCalc();
    calc->AngleCalculator_DT(cm_mc_rec.daughter1[i],cm_mc_rec.daughter2[i]);
    yieldCsTRec->Fill(calc->ct_dt);
    yieldPhiRec->Fill(calc->phi_dt);


  }

  //return;

  yieldPhiGen->Sumw2();
  yieldPhiRec->Sumw2();

  yieldCsTGen->Sumw2();
  yieldCsTRec->Sumw2();

  //yieldCsTGen->Draw("text0");

  yieldPhiAXE->Divide(yieldPhiRec,yieldPhiGen,1,1,"B");
  yieldCsTAXE->Divide(yieldCsTRec,yieldCsTGen,1,1,"B");

  //yieldPhiAXE->Sumw2();
  //yieldCsTAXE->Sumw2();

  //yieldCsTAXE->Draw("text0");

  yieldPhiCorrected->Divide(yieldPhiStore_read,yieldPhiAXE);
  yieldCsTCorrected->Divide(yieldCsTStore_read,yieldCsTAXE);
  //yieldPhiCorrected->Divide(yieldPhiRec,yieldPhiAXE);
  //yieldCsTCorrected->Divide(yieldCsTRec,yieldCsTAXE);


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




  yieldCsTCorrected->Draw("text0");
  TCanvas* canvas = new TCanvas("canvas","", 900,600);
  canvas->Divide(3,2);

  canvas->cd(1);
  yieldCsTStore_read->Draw("text0");

  canvas->cd(2);

  yieldCsTAXE->Draw("text0");


  canvas->cd(3);

  yieldCsTCorrected->Draw("text0");


  canvas->cd(4);
  yieldPhiStore_read->Draw("text0");

  canvas->cd(5);
  yieldPhiAXE->Draw("text0");

  canvas->cd(6);
  yieldPhiCorrected->Draw("text0");


  TCanvas* canvas2 = new TCanvas("canvas2","", 900,600);
  canvas2->Divide(3,2);

  canvas2->cd(1);
  yieldCsTStore_read->Draw("text0");

  canvas2->cd(2);
  yieldCsTGen->Draw("text0");
  canvas2->cd(3);
  yieldCsTRec->Draw("text0");


  canvas2->cd(4);

  yieldPhiStore_read->Draw("text0");
  canvas2->cd(5);

  yieldPhiGen->Draw("text0");

  canvas2->cd(6);
  yieldPhiRec->Draw("text0");

 Double_t Ngen1=  yieldPhiGen->GetBinContent(2);
 Double_t Nrec1=  yieldPhiRec->GetBinContent(2);

 Double_t eff = Nrec1/Ngen1;
 cout << "read value is: " << Ngen1 << ", " << Nrec1 << " , "<<eff << endl;
 Double_t Calc_Error = TMath::Sqrt((eff*(1-eff))/Ngen1);
 cout << "calculated error is : " << Calc_Error<< endl;
 cout << "error from plot is : " << yieldPhiAXE->GetBinError(2)<< endl;

 // Saving the Corrected Distributions into files
 TFile * file = new TFile("1D_Corrected.root","recreate");
 yieldCsTCorrected->Write();
 yieldPhiCorrected->Write();
























  //yieldPhiCsTStore_read->Draw("COLZtext");


  }
