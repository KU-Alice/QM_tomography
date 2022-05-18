#ifndef FitterMichal_h
#define FitterMichal_h

#include "TLegend.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include <iostream>
using namespace RooFit;



class FitterMichal{


public:
  //FitterMichal();
//  ~FitterMichal();
  void setmassthetaphilist(vector <float> mass, vector <float> costheta,vector <float> phi ){
    rmass = mass;
    rcostheta = costheta;
    rphi = phi;
  }
  void setmcmassthetaphilist(vector <float> mcmass, vector <float> mccostheta,vector <float> mcphi ){
    rmcmass = mcmass;
    rmccostheta = mccostheta;
    rmcphi = mcphi;
  }
  void setyieldfname(TString fname){
    yieldfname = fname;
  }

  void FitJPsiMass(Int_t iPoint = -1,Float_t maxPt = 0.2, Float_t minPt = 0.0,
  		Float_t minPhi = 0.0, Float_t maxPhi =  6.29,
  		Float_t minCsT = -1.0, Float_t maxCsT =  1.0,
  		Int_t nBins = 116, Float_t fitMin = 2.5, Float_t fitMax = 4.5,int type =2);





private:

  vector <float> rmass, rcostheta, rphi;
  vector <float> rmcmass, rmccostheta, rmcphi;

   void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07,int columns=2);
   void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);
   TString yieldfname = "Yield_JPsi.root";


};
#endif








void FitterMichal::FitJPsiMass(Int_t iPoint = -1,Float_t maxPt = 0.2, Float_t minPt = 0.0,
		Float_t minPhi = 0.0, Float_t maxPhi =  6.29,
		Float_t minCsT = -1.0, Float_t maxCsT =  1.0,
		Int_t nBins = 116, Float_t fitMin = 2.5, Float_t fitMax = 4.5,int type =2) {



    UInt_t nDT = rmass.size();
    cout << nDT << "#############@#@##@@#@@@@@@@@@@@@@@@@@@@@@@@@@@@"<< endl;
    RooRealVar MassDT("MassDT","MassDT",fitMin,fitMax);
    MassDT.SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    RooDataSet DT_UnBin("DT_UnBin","DT_UnBin",RooArgSet(MassDT));

    //event loop
    for(Int_t i=0; i<nDT; i++){
      //DT_Tree->GetEntry(i);

      //if(TMath::Abs(DiMuY)>0.8) continue;
      if((rmass[i] <= fitMin) || (rmass[i] >= fitMax)) continue;
      //if(DiMuPt < minPt || DiMuPt > maxPt) continue;
      if((rcostheta[i] > maxCsT) || (rcostheta[i] < minCsT))continue;
      if((rphi[i] > maxPhi) || (rphi[i] < minPhi))continue;

      MassDT = rmass[i];

      DT_UnBin.add(RooArgSet(MassDT));
      }


    UInt_t nMC = rmcmass.size();

    RooRealVar MassMC("MassMC","MassMC",fitMin,fitMax);
    MassMC.SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
    RooDataSet MC_UnBin("MC_UnBin","MC_UnBin",RooArgSet(MassMC));


    //event loop
    for(Int_t i=0; i<nMC; i++){
      //MC_Tree->GetEntry(i);

      //if(TMath::Abs(DiMuY)>0.8) continue;
      if(rmcmass[i] <= fitMin || rmcmass[i] >= fitMax) continue;
      //if(DiMuPt < minPt || DiMuPt > maxPt) continue;
      if( rmccostheta[i] > maxCsT || rmccostheta[i] < minCsT)continue;
      if(rmcphi[i] > maxPhi || rmcphi[i] < minPhi)continue;
      MassMC = rmcmass[i];


      MC_UnBin.add(RooArgSet(MassMC));
      }

//________________________________________________FIT_________________________________________________________
// Crystal Ball P.D.F. variables for J/psi
RooRealVar Mean_JPsi("Mean_JPsi","m_{J/#psi}",3.1,3.0,3.2);
RooRealVar Sigma_JPsi("Sigma_JPsi","#sigma_{J/#psi}",0.025,0.01,0.04);
RooRealVar Alpha_JPsi("Alpha_JPsi","Alpha_JPsi",0.8,0.0,10);
RooRealVar N_JPsi("N_JPsi","N_JPsi",10,0,100);

//Fit Crystal Ball to JPsi MC signal and fix the alpha and N
RooCBShape CrystalBall_JPsi_MC("CrystalBall_JPsi_MC","CrystalBall_JPsi_MC",MassMC,Mean_JPsi,Sigma_JPsi,Alpha_JPsi,N_JPsi);
CrystalBall_JPsi_MC.fitTo(MC_UnBin,Range("JPsi_Signal"));

RooPlot* plot_MC3 = MassMC.frame(Title(" "));
MC_UnBin.plotOn(plot_MC3);
CrystalBall_JPsi_MC.plotOn(plot_MC3,Name("CrystalBall_Single"),LineColor(kBlue));

TCanvas* c31 = new TCanvas("c31","c31",800,600) ;
c31->cd() ; gPad->SetLeftMargin(0.15) ;
plot_MC3->GetYaxis()->SetTitleOffset(1.3);
plot_MC3->Draw();

//RooCBShape CrystalBall_JPsi_DT("CrystalBall_JPsi_DT","CrystalBall_JPsi_DT",MassDT,Mean_JPsi,Sigma_JPsi,Alpha_JPsi,N_JPsi);
RooCBShape CrystalBall_JPsi_DT("CrystalBall_JPsi_DT","CrystalBall_JPsi_DT",MassDT,Mean_JPsi,Sigma_JPsi,RooConst(Alpha_JPsi.getVal()),RooConst(N_JPsi.getVal()));
//RooCBShape CrystalBall_JPsi_DT("CrystalBall_JPsi_DT","CrystalBall_JPsi_DT",MassDT,Mean_JPsi,RooConst(Sigma_JPsi.getVal()),Alpha_JPsi,N_JPsi);

//Exponential function for background
RooRealVar Slope("Slope","Slope",-1.103,-2,2);
RooExponential BackgroundShape_DT("BackgroundShape_DT","BackgroundShape_DT",MassDT,Slope);

// Fit a sum of the functions to the data
RooRealVar nSignal_JPsi("nSignal_JPsi","N_{J/#psi}",1,0,1e06);
//RooRealVar nSignal_Psi2s("nSignal_Psi2s","N_{#psi2S}",1,0,1e06);
RooRealVar nBackground("nBackground","nBackground",1,0,1e06);

RooAddPdf DT_FitFunction("DT_FitFunction","DT_FitFunction",RooArgList(CrystalBall_JPsi_DT,BackgroundShape_DT),RooArgList(nSignal_JPsi,nBackground));

MassDT.setRange("Fit_Range",fitMin,fitMax);
DT_FitFunction.fitTo(DT_UnBin,Range("Fit_Range"));

// Make plot of binned dataset showing Poisson error bars (RooFit default)
RooPlot* plot_DT = MassDT.frame(Bins(nBins), Title(" ")) ;
DT_UnBin.plotOn(plot_DT);
DT_FitFunction.plotOn(plot_DT,Components(BackgroundShape_DT),LineStyle(9),LineWidth(iPoint == -1 ? 2:1),LineColor(kGreen+1));
DT_FitFunction.plotOn(plot_DT,LineWidth(iPoint == -1 ? 2:1),LineColor(kBlue));
DT_FitFunction.plotOn(plot_DT,Components(CrystalBall_JPsi_DT),LineColor(kMagenta),LineWidth(iPoint == -1 ? 2:1));

DT_FitFunction.plotOn(plot_DT,LineWidth(iPoint == -1 ? 2:1),LineColor(kBlue));
DT_UnBin.plotOn(plot_DT,MarkerSize(0.7),MarkerStyle(1));

TCanvas* c1 = new TCanvas("c1","c1",285*3,750) ;
c1->cd();
gPad->SetBottomMargin(0.11);gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.02);gPad->SetTopMargin(0.01);
plot_DT->GetYaxis()->SetTitleOffset(1.4);
plot_DT->GetYaxis()->SetTitleSize(0.04);
plot_DT->GetXaxis()->SetTitleSize(0.04);
plot_DT->GetYaxis()->SetTitle(TString::Format("Counts per %.0f MeV/#it{c}^{2}",(fitMax-fitMin)*1000/nBins));
plot_DT->Draw();

TLegend *myLegend1 = new TLegend(0.375,0.27,0.860,0.51);
myLegendSetUp(myLegend1,0.035,1);
myLegend1->AddEntry((TObject*)0,TString::Format("N_{J/#psi} = %.1f #pm %.1f",nSignal_JPsi.getVal(),nSignal_JPsi.getError())," ");
//myLegend1->AddEntry((TObject*)0,TString::Format("m_{J/#psi} = %1.5f #pm %1.5f GeV/#it{c}^{2}",Mean_JPsi.getVal(),Mean_JPsi.getError())," ");
//myLegend1->AddEntry((TObject*)0,TString::Format("#sigma_{J/#psi} = %5.3f #pm %5.3f GeV/#it{c}^{2}",Sigma_JPsi.getVal(),Sigma_JPsi.getError())," ");
//myLegend1->AddEntry((TObject*)0,TString::Format("#alpha_{J/#psi} = %5.3f #pm %5.3f",Alpha_JPsi.getVal(),Alpha_JPsi.getError())," ");
//myLegend1->AddEntry((TObject*)0,TString::Format("n_{J/#psi} = %5.3f #pm %5.3f",N_JPsi.getVal(),N_JPsi.getError())," ");
//myLegend1->AddEntry((TObject*)0,TString::Format("Slope = %5.3f #pm %5.3f",Slope.getVal(),Slope.getError())," ");
myLegend1->AddEntry((TObject*)0,TString::Format("#chi^{2} = %5.3f",plot_DT->chiSquare("DT_FitFunction_Norm[MassDT]_Range[fit_nll_DT_FitFunction_DT_UnBin]_NormRange[fit_nll_DT_FitFunction_DT_UnBin]", "h_DT_UnBin", 6))," ");

//myLegend1->Draw();

Float_t lumi = 233.0;
TLatex* latex = new TLatex();
latex->SetNDC();
latex->SetTextSize(0.045);
latex->SetTextAlign(22);
latex->SetTextFont(42);
latex->DrawLatex(0.56,0.94,"ALICE, Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
latex->SetTextSize(0.045);
latex->SetTextAlign(12);
latex->DrawLatex(0.55,0.85,"J/#psi #rightarrow #mu^{#plus} #mu^{#minus}");
//latex->DrawLatex(0.55,0.77,Form("UPC, L_{#lower[-0.3]{int}} = %.0f #pm %.0f #mub^{-1}",lumi,lumi*0.027));
latex->DrawLatex(0.55,0.77,"0.2 < #it{p}_{#lower[-0.3]{T}} < 1.0 GeV/#it{c}");
latex->DrawLatex(0.18,0.77,Form("|#it{y}| < 0.8"));
latex->DrawLatex(0.55,0.69,Form("%1.1f#pi < #phi < %1.1f#pi",minPhi/TMath::Pi(),maxPhi/TMath::Pi()));
latex->DrawLatex(0.55,0.61,Form("%1.1f < cos(#theta) < %1.1f",minCsT,maxCsT));

latex->DrawLatex(0.55,0.43,Form("N_{J/#psi} = %.1f #pm %.1f",nSignal_JPsi.getVal(),nSignal_JPsi.getError()));
latex->DrawLatex(0.55,0.35,Form("#chi^{2}/#it{dof} = %.2f",plot_DT->chiSquare("DT_FitFunction_Norm[MassDT]_Range[fit_nll_DT_FitFunction_DT_UnBin]_NormRange[fit_nll_DT_FitFunction_DT_UnBin]", "h_DT_UnBin", 6)));



if(iPoint != -1){
  TFile *yieldFile = new TFile(yieldfname.Data(),"UPDATE");
  TH2F *yieldPhiCsT =  0x0;
  TH1F *yieldPhi =  0x0;
  TH1F *yieldCsT =  0x0;
  yieldPhiCsT = (TH2F*)yieldFile->Get("yieldPhiCsTStore");
  yieldCsT = (TH1F*)yieldFile->Get("yieldCsTStore");
  yieldPhi = (TH1F*)yieldFile->Get("yieldPhiStore");
  Float_t binsPhi[6] = {0.,0.4*TMath::Pi(),0.8*TMath::Pi(),1.2*TMath::Pi(),1.6*TMath::Pi(),2*TMath::Pi()};
  Float_t binsPhi2d[5] = {0.,0.5*TMath::Pi(),1*TMath::Pi(),1.5*TMath::Pi(),2*TMath::Pi()};
  //Float_t binsCsT[4] = {-1,-0.5,0.5,1};
  Float_t binsCsT[5] = {-0.6,-0.3,0.0,0.3,0.6};
  Float_t binsCsT2d[4] = {-0.6,-0.3,0.3,0.6};



  if (type==0){
    if(!yieldPhi){

      //Float_t binsPhi[5] = {0,0.5*TMath::Pi(),TMath::Pi(),1.5*TMath::Pi(),2*TMath::Pi()};
      yieldPhi = new TH1F("yieldPhi","yieldPhi",5,binsPhi);
      yieldPhi->SetTitle("J/#psi #rightarrow l^{+} l^{-}");
      yieldPhi->GetXaxis()->SetTitle("#phi");


    }

    yieldPhi->SetName("yieldPhi");
    yieldPhi->SetBinContent(iPoint+1,nSignal_JPsi.getVal());
    yieldPhi->SetBinError(iPoint+1,nSignal_JPsi.getError());


    yieldFile->Delete("yieldPhiStore;*");
    yieldPhi->SetName("yieldPhiStore");
    yieldPhi->Write();


  }

  else if (type==1){
    if(!yieldCsT){
      yieldCsT = new TH1F("yieldCsT","yieldCsT",4,binsCsT);
      yieldCsT->SetTitle("J/#psi #rightarrow l^{+} l^{-}");
      yieldCsT->GetXaxis()->SetTitle("cos(#theta)");



    }
    yieldCsT->SetName("yieldCsT");
    yieldCsT->SetBinContent(iPoint+1,nSignal_JPsi.getVal());
    yieldCsT->SetBinError(iPoint+1,nSignal_JPsi.getError());
    yieldFile->Delete("yieldCsTStore;*");
    yieldCsT->SetName("yieldCsTStore");
    yieldCsT->Write();
  }

  else{
    if(!yieldPhiCsT ){
      cout<<"Creating new histo"<<endl;

      yieldPhiCsT = new TH2F("yieldPhiCsT","yieldPhiCsT",4,binsPhi2d,3,binsCsT2d);
      yieldPhiCsT->SetTitle("J/#psi #rightarrow l^{+} l^{-}");
      yieldPhiCsT->GetXaxis()->SetTitle("#phi");
      yieldPhiCsT->GetYaxis()->SetTitle("cos(#theta)");
      }
      yieldPhiCsT->SetName("yieldPhiCsT");
      yieldPhiCsT->SetBinContent(iPoint+1,nSignal_JPsi.getVal());
      yieldPhiCsT->SetBinError(iPoint+1,nSignal_JPsi.getError());
      //cout << "###########bincontent is ######"<<yieldPhiCsT->GetBin(0.5, 0.6) <<endl;
      yieldFile->Delete("yieldPhiCsTStore;*");
      yieldPhiCsT->SetName("yieldPhiCsTStore");
      yieldPhiCsT->Write();
  }











  c1->SetName(TString::Format("FitDiLepton_%d",iPoint));
  c1->Write();
  yieldFile->Close();
  }


}


void FitterMichal::myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07,int columns=2){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  currentLegend->SetNColumns(columns);
  return;
}

void FitterMichal::myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}
