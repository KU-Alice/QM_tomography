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
#include "TLegend.h"
#include "RooAddition.h"
#include <complex>
#include <iostream>
using namespace RooFit;


void Final_Fit(){
TFile * One_D_Corrected = TFile::Open("1D_Corrected.root");
//One_D_Corrected->ls();

//TFile *	Two_D_Corrected	= TFile::Open("2D_Corrected.root");
//Two_D_Corrected->ls();
gStyle->SetOptStat(0);

yieldTildePhi->Divide(yieldTildePhiRec,yieldTildePhiAxE);
TH1D* Phi_1D = (TH1D*) One_D_Corrected->Get("yieldPhiCorrected");
TH1D* Cst_1D = (TH1D*) One_D_Corrected->Get("yieldCstCorrected");

//TH1D* Phi_2D = (TH1D*) Two_D_Corrected->Get("yieldPhiCsCorrected_px");
//TH1D* Cst_2D = (TH1D*) Two_D_Corrected->Get("yieldPhiCsCorrected_py");




TF1* wphi = new TF1("wphi","[0]*(1 + 2*[1]/(3+[2]) * TMath::Cos(2*x))",0,2*3.141616);
wphi->SetParameters(4000,0.049,1.208);
wphi->SetParLimits(2,0.8,1.5);
//wphi.FixParameter(1,0.4);
//wphi.FixParameter(2,1.5);
TF1* wtheta =new TF1("wtheta", "([0]/(3+[1]))*(1+[1]*x*x)",-0.6,0.6);
wtheta->SetParameters(7000,1.4);


//Phi_1D->Draw();
//Phi_2D->Draw();
//Cst_1D->Draw();
//Cst_2D->Draw();

RooRealVar mt1 = {"mt1", "mt1", 0};
RooRealVar mt2 = {"mt2", "mt2", 0.02,-1,1};
RooRealVar mt3= {"mt3", "mt3",0.7 -1,1};
RooRealVar mt4= {"mt4", "mt4",0};//-1.0,1.};
RooRealVar mt5= {"mt5", "mt5",0};//-1,1};
RooRealVar mt6= {"mt6", "mt6",0.5,-1,1.};
RooRealVar mt7= {"mt7", "mt7",0};//-1.0,1.};
RooRealVar mt8= {"mt8", "mt8",0};//,-1.0,1.};
RooRealVar mt9= {"mt9", "mt9",0};//,-1.0,1.};

RooRealVar t = {"t", "Cos(#theta)", -1., 1.} ;
RooRealVar p = {"p", "Phi", 0., 2*3.141616} ;


//RooRealVar t_tot("t_tot", "Cos(theta) [sig+bkg]", -1., 1.) ;
//RooRealVar p_tot("p_tot", "phi angle [sig+bkg]", 0., 2*3.141616) ;
RooGenericPdf Pdf_theta("Pdf_theta", "Pdf_theta", "(1.5708*(1+ mt3*mt3) + 1.3823*(mt2*mt5+mt7*mt8-mt6*mt9)*t + (1.5708-4.71239*mt3*mt3)*t*t)", RooArgSet(t, mt2,mt3,mt5,mt6,mt7,mt8,mt9));


RooGenericPdf Pdf_phi("Pdf_phi", "Pdf_phi", "(0.66666667 + 0.345575*mt3*mt9*cos(p)+0.33333333*(-1+2*mt2*mt2+mt3*mt3+2*mt8*mt8+2*mt9*mt9)*cos(2*p)-0.345575*mt3*mt7*sin(p)-0.66666667*(mt2*mt4+mt6*mt8+mt7*mt8)*sin(2*p))", RooArgSet(p, mt2,mt3,mt4,mt6,mt7,mt8,mt9));

RooGenericPdf Pdf_dN_dOmega("Pdf_dN_dOmega", "Pdf_dN_dOmega", "(0.25*(1+mt3*mt3) +0.25*(1-3*mt3*mt3)*t*t -0.5* mt3* mt6* 2* TMath::Sqrt(1- t*t)* t* cos(p)+0.25*(2*mt2*mt2+mt3*mt3-1)*(1-t*t)*cos(2*p))", RooArgSet(p,t, mt2,mt3,mt6));

//TF1* pdfdndomega_root = new TF1("pdfdndomega_root","0.25*(1+[3]*[3]) +0.25*(1-3*[3]*[3])*t*t -0.5* [3]* mt6* 2* TMath::Sqrt(1- t*t)* t* cos(p)+0.25*(2*[2]*[2]+[3]*[3]-1)*(1-t*t)*cos(2*p))",);

TF1* pdftheta_root = new TF1("pdftheta_root","[0]*(1.5708*(1+ [1]*[1])+ (1.5708-4.71239*[1]*[1])*x*x)",-1,1);
TF1* pdfphi_root = new TF1("pdfphi_root","[0]*(0.66666667 +0.33333333*(-1+2*[1]*[1]+[2]*[2])*cos(2*x))",0,2*TMath::Pi());
//TF1* pdfphi_root = new TF1("pdfphi_root","[0]*(0.66666667)",0,2*3.141616);
pdftheta_root->SetParameter(0,100);
pdftheta_root->SetParameter(1,0.5);

pdfphi_root->SetParameter(0,100);
pdfphi_root->SetParameter(1,0.05);
pdfphi_root->SetParameter(2,0.05);
pdfphi_root->SetParLimits(1,-1,1);
pdfphi_root->SetParLimits(2,-1,1);
//pdfphi_root->FixParameter(2,7.42211e-01);
//pdfphi_root->FixParameter(0,1.00888e+03*((4.*TMath::Pi())/(3.)));

//pdfphi_root->FixParameter(3,7.94132e-01);





RooDataHist cst1d_histo( "cst1d_histo", "costheta 1d data", t, Cst_1D );
RooDataHist phi1d_histo( "phi1d_histo", "phi 1d data", p, Phi_1D );
//RooAddition phi_cst("phi_cst","phi_cst",RooArgSet(cst1d_histo,phi1d_histo));
//RooWorkspace w(“w”,1) ;

//RooDataHist phi_cst("phi_cst","phi_cst",RooArgSet(cst1d_histo,phi1d_histo));
RooCategory c("c","c");
c.defineType("Phi");
c.defineType("Cst");

RooDataHist phi_cst("phi_cst","phi_cst",RooArgSet(t,p),Index(c),Import("Phi",phi1d_histo),Import("Cst",cst1d_histo));
//RooDataSet *dataSum = Pdf_dN_dOmega.generate(RooArgSet(t,p), 10000);
//RooDataHist phi_cst("phi_cst","phi_cst",RooArgSet(t,p),Index(c),Import("Cst",*Cst_1D),Import("Phi",*Phi_1D));
//RooDataHist()
 //RooFitResult *rfit = Pdf_dN_dOmega.fitTo(phi_cst,SumW2Error(true));
//RooFitResult *rcst_phi_1d = Pdf_dN_dOmega.fitTo(phi_cst);
RooFitResult *rcst_1d = Pdf_theta.fitTo(cst1d_histo,SumW2Error(true));
RooFitResult *rphi_1d = Pdf_phi.fitTo(phi1d_histo,SumW2Error(true));
//RooFitResult *rphi_1d = Pdf_dN_dOmega.fitTo(phi1d_histo, Range(0,2*3.14));
RooPlot *pframe = p.frame();
phi1d_histo.plotOn(pframe);
//phi_cst.plotOn(pframe);
Pdf_dN_dOmega.plotOn(pframe);
//Pdf_phi.plotOn(pframe);


RooPlot *tframe = t.frame();
//phi_cst.plotOn(tframe);
cst1d_histo.plotOn(tframe);
Pdf_dN_dOmega.plotOn(tframe);

//Pdf_theta.plotOn(tframe);


TCanvas *newcanvas = new TCanvas("newcanvas","mycanvas",1800,800);
newcanvas->Divide(2,1);
newcanvas->cd(1);
//tframe->Draw();
Cst_1D->Draw();
Cst_1D->GetYaxis()->SetRangeUser(0,10000);
Cst_1D->Fit("pdftheta_root","","IL");
//Cst_1D->Fit("wtheta","","L");
//pdftheta_root->Draw("SAME");
//newcanvas.Draw();


newcanvas->cd(2);
Phi_1D->Draw();
Phi_1D->GetYaxis()->SetRangeUser(0,10000);
pdfphi_root->FixParameter(2,pdftheta_root->GetParameter(1));
Phi_1D->Fit("pdfphi_root","","IL");
//Phi_1D->Fit("wphi");
//pdfphi_root->Draw();
//pframe->Draw();




}
