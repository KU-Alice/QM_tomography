#ifndef Fitter_h
#define Fitter_h
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
#include <complex>
#include <iostream>

class Fitter{


public:
  //Fitter();
//  ~Fitter();

  // values of m parameters to access.
  Float_t m1;
  Float_t m2;
  Float_t m3;
  Float_t m4;
  Float_t m5;
  Float_t m6;
  Float_t m7;
  Float_t m8;
  Float_t m9;
  Float_t rho[3][3];
  Float_t rho_sig[3][3];
  virtual void runFit(vector <float>costheta , vector <float> phi);
  virtual void rho_calc(Float_t rho_mix[3][3],Float_t rho_bkg[3][3]);
  virtual void Plot();


private:


  RooRealVar t = {"t", "Cos(#theta)", -1., 1.} ;
  RooRealVar p = {"p", "Phi", 0., 2*3.141616} ;
  RooDataSet* data;

  RooRealVar mt1 = {"mt1", "mt1", 0};
  RooRealVar mt2 = {"mt2", "mt2", 0.3,-1,1};
  RooRealVar mt3= {"mt3", "mt3",0.4 -1,1};
  RooRealVar mt4= {"mt4", "mt4",0,-1.0,1.};
  RooRealVar mt5= {"mt5", "mt5",0,-1,1};
  RooRealVar mt6= {"mt6", "mt6",0.5,-1.0,1.};
  RooRealVar mt7= {"mt7", "mt7",0,-1.0,1.};
  RooRealVar mt8= {"mt8", "mt8",0,-1.0,1.};
  RooRealVar mt9= {"mt9", "mt9",0,-1.0,1.};



  RooGenericPdf Pdf_theta= {"Pdf_theta", "Pdf_theta", "(1.5708*(1+ mt3*mt3) + 1.3823*(mt2*mt5+mt7*mt8-mt6*mt9)*t + (1.5708-4.71239*mt3*mt3)*t*t)", RooArgSet(t, mt2,mt3,mt5,mt6,mt7,mt8,mt9)};


  RooGenericPdf Pdf_phi= {"Pdf_phi", "Pdf_phi", "(0.66666667 + 0.345575*mt3*mt9*cos(p)+0.33333333*(-1+2*mt2*mt2+mt3*mt3+2*mt8*mt8+2*mt9*mt9)*cos(2*p)-0.345575*mt3*mt7*sin(p)-0.66666667*(mt2*mt4+mt6*mt8+mt7*mt8)*sin(2*p))", RooArgSet(p, mt2,mt3,mt4,mt6,mt7,mt8,mt9)};

  RooGenericPdf Pdf_dN_dOmega= {"Pdf_dN_dOmega", "Pdf_dN_dOmega", "(0.25*(1+mt3*mt3) +0.25*(1-3*mt3*mt3)*t*t -0.5* mt3* mt6* 2* TMath::Sqrt(1- t*t)* t* cos(p)+0.25*(2*mt2*mt2+mt3*mt3-1)*(1-t*t)*cos(2*p))", RooArgSet(p,t, mt2,mt3,mt6)};




};
#endif






void Fitter::runFit(vector <float> costheta, vector <float> phi ){

  data = new RooDataSet("data", "data", RooArgSet(t,p));
  for (int i =0 ; i<costheta.size();i++){
    t = costheta[i];
    p = phi[i];
    data->add(RooArgSet(t,p));
  }
  //data = data1->Clone();
  RooFitResult *rtot = Pdf_dN_dOmega.fitTo(*data);

  m1 = mt1.getVal();
  m2 = mt2.getVal();
  m3 = mt3.getVal();
  m4 = mt4.getVal();
  m5 = mt5.getVal();
  m6 = mt6.getVal();
  m7 = mt7.getVal();
  m8 = mt8.getVal();
  m9 = mt9.getVal();

  rho[0][0] = mt6.getVal() * mt6.getVal();
  rho[0][1] = 0;
  rho[0][2] = mt3.getVal()*mt6.getVal();
  rho[1][0] = 0;
  rho[1][1] = mt2.getVal()*mt2.getVal();
  rho[1][2] = 0;
  rho[2][0] = mt3.getVal()*mt6.getVal();
  rho[2][1] = 0;
  rho[2][2] = mt3.getVal()*mt3.getVal();

}

void Fitter::rho_calc(Float_t rho_mix[3][3],Float_t rho_bkg[3][3]){

  int rows=3;
  int cols=3;

  for(int i = 0; i < rows; ++i){
              for(int j = 0; j < cols; ++j)
              {
                  rho_mix[i][j]=0;
                  rho_bkg[i][j]=0;
              }
          }

          for(int i = 0; i < rows; ++i){
              for(int j = 0; j < cols; ++j){
                  // std::cout<<"rho_bkg["<<i<<"]["<<j<<"] ="<<rho_bkg[i][j]<<std::endl;
                  // std::cout<<"rho_tot["<<i<<"]["<<j<<"] ="<<rho_tot[i][j]<<std::endl;
                  for(int k = 0; k < cols; ++k)
                  {
                      int Factor = 0;
                      if (i==k){
                        Factor = 1;
                      }
                      rho_sig[i][j] += (1* Factor - rho_bkg[i][k]) * rho_mix[k][j];

                  }
                  //std::cout<<i<<"--"<<j<<"--"<<rho_final[i][j]<<"------------ "<< rho_sig[i][j]<<std::endl;

              }

          }

}

void Fitter::Plot(){

    RooPlot *tframe = t.frame();
    data->plotOn(tframe);
    Pdf_dN_dOmega.plotOn(tframe,RooFit::LineColor(kRed));


    RooPlot *pframe = p.frame();
    data->plotOn(pframe);
    Pdf_dN_dOmega.plotOn(pframe,RooFit::LineColor(kRed));

    TCanvas *c = new TCanvas("c","c",300,400);
    c->Divide(1,2);
    c->cd(1);
    tframe->Draw();
    c->cd(2);
    pframe->Draw();


}
