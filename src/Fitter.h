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
#include "TLegend.h"
#include <complex>
#include <iostream>


class Fitter{


public:
  //Fitter();
//  ~Fitter();

  // values of m parameters to access.
  Float_t m1;
  Float_t m1err;
  Float_t m2;
  Float_t m2err;
  Float_t m3;
  Float_t m3err;
  Float_t m4;
  Float_t m4err;
  Float_t m5;
  Float_t m5err;
  Float_t m6;
  Float_t m6err;
  Float_t m7;
  Float_t m7err;
  Float_t m8;
  Float_t m8err;
  Float_t m9;
  Float_t m9err;
  Float_t rho[3][3];
  Float_t rho_sig[3][3];
  float yield_bkg;
  float yield_sig;
  float yield_sig_err;
  float yield_bkg_err;
  virtual void runFit(vector <float>costheta , vector <float> phi);
  virtual void runFitData(RooDataSet* data);
  virtual void rho_calc(Float_t rho_mix[3][3],Float_t rho_bkg[3][3]);

  virtual void Plot();
  virtual void traditional_fit(vector <float> mass,vector <float> mass_mc,vector <float> angle,vector <float> angle_mc ,float ang_min,float ang_max,int ncanvas);

  //TF1 setsignalfunction(TF1 myfunc) {return mysigfunc};
  //TF1 setbkgfucntion(TF1 mybkgfunc){ return mybkgfunc};



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
  //RooGenericPdf Pdf_dN_dOmega= {"Pdf_dN_dOmega", "Pdf_dN_dOmega", "(0.25*(1+mt3*mt3) +0.25*(1-3*mt3*mt3)*t*t -0.5* mt3* TMath::Sqrt(1-mt2*mt2-mt3*mt3)* 2* TMath::Sqrt(1- t*t)* t* cos(p)+0.25*(2*mt2*mt2+mt3*mt3-1)*(1-t*t)*cos(2*p))", RooArgSet(p,t, mt2,mt3)};








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

  m1err= mt1.getError();
  m2err= mt2.getError();
  m3err= mt3.getError();
  m4err= mt4.getError();
  m5err= mt5.getError();
  m6err= mt6.getError();
  m7err= mt7.getError();
  m8err= mt8.getError();
  m9err= mt9.getError();


  rho[0][0] = m6 * m6;
  rho[0][1] = 0;
  rho[0][2] = m3*m6;
  rho[1][0] = 0;
  rho[1][1] = m2*m2;
  rho[1][2] = 0;
  rho[2][0] = m3*m6;
  rho[2][1] = 0;
  rho[2][2] = m3*m3;

}

void Fitter::runFitData(RooDataSet*data ){

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

  m1err= mt1.getError();
  m2err= mt2.getError();
  m3err= mt3.getError();
  m4err= mt4.getError();
  m5err= mt5.getError();
  m6err= mt6.getError();
  m7err= mt7.getError();
  m8err= mt8.getError();
  m9err= mt9.getError();


  rho[0][0] = m6 * m6;
  rho[0][1] = 0;
  rho[0][2] = m3*m6;
  rho[1][0] = 0;
  rho[1][1] = m2*m2;
  rho[1][2] = 0;
  rho[2][0] = m3*m6;
  rho[2][1] = 0;
  rho[2][2] = m3*m3;

}

void Fitter::rho_calc(Float_t rho_mix[3][3],Float_t rho_bkg[3][3]){

  int rows=3;
  int cols=3;

  for(int i = 0; i < rows; ++i){
              for(int j = 0; j < cols; ++j)
              {
                  rho_sig[i][j]=0;
                  //rho_bkg[i][j]=0;
              }
          }

          for(int i = 0; i < rows; ++i){
              for(int j = 0; j < cols; ++j){
                  // std::cout<<"rho_bkg["<<i<<"]["<<j<<"] ="<<rho_bkg[i][j]<<std::endl;
                  // std::cout<<"rho_tot["<<i<<"]["<<j<<"] ="<<rho_tot[i][j]<<std::endl;
                  for(int k = 0; k < cols; ++k)
                  {
                      int Factor = 0;RooRealVar mt3= {"mt3", "mt3",0.4 -1,1};

  RooRealVar mt6= {"mt6", "mt6",0.5,-1.0,1.};

                      if (i==k){
                        Factor = 1;
                      }

                      rho_sig[i][j] += (1* Factor - rho_bkg[i][k]) * rho_mix[k][j];

                  }
                  std::cout<<i<<"--"<<j<<"--"<<rho_sig[i][j]<<std::endl;
                  //std::cout<<i<<"--"<<j<<"--"<<rho_sig[i][j]<<std::endl;

              }RooRealVar mt3= {"mt3", "mt3",0.4 -1,1};

  RooRealVar mt6= {"mt6", "mt6",0.5,-1.0,1.};


          }

}









void Fitter::traditional_fit(vector <float> mass,vector <float> mass_mc,vector <float> angle,vector <float> angle_mc ,float ang_min,float ang_max,int ncanvas){

    float mass_max = *max_element(mass.begin(),mass.end());
    float mass_min = *min_element(mass.begin(),mass.end());


    TH1F *masshisto = new TH1F("masshisto",";mass GeV/C ; Events ",100, 2, 6);
    TH1F *masshisto_mc = new TH1F("masshisto_mc","",200, 2, 6);









    TF1* expo_tail = new TF1("expo_tail", "[0]*TMath::Exp(-[1]*x)",2,6);
    expo_tail->SetParLimits(1,0.0001,1);

    //TF1* sig_func_jpsi =  new TF1("sig_func_jpsi","gaus",2,4);

    TF1* sig_func_jpsi = new  TF1("sig_func_jpsi","crystalball",0,5);
    sig_func_jpsi->SetParameters(100,3.15,0.090,0.8,115);
    //sig_func_jpsi->SetParameters(1000,3.0969,0.09);
    sig_func_jpsi->SetParLimits(1,3.05,3.17);

    sig_func_jpsi->SetLineColor(kBlack);

    TF1* sig_func_jpsi2 = new  TF1("sig_func_jpsi2","crystalball",0,5);
    sig_func_jpsi2->SetParameters(100,3.15,0.090,0.8,115);
    sig_func_jpsi2->SetParLimits(1,3.08,3.13);

    TF1* comb_func = new TF1("comb_func","[0]*TMath::Exp(-[1]*x)+crystalball(2)",2,5);
    comb_func->SetParameters(1,0.001,100,3.15,0.090,0.8,115);
    comb_func->SetLineColor(kRed);
    for(int i =0; i< mass.size();i++){
      if((angle[i]>ang_max) || (angle[i] < ang_min)) continue;
      masshisto ->Fill(mass[i]);


    }

    for(int i =0; i< mass_mc.size();i++){

      if((angle_mc[i]>ang_max) || (angle_mc[i] < ang_min)) continue;
      masshisto_mc->Fill(mass_mc[i]);



    }



    TLegend *legend = new  TLegend(0.45,0.8,0.8,0.4);
    legend->SetFillStyle(0);
    legend->SetHeader("UPC Run2 Dataset","C");
    TString cut_editor;
    cut_editor.Form("%f< angle<%f",ang_min,ang_max);
    legend->AddEntry((TObject*)0,cut_editor.Data(),"");
    legend->SetLineColor(0);
    legend->SetTextSize(0.025);


    TCanvas* mc = new TCanvas("mc","reconstructed mass",600,300);
    mc->cd();



    masshisto_mc->Draw("e");


    masshisto_mc->Fit("sig_func_jpsi","R","",2.5,4);



    TString cname;
    cname.Form("%f< anglef",ang_min);
    TCanvas* dat = new TCanvas(cname.Data(),"",400,300);
    dat->cd();
    //mycanvas->cd(ncanvas);
    masshisto->Draw("e");
    expo_tail->SetLineColor(kGreen);


    //masshisto_corrected->Fit("expo_tail","R+Q","",4,6);
    sig_func_jpsi2->FixParameter(1,sig_func_jpsi->GetParameter(1));
    sig_func_jpsi2->FixParameter(2,sig_func_jpsi->GetParameter(2));
    sig_func_jpsi2->SetParameter(3,sig_func_jpsi->GetParameter(3));
    sig_func_jpsi2->SetParameter(4,sig_func_jpsi->GetParameter(4));

    //masshisto->Fit("sig_func_jpsi");
    //comb_func->SetParameter(1,expo_tail->GetParameter(1));
    //comb_func->FixParameter(2,expo_tail->GetParameter(2));
    sig_func_jpsi2->SetLineColor(kBlue);
    masshisto->Fit("sig_func_jpsi2","R+NQ","",2.95,3.25);

    comb_func->FixParameter(3,sig_func_jpsi2->GetParameter(1));
    comb_func->FixParameter(4,sig_func_jpsi2->GetParameter(2));
    comb_func->FixParameter(5,sig_func_jpsi2->GetParameter(3));
    comb_func->FixParameter(6,sig_func_jpsi2->GetParameter(4));
    masshisto->Fit("comb_func","R+","",2.5,3.4);


  TF1* sig_func_jpsi3 = new  TF1("sig_func_jpsi3","crystalball",0,5);
  sig_func_jpsi3->SetParameters(comb_func->GetParameter(2),comb_func->GetParameter(3),comb_func->GetParameter(4),comb_func->GetParameter(5),comb_func->GetParameter(6));
  sig_func_jpsi3->SetLineColor(kBlack);
  expo_tail->SetParameters(comb_func->GetParameter(0),comb_func->GetParameter(1));
  expo_tail->Draw("SAME");
  sig_func_jpsi3->Draw("SAME");
  float   mu = comb_func->GetParameter(3);
  float  sigma = comb_func->GetParameter(4);
  float min = mu-3*sigma;
  float max = mu +3*sigma;

  float sig = sig_func_jpsi3->Integral(min,max)/masshisto->GetBinWidth(1);
  float bkg = expo_tail->Integral(min,max)/masshisto->GetBinWidth(1);
  cout << "sig is " << sig << endl;
  cout << "bkg is " << bkg << endl;
  cout << "bin " << masshisto->GetBinWidth(1) << endl;
  cout << "bin " << masshisto->GetBinWidth(20) << endl;

  yield_sig = sig;
  yield_sig_err = (yield_sig * comb_func->GetParError(2))/comb_func->GetParameter(2);
  yield_bkg_err = (yield_bkg * comb_func->GetParError(0))/comb_func->GetParameter(0);

  TString YeildSig;
  YeildSig.Form("N_{j/#Psi}%f #pm %f",yield_sig,yield_sig_err);


  TString YeildBkg;
  YeildBkg.Form("N_{Bkg}=%f #pm %f",bkg,yield_sig_err);
  legend->AddEntry((TObject*)0,YeildSig.Data() ,"");
  legend->AddEntry((TObject*)0,YeildBkg.Data(),"");
    //legend->AddEntry(0,"#tilde{#chi}^{2} ="+chisq_ndf_str,"")

  legend->Draw("SAME");

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
