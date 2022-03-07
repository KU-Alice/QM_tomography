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
using namespace RooFit;


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
  vector <float> GetParametersPdf();


  void FitJPsiMass(Int_t iPoint = -1,Float_t maxPt = 0.2, Float_t minPt = 0.0,
  		Float_t minPhi = 0.0, Float_t maxPhi =  6.29,
  		Float_t minCsT = -1.0, Float_t maxCsT =  1.0,
  		Int_t nBins = 116, Float_t fitMin = 2.5, Float_t fitMax = 4.5);

  void setmassthetaphilist(vector <float> mass, vector <float> costheta,vector <float> phi ){
    rmass = mass;
    rcostheta = costheta;
    rphi = phi;
    //cout << rmass.size()<< "##################################"<<endl;

  }
  void setmcmassthetaphilist(vector <float> mcmass, vector <float> mccostheta,vector <float> mcphi ){
    rmcmass = mcmass;
    rmccostheta = mccostheta;
    rmcphi = mcphi;
  }

  //TF1 setsignalfunction(TF1 myfunc) {return mysigfunc};
  //TF1 setbkgfucntion(TF1 mybkgfunc){ return mybkgfunc};



private:

  vector <float> rmass, rcostheta, rphi;
  vector <float> rmcmass, rmccostheta, rmcphi;


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

  TH1F *masshisto;
  TH1F *masshisto_mc;

  TF1* expo_tail;
  TF1* sig_func_jpsi;
  TF1* sig_func_jpsi2;
  TF1* comb_func;


  TF1* funccos;
  TF1* funcphi;



  RooGenericPdf Pdf_theta= {"Pdf_theta", "Pdf_theta", "(1.5708*(1+ mt3*mt3) + 1.3823*(mt2*mt5+mt7*mt8-mt6*mt9)*t + (1.5708-4.71239*mt3*mt3)*t*t)", RooArgSet(t, mt2,mt3,mt5,mt6,mt7,mt8,mt9)};


  RooGenericPdf Pdf_phi= {"Pdf_phi", "Pdf_phi", "(0.66666667 + 0.345575*mt3*mt9*cos(p)+0.33333333*(-1+2*mt2*mt2+mt3*mt3+2*mt8*mt8+2*mt9*mt9)*cos(2*p)-0.345575*mt3*mt7*sin(p)-0.66666667*(mt2*mt4+mt6*mt8+mt7*mt8)*sin(2*p))", RooArgSet(p, mt2,mt3,mt4,mt6,mt7,mt8,mt9)};

  RooGenericPdf Pdf_dN_dOmega= {"Pdf_dN_dOmega", "Pdf_dN_dOmega", "(0.25*(1+mt3*mt3) +0.25*(1-3*mt3*mt3)*t*t -0.5* mt3* mt6* 2* TMath::Sqrt(1- t*t)* t* cos(p)+0.25*(2*mt2*mt2+mt3*mt3-1)*(1-t*t)*cos(2*p))", RooArgSet(p,t, mt2,mt3,mt6)};
  //RooGenericPdf Pdf_dN_dOmega= {"Pdf_dN_dOmega", "Pdf_dN_dOmega", "(0.25*(1+mt3*mt3) +0.25*(1-3*mt3*mt3)*t*t -0.5* mt3* TMath::Sqrt(1-mt2*mt2-mt3*mt3)* 2* TMath::Sqrt(1- t*t)* t* cos(p)+0.25*(2*mt2*mt2+mt3*mt3-1)*(1-t*t)*cos(2*p))", RooArgSet(p,t, mt2,mt3)};

  virtual void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07,int columns=2);
  virtual void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15);







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


/*void Fitter::runFit2(vector <float> costheta, vector <float> phi ){

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

}*/

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


  masshisto = new TH1F("masshisto",";mass GeV/C ; Events ",100, 2, 6);
  masshisto_mc = new TH1F("masshisto_mc","",100, 2, 6);









    expo_tail = new TF1("expo_tail", "[0]*TMath::Exp(-[1]*x)",2,6);
    expo_tail->SetParLimits(1,0.0001,1);

    //TF1* sig_func_jpsi =  new TF1("sig_func_jpsi","gaus",2,4);

    sig_func_jpsi = new  TF1("sig_func_jpsi","crystalball",0,5);
    sig_func_jpsi->SetParameters(100,3.15,0.090,0.8,115);
    //sig_func_jpsi->SetParameters(1000,3.0969,0.09);
    sig_func_jpsi->SetParLimits(1,3.05,3.17);

    sig_func_jpsi->SetLineColor(kBlack);

    sig_func_jpsi2 = new  TF1("sig_func_jpsi2","crystalball",0,5);

    sig_func_jpsi2->SetParameters(100,3.15,0.090,0.8,115);
    sig_func_jpsi2->SetParLimits(1,3.08,3.13);

    comb_func = new TF1("comb_func","[0]*TMath::Exp(-[1]*x)+crystalball(2)",2,5);
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


    masshisto_mc->Fit("sig_func_jpsi","QR","",2.5,4);



    TString cname;
    cname.Form("%f< angle ",ang_min);
    TCanvas* dat = new TCanvas(cname.Data(),"",400,300);
    dat->cd();
    //mycanvas->cd(ncanvas);
    masshisto->Draw("e");
    expo_tail->SetLineColor(kGreen);

    //masshisto_corrected->Fit("expo_tail","R+Q","",4,6);
    sig_func_jpsi2->FixParameter(1,sig_func_jpsi->GetParameter(1));
    sig_func_jpsi2->FixParameter(2,sig_func_jpsi->GetParameter(2));
    sig_func_jpsi2->FixParameter(3,sig_func_jpsi->GetParameter(3));
    sig_func_jpsi2->FixParameter(4,sig_func_jpsi->GetParameter(4));

    //masshisto->Fit("sig_func_jpsi");
    //comb_func->SetParameter(1,expo_tail->GetParameter(1));
    //comb_func->FixParameter(2,expo_tail->GetParameter(2));
    sig_func_jpsi2->SetLineColor(kBlue);
    masshisto->Fit("sig_func_jpsi2","NQR","",2.95,3.2);

    comb_func->FixParameter(3,sig_func_jpsi->GetParameter(1));
    comb_func->FixParameter(4,sig_func_jpsi->GetParameter(2));
    comb_func->FixParameter(5,sig_func_jpsi2->GetParameter(3));
    comb_func->FixParameter(6,sig_func_jpsi2->GetParameter(4));
    masshisto->Fit("comb_func","QR+","",2.5,3.4);


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
  YeildSig.Form("N_{j/#Psi} = %f #pm %f",yield_sig,yield_sig_err);


  TString YeildBkg;
  YeildBkg.Form("N_{Bkg} = %f #pm %f",bkg,yield_sig_err);
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

vector <float> Fitter::GetParametersPdf(){

  vector <float> mparameters;
  mparameters.push_back(m1);
  mparameters.push_back(m2);
  mparameters.push_back(m3);
  mparameters.push_back(m4);
  mparameters.push_back(m5);
  mparameters.push_back(m6);
  mparameters.push_back(m7);
  mparameters.push_back(m8);

  return mparameters;



}

void Fitter::FitJPsiMass(Int_t iPoint = -1,Float_t maxPt = 0.2, Float_t minPt = 0.0,
		Float_t minPhi = 0.0, Float_t maxPhi =  6.29,
		Float_t minCsT = -1.0, Float_t maxCsT =  1.0,
		Int_t nBins = 116, Float_t fitMin = 2.5, Float_t fitMax = 4.5) {



    UInt_t nDT = rmass.size();
cout << nDT << "#############@#@##@@#@@@@@@@@@@@@@@@@@@@@@@@@@@@"<< endl;
RooRealVar MassDT("MassDT","MassDT",fitMin,fitMax);
MassDT.SetTitle("#it{m}_{#mu#mu} (GeV/#it{c}^{2})");
RooDataSet DT_UnBin("DT_UnBin","DT_UnBin",RooArgSet(MassDT));

    //event loop
for(Int_t i=0; i<nDT; i++){
  //DT_Tree->GetEntry(i);

  //if(TMath::Abs(DiMuY)>0.8) continue;
  if(rmass[i] <= fitMin || rmass[i] >= fitMax) continue;
  //if(DiMuPt < minPt || DiMuPt > maxPt) continue;
  if(rcostheta[i] > maxCsT || rcostheta[i] < minCsT)continue;
  if(rphi[i] > maxPhi || rphi[i] < minPhi)continue;

  MassDT = rmass[i];

  DT_UnBin.add(RooArgSet(MassDT));
  }

//________________________________________________MONTE CARLO_________________________________________________________
//TTree *MC_Tree  = (TTree*) MC_file->Get("myTree"));
//MC_Tree->SetName("MC_Tree");
//MC_Tree->SetBranchAddress("DiMuM", &DiMuM);
//MC_Tree->SetBranchAddress("DiMuPt", &DiMuPt);
//MC_Tree->SetBranchAddress("DiMuY", &DiMuY);
//MC_Tree->SetBranchAddress("CosTheta", &CosTheta);
//MC_Tree->SetBranchAddress("SinThetaCosPhi", &SinThetaCosPhi);
//MC_Tree->SetBranchAddress("SinThetaSinPhi", &SinThetaSinPhi);
//MC_Tree->SetBranchAddress("Phi", &Phi);
//MC_Tree->SetBranchAddress("Theta", &Theta);

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
//RooCBShape CrystalBall_JPsi_DT("CrystalBall_JPsi_DT","CrystalBall_JPsi_DT",MassDT,Mean_JPsi,RooConst(Sigma_JPsi.getVal()),RooConst(Alpha_JPsi.getVal()),RooConst(N_JPsi.getVal()));
RooCBShape CrystalBall_JPsi_DT("CrystalBall_JPsi_DT","CrystalBall_JPsi_DT",MassDT,Mean_JPsi,RooConst(Sigma_JPsi.getVal()),Alpha_JPsi,N_JPsi);

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
latex->DrawLatex(0.55,0.77,"#it{p}_{#lower[-0.3]{T}} < 0.2 GeV/#it{c}");
latex->DrawLatex(0.18,0.77,Form("|#it{y}| < 0.8"));
latex->DrawLatex(0.55,0.69,Form("%1.1f#pi < #phi < %1.1f#pi",minPhi/TMath::Pi(),maxPhi/TMath::Pi()));
latex->DrawLatex(0.55,0.61,Form("%1.1f < cos(#theta) < %1.1f",minCsT,maxCsT));

latex->DrawLatex(0.55,0.43,Form("N_{J/#psi} = %.1f #pm %.1f",nSignal_JPsi.getVal(),nSignal_JPsi.getError()));
latex->DrawLatex(0.55,0.35,Form("#chi^{2}/#it{dof} = %.2f",plot_DT->chiSquare("DT_FitFunction_Norm[MassDT]_Range[fit_nll_DT_FitFunction_DT_UnBin]_NormRange[fit_nll_DT_FitFunction_DT_UnBin]", "h_DT_UnBin", 6)));


if(iPoint != -1){
  TFile *yieldFile = new TFile("Yield_JPsi.root","UPDATE");
  TH2F *yieldPhiCsT =  0x0;
  TH1F *yieldPhi =  0x0;
  TH1F *yieldCsT =  0x0;
  yieldPhiCsT = (TH2F*)yieldFile->Get("yieldPhiCsTStore");
  yieldCsT = (TH1F*)yieldFile->Get("yieldCsTStore");
  yieldPhi = (TH1F*)yieldFile->Get("yieldPhiStore");

  if(!yieldPhiCsT){
    cout<<"Creating new histo"<<endl;
    Float_t binsPhi[5] = {0,0.5*TMath::Pi(),TMath::Pi(),1.5*TMath::Pi(),2*TMath::Pi()};
    Float_t binsCsT[4] = {-1,-0.5,0.5,1};
    yieldPhi = new TH1F("yieldPhi","yieldPhi",4,binsPhi);
    yieldCsT = new TH1F("yieldCsT","yieldCsT",3,binsCsT);
    yieldPhiCsT = new TH2F("yieldPhiCsT","yieldPhiCsT",4,binsPhi,3,binsCsT);

    yieldPhiCsT->SetTitle("J/#psi #rightarrow l^{+} l^{-}");
    yieldPhi->SetTitle("J/#psi #rightarrow l^{+} l^{-}");
    yieldCsT->SetTitle("J/#psi #rightarrow l^{+} l^{-}");
    yieldPhiCsT->GetXaxis()->SetTitle("#phi");
    yieldPhiCsT->GetYaxis()->SetTitle("cos(#theta)");

    yieldPhi->GetXaxis()->SetTitle("#phi");
    yieldCsT->GetXaxis()->SetTitle("cos(#theta)");

    }
    yieldPhiCsT->SetName("yieldPhiCsT");
    yieldPhi->SetName("yieldPhi");
    yieldCsT->SetName("yieldCsT");

  yieldPhiCsT->SetBinContent(iPoint+1,nSignal_JPsi.getVal());
  yieldPhiCsT->SetBinError(iPoint+1,nSignal_JPsi.getError());

  yieldPhi->SetBinContent(iPoint+1,nSignal_JPsi.getVal());
  yieldPhi->SetBinError(iPoint+1,nSignal_JPsi.getError());

  yieldCsT->SetBinContent(iPoint+1,nSignal_JPsi.getVal());
  yieldCsT->SetBinError(iPoint+1,nSignal_JPsi.getError());

  yieldFile->Delete("yieldPhiStore;*");
  yieldFile->Delete("yieldCsTStore;*");
  yieldFile->Delete("yieldPhiCsTStore;*");
  yieldPhiCsT->SetName("yieldPhiCsTStore");
  yieldPhi->SetName("yieldPhiStore");
  yieldCsT->SetName("yieldCsTStore");
  yieldPhiCsT->Write();
  yieldPhi->Write();
  yieldCsT->Write();

  c1->SetName(TString::Format("FitDiLepton_%d",iPoint));
  c1->Write();
  yieldFile->Close();
  }

}


void Fitter::myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07,int columns=2){
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

void Fitter::myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}
