#include "Fitter.h"
#include "Correlation.h"
void runAna(){
  //gROOT->ProcessLine(".L Correlation.cxx");
  //gROOT->ProcessLine(".L Fitter.cxx");
  Fitter myfit1; //= new Fitter();
  //Fitter myfit2;
  Fitter myfit3;
  Correlation runcorrelation1;
  Correlation runcorrelation2;
  Correlation runcorrelation3;

  //TH2D* hist2 = new TH2D("hist2", "pt vs m",10,0.,10.,10,-1.,1.);
  TH1D* mass = new TH1D("mass","mass",100,0.,10.);

  TH1D* mpar = new TH1D ("mpar", "mpar",10 ,-1 ,1);


  //runcorrelation1.ReadFile("../datafiles/incohjpsitomumu.dat");
  //runcorrelation2.ReadFile("../datafiles/starlightgammatomubig.dat");
  runcorrelation2.ReadFile("../datafiles/incohjpsitomumu.dat");
  //runcorrelation3.ReadFile("../datafiles/incohjpsi_gammagamma_to_mumu.dat");

  int n_bins = 10;
  float max_mass = *max_element(runcorrelation2.DiMuMass.begin(),runcorrelation2.DiMuMass.end())-12;
  float min_mass = *min_element(runcorrelation2.DiMuMass.begin(),runcorrelation2.DiMuMass.end());
  float max_pt = *max_element(runcorrelation2.DiMuPT.begin(),runcorrelation2.DiMuPT.end())-13;
  float min_pt = *min_element(runcorrelation2.DiMuPT.begin(),runcorrelation2.DiMuPT.end());

  TH2D* hist = new TH2D("hist", "mass vs m1",n_bins,min_mass,max_mass,n_bins,-1.,1.);
  TH2D* hist2 = new TH2D("hist2", "pt vs m",n_bins,min_pt,max_pt,n_bins,-1.,1.);
  TH1D* pt = new TH1D("pt","pt",n_bins,min_pt,1);

  float width = (max_mass -min_mass)/(n_bins);
  float widthpt = (max_pt -min_pt)/(n_bins);
//  cout << "min_mass is "<< max_mass << endl;
  for (int i=0 ;i <10 ; i++){

    Fitter* myfit2 = new Fitter();

  //  myfit2->Plot();
    float min =min_mass+(i)* width;
    float max =min_mass+(i+1)*width;
    float minpt =min_pt+(i)* widthpt;
    float maxpt =min_pt+(i+1)*widthpt;
    //cout << i << ","<< myfit2->m6 << "hist" << endl;
    float cent = (min+max) *0.5;
    float centpt = (minpt+maxpt) *0.5;
    vector <float> ct;
    vector <float> ph;
  //  cout << ct.size() <<"," << ph.size() <<"hist"<< endl;
    for(int j=0 ; j < runcorrelation2.DiMuMass.size(); j++ ){

      //if (runcorrelation2.DiMuMass[j]>min && runcorrelation2.DiMuMass[j]<max){
      if (runcorrelation2.DiMuPT[j]>minpt && runcorrelation2.DiMuPT[j]<maxpt){
        //cout << runcorrelation2.DiMuMass[j] << endl;
        ct.push_back(runcorrelation2.costheta[j]);
        ph.push_back(runcorrelation2.phi[j]);
      }


    }
    // cout << ct.size() <<"," << ph.size() <<"hist"<< runcorrelation2.costheta.size()<<endl;

     myfit2->runFit(ct,ph);
     myfit2->Plot();
     //myfit2->Plot();
     mpar->Fill(myfit2->m6);

     //cout << "hist" << myfit2->m6<< "2"<< endl;

     //hist->Fill(cent,myfit2->m1);
     hist2->Fill(cent,myfit2->m1);
     //myfit2 = new Fitter();
  }
  for(int j=0 ; j < runcorrelation2.DiMuMass.size(); j++ ){
     mass->Fill(runcorrelation2.DiMuMass[j]);
     pt->Fill(runcorrelation2.DiMuPT[j]);

  }



//  myfit1.runFit(runcorrelation1.costheta,runcorrelation1.phi);

  //myfit3.runFit(runcorrelation3.costheta,runcorrelation3.phi);
  //runcorrelation.AngleCalculator(runcorrelation.d1,runcorrelation.d2);
  //cout<< myfit1.m2 <<", "<< myfit1.m6 <<", " << myfit1.m3<< endl;
  //cout<< myfit2->m2 <<", "<< myfit2->m6 <<", " << myfit2->m3<< endl;
  //cout<< myfit3.m2 <<", "<< myfit3.m6 <<", " << myfit3.m3<< endl;
  //myfit2->Plot();
 //mpar->Draw();
 hist2->Draw("e");

}
