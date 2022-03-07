#include "Fitter.h"
#include "Correlation.h"
#include "generator.h"

void traditional_fit_data(int choose = 2){


  Fitter ang;
  Correlation cm;
  Correlation cm_mc_rec;
  Correlation cm_mc_gen;
  cm.ReadTree("../../david_data/download.root",0,0);

  cm.Plot1D(cm.costheta);
  //cm.SaveTree("AngularDist_Data.root");

  cm_mc_rec.ReadTree("../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1,0);

  //cm.SaveTree("AngularDist_Data.root");
  cm_mc_gen.ReadTree("../../david_data/AnalysisResults_MC_kIncohJpsiToMu_Michal.root",1,1);


  //cout<< "saved dimumass is "<<cm.DiMuMass[1000]<<endl;
    /*TH1F *masshisto = new TH1F("masshisto",";mass GeV/C ; Events/0.13 ",300, 2, 6);
  for (int i =0 ; i < cm_mc_rec.DiMuMass.size(); i++){
    masshisto->Fill(cm_mc_rec.DiMuMass[i]);

  }
  masshisto->Draw();*/


  int n_bins = 4;

  float max_costheta = *max_element(cm.costheta.begin(),cm.costheta.end());
  float min_costheta = *min_element(cm.costheta.begin(),cm.costheta.end());
  float max_phi = *max_element(cm.phi.begin(),cm.phi.end());
  float min_phi = *min_element(cm.phi.begin(),cm.phi.end());
  float width_costheta = (max_costheta -min_costheta)/(n_bins);
  TCanvas* costcanvas0 = new TCanvas("costcanvas0","",600,600);
  //costcanvas->Divide(2,3);
  TH1F* Unfolded_Costheta = new TH1F("Unfolded_Costheta",";costheta;",n_bins,min_costheta,max_costheta);
  TH1F* Unfolded_Phi = new TH1F("Unfolded_Phi",";phi;",n_bins,min_phi,max_phi);

  TH1F* Unfolded_Costheta_Gen = new TH1F("Unfolded_Costheta_Gen",";costheta_{Gen};",n_bins,min_costheta,max_costheta);
  TH1F* Unfolded_Phi_Gen = new TH1F("Unfolded_Phi_Gen",";phi_{Gen};",n_bins,min_phi,max_phi);


  TH1F* Unfolded_Costheta_Rec = new TH1F("Unfolded_Costheta_Rec",";costheta_{Rec};",n_bins,min_costheta,max_costheta);
  TH1F* Unfolded_Phi_Rec = new TH1F("Unfolded_Phi_Rec",";phi_{Rec};",n_bins,min_phi,max_phi);
//  cout << "#######################"<<cm_mc_gen.phi.size() << endl;
  for (int i =0 ;i<cm_mc_gen.phi.size(); i++){
  //  cout <<"$$$$$$$$$$$$$$$$$$$$$$$$$$" <<cm_mc_gen.phi[i] << endl;
    Unfolded_Phi_Gen->Fill(cm_mc_gen.phi[i]);
    Unfolded_Costheta_Gen->Fill(cm_mc_gen.costheta[i]);
  }

  for (int i =0 ;i<cm_mc_rec.phi.size(); i++){
  //  cout <<"$$$$$$$$$$$$$$$$$$$$$$$$$$" <<cm_mc_gen.phi[i] << endl;
    Unfolded_Phi_Rec->Fill(cm_mc_rec.phi[i]);
    Unfolded_Costheta_Rec->Fill(cm_mc_rec.costheta[i]);
  }



  float widthphi = (max_phi -min_phi)/(n_bins);



  for (int i=0 ;i <n_bins; i++){

      Fitter* mix = new Fitter();
      Fitter* mix2 = new Fitter();




  //  myfit2->Plot();




    if(choose==1){
      float min_cost =min_costheta+(i)* width_costheta;
      float max_cost =min_costheta+(i+1)*width_costheta;
      mix->traditional_fit(cm.DiMuMass,cm_mc_rec.DiMuMass,cm.costheta,cm_mc_rec.costheta,min_cost,max_cost,i);
      float centcost = (min_cost+max_cost) *0.5;

      Unfolded_Costheta->Fill(centcost,mix->yield_sig);
      Unfolded_Costheta->SetBinError(i+1,mix->yield_sig_err);
    }
    else{
      float minph =min_phi+(i)* widthphi;
      float maxph =min_phi+(i+1)*widthphi;
      float centph = (minph+maxph) *0.5;
      mix2->traditional_fit(cm.DiMuMass,cm_mc_rec.DiMuMass,cm.phi,cm_mc_rec.phi,minph,maxph, i);
      Unfolded_Phi->Fill(centph,mix2->yield_sig);
      Unfolded_Phi->SetBinError(i+1,mix2->yield_sig_err);
    }

    //mix2->traditional_fit(cm_mc_gen.DiMuMass,cm_mc_rec.DiMuMass,cm.phi,cm_mc_rec.phi,minph,maxph, i);






  }
  Unfolded_Phi_Gen->Sumw2();
  Unfolded_Phi_Rec->Sumw2();
  Unfolded_Phi->Sumw2();
  costcanvas0->Divide(2,2);
  if(choose==1){
    Unfolded_Costheta_Gen->Sumw2();
    Unfolded_Costheta_Rec->Sumw2();
    Unfolded_Costheta->Sumw2();
    costcanvas0->cd(1);
    Unfolded_Costheta->Draw("e");
    costcanvas0->cd(2);
    Unfolded_Costheta_Gen->Draw("e");
    TH1D * axecostheta = ( TH1D* ) Unfolded_Costheta_Rec->Clone("axecostheta");
    costcanvas0->cd(3);
    axecostheta->Draw("e");
    TH1D * correctedcostheta = ( TH1D* ) Unfolded_Costheta_Rec->Clone("correctedcostheta");
    axecostheta->Divide(Unfolded_Costheta_Rec,Unfolded_Costheta_Gen);
    correctedcostheta->Divide(Unfolded_Costheta,axecostheta);

    costcanvas0->cd(4);
    correctedcostheta->Draw("e");




  }
  else{

    costcanvas0->cd(1);
    Unfolded_Phi->Draw("e");

    costcanvas0->cd(2);
    Unfolded_Phi_Gen->Draw("e");
    Unfolded_Phi_Gen->Sumw2();

    TH1D * axephi =( TH1D* )Unfolded_Phi->Clone("axephi");
    axephi->Divide(Unfolded_Phi_Rec,Unfolded_Phi_Gen);

    costcanvas0->cd(3);
    axephi->Draw("e");

    TH1D * correctedphi =( TH1D* )Unfolded_Phi->Clone("correctedphi");

    correctedphi->Divide(Unfolded_Phi,axephi);

    costcanvas0->cd(4);
    correctedphi->Draw("e");

  }

















}
