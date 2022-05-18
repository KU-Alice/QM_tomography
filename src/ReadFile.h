#ifndef ReadFile_h
#define ReadFile_h
#include "TLorentzVector.h"
#include "TMath.h"
#include "AngleCalc.h"
#include "Plotter.h"
#include "iostream"

class ReadFile{

public:
  vector <TLorentzVector> daughter1;
  vector <TLorentzVector> daughter2;
  vector <TLorentzVector> daughter_pos;
  vector <TLorentzVector> daughter_neg;
  vector <TLorentzVector> parent;
  void ReadFileAsci(TString inFileName);
  void ReadFileTreeData(TString inFileNameRoot,bool ismc);
  void ReadFileTreeGen(TString inFileNameRoot);

private:

  Float_t mmuon=0.105658; //mass of muon




};
#endif


void ReadFile::ReadFileTreeGen(TString inFileNameRoot){
  TFile* file = TFile::Open(inFileNameRoot);
  TTree* tree;
  //return;
  cout <<"is mc gen"<<endl;
  tree = (TTree*) file->Get("AnalysisOutput/fTreeJPsiMCGen");
  //tree->Print();
  //break;

  //tree->Print();
  int entries = tree->GetEntries();
  double fQ1Gen;
  double fQ2Gen;
  double fMGen;
  double fPtGen;
  double fYGen;
  double fPt1Gen;
  double fEta1Gen;
  double fPhi1Gen;
  double fPt2Gen;
  double fEta2Gen;
  double fPhi2Gen;

  tree->SetBranchAddress("fMGen",&fMGen);
  tree->SetBranchAddress("fPtGen",&fPtGen);
  tree->SetBranchAddress("fYGen",&fYGen);
  tree->SetBranchAddress("fPt2Gen",&fPt2Gen);
  tree->SetBranchAddress("fEta2Gen",&fEta2Gen);
  tree->SetBranchAddress("fPhi2Gen",&fPhi2Gen);
  tree->SetBranchAddress("fPt1Gen",&fPt1Gen);
  tree->SetBranchAddress("fEta1Gen",&fEta1Gen);
  tree->SetBranchAddress("fPhi1Gen",&fPhi1Gen);
  tree->SetBranchAddress("fQ1Gen",&fQ1Gen);
  tree->SetBranchAddress("fQ2Gen",&fQ2Gen);

  for(int i =0; i<entries ; i++){
    tree->GetEntry(i);

    if ((fPtGen<0.2)) continue;
    if ((fPtGen>1)) continue;
    TLorentzVector dneg;
    TLorentzVector dpos;
    TLorentzVector d1;
    TLorentzVector d2;

    d1.SetPtEtaPhiM(fPt1Gen,fEta1Gen,fPhi1Gen,mmuon);
    d2.SetPtEtaPhiM(fPt2Gen,fEta2Gen,fPhi2Gen,mmuon);
    if (fQ1Gen<0){
      dneg.SetPtEtaPhiM(fPt1Gen,fEta1Gen,fPhi1Gen,mmuon);
      dpos.SetPtEtaPhiM(fPt2Gen,fEta2Gen,fPhi2Gen,mmuon);
    }
    else{
      dpos.SetPtEtaPhiM(fPt1Gen,fEta1Gen,fPhi1Gen,mmuon);
      dneg.SetPtEtaPhiM(fPt2Gen,fEta2Gen,fPhi2Gen,mmuon);
    }

    TLorentzVector p = d1 +d2;
    daughter1.push_back(d1);
    daughter2.push_back(d2);
    parent.push_back(p);
    daughter_neg.push_back(dneg);
    daughter_pos.push_back(dpos);

    }

  }


//Read data
void ReadFile::ReadFileTreeData(TString inFileNameRoot,bool ismc =0){
  TFile* file = TFile::Open(inFileNameRoot);
  //file->cd("AnalysisOutput");

  TTree* tree;


  if(ismc){
    cout <<"is mc rec"<<endl;
    tree = (TTree*) file->Get("AnalysisOutput/fTreeJPsiMCRec");
  }

  else{
    // reading tree file is inside directrory AnalysisOutput/Tree(fTreeJpsi)
    cout <<"is data"<<endl;
    tree = (TTree*) file->Get("AnalysisOutput/fTreeJPsi");
  }


  int entries = tree->GetEntries();
  double fM;
  double fPt;
  double fY;
  double fPt1;
  double fEta1;
  double fPhi1;
  double fPt2;
  double fEta2;
  double fPhi2;
  double fQ1;
  double fQ2;
  double fTrk1SigIfMu;
  double fTrk2SigIfMu;
  double fTrk1SigIfEl;
  double fTrk2SigIfEl;
  // detector info
  int fV0A_dec;
  int fV0C_dec;
  int fADA_dec;
  int fADC_dec;
  bool fMatchingSPD;
  int fRunNumber;



  tree->SetBranchAddress("fM",&fM);
  tree->SetBranchAddress("fPt",&fPt);
  tree->SetBranchAddress("fY",&fY);
  tree->SetBranchAddress("fPt2",&fPt2);
  tree->SetBranchAddress("fEta2",&fEta2);
  tree->SetBranchAddress("fPhi2",&fPhi2);
  tree->SetBranchAddress("fPt1",&fPt1);
  tree->SetBranchAddress("fEta1",&fEta1);
  tree->SetBranchAddress("fPhi1",&fPhi1);
  tree->SetBranchAddress("fRunNumber",&fRunNumber);
  tree->SetBranchAddress("fTrk1SigIfMu",&fTrk1SigIfMu);
  tree->SetBranchAddress("fTrk2SigIfMu",&fTrk2SigIfMu);
  tree->SetBranchAddress("fTrk1SigIfEl",&fTrk1SigIfEl);
  tree->SetBranchAddress("fTrk2SigIfEl",&fTrk2SigIfEl);
  tree->SetBranchAddress("fV0A_dec",&fV0A_dec);
  tree->SetBranchAddress("fV0C_dec",&fV0C_dec);
  tree->SetBranchAddress("fADA_dec",&fADA_dec);
  tree->SetBranchAddress("fADC_dec",&fADC_dec);
  tree->SetBranchAddress("fMatchingSPD",&fMatchingSPD);
  tree->SetBranchAddress("fQ1",&fQ1);
  tree->SetBranchAddress("fQ2",&fQ2);

  for(int i =0; i<entries ; i++){
    tree->GetEntry(i);



    if (fMatchingSPD==0) continue;

    if (fADA_dec!=0) continue;
    if (fADC_dec!=0) continue;
    if (fV0A_dec!=0) continue;
    if (fV0C_dec!=0) continue;

    if ((fPt<0.2)) continue;
    if ((fPt>1)) continue;

    if ((TMath::Abs(fEta1))>0.8) continue;
    if ((TMath::Abs(fEta2))>0.8) continue;
    if ((TMath::Abs(fY)>0.8)) continue;
    if (fQ1*fQ2>=0) continue;


    if (TMath::Sqrt(fTrk1SigIfMu * fTrk1SigIfMu + fTrk2SigIfMu *fTrk2SigIfMu) > TMath::Sqrt(fTrk1SigIfEl * fTrk1SigIfEl + fTrk2SigIfEl *fTrk2SigIfEl)) continue;






    TLorentzVector dneg;
    TLorentzVector dpos;
    TLorentzVector d1;
    TLorentzVector d2;

    d1.SetPtEtaPhiM(fPt1,fEta1,fPhi1,mmuon);
    d2.SetPtEtaPhiM(fPt2,fEta2,fPhi2,mmuon);
    if (fQ1<0){
      dneg.SetPtEtaPhiM(fPt1,fEta1,fPhi1,mmuon);
      dpos.SetPtEtaPhiM(fPt2,fEta2,fPhi2,mmuon);
    }

    else{
      dpos.SetPtEtaPhiM(fPt1,fEta1,fPhi1,mmuon);
      dneg.SetPtEtaPhiM(fPt2,fEta2,fPhi2,mmuon);
    }



    TLorentzVector p = d1 +d2;

    daughter1.push_back(d1);
    daughter2.push_back(d2);
    parent.push_back(p);
    daughter_neg.push_back(dneg);
    daughter_pos.push_back(dpos);

  }

}


// Read Rec
