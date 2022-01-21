  TFile *myFile;
  TTree *myTree;


// tree branches
  Float_t MuPt1;
  Float_t MuPt2;
  Float_t DiMuM;
  Float_t DiMuPt;
  Float_t DiMuY;
  Float_t cosThetaCS;
  Float_t SinThetaCosPhiCS;
  Float_t SinThetaSinPhiCS;
  Float_t accoplCut;
  Float_t Phi;
  Float_t Theta;

  //TLorentzVector Pos1, Neg1, Sum;
  //Float_t muMass = 0.1566;
  Float_t ebeam=3500*82/208;// not used
  Float_t mp=0.93827231; // not used
  Float_t mmuon=0.105658; //mass of muon
  Float_t pbeam=TMath::Sqrt(ebeam*ebeam + mp*mp);
