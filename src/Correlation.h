#ifndef Correlation_h
#define Correlation_h

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
class Correlation{

// tree branches


  //TLorentzVector Pos1, Neg1, Sum;
  //Float_t muMass = 0.1566;

public:
  //Correlation();
  //~Correlation();

  //RooRealVar theta();
  //RooRealVar phi();
  //RooDataSet* data ;//  RooDataSet d("d","d",RooArgSet(x,y,c)) ;
  vector <float> costheta;
  vector <float> phi;
  vector <float> DiMuMass;
  vector <float> DiMuPT;

  TLorentzVector d1;
  TLorentzVector d2;
  TLorentzVector p;


  void ReadFile(TString inFileName);
  virtual void SaveTree(const char* rootfile);

  //TFile *myFile;
  //TTree *myTree;
private:

  Float_t ebeam=3500*82/208;// not used
  Float_t mp=0.93827231; // not used
  Float_t mmuon=0.105658; //mass of muon
  Float_t pbeam=TMath::Sqrt(ebeam*ebeam + mp*mp);
  void AngleCalculator(TLorentzVector v1, TLorentzVector v2);


  Float_t DiMuM;
  Float_t DiMuPt;
  Float_t MuPt1;
  Float_t MuPt2;
  Float_t DiMuY;
  Float_t cosThetaCS;
  Float_t SinThetaCosPhiCS;
  Float_t SinThetaSinPhiCS;
  Float_t accoplCut;
  Float_t Phi;
  Float_t Theta;
  //RooRealVar cost = {"cost","costheta", -1,1};
  //RooRealVar ph = {"ph","phi", 2*3.141616};

};
#endif


void Correlation::AngleCalculator(TLorentzVector v4, TLorentzVector v5)

{

  //TLorentzVector v4;
  //TLorentzVector v5;
  TLorentzVector pa(1.,0.,0,1); // projectile
  TLorentzVector pb(1.,0., 0,-1); // target

  Float_t px1 = v4.Px();
  Float_t py1 = v4.Py();
  Float_t pz1 = v4.Pz();



  Float_t px2 = v5.Px();
  Float_t py2 = v5.Py();
  Float_t pz2 = v5.Pz();



  TVector3 q1(px1,py1,pz1);
  TVector3 q2(px2,py2,pz2);
  Float_t mp1 = q1.Mag();
  Float_t mp2 = q2.Mag();



  Float_t P1 = TMath::Sqrt(px1*(px1) + py1*(py1) + pz1*(pz1));
  Float_t P2 = TMath::Sqrt(px2*(px2) + py2*(py2) + pz2*(pz2));

  Double_t muonE1 = TMath::Sqrt(mmuon*mmuon+ P1*P1);
  Double_t muonE2 = TMath::Sqrt(mmuon*mmuon + P2*P2);

  v4.SetPxPyPzE(px1,py1,pz1,muonE1);
  v5.SetPxPyPzE(px2,py2,pz2,muonE2);
  //cout<< px2<< py2 << pz2 << muonE2<< endl;

   MuPt1 = v4.Pt();

   MuPt2 = v5.Pt();
   //cout<<"mupt2 is "<< MuPt2 << endl;
  //Unnormalized Q
  TLorentzVector Q;
  Q = v4+v5;

  //Compute mass, Pt and rapidity of the dilepton system
   DiMuM = Q.M();
   DiMuPt = Q.Pt();
   DiMuY = Q.Rapidity();



  //calculate Accoplanarity Angle

  Float_t deltaPhi;

  Float_t fPhi1 = v4.Phi();
  Float_t fPhi2 = v5.Phi();

  deltaPhi = fPhi1 - fPhi2;
  accoplCut = 1. - TMath::Abs(deltaPhi)/TMath::Pi();

  //calculate z;
  TLorentzVector z;
  Float_t part1, part2;


  //Dot product: v1*v2 = t1*t2-x1*x2-y1*y2-z1*z2

  part1 = Q.Dot(pb);
  part2 = Q.Dot(pa);

  Float_t part3x = pa.X()*part1;
  Float_t part3y = pa.Y()*part1;
  Float_t part3z = pa.Z()*part1;
  Float_t part3e = pa.T()*part1;

  Float_t part4x = pb.X()*part2;
  Float_t part4y = pb.Y()*part2;
  Float_t part4z = pb.Z()*part2;
  Float_t part4e = pb.T()*part2;

  TLorentzVector part3(TVector3(part3x,part3y,part3z), part3e);
  TLorentzVector part4(TVector3(part4x,part4y,part4z),part4e);

  // Q=Q; pb=Pbar; pa=P; from paper

  //Un-normalized Z
  z = part3 - part4;

     //Normalized z
     Float_t normz = TMath::Sqrt(-z*z);
     Float_t znx = z.X()/normz;
     Float_t zny = z.Y()/normz;
     Float_t znz = z.Z()/normz;
     Float_t zne = z.E()/normz;

     //Normalized z
     TLorentzVector zhat(TVector3(znx,zny,znz),zne);

  // calculate x
  TLorentzVector x;

  Float_t constant1 = (Q.Dot(Q))/(2*(Q.Dot(pa)));
  Float_t constant2 = (Q.Dot(Q))/(2*(Q.Dot(pb)));

  Float_t comp1x = pa.X()*constant1;
  Float_t comp1y = pa.Y()*constant1;
  Float_t comp1z = pa.Z()*constant1;
  Float_t comp1e = pa.T()*constant1;

  TLorentzVector comp1(TVector3(comp1x,comp1y,comp1z),comp1e);

  Float_t comp2x = pb.X()*constant2;
  Float_t comp2y = pb.Y()*constant2;
  Float_t comp2z = pb.Z()*constant2;
  Float_t comp2e = pb.T()*constant2;

  TLorentzVector comp2(TVector3(comp2x,comp2y, comp2z),comp2e);

  //Un-normalized x
  x = Q - comp1 - comp2;


    //normalize x
    Float_t normx = TMath::Sqrt(-x*x);
    Float_t xnx = x.X()/normx;
    Float_t xny = x.Y()/normx;
    Float_t xnz = x.Z()/normx;
    Float_t xne = x.E()/normx;

     //Normalized x
    TLorentzVector xhat(TVector3(xnx,xny,xnz),xne);


  // calculate y
  //TLorentzVector y;
  Float_t yone = pa.Y()*pb.Z()*Q.E() - pa.Z()*pb.Y()*Q.E() + pa.Z()*pb.E()*Q.Y() + pa.E()*pb.Y()*Q.Z() - pa.Y()*pb.E()*Q.Z() - pa.E()*pb.Z()*Q.Y();
  Float_t ytwo = -pa.Z()*pb.E()*Q.X() + pa.Z()*pb.X()*Q.E() - pa.X()*pb.Z()*Q.E() + pa.X()*pb.E()*Q.Z() - pa.E()*pb.X()*Q.Z() + pa.E()*pb.Z()*Q.X();
  Float_t ythree = pa.X()*pb.Y()*Q.E() - pa.Y()*pb.X()*Q.E() + pa.Y()*pb.E()*Q.X() - pa.X()*pb.E()*Q.Y() + pa.E()*pb.X()*Q.Y() - pa.E()*pb.Y()*Q.X();
  Float_t yfour = -pa.X()*pb.Y()*Q.Z() + pa.X()*pb.Z()*Q.Y() - pa.Z()*pb.X()*Q.Y() + pa.Z()*pb.Y()*Q.X() - pa.Y()*pb.Z()*Q.X() + pa.Y()*pb.X()*Q.Z();

  //Un-normalized y
  TLorentzVector y(TVector3(yone,ytwo,ythree),yfour);

   //normalize y
   Float_t normy = TMath::Sqrt(-y*y);
   Float_t ynx = y.X()/normy;
   Float_t yny = y.Y()/normy;
   Float_t ynz = y.Z()/normy;
   Float_t yne = y.E()/normy;

  //normalized y
   TLorentzVector yhat(TVector3(ynx,yny,ynz),yne);

     //Lepton momentum difference
     TLorentzVector diff;
     diff = (v4 - v5);
     Float_t diff2x = diff.X()/2.;
     Float_t diff2y = diff.Y()/2.;
     Float_t diff2z = diff.Z()/2.;
     Float_t diff2e = diff.E()/2.;
     TLorentzVector diff2(TVector3(diff2x,diff2y,diff2z),diff2e);

     //Normalize diff2
     Float_t norm2 = TMath::Sqrt(-diff2*diff2);
     Float_t diff3x = diff2.X()/norm2;
     Float_t diff3y = diff2.Y()/norm2;
     Float_t diff3z = diff2.Z()/norm2;
     Float_t diff3e = diff2.E()/norm2;

     TLorentzVector diff3(TVector3(diff3x,diff3y,diff3z),diff3e);

     //computing the angles
      cosThetaCS =  zhat*diff3;
      SinThetaCosPhiCS = xhat*diff3;
      SinThetaSinPhiCS = yhat*diff3;
    //**************************************

  //Check that the CS frame was built correctly
  //   cout << "for Q "<< endl;
  //cout << Q.Dot(x) << "--" << Q.Dot(y)  << "--" << Q.Dot(z) << endl;
  //cout << "for x y and z"<< endl;
  //cout << x.Dot(y) <<"--" << y.Dot(z)  << "--"  << x.Dot(z) << endl;


  float phi   = atan2(SinThetaSinPhiCS,SinThetaCosPhiCS);
  if (phi>=0) phi = phi;
  if (phi<0) phi = phi + 2*TMath::Pi();

  //Make an ASCII fi  le as output

  if (isnan(cosThetaCS)){
  cosThetaCS = 0.0;
  }

  if (isnan(phi)){
  phi = 0.0;
  }
  Phi = phi; //for tree
  Theta = TMath::ACos(cosThetaCS);


}

void Correlation::ReadFile(TString inFileName = "file.dat"){
  ifstream file(inFileName);
  if (!file.is_open()){
    cout << "unable to open file"<< endl;
  }
  //data = new RooDataSet("data", "data", RooArgSet(cost,ph));
  string line;
  while (getline(file,line)){

    vector<float> ret;
    string num;
    while(line.find(" ") != string::npos){
      size_t p = line.find(" ");
      //cout << p<< endl;
      num = line.substr(0,p);
      //cout << num << endl;
      line.erase(0,p+1);
      while(line[0] == string(" ")[0])
        line.erase(0,1);

      ret.push_back(stof(num));
    } // end of reading file
    float px1 = ret[0];
    float py1 = ret[1];
    float pz1 = ret[2];
    float e1 = ret[3];
    float px2 = ret[4];
    float py2 = ret[5];
    float pz2 = ret[6];
    float e2 = ret[7];
    d1.SetPxPyPzE(px1,py1,pz1,e1);
    d2.SetPxPyPzE(px2,py2,pz2,e2);
    p = d1 +d2;
    AngleCalculator(d1,d2);

    //if(DiMuM>mass_min && DiMuM<mass_max && DiMuPt>Pt_min && DiMuPt>Pt_max ){
    costheta.push_back(cosThetaCS);
    phi.push_back(Phi);

    DiMuMass.push_back(DiMuM);
    DiMuPT.push_back(DiMuPt);


  }
}


void Correlation::SaveTree(const char* rootfile = "myTree.root"){


TFile*  myFile = new TFile(rootfile,"RECREATE");
TTree*  myTree = new TTree("MyTree", "MyTree");

myTree->Branch("MuPt1", &MuPt1);
myTree->Branch("MuPt2", &MuPt2);
myTree->Branch("DiMuM", &DiMuM);
myTree->Branch("DiMuPt", &DiMuPt);
myTree->Branch("DiMuY", &DiMuY);
myTree->Branch("AccoplAngle", &accoplCut);
myTree->Branch("CosTheta", &cosThetaCS);
myTree->Branch("SinThetaCosPhi", &SinThetaCosPhiCS);
myTree->Branch("SinThetaSinPhi", &SinThetaSinPhiCS);
myTree->Branch("Phi", &Phi);
myTree->Branch("Theta", &Theta);
}
