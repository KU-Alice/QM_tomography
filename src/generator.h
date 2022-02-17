#ifndef generator_h
#define generator_h

class generator{
public:
  vector <float> costheta;
  vector <float> phi;
  RooDataSet* data;
  //virtual void Generate(float m2 , float m3, float m6, int ngen);
  virtual void Generate(float m2 , float m3, float m6,RooRealVar mytheta, RooRealVar myPhi, int ngen);

private:
  RooRealVar mytheta;
  RooRealVar myPhi;

};


#endif

/*void generator::Generate(float mp2, float mp3, float mp6, int ngen =1000){

  RooRealVar md2("md2", "md2", mp2);
  RooRealVar md3("md3", "md3", mp3);
  RooRealVar md6("md6", "md6", mp6);
  RooRealVar mytheta = {"mytheta", "Cos(#theta)", -1., 1.} ;
  RooRealVar myPhi = {"myPhi", "Phi", 0., 2*3.141616} ;

  RooGenericPdf pdf("pdf", "Pdf_dN_dOmega", "(0.25*(1+md3*md3) +0.25*(1-3*md3*md3)*mytheta*mytheta -0.5* md3* md6* 2* TMath::Sqrt(1- mytheta*mytheta)* mytheta* cos(myPhi)+0.25*(2*md2*md2+md3*md3-1)*(1-mytheta*mytheta)*cos(2*myPhi))", RooArgSet(myPhi,mytheta, md2,md3,md6));
  RooDataSet* dummydata =pdf.generate(RooArgList(mytheta,myPhi),ngen);
  data = (RooDataSet*) dummydata->Clone();
  //data->get(0)->Print("v");
  //cout<<data->sumEntries()<<endl;
  //dummydata->get(0)->Print("V");
  //dummydata->get(1)->Print("V");
  auto m = mytheta.getVal();
  auto n = myPhi.getVal();
  for (int i=0 ; i< dummydata->sumEntries();i++){

  const RooArgSet* cost = dummydata->get(i);//->Print("V");
  dummydata->GetName();

  }

  //data = (RooDataSet*) dummydata->Clone();

}*/

void generator::Generate(float mp2, float mp3, float mp6,RooRealVar mytheta,RooRealVar myPhi, int ngen =1000){

  RooRealVar md2("md2", "md2", mp2);
  RooRealVar md3("md3", "md3", mp3);
  RooRealVar md6("md6", "md6", mp6);
  RooRealVar
  //RooRealVar mytheta = {"mytheta", "Cos(#theta)", -1., 1.} ;
  //RooRealVar myPhi = {"myPhi", "Phi", 0., 2*3.141616} ;

  RooGenericPdf pdf("pdf", "Pdf_dN_dOmega", "(0.25*(1+md3*md3) +0.25*(1-3*md3*md3)*mytheta*mytheta -0.5* md3* md6* 2* TMath::Sqrt(1- mytheta*mytheta)* mytheta* cos(myPhi)+0.25*(2*md2*md2+md3*md3-1)*(1-mytheta*mytheta)*cos(2*myPhi))", RooArgSet(myPhi,mytheta, md2,md3,md6));
  RooDataSet* dummydata =pdf.generate(RooArgList(mytheta,myPhi),ngen);
  data = (RooDataSet*) dummydata->Clone();
}
