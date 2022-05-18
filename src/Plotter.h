#ifndef Plotter_H
#define Plotter_H


class Plotter{
public:

  void Plot1D(vector <float> item,TString title , TString option,TString canvasname ="canvas");
  void Plot2D(vector <float> item1,vector <float> item2,TString title , TString option ,TString canvasname);
  void SetBinX(int bin){nbinshistx = bin;}
  void SetBinY(int bin){nbinshisty = bin;}
private:
  int nbinshistx=100;
  int nbinshisty=100;


};
#endif

void Plotter::Plot1D(vector <float> item,TString title ="; ; Events" , TString option ="",TString canvasname ="canvas"){

  float max = *max_element(item.begin(),item.end());
  float min = *min_element(item.begin(),item.end());
  float bin = nbinshistx;

  TH1F* histo;
  histo = new TH1F("histo",title.Data(),bin,min,max);

  for(int i=0 ; i< item.size();i++){

    histo->Fill(item[i]);
  }
  TCanvas* canvas = new TCanvas(canvasname.Data(), "",800,800);
  canvas->cd();
  histo->Draw(option.Data());
}

void Plotter::Plot2D(vector <float> item1,vector <float> item2,TString title ="; ; Events" , TString option ="",TString canvasname ="canvas" ){

  float max1 = *max_element(item1.begin(),item1.end());
  float min1 = *min_element(item1.begin(),item1.end());

  float max2 = *max_element(item2.begin(),item2.end());
  float min2 = *min_element(item2.begin(),item2.end());
  float binx = nbinshistx;
  float biny = nbinshisty;

  //float bin1 = item1.size()/1000;
  //float bin2 = item2.size()/1000;
  TH2F* histo = new TH2F("histo","; ; Events ",binx,min1,max1,biny,min2,max2);
  for(int i=0 ; i< item1.size();i++){

    histo->Fill(item1[i],item2[i]);
  }
  TCanvas* canvas = new TCanvas("canvas", "",800,800);
  canvas->cd();
  histo->Draw(option.Data());
}
