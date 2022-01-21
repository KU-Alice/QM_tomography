#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include <complex>

using namespace RooFit;

void myfit3(const char* infilename = "angle.txt", bool dMatrix = false )
{

    //roofit
    RooRealVar t("t", "Cos(theta) ", -1., 1) ;
    RooRealVar p("p", "phi angle ", 0., 2*3.141616) ;

    RooDataSet* dataAngl = RooDataSet::read(infilename, RooArgList(t,p), "Q") ;
    //Int_t nEvts = dataAngl->numEntries();

   // ----------------------------------------------------
   // G e n e r i c   i n t e r p r e t e d   p . d . f .
   // ====================================================

   // Declare observables
    RooRealVar theta("theta", "theta", 0, 3.141516);
    RooRealVar phi("phi", "phi", 0, 2*3.141516);

   // C o n s t r u c t   g e n e r i c   p d f   f r o m   i n t e r p r e t e d   e x p r e s s i o n
   // -------------------------------------------------------------------------------------------------

   // To construct a proper pdf, the formula expression is explicitly normalized internally by dividing
   // it by a numeric integral of the expression over x in the range [-20,20]
   //
    //The m parameters are defined between -1 and 1
    RooRealVar m1("m1", "m1", 0);
    RooRealVar m2("m2", "m2", 0.5,-1,1); //
    RooRealVar m3("m3", "m3", 0.5, -1,1); //
    RooRealVar m4("m4", "m4", -1.0,1.);
    RooRealVar m5("m5", "m5", 0.5,-1,1); //
    RooRealVar m6("m6", "m6", -1.0,1.);
    RooRealVar m7("m7", "m7", -1.0,1.);
    RooRealVar m8("m8", "m8",  -1.0,1.);
    RooRealVar m9("m9", "m9",  -1.0,1.);
    // RooRealVar m1("m1", "m1", 0);
    // RooRealVar m2("m2", "m2", 1,-1,1); //
    // RooRealVar m3("m3", "m3", 1, -1,1); //
    // RooRealVar m4("m4", "m4", 0);
    // RooRealVar m5("m5", "m5", 0); //
    // RooRealVar m6("m6", "m6", 0);
    // RooRealVar m7("m7", "m7", 0);
    // RooRealVar m8("m8", "m8", 0);
    // RooRealVar m9("m9", "m9", 0);

    //PDF for dN/dCosTheta
    RooGenericPdf genpdf1("genpdf1", "genpdf1", "(1.5708*(1+ m3*m3) + 1.3823*(m2*m5+m7*m8-m6*m9)*cos(theta) + (1.5708-4.71239*m3*m3)*cos(theta)*cos(theta))", RooArgSet(theta, m2,m3,m5,m6,m7,m8,m9));

    //PDF for dN/dPhi
    RooGenericPdf genpdf2("genpdf2", "genpdf2", "(0.66666667 + 0.345575*m3*m9*cos(phi)+0.33333333*(-1+2*m2*m2+m3*m3+2*m8*m8+2*m9*m9)*cos(2*phi)-0.345575*m3*m7*sin(phi)-0.66666667*(m2*m4+m6*m8+m7*m8)*sin(2*phi))", RooArgSet(phi, m2,m3,m4,m6,m7,m8,m9));

    RooGenericPdf genpdf3("genpdf3", "genpdf3", "(1.5708*(1+ m3*m3) + 1.3823*(m2*m5+m7*m8-m6*m9)*t + (1.5708-4.71239*m3*m3)*t*t)", RooArgSet(t, m2,m3,m5,m6,m7,m8,m9));

    RooGenericPdf genpdf4("genpdf4", "genpdf4", "(0.66666667 + 0.345575*m3*m9*cos(p)+0.33333333*(-1+2*m2*m2+m3*m3+2*m8*m8+2*m9*m9)*cos(2*p)-0.345575*m3*m7*sin(p)-0.66666667*(m2*m4+m6*m8+m7*m8)*sin(2*p))", RooArgSet(p, m2,m3,m4,m6,m7,m8,m9));


    //RooGenericPdf genpdf4("genpdf4", "genpdf4", "(0.66666667)", RooArgSet(p, m2,m3));

    //Add the two PDF components
    RooRealVar g1frac("g1frac","fraction of pdf 1",1);
    RooRealVar g2frac("g2frac","fraction of pdf 2",0.2);
    RooAddPdf sum("sum","genpdf3+genpdf4",RooArgList(genpdf3,genpdf4),RooArgList(g1frac,g2frac));

   // S a m p l e ,   f i t   a n d   p l o t   g e n e r i c   p d f
   // ---------------------------------------------------------------

   //Generate a toy dataset from the interpreted pdf
  // RooDataSet *data1 = genpdf1.generate(theta, 10000);
  // RooDataSet *data2 = genpdf2.generate(phi, 10000);

  // RooDataSet *dataSum = sum.generate(RooArgSet(theta,phi), 10000);

   // Fit the interpreted pdf to the generated data
  // genpdf1.fitTo(*data1);
  // RooFitResult *r1 = genpdf1.fitTo(*data1, Save());
   //r1->Print();

  //RooFitResult *r2 = sum.fitTo(*data2, Save());
    //r2->Print();

  //  RooFitResult *rSum = sum.fitTo(*dataSum, Save());
  //  rSum->Print();


    RooFitResult *r3 = sum.fitTo(*dataAngl, Save());
    Float_t param[8]; 
    r3->Print();

  if(dMatrix == true){
    const RooArgList & fitParams = r3->floatParsFinal();
    for ( int i = 0; i < fitParams.getSize(); ++i)
    {                   
      auto & fitPar = (RooRealVar &) fitParams[i];
      param[i] = fitPar.getVal(); 
      std::cout << fitPar.GetName() << " " << fitPar.getVal() << std::endl;
    }

   std::ofstream outfile ("matrix_signal.dat");
    // Building the density matrix (with 3 parameters)
    Float_t rho[3][3]; 
    rho[0][0] = param[5]*param[5]; 
    rho[0][1] = 0; 
    rho[0][2] = param[3]*param[6]; 
    rho[1][0] = 0; 
    rho[1][1] = param[1]*param[1]; 
    rho[1][2] = 0; 
    rho[2][0] = param[3]*param[6]; 
    rho[2][1] = 0;
    rho[2][2] = param[2]*param[2];

    int rows =  sizeof(rho) / sizeof(rho[0]);   
    int cols = sizeof(rho[0]) / sizeof(rho[0][0]); 
    
    if (outfile.is_open()){
      for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
          outfile <<rho[i][j]<<" ";
        }
        outfile<<std::endl; 
      }
    }
    else{std:cout << "Unable to open file";}
    

    outfile.close();
  } 


   // Make a plot of the data and the pdf overlaid
 //  RooPlot *tframe = theta.frame(Title("dN/dCosTheta expression pdf"));
 //  data1->plotOn(tframe);
 //  genpdf1.plotOn(tframe);

  //  RooPlot *pframe = theta.frame(Title("dN/dPhi expression pdf"));
  //  data2->plotOn(pframe);
   // sum.plotOn(pframe);

    //Make a plot of the actual angles
    RooPlot *t2frame = t.frame(Title("dN/dCosTheta expression pdf"));
    dataAngl->plotOn(t2frame);
    genpdf3.plotOn(t2frame);

    RooPlot *p2frame = p.frame(Title("dN/dPhi expression pdf"));
    dataAngl->plotOn(p2frame);
    genpdf4.plotOn(p2frame);



   // Draw all frames on a canvas
   TCanvas *c = new TCanvas("Angular dist", "Angular dist", 800, 400);
   c->Divide(2);

    /*
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   tframe->GetYaxis()->SetTitleOffset(1.4);
   tframe->Draw();

    c->cd(2);
    gPad->SetLeftMargin(0.15);
    pframe->GetYaxis()->SetTitleOffset(1.4);
    pframe->Draw();
    */

    c->cd(1);
   gPad->SetLeftMargin(0.15);
    t2frame->GetYaxis()->SetTitleOffset(1.4);
    t2frame->Draw();


    c->cd(2);
    gPad->SetLeftMargin(0.15);
    p2frame->GetYaxis()->SetTitleOffset(1.4);
    p2frame->Draw();


}
