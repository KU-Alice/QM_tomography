#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include <complex>
#include <iostream>

using namespace RooFit;

void totFit(const char* infilename_tot = "ang_tot.txt", const char* infilename_sig = "ang_sig.txt",  const char* infilename_bkg = "ang_bkg_low.txt")
{
    //roofit
    RooRealVar t_tot("t_tot", "Cos(theta) [sig+bkg]", -1., 1) ;
    RooRealVar p_tot("p_tot", "phi angle [sig+bkg]", 0., 2*3.141616) ;

    RooRealVar t_sig("t_sig", "Cos(theta) [sig]", -1., 1) ;
    RooRealVar p_sig("p_sig", "phi angle [sig]", 0., 2*3.141616) ;

    RooRealVar t_bkg("t_bkg", "Cos(theta) [bkg]", -1., 1) ;
    RooRealVar p_bkg("p_bkg", "phi angle [bkg]", 0., 2*3.141616) ;

    RooDataSet* dataAngl_tot = RooDataSet::read(infilename_tot, RooArgList(t_tot,p_tot), "Q") ;
    RooDataSet* dataAngl_sig = RooDataSet::read(infilename_sig, RooArgList(t_sig,p_sig), "Q") ;
    RooDataSet* dataAngl_bkg = RooDataSet::read(infilename_bkg, RooArgList(t_bkg,p_bkg), "Q") ;
    //Int_t nEvts = dataAngl->numEntries();

   // ----------------------------------------------------
   // G e n e r i c   i n t e r p r e t e d   p . d . f .
   // ====================================================

   // Declare observables
    RooRealVar theta_tot("theta_tot", "theta_tot", -1, 1);
    RooRealVar phi_tot("phi_tot", "phi_tot", 0, 2*3.141516);

    RooRealVar theta_sig("theta_sig", "theta_sig", -1, 1);
    RooRealVar phi_sig("phi_sig", "phi_sig", 0, 2*3.141516);

    RooRealVar theta_bkg("theta_bkg", "theta_bkg", -1, 1);
    RooRealVar phi_bkg("phi_bkg", "phi_bkg", 0, 2*3.141516);

   // C o n s t r u c t   g e n e r i c   p d f   f r o m   i n t e r p r e t e d   e x p r e s s i o n

    //The m parameters are defined between -1 and 1
    //signal + background
    RooRealVar mt1("mt1", "mt1", 0);
    RooRealVar mt2("mt2", "mt2", 0.5,-1,1);
    RooRealVar mt3("mt3", "mt3", 0.5, -1,1);
    RooRealVar mt4("mt4", "mt4",0);//-1.0,1.);
    RooRealVar mt5("mt5", "mt5",0);// 0.5,-1,1);
    RooRealVar mt6("mt6", "mt6", -1.0,1.);
    RooRealVar mt7("mt7", "mt7",0);//0.5,-1.0,1.);
    RooRealVar mt8("mt8", "mt8", -1.0,1.);
    RooRealVar mt9("mt9", "mt9",  -1.0,1.);
    //signal
    RooRealVar ms1("ms1", "ms1", 0);
    RooRealVar ms2("ms2", "ms2", 0.5,-1,1);
    RooRealVar ms3("ms3", "ms3", 0.5, -1,1);
    RooRealVar ms4("ms4", "ms4", 0);//-1.0,1.);
    RooRealVar ms5("ms5", "ms5", 0);//0.5,-1,1);
    RooRealVar ms6("ms6", "ms6", -1.0,1.);
    RooRealVar ms7("ms7", "ms7",0);//0.5,-1.0,1.);
    RooRealVar ms8("ms8", "ms8",-1.0,1.);
    RooRealVar ms9("ms9", "ms9", -1.0,1.);
    //background
    RooRealVar mb1("mb1", "mb1", 0);
    RooRealVar mb2("mb2", "mb2", 0.5,-1,1);
    RooRealVar mb3("mb3", "mb3", 0.5, -1,1);
    RooRealVar mb4("mb4", "mb4",0);//-1.0,1.);
    RooRealVar mb5("mb5", "mb5",0);//0.5,-1,1);
    RooRealVar mb6("mb6", "mb6", -1.0,1.);
    RooRealVar mb7("mb7", "mb7",0);// -1.0,1.);
    RooRealVar mb8("mb8", "mb8",-1.0,1.);
    RooRealVar mb9("mb9", "mb9", -1.0,1.);

    //PDF for Fit
    //SIGNAL + BACKGROUND -------------------------------------------------------------------------------------------------------------------------
    RooGenericPdf Pdf_theta_tot("Pdf_theta_tot", "Pdf_theta_tot", "(1.5708*(1+ mt3*mt3) + 1.3823*(mt2*mt5+mt7*mt8-mt6*mt9)*t_tot + (1.5708-4.71239*mt3*mt3)*t_tot*t_tot)", RooArgSet(t_tot, mt2,mt3,mt5,mt6,mt7,mt8,mt9));
    RooGenericPdf Pdf_phi_tot("Pdf_phi_tot", "Pdf_phi_tot", "(0.66666667 + 0.345575*mt3*mt9*cos(p_tot)+0.33333333*(-1+2*mt2*mt2+mt3*mt3+2*mt8*mt8+2*mt9*mt9)*cos(2*p_tot)-0.345575*mt3*mt7*sin(p_tot)-0.66666667*(mt2*mt4+mt6*mt8+mt7*mt8)*sin(2*p_tot))", RooArgSet(p_tot, mt2,mt3,mt4,mt6,mt7,mt8,mt9));

    //SIGNAL ----------------------------------------------------------------------------------------------------------------------------------------
    RooGenericPdf Pdf_theta_sig("Pdf_theta_sig", "Pdf_theta_sig", "(1.5708*(1+ ms3*ms3) + 1.3823*(ms2*ms5+ms7*ms8-ms6*ms9)*t_sig + (1.5708-4.71239*ms3*ms3)*t_sig*t_sig)", RooArgSet(t_sig, ms2,ms3,ms5,ms6,ms7,ms8,ms9));
    RooGenericPdf Pdf_phi_sig("Pdf_phi_sig", "Pdf_phi_sig", "(0.66666667 + 0.345575*ms3*ms9*cos(p_sig)+0.33333333*(-1+2*ms2*ms2+ms3*ms3+2*ms8*ms8+2*ms9*ms9)*cos(2*p_sig)-0.345575*ms3*ms7*sin(p_sig)-0.66666667*(ms2*ms4+ms6*ms8+ms7*ms8)*sin(2*p_sig))", RooArgSet(p_sig, ms2,ms3,ms4,ms6,ms7,ms8,ms9));

    //BACKGROUND ------------------------------------------------------------------------------------------------------------------------------------
    RooGenericPdf Pdf_theta_bkg("Pdf_phi_bkg", "Pdf_phi_bkg", "(1.5708*(1+ mb3*mb3) + 1.3823*(mb2*mb5+mb7*mb8-mb6*mb9)*t_bkg + (1.5708-4.71239*mb3*mb3)*t_bkg*t_bkg)", RooArgSet(t_bkg, mb2,mb3,mb5,mb6,mb7,mb8,mb9));
    RooGenericPdf Pdf_phi_bkg("Pdf_phi_bkg", "Pdf_phi_bkg", "(0.66666667 + 0.345575*mb3*mb9*cos(p_bkg)+0.33333333*(-1+2*mb2*mb2+mb3*mb3+2*mb8*mb8+2*mb9*mb9)*cos(2*p_bkg)-0.345575*mb3*mb7*sin(p_bkg)-0.66666667*(mb2*mb4+mb6*mb8+mb7*mb8)*sin(2*p_bkg))", RooArgSet(p_bkg, mb2,mb3,mb4,mb6,mb7,mb8,mb9));
    //-----------------------------------------------------------------------------------------------------------------------------------------------


    // //PDF for data generation
    // //SIGNAL + BACKGROUND
    //     //-- data gen // -----------------------------------------------------------------------
    // RooGenericPdf Gen_theta_tot("Gen_theta_tot", "Gen_theta_tot", "(1.5708*(1+ mt3*mt3) + 1.3823*(mt2*mt5+mt7*mt8-mt6*mt9)*cos(theta_tot) + (1.5708-4.71239*mt3*mt3)*cos(theta_tot)*cos(theta_tot))", RooArgSet(theta_tot, mt2,mt3,mt5,mt6,mt7,mt8,mt9));
    // RooGenericPdf Gen_phi_tot("Gen_phi_tot", "Gen_phi_tot", "(0.66666667 + 0.345575*mt3*mt9*cos(phi_tot)+0.33333333*(-1+2*mt2*mt2+mt3*mt3+2*mt8*mt8+2*mt9*mt9)*cos(2*phi_tot)-0.345575*mt3*mt7*sin(phi_tot)-0.66666667*(mt2*mt4+mt6*mt8+mt7*mt8)*sin(2*phi_tot))", RooArgSet(phi_tot, mt2,mt3,mt4,mt6,mt7,mt8,mt9));

    // //SIGNAL
    //     //-- data gen --//-----------------------------------------------------------------------
    // RooGenericPdf Gen_theta_sig("Gen_theta_sig", "Gen_theta_sig", "(1.5708*(1+ ms3*ms3) + 1.3823*(ms2*ms5+ms7*ms8-ms6*ms9)*cos(theta_sig) + (1.5708-4.71239*ms3*ms3)*cos(theta_sig)*cos(theta_sig))", RooArgSet(theta_sig, ms2,ms3,ms5,ms6,ms7,ms8,ms9));
    // RooGenericPdf Gen_phi_sig("Gen_phi_sig", "Gen_phi_sig", "(0.66666667 + 0.345575*ms3*ms9*cos(phi_sig)+0.33333333*(-1+2*ms2*ms2+ms3*ms3+2*ms8*ms8+2*ms9*ms9)*cos(2*phi_sig)-0.345575*ms3*ms7*sin(phi_sig)-0.66666667*(ms2*ms4+ms6*ms8+ms7*ms8)*sin(2*phi_sig))", RooArgSet(phi_sig, ms2,ms3,ms4,ms6,ms7,ms8,ms9));

    // //BACKGROUND ------------------------------------------------------------------------------------------------------------------------------------
    //     //-- data gen --//-----------------------------------------------------------------------
    // RooGenericPdf Gen_theta_bkg("Gen_theta_bkg", "Gen_theta_bkg", "(1.5708*(1+ mb3*mb3) + 1.3823*(mb2*mb5+mb7*mb8-mb6*mb9)*cos(theta_bkg) + (1.5708-4.71239*mb3*mb3)*cos(theta_bkg)*cos(theta_bkg))", RooArgSet(theta_bkg, mb2,mb3,mb5,mb6,mb7,mb8,mb9));
    // RooGenericPdf Gen_phi_bkg("Gen_phi_bkg", "Gen_phi_bkg", "(0.66666667 + 0.345575*mb3*mb9*cos(phi_bkg)+0.33333333*(-1+2*mb2*mb2+mb3*mb3+2*mb8*mb8+2*mb9*mb9)*cos(2*phi_bkg)-0.345575*mb3*mb7*sin(phi_bkg)-0.66666667*(mb2*mb4+mb6*mb8+mb7*mb8)*sin(2*phi_bkg))", RooArgSet(phi_bkg, mb2,mb3,mb4,mb6,mb7,mb8,mb9));
    // //-----------------------------------------------------------------------------------------------------------------------------------------------


    //Add the two PDF components -- signal + background
    RooRealVar g1frac_tot("g1frac_tot","fraction of pdf 1",0.2);
    RooRealVar g2frac_tot("g2frac_tot","fraction of pdf 2",0.2);
    RooAddPdf sum_tot("sum_tot","Pdf_theta_tot+Pdf_phi_tot",RooArgList(Pdf_theta_tot,Pdf_phi_tot),RooArgList(g1frac_tot,g2frac_tot));
    //signal
    RooRealVar g1frac_sig("g1frac_sig","fraction of pdf 1",0.2);
    RooRealVar g2frac_sig("g2frac_sig","fraction of pdf 2",0.2);
    RooAddPdf sum_sig("sum_sig","Pdf_theta_sig+Pdf_phi_sig",RooArgList(Pdf_theta_sig,Pdf_phi_sig),RooArgList(g1frac_sig,g2frac_sig));
    //background
    RooRealVar g1frac_bkg("g1frac_bkg","fraction of pdf_bkg 1",0.2);
    RooRealVar g2frac_bkg("g2frac_bkg","fraction of pdf_bkg 2",0.2);
    RooAddPdf sum_bkg("sum_bkg","Pdf_theta_bkg+Pdf_phi_bkg",RooArgList(Pdf_theta_bkg,Pdf_phi_bkg),RooArgList(g1frac_bkg,g2frac_bkg));

    RooFitResult *rtot = sum_tot.fitTo(*dataAngl_tot, Save(1));
    // rtot->Print();

    RooFitResult *rsig = sum_sig.fitTo(*dataAngl_sig, Save(1));
    // rsig->Print();

    RooFitResult *rbkg = sum_bkg.fitTo(*dataAngl_bkg, Save(1));
    // rbkg->Print();


    Float_t rho_final[3][3];
// if dMatrix==True --> print desnity matrix in file
// if(dMatrix == true){
    const RooArgList & fitParams_tot = rtot->floatParsFinal();
    const RooArgList & fitParams_sig = rsig->floatParsFinal();
    const RooArgList & fitParams_bkg = rbkg ->floatParsFinal();
    int size = fitParams_tot.getSize();
    cout<< "size of the paramters is "<< size<< endl;
    Float_t param_tot[size];
    Float_t param_sig[size];
    Float_t param_bkg[size];
    cout <<"number of fit parameters are"<< endl;
    for ( int i = 0; i < fitParams_tot.getSize(); ++i)
    {
    auto & fitPar_tot = (RooRealVar &) fitParams_tot[i];
    auto & fitPar_sig = (RooRealVar &) fitParams_sig[i];
    auto & fitPar_bkg = (RooRealVar &) fitParams_bkg[i];

    param_tot[i] = fitPar_tot.getVal();
    param_sig[i] = fitPar_sig.getVal();
    param_bkg[i] = fitPar_bkg.getVal();
  }
    // const char* outfilename = "dmatrix.dat";
    // std::ofstream outfile (outfilename);
    for(int i = 0 ;i<size; i++)
    {
      cout << "fit parameters:: "<< param_tot[i] << endl;
    }


    // Building the density matrix (with 3 parameters)
    Float_t rho_tot[3][3];

    rho_tot[0][0] = param_tot[3]*param_tot[3];//param_tot[4]*param_tot[4];
    rho_tot[0][1] = 0;
    rho_tot[0][2] = param_tot[1]*param_tot[3];
    rho_tot[1][0] = 0;
    rho_tot[1][1] = param_tot[0]*param_tot[0];
    rho_tot[1][2] = 0;
    rho_tot[2][0] = param_tot[1]*param_tot[3];
    rho_tot[2][1] = 0;
    rho_tot[2][2] = param_tot[1]*param_tot[1];

    Float_t rho_bkg[3][3];
    rho_bkg[0][0] = 1 - param_bkg[3]*param_bkg[3];
    rho_bkg[0][1] = 0;
    rho_bkg[0][2] = param_bkg[1]*param_bkg[3];
    rho_bkg[1][0] = 0;
    rho_bkg[1][1] = 1 - param_bkg[0]*param_bkg[0];
    rho_bkg[1][2] = 0;
    rho_bkg[2][0] =  param_bkg[1]*param_bkg[3];
    rho_bkg[2][1] = 0;
    rho_bkg[2][2] = 1 - param_bkg[1]*param_bkg[1];


    Float_t rho_sig[3][3];
    rho_sig[0][0] = param_sig[3]*param_sig[3];
    rho_sig[0][1] = 0;
    rho_sig[0][2] = param_sig[1]*param_sig[3];
    rho_sig[1][0] = 0;
    rho_sig[1][1] = param_sig[0]*param_sig[0];
    rho_sig[1][2] = 0;
    rho_sig[2][0] = param_sig[1]*param_sig[3];
    rho_sig[2][1] = 0;
    rho_sig[2][2] = param_sig[1]*param_sig[1];

    int rows =  sizeof(rho_tot) / sizeof(rho_tot[0]);
    int cols = sizeof(rho_tot[0]) / sizeof(rho_tot[0][0]);


    // if (outfile.is_open()){
    //     for(int i = 0; i < rows; ++i){
    //         for(int j = 0; j < cols; ++j)
    //         {
    //             rho_final[i][j]=0;
    //         }
    //     }

    //     for(int i = 0; i < rows; ++i){
    //         for(int j = 0; j < cols; ++j){
    //             for(int k = 0; k < cols; ++k)
    //             {
    //                 rho_final[i][j] += rho_bkg[i][k] * rho_tot[k][j];
    //             }
    //             outfile <<rho_final[i][j]<<" ";
    //         }
    //         outfile<<std::endl;
    //     }
    //     outfile.close();
    // }
    // else if(!outfile.is_open()){
    //     std:cout << "Unable to open file";
    //     }
    // }


    std::cout<<"************"<<std::endl;
    std::cout<<"***************"<<std::endl;
    std::cout<<"*****************"<<std::endl;
    std::cout<<"*******************"<<std::endl;

    for(int i = 0; i < rows; ++i){
                for(int j = 0; j < cols; ++j)
                {
                    rho_final[i][j]=0;
                }
            }

            for(int i = 0; i < rows; ++i){
                for(int j = 0; j < cols; ++j){
                    // std::cout<<"rho_bkg["<<i<<"]["<<j<<"] ="<<rho_bkg[i][j]<<std::endl;
                    // std::cout<<"rho_tot["<<i<<"]["<<j<<"] ="<<rho_tot[i][j]<<std::endl;
                    for(int k = 0; k < cols; ++k)
                    {

                        rho_final[i][j] += rho_bkg[i][k] * rho_tot[k][j];

                    }
                    std::cout<<i<<"--"<<j<<"--"<<rho_final[i][j]<<"------------ "<< rho_sig[i][j]<<std::endl;

                }

            }
    // }


    std::cout<<"*******************"<<std::endl;
    std::cout<<"*****************"<<std::endl;
    std::cout<<"***************"<<std::endl;
    std::cout<<"************"<<std::endl;


  Float_t m6_ratio = TMath::Sqrt(rho_sig[0][0])/TMath::Sqrt(rho_final[0][0]);
  Float_t m2_ratio = TMath::Sqrt(rho_sig[1][1])/TMath::Sqrt(rho_final[1][1]);
  Float_t m3_ratio = TMath::Sqrt(rho_sig[2][2])/TMath::Sqrt(rho_final[2][2]);
  std::cout<< "#####################################################"<< m6_ratio<<" ...." <<m2_ratio <<" ...." <<m3_ratio << std::endl;




    RooRealVar theta("theta", "theta", -1, 1);
    RooRealVar phi("phi", "phi", 0, 2*3.141516);

    RooRealVar theta2("theta2", "theta2", -1, 1);
    RooRealVar phi2("phi2", "phi2", 0, 2*3.141516);


    RooRealVar m1("m1", "m1",  0);//rho_final[0][0]);
    RooRealVar m2("m2", "m2", rho_final[0][1]);
    RooRealVar m3("m3", "m3", rho_final[0][2]);
    RooRealVar m4("m4", "m4",  0);//rho_final[1][0]);
    RooRealVar m5("m5", "m5", 0);// rho_final[1][1]);
    RooRealVar m6("m6", "m6", rho_final[1][2]);
    RooRealVar m7("m7", "m7",  0);//rho_final[2][0]);
    RooRealVar m8("m8", "mb",  0);//rho_final[2][1]);
    RooRealVar m9("m9", "m9",  0);//rho_final[2][2]);


    RooRealVar mm1("mm1", "mm1",  0);//rho_final[0][0]);
    RooRealVar mm2("mm2", "mm2", 0.2,-1,1);
    RooRealVar mm3("mm3", "mm3", 0.2,-1,1);
    RooRealVar mm4("mm4", "mm4",  0);//rho_final[1][0]);
    RooRealVar mm5("mm5", "mm5", 0);// rho_final[1][1]);
    RooRealVar mm6("mm6", "mm6", 0.2,-1,1);
    RooRealVar mm7("mm7", "mm7",  0);//rho_final[2][0]);
    RooRealVar mm8("mm8", "mm8",  0);//rho_final[2][1]);
    RooRealVar mm9("mm9", "mm9",  0);//rho_final[2][2]);

    //PDF for dN/dCosTheta
    RooGenericPdf dNdcosTheta("dNdcosTheta", "dNdcosTheta", "(1.5708*(1+ m3*m3) + 1.3823*(m2*m5+m7*m8-m6*m9)*theta + (1.5708-4.71239*m3*m3)*theta*theta)", RooArgSet(theta, m2,m3,m5,m6,m7,m8,m9));
    //PDF for dN/dPhi
    RooGenericPdf dndPhi("dndPhi", "dndPhi", "(0.66666667 + 0.345575*m3*m9*cos(phi)+0.33333333*(-1+2*m2*m2+m3*m3+2*m8*m8+2*m9*m9)*cos(2*phi)-0.345575*m3*m7*sin(phi)-0.66666667*(m2*m4+m6*m8+m7*m8)*sin(2*phi))", RooArgSet(phi, m2,m3,m4,m6,m7,m8,m9));

    //PDF for dN/dCosTheta
    RooGenericPdf dNdcostheta2("dNdcostheta2", "dNdcostheta2", "(1.5708*(1+ mm3*mm3) + 1.3823*(mm2*mm5+mm7*mm8-mm6*mm9)*theta + (1.5708-4.71239*mm3*mm3)*theta*theta)", RooArgSet(theta, mm2,mm3,mm5,mm6,mm7,mm8,mm9));
    //PDF for dN/dPhi
    RooGenericPdf dndphi2("dndphi2", "dndphi2", "(0.66666667 + 0.345575*mm3*mm9*cos(phi)+0.33333333*(-1+2*mm2*mm2+mm3*mm3+2*mm8*mm8+2*mm9*mm9)*cos(2*phi)-0.345575*mm3*mm7*sin(phi)-0.66666667*(mm2*mm4+mm6*mm8+mm7*mm8)*sin(2*phi))", RooArgSet(phi, mm2,mm3,mm4,mm6,mm7,mm8,mm9));



    RooRealVar gg1frac("gg1frac","fraction of pdf 1",1);
    RooRealVar gg2frac("gg2frac","fraction of pdf 2",0.2);
    RooAddPdf sum("sum","dNdcosTheta+dndPhi",RooArgList(dNdcosTheta,dndPhi),RooArgList(gg1frac,gg2frac));


    RooRealVar gg3frac("gg3frac","fraction of pdf 1",1);
    RooRealVar gg4frac("gg4frac","fraction of pdf 2",0.2);
    // RooGenericPdf dNdOmega("dNdOmega", "dNdOmega", "(1.5708*(1+ m3*m3) + (1.5708-4.71239*m3*m3)*theta*theta)+ 0.2*((0.33333333*(-1+2*m2*m2+m3*m3)*cos(2*phi)))", RooArgSet(theta,phi,m2,m3,m6));
    RooAddPdf sum2("sum2","dNdcosTheta2+dndPhi2",RooArgList(dNdcostheta2,dndphi2),RooArgList(gg3frac,gg4frac));


    // //Generate a toy dataset from the interpreted
    RooDataSet *data1 = dNdcosTheta.generate(theta, 10000);
    RooDataSet *data2 = dndPhi.generate(phi, 10000);
    RooDataSet *dataSum = sum.generate(RooArgSet(theta,phi), 10000);

    //RooGenericPdf dratio("dratio", "dratio", theta/t_sig, RooArgSet(theta,t_sig));

    //RooDataSet *dataSum = sum.generate(RooArgSet(theta,phi), 10000);




    // Fit the interpreted pdf to the generated data
    //dNdcosTheta.fitTo(*data1);

    //RooFitResult *t_final = dNdcosTheta.fitTo(*data1, Save(1));

    //dndPhi.fitTo(*data2);
    //RooFitResult *p_final = dndPhi.fitTo(*data2, Save(1));
    //p_final->Print(0);

    //sum.fitTo(*dataSum);
    //RooFitResult *rSum = sum.fitTo(*dataSum,Save());

    //rSum->Print();


  //  RooRealVar ratiotheta = theta/t_sig;


    //Make a plot of the actual angles -- signal + background
    RooPlot *tframe_tot = t_tot.frame(Title("dN/dCosTheta [signal + bkg] expression pdf"));
  //RooPlot *tframe_tot = theta.frame(Title("dN/dCosTheta [signal + bkg] expression pdf"));
    dataAngl_tot->plotOn(tframe_tot);
    Pdf_theta_tot.plotOn(tframe_tot);
    //dratio->plotOn(tframe_tot);


    RooPlot *pframe_tot = p_tot.frame(Title("dN/dPhi [signal + bkg] expression pdf"));
    dataAngl_tot->plotOn(pframe_tot);
    Pdf_phi_tot.plotOn(pframe_tot);

    //Make a plot of the actual angles -- background
    RooPlot *tframe_bkg= t_bkg.frame(Title("dN/dCosTheta [bkg] expression pdf"));
    dataAngl_bkg->plotOn(tframe_bkg);
    Pdf_theta_bkg.plotOn(tframe_bkg);

    RooPlot *pframe_bkg = p_bkg.frame(Title("dN/dPhi [bkg] expression pdf"));
    dataAngl_bkg->plotOn(pframe_bkg);
    Pdf_phi_bkg.plotOn(pframe_bkg);

    // Make a plot of the data and the pdf overlaid
    RooPlot *tframe_final = theta.frame(Title("dN/dCosTheta expression pdf"));
    data1->plotOn(tframe_final);
    dNdcosTheta.plotOn(tframe_final);

    RooPlot *pframe_final = phi.frame(Title("dN/dPhi expression pdf"));
    data2->plotOn(pframe_final);
    dndPhi.plotOn(pframe_final);

    //std::cout<< "parameter" << mm1 << std::endl;
    // RooPlot *sum_final = sum.frame(Title("dN/dPhi expression pdf"));
    // dataSum->plotOn(sum_final);
    // sum.plotOn(sum_final);


    // std::cout<<"HEREEEEEEEE"<<std::endl;
    // std::cout<<"************"<<std::endl;
    // std::cout<<"***************"<<std::endl;
    // std::cout<<"*****************"<<std::endl;
    // std::cout<<"*******************"<<std::endl;
    /*const RooArgList & fitParams_final =rSum->floatParsFinal();
    int size_fin = fitParams_final.getSize();
    Float_t param_final[size_fin];
    Float_t ratio[size_fin];
    for ( int i = 0; i < fitParams_final.getSize(); ++i)
    {
    auto & fitPar_final = (RooRealVar &) fitParams_tot[i];
    param_final[i] = fitPar_final.getVal();
    ratio[i] = param_final[i]/param_sig[i];

    std::cout<<ratio<<std::endl;
  }*/

    //Float_t ratio1 = mt1/mm1

    //std::cout<< mm2<< std::endl;
    // std::cout<<"HEREEEEEEEE"<<std::endl;
    // std::cout<<"************"<<std::endl;
    // std::cout<<"***************"<<std::endl;
    // std::cout<<"*****************"<<std::endl;
    // std::cout<<"*******************"<<std::endl;

    // Draw all frames on a canvas
    TCanvas *c = new TCanvas("Angular dist", "Angular dist", 800, 400);
    c->Divide(2,2);

    c->cd(1);
    gPad->SetLeftMargin(0.15);
    tframe_tot->GetYaxis()->SetTitleOffset(1.4);
    tframe_tot->Draw();


    c->cd(2);
    gPad->SetLeftMargin(0.15);
    pframe_tot->GetYaxis()->SetTitleOffset(1.4);
    pframe_tot->Draw();


    c->cd(3);
    gPad->SetLeftMargin(0.15);
    tframe_bkg->GetYaxis()->SetTitleOffset(1.4);
    tframe_bkg->Draw();


    c->cd(4);
    gPad->SetLeftMargin(0.15);
    pframe_bkg->GetYaxis()->SetTitleOffset(1.4);
    pframe_bkg->Draw();


    // Draw all frames on a canvas
    TCanvas *c1 = new TCanvas("Angular dist final", "Angular dist final", 800, 400);
    c1->Divide(2);

    //Make a plot of the final angle distributions
    c1->cd(1);
    gPad->SetLeftMargin(0.15);
    tframe_final->GetYaxis()->SetTitleOffset(1.4);
    tframe_final->Draw();


    c1->cd(2);
    gPad->SetLeftMargin(0.15);
    pframe_final->GetYaxis()->SetTitleOffset(1.4);
    pframe_final->Draw();

    // c1->cd(3);
    // gPad->SetLeftMargin(0.15);
    // sum_final->GetYaxis()->SetTitleOffset(1.4);
    // sum_final->Draw();

}
