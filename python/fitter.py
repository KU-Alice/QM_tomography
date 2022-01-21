#!/usr/bin/python
import code
import math
#import histograms as h1

import ROOT as rt
import scipy

from ROOT import TFile, TTree , TCanvas, TH1F, TList , TH2F ,TH3F ,TMath ,TF1, TStyle ,gStyle , TRefArray, TClonesArray, TObjArray, gPad, TPaveText, TLegend ,TString, TObject,gROOT,TFormula, TEllipse
from ROOT import TMath as mt
#Opening the file>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#file = TFile.Open("Incohernt_Phi_PID1.root")
#list = TList()
#list = file.Get("UPCPhiTaskTest/UPCPhiWithTrigger")
#tree = list.FindObject("scatterplot")
#tree2 = file.Get("tree2")
#tree.SetMarkerStyle(3)
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

##Defining crystalball function>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

func = TF1("func","crystalball",0,8)
wphi = TF1("wphi","[0]*(1 + 2*[1]/(3+[2]) * TMath::Cos(2*x))",-3.14,3.14)
wphi.SetParameters(4000,0.049,1.208)
wphi.SetParLimits(2,0.8,1.5)
#wphi.FixParameter(1,0.4)
#wphi.FixParameter(2,1.5)
wtheta =TF1("wtheta", "([0]/(3+[1]))*(1+[1]*x*x)",-0.6,0.6)
wtheta.SetParameters(7000,1.4)

wtildetheta = TF1("wtildetheta","[0]*(1+(TMath::Sqrt(2)*[1]/(3+[2])*TMath::Cos(x)))",0,6.280)
wtildetheta.SetParameters(1,1,1)
wphi.SetParLimits(2,0.8,1.5)







sig_func_psi =  TF1("sig_func_psi","crystalball",0,10)

sig_func_psi.SetParameters(1,3.67,0.070,1.08,20)

sig_func_psi.SetParLimits(3,0.5,1.5)
sig_func_psi.SetParLimits(1,3.64,3.74)
#sig_func_psi.SetParLimits(4,1.5,100)



sig_func_psi_data =  TF1("sig_func_psi_data","crystalball",0,5)
sig_func_psi_data.SetParameters(1,3.67,0.070,1.08,20)

sig_func_psi_data.SetParLimits(3,0.5,1.5)
sig_func_psi_data.SetParLimits(1,3.64,3.74)
#sig_func_psi_data.SetParLimits(4,1.5,100)

#sig_func2 =  TF1("sig_func2",cball,4,10,4)
sig_func_jpsi =  TF1("sig_func_jpsi","crystalball",0,5)
sig_func_jpsi.SetParameters(1,3.15,0.090,1,115)
#sig_func_jpsi.FixParameter(3,1.08)
sig_func_jpsi.SetParLimits(1,3.05,3.17)
#sig_func_jpsi.SetParLimits(4,110,120)



sig_func_jpsi_data =  TF1("sig_func_jpsi_data","crystalball",2,4)

sig_func_jpsi_data.SetParameters(1,3.115,0.090,1.08,115)
sig_func_jpsi_data.SetParLimits(1,3.05,3.17)
#sig_func_jpsi_data.FixParameter(3,1.08)
#sig_func_jpsi_data.SetParLimits(1,3.09,3.17)
#sig_func_jpsi_data.SetParLimits(4,110,120)





#comb_my_func = TF1("comb_my_func",sig_func(0)+sig_func2(5)+tail_func(9),0,10)




#function implemented from simone's analysis code

#GammaGammaFit   = TF1( "GammaGammaFit","[0]*TMath::Exp(-[1]*x)*( (x < 4) ? 1 : 1 + [2]*(x-4)*(x-4) + [3]*(x-4)*(x-4)*(x-4) +[4]*(x-4)*(x-4)*(x-4)*(x-4) )",0,10)
GammaGammaFit   = TF1( "GammaGammaFit","[0]*TMath::Exp(-[1]*x)*( (x < [2]) ? 1 : 1 + [3]*(x-[2])*(x-[2]) + [4]*(x-[2])*(x-[2])*(x-[2]) +[5]*(x-[2])*(x-[2])*(x-[2])*(x-[2]))",2,10)
#GammaGammaFit.SetParameters(0.025,0.9,0.0001,0.0001,0.0001,0.0001)
GammaGammaFit.SetParameters(19,1.73,1.51,5881,255,-170)
GammaGammaFit.SetParLimits(0,0.0001,20)
GammaGammaFit.SetParLimits(1,0.001,10)
GammaGammaFit.SetParLimits(2,0.5,10)
#GammaGammaFit.FixParameter(4,4)
GammaGammaFit.SetParLimits(3,0.00000001,10000)
GammaGammaFit.SetParLimits(4,0.0000001,1000)


GammaGammaFit2   = TF1( "GammaGammaFit2","[0]*TMath::Exp(-[1]*x)*( (x < [2]) ? 1 : 1 + [3]*(x-[2])*(x-[2]) + [4]*(x-[2])*(x-[2])*(x-[2]) +[5]*(x-[2])*(x-[2])*(x-[2])*(x-[2]))",2,8)
GammaGammaFit2.SetParameters(19,1.73,1.51,5881,255,-170)
GammaGammaFit2.SetParLimits(0,0.0001,30)
GammaGammaFit2.SetParLimits(1,0.001,10)
GammaGammaFit2.SetParLimits(2,0.5,10)
GammaGammaFit2.SetParLimits(3,0.00000001,10000)
GammaGammaFit2.SetParLimits(4,0.0000001,1000)



#GammaGammaFit  = TF1( "GammaGammaFit","[0]*TMath::Exp(-[1]*(x))*(1  + [2]*(x-4)*(x-4) + [3]*(x-4)*(x-4)*(x-4) +[4]*(x-4)*(x-4)*(x-4)*(x-4))",3,6)
#GammaGammaFit   = TF1( "GammaGammaFit","[0]*TMath::Exp(-[1]*x)*( (x < 3) ? 1 : 1 + [2]*(x-3)*(x-3) + [3]*(x-3)*(x-3)*(x-3) +[4]*(x-3)*(x-3)*(x-3)*(x-3) )",0,7)
#GammaGammaFit   = TF1( "GammaGammaFit","[0]*TMath::Exp(-[1]*(x))*(1*(x>=4)+(x<4)*(1  + [2]*(x-4)*(x-4) + [3]*(x-4)*(x-4)*(x-4) +[4]*(x-4)*(x-4)*(x-4)*(x-4) ))",0,4)

#GammaGammaFit = TF1("GammaGammaFit","[0]*TMath::Exp(-[1]*(x))",7,10)
#myfunc_test = TF1("myfunc_test",)
#low_mass_fit = TF1("low_mass_fit","(1+ [0]*(x-4)**2 +[1]*(x-4)**3+[2]*(x-4)**4)",0,4)
expo_tail = TF1("expo_tail", "[0]*TMath::Exp(-[1]*x)",3.5,6)
expo_tail.SetParLimits(1,0.0001,1)

#low_mass_fit = TF1("low_mass_fit","(1+ [0]*(x-4)*(x-4) +[1]*(x-4)*(x-4)*(x-4)+[2]*(x-4)*(x-4)*(x-4)*(x-4))",3,4)
#GammaGammaFit = TF1("GammaGammaFit",tail,4,10,5)
#low_mass_fit = TF1("low_mass_fit","([0]+ [1]*(x-4)*(x-4) +[2]*(x-4)*(x-4)*(x-4)+[3]*(x-4)*(x-4)*(x-4)*(x-4))",3,4)


#low_mass_fit = TF1("low_mass_fit","[0]*TMath::Exp(-[1]*(x))*(1+ [2]*(x-4)*(x-4) +[3]*(x-4)*(x-4)*(x-4)+[4]*(x-4)*(x-4)*(x-4)*(x-4))",2,3)
#low_mass_fit = TF1("low_mass_fit","TMath::Exp(-[0]*(x))*(1+ [1]*(x-4)*(x-4) +[2]*(x-4)*(x-4)*(x-4)+[3]*(x-4)*(x-4)*(x-4)*(x-4))",3,4)




#GammaGammaFit.SetParameters(2418,0.5,-17.5,-35.6,-19.39)
#GammaGammaFit.SetParameters(1,1,1,1,1,1)
#GammaGammaFit.SetParameters(1,1,1,1,1)
#GammaGammaFit.Draw()
#expo_tail.SetParameters(1,1)
#low_mass_fit.SetParameters(1,1,1,1,1)
#GammaGammaFit.SetParameter(1,1)
#GammaGammaFit.SetParameters(1,1,1,1)

#[10]*TMath::Exp(-[11]*x)*( (x < [12]) ? 1 : 1 + [13]*(x-[12])*(x-[12]) + [14]*(x-[12])*(x-[12])*(x-[12]) +[15]*(x-[12])*(x-[12])*(x-[12])*(x-[12]))
comb_func = TF1("comb_func","crystalball(0)+crystalball(5)+([10]*TMath::Exp(-[11]*x)*( (x < [12]) ? 1 : 1 + [13]*(x-[12])*(x-[12]) + [14]*(x-[12])*(x-[12])*(x-[12]) +[15]*(x-[12])*(x-[12])*(x-[12])*(x-[12])))",2,4)

#comb_func = TF1("comb_func",sig_bkg,0,10,9)
#comb_func.SetParameter(9,1)
#tail_func_mc_psi = TF1("tail_func_mc_psi",tail,4,10,4)

#comb_func_mc_psi = TF1("comb_func_mc_psi",sig_bkg,0,10,9)

#final_sig_bkg_func = TF1("final_sig_func",sig_func(0)+sign_func2(4)+tail_func(8))

#ft.comb_func.SetParameters(1,3.1,0.090,1,115,1,3.67,0.070,1.08,20,19,1.73,1.51,5881,255,-170)
#ft.comb_func.SetParameters(0,1,3.1,0.090,1,115,1,3.67,0.070,1.08,20,19,1.73,1.51,5881,255,-170)

comb_func.SetParName(0,"constjp")
comb_func.SetParName(1,"meanjp")
comb_func.SetParName(2,"sigmajp")
comb_func.SetParName(3,"alphajp")
comb_func.SetParName(4,"Njp")
#comb_func.SetParName(17,"scaleconstant2")
comb_func.SetParName(5,"constpsi")
comb_func.SetParName(6,"meanpsi")
comb_func.SetParName(7,"sigmapsi")
comb_func.SetParName(8,"alphapsi")
comb_func.SetParName(9,"Npsi")
#comb_func.SetParName(18,"scaleconstant3")
comb_func.SetParName(10,"gammaconst")
comb_func.SetParName(11,"b")
comb_func.SetParName(12,"mth")
comb_func.SetParName(13,"a1")
comb_func.SetParName(14,"a2")
comb_func.SetParName(15,"a3")


comb_func.SetParameter(0,100)
comb_func.SetParLimits(0,40,700)

comb_func.SetParameter(1,3.1)
comb_func.SetParLimits(1,3.,3.17)
comb_func.SetParameter(2,0.090)
comb_func.SetParLimits(2,0.05,0.15)
comb_func.SetParameter(3,1)
comb_func.SetParameter(4,1)

comb_func.SetParameter(5,1)
comb_func.SetParLimits(5,0.00000001,10)
comb_func.SetParameter(6,3.67)
comb_func.SetParLimits(6,3.64,3.74)
comb_func.SetParameter(7,0.070)
comb_func.SetParLimits(7,0.01,0.15)

#comb_func.SetParLimits(8,1,1.5)
comb_func.SetParameter(8,1)
comb_func.SetParameter(9,1)
comb_func.SetParameter(10,19)
comb_func.SetParLimits(10,0.001,100)
comb_func.SetParameter(11,0.9)
comb_func.SetParLimits(11,0.8,10)
comb_func.SetParameter(12,1.3)
comb_func.SetParLimits(12,0.1,5)
comb_func.SetParameter(13,50)
comb_func.SetParameter(14,255)
comb_func.SetParameter(15,-170)
#comb_func.SetParLimits(0,1E-9,1E9)
#comb_func.SetParLimits(5,1E-9,1E9)
#comb_func.SetParLimits(10,1E-9,1E9)
#comb_func.SetParLimits(1,3.09,3.2)
#comb_func.SetParLimits(8,0.5,1.5)
partest = []
k  = 0
while k<16:
  #print comb_func.GetParameter(k)
  partest.append(comb_func.GetParameter(k))
  #print partest[k]
  k = k + 1
