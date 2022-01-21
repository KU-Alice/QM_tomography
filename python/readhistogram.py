#!/usr/bin/python
#include<stdio.h>
import code
import math
import ROOT as rt
import fitter as ft
from ROOT import TFile, TTree, TCanvas, TH1F, TList, TH2F,TH3F ,TMath, TF1, TStyle, gStyle, TRefArray, TClonesArray, TObjArray, gPad, TPaveText, TLegend ,TString, TObject,gROOT,TFormula, TEllipse, TDirectory,TLorentzVector
from ROOT import TMath as mt
from datetime import date

mixedfilename = "../root/anglulardist_big.root"
mixedfile  = TFile.Open(mixedfilename)
tree = TTree()
tree = mixedfile.Get("MyTree")
#nentries = tree.GetEntries()

#bkg_hist = TH1F("bkg_hist","histogram for background",500,3.5,6)
sig_hist = TH1F("sig_hist","histogram for signal",500,3.095,3.11);
mixed_hist = TH1F("mixed_hist", "histogram for mixed data",1000,0,3.2)
#canvas_bkg = TCanvas("canvas_bkg","canvas background",600,600)
canvas_sig = TCanvas("canvas_sig")
canvas_mix = TCanvas("canvas_mix")
expo_tail = TF1("expo_tail", "[0]*TMath::Exp(-[1]*x)",3.5,6)
expo_tail.SetParLimits(1,0.0001,1)

sig_func_jpsi =  TF1("sig_func_jpsi","gaus",2,4)

sig_func_jpsi.SetParameters(1000,3.0969,0.002)
sig_func_jpsi.SetParLimits(1,3.09689,3.09693)
sig_func_jpsi.SetLineColor(rt.kBlack)

sig_func_jpsi2 =  TF1("sig_func_jpsi2","crystalball",0,5)
sig_func_jpsi2.SetParameters(1,3.15,0.090,1,115)

comb_func = TF1("comb_func","[0]*TMath::Exp(-[1]*x)+crystalball(2)",2,5)



canvas_sig.cd()
tree.Draw("DiMuM>>sig_hist","","e")
sig_hist.SetDirectory(0)
sig_hist.Fit("sig_func_jpsi","R","",3.09685,3.097)
sig_func_jpsi2.FixParameter(1,sig_func_jpsi.GetParameter(1))
sig_func_jpsi2.FixParameter(2,sig_func_jpsi.GetParameter(2))
sig_func_jpsi2.SetLineColor(rt.kGreen)
sig_hist.Fit("sig_func_jpsi2","R+","",3.0966,3.097)
canvas_mix.cd()
tree.Draw("DiMuM>>mixed_hist","","e")
mixed_hist.SetDirectory(0)
#mixed_hist.Draw()
mixed_hist.Fit("expo_tail","R")
comb_func.SetParameters(100,1,1,100,1,1,1,1)
comb_func.FixParameter(1,expo_tail.GetParameter(1))
#comb_func.FixParameter(2,expo_tail.GetParameter(2))
comb_func.FixParameter(3,sig_func_jpsi2.GetParameter(1))
comb_func.FixParameter(4,sig_func_jpsi2.GetParameter(2))
comb_func.FixParameter(5,sig_func_jpsi2.GetParameter(3))
comb_func.FixParameter(6,sig_func_jpsi2.GetParameter(4))
comb_func.SetLineColor(rt.kGreen)
mixed_hist.Fit("comb_func","R+")

mu = comb_func.GetParameter(3)
sigma = comb_func.GetParameter(4)

expo2 = TF1("expo2","[0]*TMath::Exp(-[1]*x)",2,5)
sig_funct = TF1("sig_funct","crystalball",2,5)
expo2.FixParameter(0,comb_func.GetParameter(0))
expo2.FixParameter(1,comb_func.GetParameter(1))
sig_funct.FixParameter(0,comb_func.GetParameter(2))
sig_funct.FixParameter(1,comb_func.GetParameter(3))
sig_funct.FixParameter(2,comb_func.GetParameter(4))
sig_funct.FixParameter(3,comb_func.GetParameter(5))
sig_funct.FixParameter(4,comb_func.GetParameter(6))
sig_funct.SetLineColor(rt.kBlack)
sig_funct.Draw("SAME")

nsignalpb = comb_func.Integral(mu+2*sigma,mu-2*sigma)
nbkg = expo2.Integral(mu+4*sigma,mu-4*sigma)
nsignal = nsignalpb-nbkg
signaltobkg = nsignal/nbkg* mixed_hist.GetBinWidth(100)
print ("signal to background ratio is : ", signaltobkg)
#canvas_sig.cd()
#tree.Draw("DiMuM>>sig_hist","","e")
#sig_hist.SetDirectory(0)
#sig_hist.Draw()


#anvas_mix.cd()
#tree.Draw("DiMuM>>mix","","e")
#mixed_hist.SetDirectory(0)
#mixed_hist.Draw()
#The lines below are so that pyroot do not exit on completing the analysis
vars = globals()
vars.update(locals())
shell = code.InteractiveConsole(vars)
shell.interact()
