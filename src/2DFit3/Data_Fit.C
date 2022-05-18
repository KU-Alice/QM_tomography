void Data_Fit(){
gStyle->SetOptStat(0);
gStyle->SetLineWidth(2);
TFile * One_D_Corrected = TFile::Open("2D_Corrected.root");
//yieldTildePhi->Divide(yieldTildePhiRec,yieldTildePhiAxE);
TH1D* yieldPhi=(TH1D*) One_D_Corrected->Get("yieldPhiCsCorrected_px");
TH1D* yieldCsT=(TH1D*) One_D_Corrected->Get("yieldPhiCsCorrected_py");






TF1* pdftheta_root = new TF1("pdftheta_root","[0]*(1.5708*(1+ [1]*[1])+ (1.5708-4.71239*[1]*[1])*x*x)",-1,1);
TF1* pdfphi_root = new TF1("pdfphi_root","[0]*(0.66666667 +0.33333333*(-1+2*[1]*[1]+[2]*[2])*cos(2*x))",0,2*TMath::Pi());
pdftheta_root->SetLineWidth(3);
pdfphi_root->SetLineWidth(3);

TF1* wphi = new TF1("wphi","[0]*(1 + 2*[1]/(3+[2]) * TMath::Cos(2*x))",0,2*TMath::Pi());
wphi->SetParameters(1000,0.049,1.208);
wphi->SetParLimits(2,0.34,1.5);
wphi->SetLineColor(kBlue);
//wphi.FixParameter(1,0.4);
//wphi.FixParameter(2,1.5);
TF1* wtheta =new TF1("wtheta", "([0]/(3+[1]))*(1+[1]*x*x)",-1,1);
wtheta->SetParameters(1000,1.4);
wtheta->SetLineColor(kBlue);

TF1* wtheta2 =new TF1("wtheta2", "([0]/(3+[1]))*(1+[1]*x*x)",-1,1);
wtheta2->SetParameters(1000,1.4);
wtheta2->SetLineColor(kBlue);


///TF1* wtildephi =new TF1("wtildephi","[0]*(1+(TMath::Sqrt(2)*[1]/(3+[2])*TMath::Cos(x)))",0,2*TMath::Pi());
//wtildephi->SetParameters(1,1,1);
//TF1* pdfphi_root = new TF1("pdfphi_root","[0]*(0.66666667)",0,2*3.141616);

pdftheta_root->SetParName(0,"Constant");
pdftheta_root->SetParName(1,"m3");

pdftheta_root->SetParameter(0,100);
pdftheta_root->SetParameter(1,0.5);

pdfphi_root->SetParameter(0,100);
pdfphi_root->SetParName(0,"Constant");
pdfphi_root->SetParName(1,"m2");
pdfphi_root->SetParName(2,"m3");



pdfphi_root->SetParameter(1,0.05);
pdfphi_root->SetParameter(2,0.05);
pdfphi_root->SetParLimits(1,-1,1);
pdfphi_root->SetParLimits(2,-1,1);
//pdfphi_root->FixParameter(2,7.42211e-01);
//pdfphi_root->FixParameter(0,1.00888e+03*((4.*TMath::Pi())/(3.)));

TCanvas *newcanvas = new TCanvas("newcanvas","mycanvas",1800,800);
TLegend *legend = new  TLegend(0.6,0.85,0.7,0.55);
legend->SetFillStyle(0);
legend->SetHeader("Corrected J/#Psi","");

legend->SetLineColor(0);
legend->SetTextSize(0.025);

TLegend *legend2 = new  TLegend(0.6,0.85,0.7,0.55);
legend2->SetFillStyle(0);
legend2->SetHeader("Corrected J/#Psi","");

legend2->SetLineColor(0);
legend2->SetTextSize(0.025);

TLegend *legend3 = new  TLegend(0.6,0.85,0.7,0.55);
legend3->SetFillStyle(0);
legend3->SetHeader("Corrected J/#Psi","");

legend3->SetLineColor(0);
legend3->SetTextSize(0.05);


TLegend *legend4 = new  TLegend(0.6,0.85,0.7,0.55);
legend4->SetFillStyle(0);
legend4->SetHeader("Corrected J/#Psi","");

legend4->SetLineColor(0);
legend4->SetTextSize(0.015);





newcanvas->Divide(2,1);

newcanvas->cd(1);
yieldCsT->SetMarkerSize(2);
yieldCsT->SetLineWidth(2);
//tframe->Draw();
yieldCsT->Draw("text45e");
yieldCsT->GetYaxis()->SetRangeUser(0,3000);

yieldCsT->Fit("pdftheta_root","","IL");
yieldCsT->Fit("wtheta","","IL");


TString m_3;
m_3.Form("m3 = %f #pm %f",pdftheta_root->GetParameter(1),pdftheta_root->GetParError(1));

TString lambda_theta_2;
lambda_theta_2.Form("#lambda_{#theta} = %f #pm %f",wtheta->GetParameter(1),wtheta->GetParError(1));

legend2->AddEntry((TObject*)0,m_3.Data() ,"");
legend2->AddEntry((TObject*)0,lambda_theta_2.Data() ,"");
legend2->Draw("SAME");

newcanvas->cd(2);
yieldPhi->SetMarkerSize(2);
yieldPhi->SetLineWidth(2);

yieldPhi->Draw("text45e");
//yieldPhir->Draw("e");
yieldPhi->GetYaxis()->SetRangeUser(0,3000);
pdfphi_root->FixParameter(2,pdftheta_root->GetParameter(1));
yieldPhi->Fit("pdfphi_root","","IL");
wphi->FixParameter(2,wtheta->GetParameter(1));
yieldPhi->Fit("wphi","","IL");


TString m_3_2;
m_3_2.Form("m3 = %f #pm %f",pdfphi_root->GetParameter(2),pdftheta_root->GetParError(1));

TString m_2;
m_2.Form("m2 = %f #pm %f",pdfphi_root->GetParameter(1),pdfphi_root->GetParError(1));

TString lambda_phi;
lambda_phi.Form("#lambda_{#phi} = %f #pm %f",wphi->GetParameter(1),wphi->GetParError(1));

TString lambda_theta;
lambda_theta.Form("#lambda_{#theta} = %f #pm %f",wphi->GetParameter(2),wtheta->GetParError(1));

legend->AddEntry((TObject*)0,m_2.Data() ,"");
legend->AddEntry((TObject*)0,m_3_2.Data() ,"");
legend->AddEntry((TObject*)0,lambda_theta.Data() ,"");
legend->AddEntry((TObject*)0,lambda_phi.Data() ,"");
legend->Draw("SAME");


}
