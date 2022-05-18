#ifndef AngleCalc_H
#define AngleCalc_H

class AngleCalc{
public:
  /*vector <float> costheta_dt;
  vector <float> phi_dt;
  vector <Double_t> costheta_cs;
  vector <Double_t> costheta_he;
  vector <float> phi_cs;
  vector <float> phi_he;*/
  float phi_he;
  float phi_cs;
  float ct_he;
  float ct_cs;
  float ct_dt;
  float phi_dt;


  void CosThetaHelicityFrame( TLorentzVector muonPositive,TLorentzVector muonNegative);
  void PhiHelicityFrame(  TLorentzVector muonPositive,TLorentzVector muonNegative );
  void CosThetaCollinsSoper( TLorentzVector muonPositive,TLorentzVector muonNegative);
  void PhiCollinsSoper(  TLorentzVector muonPositive,TLorentzVector muonNegative);

  void AngleCalculator_DT(TLorentzVector v1, TLorentzVector v2); // both phi and costheta

private:

};
#endif

void AngleCalc::PhiHelicityFrame(  TLorentzVector muonPositive,TLorentzVector muonNegative)
{
  /* - This function computes the helicity phi for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
  */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TLorentzVector possibleJPsi = muonPositive +muonNegative;
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle
  // (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis = (possibleJPsi.Vect()).Unit();
  TVector3 yaxis = ((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
  TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
  //
  // --- Calculation of the azimuthal angle (Helicity)
  //
  phi_he = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
  if (isnan(phi_he)){
  phi_he = 0.0;
  }
  if (phi_he<0){
    phi_he = phi_he+2*TMath::Pi();
  }
}




void AngleCalc::CosThetaHelicityFrame( TLorentzVector muonPositive,TLorentzVector muonNegative)
  {
    /* - This function computes the Helicity cos(theta) for the
       - helicity of the J/Psi.
       - The idea should be to get back to a reference frame where it
       - is easier to compute and to define the proper z-axis.
       -
     */

    /* - Half of the energy per pair of the colliding nucleons.
       -
     */
    Double_t HalfSqrtSnn   = 2510.;
    Double_t MassOfLead208 = 193.6823;
    Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
    /* - Fill the Lorentz vector for projectile and target.
       - For the moment we do not consider the crossing angle.
       - Projectile runs towards the MUON arm.
       -
     */
    TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
    TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
    /* - Translate the dimuon parameters in the dimuon rest frame
       -
     */
    TLorentzVector possibleJPsi = muonPositive +muonNegative;
    TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
    TLorentzVector pMu1Dimu  = muonPositive;
    TLorentzVector pMu2Dimu  = muonNegative;
    TLorentzVector pProjDimu = pProjCM;
    TLorentzVector pTargDimu = pTargCM;
    pMu1Dimu.Boost(beta);
    pMu2Dimu.Boost(beta);
    pProjDimu.Boost(beta);
    pTargDimu.Boost(beta);
    //
    // --- Determine the z axis for the calculation of the polarization angle
    // (i.e. the direction of the dimuon in the CM system)
    //
    TVector3 zaxis = (possibleJPsi.Vect()).Unit();
    /* - Determine the He angle (angle between mu+ and the z axis defined above)
       -
     */
    ct_he = zaxis.Dot((pMu1Dimu.Vect()).Unit());
    if (isnan(ct_he)){
      ct_he = 0.0;
    }

}

void AngleCalc::CosThetaCollinsSoper( TLorentzVector muonPositive,TLorentzVector muonNegative)
{
  /* - This function computes the Collins-Soper cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TLorentzVector possibleJPsi = muonPositive + muonNegative;
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */


  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  /* - Determine the CS angle (angle between mu+ and the z axis defined above)
     -
   */
  //cout << "%%%%%%%%%%%%%%%%%%%"<< (pMu1Dimu.Vect()).Dot(yaxisCS)<< "&&&&&&&&&&&&&&&&&"<<((pMu1Dimu.Vect()).Dot(xaxisCS))<< endl;
  ct_cs = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
  if (isnan(ct_cs)){
  ct_cs = 0.0;
  }


}

void AngleCalc::PhiCollinsSoper( TLorentzVector muonPositive,TLorentzVector muonNegative)
{
  /* - This function computes the Collins-Soper PHI for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  //TLorentzVector pProjCM(0.,0.,-1., 1.); // projectile
  //TLorentzVector pTargCM(0.,0., 1., 1.); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TLorentzVector possibleJPsi =  muonPositive + muonNegative ;
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
  TVector3 yaxisCS=(((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
  TVector3 xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit();


  phi_cs = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
  if (isnan(phi_cs)){
  phi_cs = 0.0;
  }
  if (phi_cs<0){
    phi_cs = phi_cs +2* TMath::Pi();
  }
}

void AngleCalc::AngleCalculator_DT(TLorentzVector v4, TLorentzVector v5)

{

  //TLorentzVector v4;
  //TLorentzVector v5;
  //TLorentzVector pa(1.,0.,0,1); // projectile
  //TLorentzVector pb(1.,0., 0,-1); // target
  TLorentzVector pa(0.,0.,-1.,1.); // projectile
  TLorentzVector pb(0.,0., 1.,1); // target
  Float_t px1 = v4.Px();
  Float_t py1 = v4.Py();
  Float_t pz1 = v4.Pz();



  Float_t px2 = v5.Px();
  Float_t py2 = v5.Py();
  Float_t pz2 = v5.Pz();



  //Unnormalized Q
  TLorentzVector Q;
  Q = v4+v5;

  //Compute mass, Pt and rapidity of the dilepton system




  //calculate Accoplanarity Angle

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
      float cosThetaCS =  zhat*diff3;
      float SinThetaCosPhiCS = xhat*diff3;
      float SinThetaSinPhiCS = yhat*diff3;
    //**************************************

  //Check that the CS frame was built correctly
  //   cout << "for Q "<< endl;
  //cout << Q.Dot(x) << "--" << Q.Dot(y)  << "--" << Q.Dot(z) << endl;
  //cout << "for x y and z"<< endl;
  //cout << x.Dot(y) <<"--" << y.Dot(z)  << "--"  << x.Dot(z) << endl;


  float phi   = atan2(SinThetaSinPhiCS,SinThetaCosPhiCS);
  if (phi>=0) {phi = phi;}
  if (phi<0) {phi = phi + 2*TMath::Pi();}

  //Make an ASCII fi  le as output

  if (isnan(cosThetaCS)){
  cosThetaCS = 0.0;
  }

  if (isnan(phi)){
  phi = 0.0;
  }
  ct_dt = cosThetaCS;
  phi_dt = phi;

}
