struct trkState
{
  // vector and covariance matrix of track paramters 
  // (see https://www-jlc.kek.jp/subg/offl/kaltest/doc/ReferenceManual.pdf)
  // 0 - d_rho
  // 1 - phi_0
  // 2 - k
  // 3 - d_x
  // 4 - tanl
  ROOT::Math::SVector>double,5> parVec;    
  ROOT::Math::SMatrix<double,5,5> covMat;
  
}

struct procState
{
  trkState Pd; // Prediction
  trkState Ft; // Filtering
  trkState Sm; // Smoothing
}

struct trkKal
{
  std::vector<procState> KF; 
}

void propagateParVec(ROOT::Math::SVector>double,5>& currParVec, ROOT::Math::SVector>double,5>& prevParVec, TVector3& p0, TVector3& p1, double alpha)
{
  double zc = p0.Z() + (prevParVec[0] + alpha/prevParVec[2]) * TMath::Cos(prevParVec[1]);
  double yc = p0.Y() + (prevParVec[0] + alpha/prevParVec[2]) * TMath::Sin(prevParVec[1]);
  
  currParVec[1] = (prevParVec[2] > 0) ? TMath::ATan2(yc - p1.Y(),zc - p1.Z()) : TMath::ATan2(p1.Y() - yc,p1.Z() - zc);  
  currParVec[0] = (zc - p1.Z()) * TMath::Cos(currParVec[1]) + (yc - p1.Y()) * TMath::Sin(currParVec[1]) - alpha/prevParVec[2];
  currParVec[2] = prevParVec[2];
  currParVec[3] = p0.X() - p1.X() + prevParVec[3] - alpha/prevParVec[2] * (currParVec[1] - prevParVec[1]) * prevParVec[4];
  currParVec[4] = prevParVec[4]; 
}

void propagatorMtx(ROOT::Math::SVector>double,5>& currParVec, ROOT::Math::SVector>double,5>& prevParVec, TVector3& p0, TVector3& p1, double alpha)
{
  
}

void propagateHelix(TVector3& p1, TVector3& p0, ROOT::Math::SVector>double,5>& parVec, double phi, double alpha)
{
  p1.SetX(p0.X() + parVec[3] - alpha/parVec[2] * parVec[4] * phi);
  p1.SetZ(p0.Z() + parVec[0] * TMath::Cos(parVec[1]) + alpha/parVec[2] * (TMath::Cos(parVec[1]) - TMath::Cos(parVec[1] + phi)));
  p1.SetY(p0.Y() + parVec[0] * TMath::Sin(parVec[1]) + alpha/parVec[2] * (TMath::Sin(parVec[1]) - TMath::Sin(parVec[1] + phi)));
}

void processNoise(ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepSym<double,5> >& Q, ROOT::Math::SVector>double,5>& parVec, double Xl)
{
  double p = 1./parVec[2] * TMath::Sqrt(1 + parVec[4]*parVec[4]);
  double sigma_MS = 0.0141/p * TMath::Sqrt(Xl) * (1 + TMath::Log10(Xl)/9.);
  
  Q[0][0] = Q[0][1] = Q[0][2] = Q[0][3] = Q[0][4] = 0.;
  Q[1][2] = Q[1][3] = Q[1][4] = 0.;
  Q[2][3] = 0.;
  Q[3][3] = Q[3][4] = 0.;
  
  Q[1][1] = 1 + parVec[4] * parVec[4];
  Q[1][1] = 1 + parVec[4] * parVec[4];
  Q[2][2] = (parVec[2]*parVec[4])*(parVec[2]*parVec[4]);
  Q[0][2] = (parVec[2]*parVec[4])*(1+parVec[4]*parVec[4]);
  Q[4][4] = (1 + parVec[4] * parVec[4])*(1 + parVec[4] * parVec[4]);
  
  
}

