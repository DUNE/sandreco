
void ParToMeasure(double* next_meas, double* prev_meas, double* par, double alpha, double phi)
{
  next_meas[0] = prev_meas[0] + par[0] * TMath::Cos(par[1]) + alpha/par[2] * (TMath::Cos(par[1]) - TMath::Cos(par[1]+phi));
  next_meas[1] = prev_meas[1] + par[0] * TMath::Sin(par[1]) + alpha/par[2] * (TMath::Sin(par[1]) - TMath::Sin(par[1]+phi));
  next_meas[2] = prev_meas[2] + par[3] - alpha/par[2] * par[4] * phi;
}

void PropagateState(double* prev_meas, double* this_meas, double* this_par, double* prev_par, double alpha)
{
  double xc = prev_meas[0] + (prev_par[0] + alpha/prev_par[2]) * TMath::Cos(prev_par[1]);
  double yc = prev_meas[1] + (prev_par[0] + alpha/prev_par[2]) * TMath::Sin(prev_par[1]);
  
  this_par[0] = (xc - this_meas[0]) * TMath::Cos(prev_par[1]) + (yc - this_meas[1]) * TMath::Sin(prev_par[1]) - alpha/prev_par[2];
  if(par[2] > 0)
  {
    this_par[1] = TMath::ATan((yc - this_meas[1])/(xc - this_meas[0]));
  }
  else
  {
    this_par[1] = TMath::ATan((this_meas[1] - yc)/(this_meas[0] - xc));
  }
  this_par[2] = prev_par[2];
  this_par[3] = prev_meas[2] - this_meas[2] + prev_par[3] - alpha/prev_par[2] * (this_par[1] - prev_par[1]) * prev_par[4];
  this_par[4] = prev_par[4]; 
}

