#include <TCanvas.h>
#include <TF1.h>
#include <TTree.h>
#include <TGraph.h>
#include <TFile.h>
#include <TEllipse.h>
#include <TLorentzVector.h>
#include <TLine.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TBox.h>
#include <TVector3.h>

#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>

const int max_hits = 2000;
const int max_hits_for_y_fit = 2000;
const int max_hits_for_x_fit = 2000;

TFile* f = 0;
TTree* t = 0;

int hit_n = 0;
int eventid = 0;
double hit_x[max_hits];
double hit_y[max_hits];
double hit_z[max_hits];
double hit_t[max_hits];
double hit_e[max_hits];
bool hit_hor[max_hits];
std::vector<string>* hit_det = new  std::vector<string>(max_hits);
TLorentzVector* mom = 0;
TLorentzVector* vtx = 0;
TVector3* mom_last_in = 0;
TLorentzVector* pos_last_in = 0;

bool exiting;

bool initialized = false;

int fitLinear(int nPoints, double *x, double *y, bool withres, double res, 
        double &a, double &b, double** cov, double &chi2)
{
    a = -999;
    b = -999;
    chi2 = -999;
    
    cov[0][0] = -999;;
    cov[0][1] = -999;
    cov[1][0] = -999;
    cov[1][1] = -999;
    
    if (nPoints < 2)
    {
        std::cout << "fitLinear    - condition: \"nPoints < 2 not satisfied" << std::endl;
        return 1;		
    }
    
    double var = res*res;
    
    double S1  = 0;
    double SX  = 0;
    double SY  = 0;
    double SYX = 0;
    double SXX = 0;
    
    for(int i = 0; i < nPoints; i++)
    {
        S1 += 1;
        SX += x[i];
        SY += y[i];
        SYX += x[i]*y[i];
        SXX += x[i]*x[i];
    }
    
    if(withres) 
    {
        S1 /= var;
        SX /= var;
        SY /= var;
        SYX /= var;
        SXX /= var;
    }
    
    double D = S1*SXX - SX*SX;
    
    a = (SY*SXX - SX*SYX)/D;  // i.p. at x = 0
    b = (S1*SYX - SX*SY)/D;   // tg(theta)
    
    cov[0][0] = SXX/D;
    cov[0][1] = -SX/D;
    cov[1][0] = -SX/D;
    cov[1][1] = S1/D;
    
    chi2 = 0.0;
    
    for(int i = 0; i < nPoints; i++)
    {
        chi2 += (y[i] - a - b * x[i]) * (y[i] - a - b * x[i]);
    }
    
    if(withres) 
        chi2 /= var;
    
    chi2 /= nPoints;
    
    return 0;
}

int fitCircle(int nPoints, double *x, double *y, double &xc, double &yc, double &r, double &errr, double &chi2)
{
    xc = -999;
    yc = -999;
    r = -999;
    errr = -999;
    chi2 = -999;
    
    if (!x || !y)
    {
        std::cout << "fitCircle    - condition: \"!x || !y not satisfied" << std::endl;
        return 1; // error code or exception to be raised;		
    }
    
    if (nPoints < 3)
    {
        std::cout << "fitCircle    - condition: \"nPoints >= 3\" not satisfied" << std::endl;
        return 2; // error code or exception to be raised;		
    }
    
    double sumx = 0, sumy = 0;                            // linear    terms
    double sumx2 = 0, sumy2 = 0, sumxy = 0;               // quadratic terms
    double sumxy2 = 0, sumx2y = 0, sumx3 = 0, sumy3 = 0;  // cubic     terms
    for (int iPoints = 0; iPoints < nPoints; iPoints++)
    {
            double xp = x[iPoints];
            double yp = y[iPoints];
            sumx   += xp;       sumy   += yp;
            sumx2  += xp*xp;    sumy2  += yp*yp;    sumxy += xp*yp;
            sumxy2 += xp*yp*yp; sumx2y += xp*xp*yp; sumx3 += xp*xp*xp; sumy3 += yp*yp*yp;
    }

    double a = nPoints*sumx2 - sumx*sumx;
    double b = nPoints*sumxy - sumx*sumy;
    double c = nPoints*sumy2 - sumy*sumy;
    double d = 0.5*(nPoints*sumxy2 - sumx*sumy2 + nPoints*sumx3 - sumx*sumx2);
    double e = 0.5*(nPoints*sumx2y - sumy*sumx2 + nPoints*sumy3 - sumy*sumy2);
    
    if(a*c - b*b == 0.)
    {
        std::cout << "fitCircle    - condition: \"a*c - b*b == 0.\" not satisfied" << std::endl;
        return 3;
    }

    xc = (d*c - b*e) / (a*c - b*b);
    yc = (a*e - b*d) / (a*c - b*b);

    double rMean = 0;
    double rrms = 0;

    for (int iPoints = 0; iPoints < nPoints; iPoints++)
    {
        double xp = x[iPoints];
        double yp = y[iPoints];
        double r2 = (xp - xc)*(xp - xc) + (yp - yc)*(yp - yc);

        rMean += sqrt(r2);
        rrms += r2;
    }


    rMean /= nPoints;
    rrms /= nPoints;
    r = rMean;

    if(rrms - rMean*rMean > 0.)
        errr = sqrt(rrms - rMean*rMean);
    else
    {
        std::cout << "fitCircle    - condition: \"rrms - rMean*rMean > 0.\" not satisfied" << std::endl;
        return 4;
    }
    
    chi2 = 0.0;
    
    for(int i = 0; i < nPoints; i++)
    {
        chi2 += TMath::Abs((y[i] - yc)*(y[i] - yc) + (x[i] - xc) * (x[i] - xc) - r*r);
    }
    
    chi2 /= nPoints;

    return 0;
}

int fitQuadratic(int nPoints, double *x, double *y, bool withres, double res, 
        double &a, double &b, double &c, double** cov, double& chi2)
{
    a = -999;
    b = -999;
    c = -999;
    chi2 = -999;
    
    cov[0][0] = -999;
    cov[0][1] = -999;
    cov[0][2] = -999;
    cov[1][0] = -999;
    cov[1][1] = -999;
    cov[1][2] = -999;
    cov[2][0] = -999;
    cov[2][1] = -999;
    cov[2][2] = -999;
    
    if (!x || !y)
    {
        std::cout << "fitQuadratic - condition: \"!x || !y not satisfied" << std::endl;
        return 1; // error code or exception to be raised;		
    }
    
    if (nPoints < 3)
    {
        std::cout << "fitQuadratic - condition: \"nPoints >= 3\" not satisfied" << std::endl;
        return 2; // error code or exception to be raised;		
    }
    
    double F0 = 0.;
    double F1 = 0.;
    double F2 = 0.;
    double F3 = 0.;
    double F4 = 0.;
    
    double Y1 = 0.;
    double Y2 = 0.;
    double Y3 = 0.;
    
    double var = res*res;
    
    for(int i = 0; i < nPoints; i++)
    {
        F0 += 1.;
        F1 += x[i];
        F2 += x[i]*x[i];
        F3 += pow(x[i],3);
        F4 += pow(x[i],4);
        
        Y1 += y[i];
        Y2 += y[i]*x[i];
        Y3 += y[i]*x[i]*x[i];
    }
    
    if(withres) 
    {
        F0 /= var;
        F1 /= var;
        F2 /= var;
        F3 /= var;
        F4 /= var;
    
        Y1 /= var;
        Y2 /= var;
        Y3 /= var;
    }
    
    double F11 = F2*F4 - F3*F3;
    double F12 = F1*F4 - F2*F3;
    double F13 = F1*F3 - F2*F2;
    double F21 = F1*F4 - F2*F3;
    double F22 = F0*F4 - F2*F2;
    double F23 = F0*F3 - F1*F2;
    double F31 = F1*F3 - F2*F2;
    double F32 = F0*F3 - F1*F2;
    double F33 = F0*F2 - F1*F1;
    
    double det = F0 * F11 - F1 * F12 + F2 * F13;
    
    if(det == 0.)
    {
        std::cout << "fitQuadratic - condition: \"det == 0.\" not satisfied" << std::endl;
        return 3;
    }
    
    a =  F11 * Y1 - F12 * Y2 + F13 * Y3;
    b = -F21 * Y1 + F22 * Y2 - F23 * Y3;
    c =  F31 * Y1 - F32 * Y2 + F33 * Y3;
    
    a /= det;
    b /= det;
    c /= det;
    
    cov[0][0] =  F11/det;
    cov[0][1] = -F12/det;
    cov[0][2] =  F13/det;
    cov[1][0] = -F21/det;
    cov[1][1] =  F22/det;
    cov[1][2] = -F23/det;
    cov[2][0] =  F31/det;
    cov[2][1] = -F32/det;
    cov[2][2] =  F33/det;
    
    chi2 = 0.0;
    
    for(int i = 0; i < nPoints; i++)
    {
        chi2 += (y[i] - a - b * x[i] - c * x[i] * x[i]) * (y[i] - a - b * x[i] - c * x[i] * x[i]);    
    }
    
    if(withres) 
        chi2 /= var;
    
    chi2 /= nPoints;
    
    return 0;
}

void evalpar(double a, double b, double c, double** cov1, 
        double &xc, double &yc, double &r, double** cov2)
{
    r = 0.5/c;
    
    xc = -0.5 * b / c;
    
    yc = a - b * c + 0.5/c;
    
}

double evalH(double z0, double y0, double z1, double yc)
{
    if((z1-z0) * (yc-y0) > 0)
        return -1.;
    else
        return 1.;
}

double evalZ(double* z)
{
    if(z[1] - z[0] > 0)
        return 1.;
    else
        return -1.;
}

int evalY(double h, int n, double* z, double* y, double zc, double yc, double R)
{
    for(int i = 0; i < n; i++)
    {
        if(TMath::Abs(z[i]-zc) > R)
        {
            std::cout << "evalY        - condition: \"TMath::Abs(z[i]-zc) > R.\" not satisfied" << std::endl;
            return 1;
        }
        y[i] = yc + h * sqrt(R*R - (z[i]-zc)*(z[i]-zc));
    }
    return 0;
}

void evalRho(double h, int n, double* z, double* y, double zc, double yc, double R, double* rho)
{
    if(n == 0)
        return;
    
    double cos =  h * (y[0] - yc)/R;
    double sin = -h * (z[0] - zc)/R;
    
    for(int i = 0; i < n; i++)
    {
        rho[i] = z[i] * cos + y[i] * sin;
    }
}

int fitP(int nPoints, double *x, double *y, bool withres, double res, 
        double &r, double &xc, double &yc, double** cov, double& chi2)
{
    double* cov_abc[3];
    cov_abc[0] = new double[3];
    cov_abc[1] = new double[3];
    cov_abc[2] = new double[3];
    double a, b, c;
    
    int retval  = fitQuadratic(nPoints, x, y, withres, res, a, b, c, cov_abc, chi2);
    
    evalpar(a,b,c,cov_abc,xc,yc,r,cov);
    
    delete cov_abc[0];
    delete cov_abc[1];
    delete cov_abc[2];
    
    return retval;
}

int fitC(int nPoints, double *x, double *y, bool withres, double res, 
        double &r, double &xc, double &yc, double** cov, double& chi2)
{
    double errr;
    
    int retval = fitCircle(nPoints, x, y, xc, yc, r, errr, chi2);
    
    cov[0][0] = errr;
    cov[0][1] = -999;
    cov[0][2] = -999;
    cov[1][0] = -999;
    cov[1][1] = -999;
    cov[1][2] = -999;
    cov[2][0] = -999;
    cov[2][1] = -999;
    cov[2][2] = -999;
    
    return retval;
}

int fit(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z_x, std::vector<double>& z_y, bool withres, double res, 
        int& nPointsx, int& nPointsy, double &r, double &zc, double &yc, double** cov, double &a, double &b, double** cov_ab, double &chi2)
{    
    double chi2_1 = 0.;
    double chi2_2 = 0.;

    nPointsy = std::min(int(y.size()),max_hits_for_y_fit);

    nPointsx = std::min(int(x.size()),max_hits_for_x_fit);
        
    if(nPointsx >= 2)
    {
        nPointsx = 2;

        while(nPointsx < std::min(int(x.size()),max_hits_for_x_fit) && (z_x[nPointsx-1]-z_x[nPointsx-2])*(z_x[1]-z_x[0]) > 0.)
            nPointsx++;
    }

    //fitP(nPoints, &(z_y[0]), &(y[0]), false, -1, r, zc, yc, cov, chi2_1);
    int retval = fitC(nPointsy, &(z_y[0]), &(y[0]), withres, res, r, zc, yc, cov, chi2_1);

    std::vector<double> newy(x.size());
    std::vector<double> rho(x.size());

    double h = evalH(z_y[0], y[0], z_y[1], yc);
    
    retval += 10 * evalY(h, nPointsx, &(z_x[0]), &(newy[0]), zc, yc, r);

    evalRho(h, nPointsx, &(z_x[0]), &(newy[0]), zc, yc, r, &(rho[0]));

    retval += 100 * fitLinear(nPointsx, &(rho[0]), &(x[0]), withres, res, a, b, cov_ab, chi2_2);
    
    chi2 = chi2_1 + chi2_2;
    
    return retval;
}

void select_hits(int hit_n, double* hit_x, double* hit_y, double* hit_z, bool* hit_hor,
        std::vector<double>& x, std::vector<double>& y, 
        std::vector<double>& z_x, std::vector<double>& z_y)
{
    x.clear();
    y.clear();
    z_x.clear();
    z_y.clear();
    
    for(int i = 0; i < hit_n; i++)
    {
        if(hit_hor[i])
        {
            y.push_back(hit_y[i]);
            z_y.push_back(hit_z[i]);
        }
        else
        {
            x.push_back(hit_x[i]);
            z_x.push_back(hit_z[i]);
        }
    }
}

void init(const char* fname)
{    
    f = new TFile(fname);
    
    t = (TTree*) f->Get("hits");
    
    mom = new TLorentzVector();
    vtx = new TLorentzVector();
    
    mom_last_in = new TVector3();
    pos_last_in = new TLorentzVector();
    
    t->SetBranchAddress("hit_n",&hit_n);
    t->SetBranchAddress("hit_x",hit_x);
    t->SetBranchAddress("hit_y",hit_y);
    t->SetBranchAddress("hit_z",hit_z);
    t->SetBranchAddress("hit_t",hit_t);
    t->SetBranchAddress("hit_e",hit_e);
    t->SetBranchAddress("hit_hor",hit_hor);
    t->SetBranchAddress("hit_det",&hit_det);
    t->SetBranchAddress("EventId",&eventid);
    t->SetBranchAddress("mom",&mom);
    t->SetBranchAddress("vtx",&vtx);
    t->SetBranchAddress("mom_last_in",&mom_last_in);
    t->SetBranchAddress("pos_last_in",&pos_last_in);
    t->SetBranchAddress("exiting",&exiting);
    
    initialized = true;
}

void reco(const char* ftxtname)
{
    if(!initialized)
    {
        std::cout << "not initialized" << std::endl;
        return;
    }
    
    double zc, yc, r;
    double a, b;
    double chi2;
    double* cov[3];
    double* cov_ab[3];
    int used_x_hits;
    int used_y_hits;
    int retval;
    
    cov[0] = new double[3];
    cov[1] = new double[3];
    cov[2] = new double[3];

    cov_ab[0] = new double[2];
    cov_ab[1] = new double[2];
            
    //std::ofstream ftxt("fitpar.txt");  
    std::ofstream ftxt(ftxtname);    
        
    ftxt << "id/I:retval:nhits_x:nhits_y:n_used_hits_x:n_used_hits_y:" << 
            "r/D:ptreco:zc:yc:a:b:chi2:cov_aa:cov_ab:cov_ba:cov_bb:" <<
            "cov_00:cov_01:cov_02:cov_10:cov_11:cov_12:cov_20:cov_21:cov_22:" <<
            "x:y:z:px:py:pz:tantheta:pt:p:x_last_in:y_last_in:z_last_in:px_last_in:py_last_in:pz_last_in:exiting" << std::endl;
    
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z_x;
    std::vector<double> z_y;
    
    bool withres = false;
    double res = 0.2; // mm
    
    for(int i = 0; i < t->GetEntries(); i++)
    {
        t->GetEntry(i);
        
        select_hits(hit_n, hit_x, hit_y, hit_z, hit_hor, x, y, z_x, z_y);
        
        retval = fit(x, y, z_x, z_y, withres, res, used_x_hits, used_y_hits, r, zc, yc, cov, a, b, cov_ab,chi2);
        
        ftxt << eventid << " "
                << retval << " "
                << x.size() << " "
                << y.size() << " "
                << used_x_hits << " "
                << used_y_hits << " "
                << r << " "
                << 0.6*0.2998*r << " "
                << zc << " "
                << yc << " "
                << a << " "
                << b*evalZ(&(z_x[0])) << " "
                << chi2 << " "
                << cov_ab[0][0] << " "
                << cov_ab[0][1] << " "
                << cov_ab[1][0] << " "
                << cov_ab[1][1] << " "
                << cov[0][0] << " "
                << cov[0][1] << " "
                << cov[0][2] << " "
                << cov[1][0] << " "
                << cov[1][1] << " "
                << cov[1][2] << " "
                << cov[2][0] << " "
                << cov[2][1] << " "
                << cov[2][2] << " "
                << vtx->X() << " "
                << vtx->Y() << " "
                << vtx->Z() << " "
                << mom->X() << " "
                << mom->Y() << " "
                << mom->Z() << " "
                << mom->X()/sqrt(mom->Z()*mom->Z() + mom->Y()*mom->Y()) << " "
                << sqrt(mom->Z()*mom->Z() + mom->Y()*mom->Y()) << " "
                << sqrt(mom->Z()*mom->Z() + mom->Y()*mom->Y() + mom->X()*mom->X()) << " "
                << pos_last_in->X() << " "
                << pos_last_in->Y() << " "
                << pos_last_in->Z() << " "
                << mom_last_in->X() << " "
                << mom_last_in->Y() << " "
                << mom_last_in->Z() << " "
                << exiting << " "
                << std::endl;
    }
    
    delete cov[0];
    delete cov[1];
    delete cov[2];
    delete cov_ab[0];
    delete cov_ab[1];
        
    ftxt.close();
}

void draw(int i, bool zoom = false)
{
    if(!initialized)
    {
        std::cout << "not initialized" << std::endl;
        return;
    }
    
    double zc, yc, r;
    double a, b;
    double chi2;
    double* cov[3];
    double* cov_ab[3];
    int used_x_hits;
    int used_y_hits;
    int retval;
    
    cov[0] = new double[3];
    cov[1] = new double[3];
    cov[2] = new double[3];

    cov_ab[0] = new double[2];
    cov_ab[1] = new double[2];
    
    bool withres = false;
    double res = 0.2; // mm
    
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z_x;
    std::vector<double> z_y;

    t->GetEntry(i);
    
    select_hits(hit_n, hit_x, hit_y, hit_z, hit_hor, x, y, z_x, z_y);

    retval = fit(x, y, z_x, z_y, withres, res, used_x_hits, used_y_hits, r, zc, yc, cov, a, b, cov_ab,chi2);
    
    
        
    std::cout << "id/I:retval:nhits_x:nhits_y:n_used_hits_x:n_used_hits_y:" << 
            "r/D:ptreco:zc:yc:a:b:chi2:cov_aa:cov_ab:cov_ba:cov_bb:" <<
            "cov_00:cov_01:cov_02:cov_10:cov_11:cov_12:cov_20:cov_21:cov_22:" <<
            "x:y:z:px:py:pz:tantheta:pt:p:x_last_in:y_last_in:z_last_in:px_last_in:py_last_in:pz_last_in:exiting" << std::endl;    
        
    std::cout << eventid << " "
            << retval << " "
            << x.size() << " "
            << y.size() << " "
            << used_x_hits << " "
            << used_y_hits << " "
            << r << " "
            << 0.6*0.2998*r << " "
            << zc << " "
            << yc << " "
            << a << " "
            << b*evalZ(&(z_x[0])) << " "
            << chi2 << " "
            << cov_ab[0][0] << " "
            << cov_ab[0][1] << " "
            << cov_ab[1][0] << " "
            << cov_ab[1][1] << " "
            << cov[0][0] << " "
            << cov[0][1] << " "
            << cov[0][2] << " "
            << cov[1][0] << " "
            << cov[1][1] << " "
            << cov[1][2] << " "
            << cov[2][0] << " "
            << cov[2][1] << " "
            << cov[2][2] << " "
            << vtx->X() << " "
            << vtx->Y() << " "
            << vtx->Z() << " "
            << mom->X() << " "
            << mom->Y() << " "
            << mom->Z() << " "
            << mom->X()/sqrt(mom->Z()*mom->Z() + mom->Y()*mom->Y()) << " "
            << sqrt(mom->Z()*mom->Z() + mom->Y()*mom->Y()) << " "
            << sqrt(mom->Z()*mom->Z() + mom->Y()*mom->Y() + mom->X()*mom->X()) << " "
            << pos_last_in->X() << " "
            << pos_last_in->Y() << " "
            << pos_last_in->Z() << " "
            << mom_last_in->X() << " "
            << mom_last_in->Y() << " "
            << mom_last_in->Z() << " "
            << exiting << " "
            << std::endl;
    
    double kloe_x = 0.;
    double kloe_y = 4330.96000 - 2088.00000;
    double kloe_z = 7396.16250 + 3244.00000;
    
    double dkloe_x = 2500.;
    double dkloe_y = 2500.;
    double dkloe_z = 2500.;
    
    double kloe_rad = 2000.;
    double kloe_hdx = 1690.;

    TCanvas* can = new TCanvas("c","",1000,500);
    can->Divide(2,1);

    if(zoom)
    {
        double xmin1 = *std::min_element(z_y.begin(), z_y.begin() + used_y_hits);
        double xmax1 = *std::max_element(z_y.begin(), z_y.begin() + used_y_hits);
        double ymin1 = *std::min_element(y.begin(), y.begin() + used_y_hits);
        double ymax1 = *std::max_element(y.begin(), y.begin() + used_y_hits);
        
        double xmin2 = *std::min_element(z_x.begin(), z_x.begin() + used_x_hits);
        double xmax2 = *std::max_element(z_x.begin(), z_x.begin() + used_x_hits);
        double ymin2 = *std::min_element(x.begin(), x.begin() + used_x_hits);
        double ymax2 = *std::max_element(x.begin(), x.begin() + used_x_hits);

        double dx1 = xmax1 - xmin1;
        double dy1 = ymax1 - ymin1;
        double xm1 = 0.5* (xmin1 + xmax1);
        double ym1 = 0.5* (ymin1 + ymax1);

        double dx2 = xmax2 - xmin2;
        double dy2 = ymax2 - ymin2;
        double xm2 = 0.5* (xmin2 + xmax2);
        double ym2 = 0.5* (ymin2 + ymax2);

        if(dx1 > dy1)
        {
            dy1 = dx1;
        }
        else
        {
            dx1 = dy1;
        }

        if(dx2 > dy2)
        {
            dy2 = dx2;
        }
        else
        {
            dx2 = dy2;
        }
        
        can->cd(1)->DrawFrame(xm1 - 0.55 * dx1, ym1 - 0.55 * dy1, xm1 + 0.55 * dx1, ym1 + 0.55 * dy1,";Z (mm); Y(mm)");
        can->cd(2)->DrawFrame(xm2 - 0.55 * dx2, ym2 - 0.55 * dy2, xm2 + 0.55 * dx2, ym2 + 0.55 * dy2,";Z (mm); X(mm)");
    }
    else
    {
        can->cd(1)->DrawFrame(kloe_z - dkloe_z, kloe_y - dkloe_y, kloe_z + dkloe_z, kloe_y + dkloe_y);
        can->cd(2)->DrawFrame(kloe_z - dkloe_z, kloe_x - dkloe_x, kloe_z + dkloe_z, kloe_x + dkloe_x);
    }

    can->cd(1);
    TEllipse* kloe_yz = new TEllipse(kloe_z, kloe_y, kloe_rad);
    kloe_yz->SetLineColor(kBlue);
    kloe_yz->SetFillStyle(0);
    kloe_yz->Draw();

    can->cd(2);
    TBox* kloe_xz = new TBox(kloe_z - kloe_rad, -kloe_hdx, kloe_z + kloe_rad, kloe_hdx);
    kloe_xz->SetLineColor(kBlue);
    kloe_xz->SetFillStyle(0);
    kloe_xz->Draw();

    can->cd(1);
    TGraph* gr_yz = new TGraph(used_y_hits, &(z_y[0]), &(y[0]));
    gr_yz->SetMarkerStyle(2);
    gr_yz->Draw("P");

    TEllipse* fit_yz = new TEllipse(zc, yc, r);
    fit_yz->SetLineColor(kRed);
    fit_yz->SetFillStyle(0);
    fit_yz->Draw();

    can->cd(2);
    TGraph* gr_xz = new TGraph(used_x_hits, &(z_x[0]), &(x[0]));
    gr_xz->SetMarkerStyle(2);
    gr_xz->Draw("P");
    
    TLine* fit_xz = new TLine(vtx->Z(),vtx->X(),kloe_z + dkloe_z, vtx->X() + b * (kloe_z + dkloe_z - vtx->Z()));
    fit_xz->SetLineColor(kRed);
    fit_xz->Draw("same");

    double pt = sqrt(mom->Z() * mom->Z() + mom->Y() * mom->Y());
    double R = pt/(0.2998*0.6);
    double s = R/pt;

    double zt =  mom->Y() * s + vtx->Z();
    double yt = -mom->Z() * s + vtx->Y();  

    can->cd(1);
    TEllipse* path_yz = new TEllipse(zt,yt,R);
    path_yz->SetLineColor(6);
    path_yz->SetFillStyle(0);
    path_yz->SetLineStyle(5);
    path_yz->Draw("same"); 

    can->cd(2);
    TLine* path_xz = new TLine(vtx->Z(),vtx->X(),kloe_z + dkloe_z, vtx->X() + mom->X()/mom->Z() * (kloe_z + dkloe_z - vtx->Z()));
    path_xz->SetLineStyle(5);
    path_xz->SetLineColor(6);
    path_xz->Draw("same");

    can->SaveAs(TString::Format("displays/Ev_%d.png",eventid));
    
    delete cov[0];
    delete cov[1];
    delete cov[2];
    delete cov_ab[0];
    delete cov_ab[1];
}