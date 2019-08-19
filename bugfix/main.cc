#include <TEveManager.h>
#include <TApplication.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <TMath.h>

#include <vector>
#include <iomanip>
#include <iostream>
#include <cmath>

#include "struct.h"

void dump_dg(digit* dg)
{
  std::cout << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::setw(15) << "t" << std::setw(15) << "de" << std::setw(15) << std::endl;
  
  std::cout << std::setw(15) << dg->x << std::setw(15) << dg->y << std::setw(15) << dg->z << std::setw(15) << dg->t << std::setw(15) << dg->de << std::endl;
  
  std::cout << std::setw(15) << "hor" << std::setw(15) << "n_hits" << std::setw(65) << "det" << std::endl;
  
  std::cout << std::setw(15) << dg->hor << std::setw(15) << dg->hindex.size() << std::setw(65) << dg->det.data() << std::endl<< std::endl;
}

void dump_tr(track* tr)
{
  std::cout << "****** track parameter *******" << std::endl;
  
  std::cout << std::setw(15) << "tid" << std::setw(15) << "yC" << std::setw(15) << "zC" << std::setw(15) << "r" << std::setw(15) << "a" << std::setw(15) << "b" << std::setw(15) << "h" << std::endl;
  
  std::cout << std::setw(15) << tr->tid << std::setw(15) << tr->yc << std::setw(15) << tr->zc << std::setw(15) << tr->r << std::setw(15) << tr->a << std::setw(15) << tr->b << std::setw(15) << tr->h << std::endl << std::endl;
  
  std::cout << std::setw(15) << "x0" << std::setw(15) << "y0" << std::setw(15) << "z0" << std::setw(15) << "t0" << std::setw(15) << "n_digits" << std::endl;
  
    std::cout << std::setw(15) << tr->x0 << std::setw(15) << tr->y0 << std::setw(15) << tr->z0 << std::setw(15) << tr->t0 << std::setw(15) << tr->digits.size() << std::endl << std::endl;
    
  std::cout << std::setw(15) << "ret_ln" << std::setw(15) << "chi2_ln" << std::setw(15) << "ret_cr" << std::setw(15) << "chi2_cr" << std::endl;
  
  std::cout << std::setw(15) << tr->ret_ln << std::setw(15) << tr->chi2_ln << std::setw(15) << tr->ret_cr << std::setw(15) << tr->chi2_cr << std::endl;
  
  std::cout << "*************************" << std::endl;
  
  for(unsigned int i = 0; i < tr->digits.size(); i++)
  {
    dump_dg(&(tr->digits[i]));
  }
}


int main(int argc, char* argv[]) {
  
  int count = 0;
  
  TFile f("tmp.root");

  TTree* tev = (TTree*) f.Get("tEvent");
  
  event* evt = new event;
  
  tev->SetBranchAddress("event",&evt);
  
  TDatabasePDG pdgdb;
  pdgdb.ReadPDGTable();
  
  const int nev = tev->GetEntries();
  
  for(int i = 0; i < nev; i++)
  {  
    tev->GetEntry(i);      
    
    for(int k = 0; k < evt->particles.size(); k++)
    {    
      if(evt->particles[k].primary == 1 && evt->particles[k].pdg == 13)
      {
        if(evt->particles[k].tr.ret_ln > 3 || evt->particles[k].tr.ret_cr > 3 || evt->particles[k].tr.ret_ln < -3 || evt->particles[k].tr.ret_cr < -2 || (evt->particles[k].tr.ret_ln == -1 && evt->particles[k].tr.ret_cr == -1))
        {
          std::cout << "entry: " << i << "\tret_ln or ret_cr error -> " << evt->particles[k].tr.ret_ln << " " << evt->particles[k].tr.ret_cr << std::endl;
          count++;
        }
      
        if(evt->particles[k].has_track == true && evt->particles[k].tr.ret_ln == 0 && evt->particles[k].tr.ret_cr == 0)
        {
        
          // problematic events. Check why
          if(evt->particles[k].tr.digits.size() == 0)
          {
            std::cout << "entry: " << i << "\tdigit vector has zero size" << std::endl;
            count++;
            //dump_tr(&(evt->particles[k].tr));            
          }
          else if(isnan(evt->particles[k].tr.x0))
          {
            std::cout << "entry: " << i << "\tx0 is nan -> x0 = " << evt->particles[k].tr.x0 << std::endl;
            count++;
            //dump_tr(&(evt->particles[k].tr));
          }
          else if(evt->particles[k].tr.t0 < 1.)
          {
            std::cout << "entry: " << i << "\tt0 is less then 1 ns -> t0 = " << evt->particles[k].tr.t0 << std::endl;
            count++;
            //dump_tr(&(evt->particles[k].tr));
          }
          continue;
        }
      }
    }
  }
  std::cout << "Total: " << count << " errors" << std::endl;
}


