#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TParameter.h>

#include <iostream>
#include <vector>
#include <map>

#include "/wd/dune-it/ext_bkg/kloe-simu/src/display.cpp"

void erase_element(int val, std::vector<int>& vec)
{
  std::vector<int>::iterator ite = std::find(vec.begin(),vec.end(),val);
  if(ite == vec.end())
  {
    std::cout << "problem: " << val << " not found in vec (size: " << vec.size() << ")" << std::endl;
    for(std::vector<int>::iterator it = vec.begin(); it != vec.end(); ++it)
    {
      std::cout << distance(vec.begin(), it) << " " << *it << std::endl;
    }
    exit(1);
  } 
  else
  {
    vec.erase(ite);
  }
}

void processMtrVct(std::vector<int>& regMtrX, std::vector<int>& reg0trX, std::map<int,int>& regionsX, TH1D& hmultX, const int minreg)
{
    while(regMtrX.size() != 0)
    {
      std::map<int,int>::iterator this_el = regionsX.find(regMtrX.front());
      std::map<int,int>::iterator next_el = std::next(this_el);
      std::map<int,int>::iterator prev_el = std::prev(this_el);
      std::map<int,int>::iterator next2_el = std::next(this_el,2);
      std::map<int,int>::iterator prev2_el = std::prev(this_el,2);
      bool mantain = false;
      
      if(next_el != regionsX.end() && this_el != regionsX.begin())
      {
        if(this_el->second < minreg && hmultX.GetBinContent(prev_el->first) == 1 && hmultX.GetBinContent(next_el->first) == 1)
        {
          for(int k = 0; k < this_el->second; k++)
          {
            hmultX.SetBinContent(this_el->first+k,1);
          }
          prev_el->second += this_el->second + next_el->second;
          regionsX.erase(this_el);
          regionsX.erase(next_el);
          regMtrX.erase(regMtrX.begin());        
          
          continue;
        }
      }
      
      if(next_el != regionsX.end())
      {        
        if(next_el->second < minreg)
        {
          mantain = true;
          
          if(next2_el != regionsX.end())
          {
            if(hmultX.GetBinContent(next2_el->first) == 2)
            {
              if(hmultX.GetBinContent(next_el->first) == 0)
                erase_element(next_el->first, reg0trX);
              for(int k = 0; k < next_el->second; k++)
              {
                hmultX.SetBinContent(next_el->first+k,2);
              }
              this_el->second += next_el->second + next2_el->second;
              erase_element(next2_el->first, regMtrX);
              regionsX.erase(next_el);
              regionsX.erase(next2_el);
            }
            else
            {
              if(hmultX.GetBinContent(next_el->first) == 0)
                erase_element(next_el->first, reg0trX);
              for(int k = 0; k < next_el->second; k++)
              {
                hmultX.SetBinContent(next_el->first+k,2);
              }
              this_el->second += next_el->second;
              regionsX.erase(next_el);
            }
          }
          else
          {
            if(hmultX.GetBinContent(next_el->first) == 0)
              erase_element(next_el->first, reg0trX);
            for(int k = 0; k < next_el->second; k++)
            {
              hmultX.SetBinContent(next_el->first+k,2);
            }
            this_el->second += next_el->second;
            regionsX.erase(next_el);
          }
        }
      }
      
      if(this_el != regionsX.begin())
      {
        if(prev_el->second < minreg)
        {
          mantain = true;
          
          if(prev_el != regionsX.begin())
          {
            if(hmultX.GetBinContent(prev2_el->first) == 2)
            {
              if(hmultX.GetBinContent(prev_el->first) == 0)
                erase_element(prev_el->first, reg0trX);
              for(int k = 0; k < prev_el->second; k++)
              {
                hmultX.SetBinContent(prev_el->first+k,2);
              }
              prev2_el->second += this_el->second + prev_el->second;
              regionsX.erase(prev_el);
              regionsX.erase(this_el);
              regMtrX.front() = prev2_el->first;
              
            }
            else
            {
              if(hmultX.GetBinContent(prev_el->first) == 0)
                erase_element(prev_el->first, reg0trX);
              for(int k = 0; k < prev_el->second; k++)
              {
                hmultX.SetBinContent(prev_el->first+k,2);
              }
              prev_el->second += this_el->second;
              regionsX.erase(this_el);
              regMtrX.front() = prev_el->first;
            }
          }
          else
          {
            if(hmultX.GetBinContent(prev_el->first) == 0)
              erase_element(prev_el->first, reg0trX);
            for(int k = 0; k < prev_el->second; k++)
            {
              hmultX.SetBinContent(prev_el->first+k,2);
            }
            prev_el->second += this_el->second;
            regionsX.erase(this_el);
            regMtrX.front() = prev_el->first;
          }
        }
      }
      if(mantain == false)
      {
        regMtrX.erase(regMtrX.begin());
      }
    }
}

void process0trVct(std::vector<int>& reg0trX, std::map<int,int>& regionsX, TH1D& hmultX, const int minreg)
{
    while(reg0trX.size() != 0)
    {
      std::map<int,int>::iterator this_el = regionsX.find(reg0trX.front());
      std::map<int,int>::iterator next_el = std::next(this_el);
      std::map<int,int>::iterator prev_el = std::prev(this_el);
      std::map<int,int>::iterator next2_el = std::next(this_el,2);
      std::map<int,int>::iterator prev2_el = std::prev(this_el,2);
      bool mantain = false;
      
      
      if(next_el != regionsX.end() && this_el != regionsX.begin())
      {
        if(this_el->second < minreg && hmultX.GetBinContent(prev_el->first) == 1 && hmultX.GetBinContent(next_el->first) == 1)
        {
          for(int k = 0; k < this_el->second; k++)
          {
            hmultX.SetBinContent(this_el->first+k,1);
          }
          prev_el->second += this_el->second + next_el->second;
          regionsX.erase(this_el);
          regionsX.erase(next_el); 
          reg0trX.erase(reg0trX.begin());
          continue;
        }
      }
      
      if(next_el != regionsX.end())
      {        
        if(next_el->second < minreg)
        {
          mantain = true;
          
          if(next2_el != regionsX.end())
          {
            if(hmultX.GetBinContent(next2_el->first) == 0)
            {
              for(int k = 0; k < next_el->second; k++)
              {
                hmultX.SetBinContent(next_el->first+k,0);
              }
              this_el->second += next_el->second + next2_el->second;
              erase_element(next2_el->first, reg0trX);
              regionsX.erase(next_el);
              regionsX.erase(next2_el);
              
            }
            else
            {
              for(int k = 0; k < next_el->second; k++)
              {
                hmultX.SetBinContent(next_el->first+k,hmultX.GetBinContent(next2_el->first));
              }
              next_el->second += next2_el->second;
              regionsX.erase(next2_el);   
            }
          }
          else
          {
            for(int k = 0; k < next_el->second; k++)
            {
              hmultX.SetBinContent(next_el->first+k,0);
            }
            this_el->second += next_el->second;
            regionsX.erase(next_el);
          }
        }
      }
      
      if(this_el != regionsX.begin())
      {        
        if(prev_el->second < minreg)
        {
          mantain = true;
          
          if(prev_el != regionsX.begin())
          {
            if(hmultX.GetBinContent(prev2_el->first) == 0)
            {
              for(int k = 0; k < prev_el->second; k++)
              {
                hmultX.SetBinContent(prev_el->first+k,0);
              }
              prev2_el->second += prev_el->second + this_el->second;
              regionsX.erase(prev_el);
              regionsX.erase(this_el);  
              reg0trX.front() = prev2_el->first;            
            }
            else
            {
              for(int k = 0; k < prev_el->second; k++)
              {
                hmultX.SetBinContent(prev_el->first+k,hmultX.GetBinContent(prev2_el->first));
              }
              prev2_el->second += prev_el->second;
              regionsX.erase(prev_el);
              reg0trX.front() = prev_el->first;
            }
          }
          else
          {
            for(int k = 0; k < prev_el->second; k++)
            {
              hmultX.SetBinContent(prev_el->first+k,0);
            }
            prev_el->second += this_el->second;
            regionsX.erase(this_el);
            reg0trX.front() = prev_el->first;
          }
        }
      }
      if(mantain == false)
      {
        reg0trX.erase(reg0trX.begin());
      }
    }
}

void findvtx()
{
  gROOT->SetBatch();
  TCanvas c;

// root [0] TGeoManager::Import("../geo/nd_hall_kloe_empty.gdml")
// (TGeoManager *) 0x3037200
// root [12] gGeoManager->cd("volWorld/rockBox_lv_0/volDetEnclosure_0/volKLOE_0")
// (bool) true
// root [16] double vorigin[3]
// (double [3]) { 0.0000000, 0.0000000, 0.0000000 }
// root [17] double vmaster[3]
// (double [3]) { 0.0000000, 0.0000000, 0.0000000 }
// root [18] gGeoManager->LocalToMaster(vorigin,vmaster)
// root [19] vmaster
// (double [3]) { 0.0000000, -238.47300, 2391.0000 }
// root [7] v = gGeoManager->GetVolume("volSTTFULL")
// (TGeoVolume *) 0x17fc5b10
// root [9] TGeoTube* tb = (TGeoTube*) v->GetShape()
// root [11] tb->GetRmin()
// (double) 0.0000000
// root [12] tb->GetRmax()
// (double) 200.00000
// root [13] tb->GetDz()
// (double) 169.00000

  double kloe_center[] = { 0.0000000, -2384.7300, 23910.000 };
  double kloe_size[] = { 2 * 1690.00000, 2 * 2000.00000, 2 * 2000.00000};

  TFile f("../files/reco/numu_internal_10k.0.reco.root");
  
  TTree* tDigit = (TTree*) f.Get("tDigit");
  TTree* tMC = (TTree*) f.Get("EDepSimEvents");
  TGeoManager* g = (TGeoManager*) f.Get("EDepSimGeometry");
  
  TString path_prefix = "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volKLOE_PV_0/MagIntVol_volume_PV_0/volSTTFULL_PV_0/";
  TGeoVolume* v = g->FindVolumeFast("volSTTFULL_PV");
  
  double origin[3];
  double master[3];
  double last_z = 0.;
  double dz;
  std::vector<double> binning;
  bool is_first = true;
  
  // assuming they are order by Z position
  for(int i = 0; i < v->GetNdaughters(); i++)
  {
    TString name = v->GetNode(i)->GetName();
    
    if(name.Contains("sttmod") || name.Contains("volfrontST"))
    {
      TString path = path_prefix + name;
      g->cd(path.Data());
      g->LocalToMaster(origin,master);
      TGeoBBox* b = (TGeoBBox*) v->GetNode(i)->GetVolume()->GetShape();
      dz = b->GetDX();
    
      if(is_first)
      {
        is_first = false;
        binning.push_back(master[2] - dz);
      }
      else
      {
        binning.push_back(0.5 * (last_z + master[2] - dz));
      }
      last_z = master[2] + dz;
    }
  }
  binning.push_back(last_z);
  /*
  for(unsigned int i = 0; i < binning.size(); i++)
  {
    std::cout << i << " " << binning.data()[i] << std::endl;
  }*/
    
  std::vector<digit>* digits = new std::vector<digit>;
  TG4Event* ev = new TG4Event;
  
  tDigit->SetBranchAddress("Stt",&digits);
  tMC->SetBranchAddress("Event",&ev);
  
  TFile fout("vtx.root","RECREATE");
    
  TH1D hrmsX("hrmsX","rmsX;Z (mm); rmsX (mm)",binning.size()-1,binning.data());
  TH1D hrmsY("hrmsY","rmsY;Z (mm); rmsY (mm)",binning.size()-1,binning.data());
  hrmsY.SetLineColor(kRed);
    
  TH1D hmeanX("hmeanX","meanX;Z (mm); meanX (mm)",binning.size()-1,binning.data());
  TH1D hmeanY("hmeanY","meanY;Z (mm); meanY (mm)",binning.size()-1,binning.data());
  hmeanY.SetLineColor(kRed);
    
  TH1I hnX("hnX","nX;Z (mm); nX",binning.size()-1,binning.data());
  TH1I hnY("hnY","nY;Z (mm); nY",binning.size()-1,binning.data());
  hnY.SetLineColor(kRed);
  
  TH1D h0trX("h0trX","h0trX;Z (mm); rmsX (mm)",binning.size()-1,binning.data());
  TH1D h1trX("h1trX","h1trX;Z (mm); rmsX (mm)",binning.size()-1,binning.data());
  TH1D hmtrX("hmtrX","hmtrX;Z (mm); rmsX (mm)",binning.size()-1,binning.data());
  
  TH1D h0trY("h0trY","h0trY;Z (mm); rmsY (mm)",binning.size()-1,binning.data());
  TH1D h1trY("h1trY","h1trY;Z (mm); rmsY (mm)",binning.size()-1,binning.data());
  TH1D hmtrY("hmtrY","hmtrY;Z (mm); rmsY (mm)",binning.size()-1,binning.data());
  
  TH1D hmultX("hmultX","hmultX;Z (mm); multipliciy",binning.size()-1,binning.data());
  TH1D hmultY("hmultY","hmultY;Z (mm); multipliciy",binning.size()-1,binning.data());
  
  h0trX.SetFillColorAlpha(kRed, 0.15);
  h1trX.SetFillColorAlpha(kBlue, 0.15);
  hmtrX.SetFillColorAlpha(kGreen, 0.15);
  
  h0trY.SetFillColorAlpha(kRed, 0.15);
  h1trY.SetFillColorAlpha(kBlue, 0.15);
  hmtrY.SetFillColorAlpha(kGreen, 0.15);
  
  h0trX.SetLineColor(0);
  h1trX.SetLineColor(0);
  hmtrX.SetLineColor(0);
  
  h0trX.SetLineWidth(0);
  h1trX.SetLineWidth(0);
  hmtrX.SetLineWidth(0);
  
  h0trY.SetLineColor(0);
  h1trY.SetLineColor(0);
  hmtrY.SetLineColor(0);
  
  h0trY.SetLineWidth(0);
  h1trY.SetLineWidth(0);
  hmtrY.SetLineWidth(0);
  
  hrmsX.SetStats(false);
  hrmsY.SetStats(false);
  hmeanX.SetStats(false);
  hmeanY.SetStats(false);
  hnX.SetStats(false);
  hnY.SetStats(false);
  h0trX.SetStats(false);
  h1trX.SetStats(false);
  hmtrX.SetStats(false);
  h0trY.SetStats(false);
  h1trY.SetStats(false);
  hmtrY.SetStats(false);
  hmultX.SetStats(false);
  hmultY.SetStats(false);
  
  TParameter<double> xv("xv",0);
  TParameter<double> yv("yv",0);
  TParameter<double> zv("zv",0);
  
  double mean, mean2, var, rms;
  int n;
  
  init("../files/reco/numu_internal_10k.0.reco.root");
  
  c.Divide(1,2);
  
  int first = 0;
  int last = 100;/*tDigit->GetEntries();
  
  for(int i = first; i < last+1; i++)
  {
    tDigit->GetEntry(i);
    tMC->GetEntry(i);
    
    TString reaction = ev->Primaries.at(0).GetReaction();
    
    //if(reaction.Contains("CC") == false)
    //  continue;
    
    xv.SetVal(ev->Primaries.at(0).GetPosition().X());
    yv.SetVal(ev->Primaries.at(0).GetPosition().Y());
    zv.SetVal(ev->Primaries.at(0).GetPosition().Z());
    
    TDirectoryFile fd(TString::Format("ev_%d",i).Data(),TString::Format("ev_%d",i).Data(),"",&fout);
    
    hrmsX.Reset("ICESM");
    hrmsY.Reset("ICESM");
    hmeanX.Reset("ICESM");
    hmeanY.Reset("ICESM");
    hnX.Reset("ICESM");
    hnY.Reset("ICESM");
    h0trX.Reset("ICESM");
    h1trX.Reset("ICESM");
    hmtrX.Reset("ICESM");
    h0trY.Reset("ICESM");
    h1trY.Reset("ICESM");
    hmtrY.Reset("ICESM");
    
    hrmsX.SetTitle(TString::Format("event: %d",i).Data());
    hrmsY.SetTitle(TString::Format("event: %d",i).Data());
    hmeanX.SetTitle(TString::Format("event: %d",i).Data());
    hmeanY.SetTitle(TString::Format("event: %d",i).Data());
    hnX.SetTitle(TString::Format("event: %d",i).Data());
    hnY.SetTitle(TString::Format("event: %d",i).Data());
    h0trX.SetTitle(TString::Format("event: %d",i).Data());
    h1trX.SetTitle(TString::Format("event: %d",i).Data());
    hmtrX.SetTitle(TString::Format("event: %d",i).Data());
    h0trY.SetTitle(TString::Format("event: %d",i).Data());
    h1trY.SetTitle(TString::Format("event: %d",i).Data());
    hmtrY.SetTitle(TString::Format("event: %d",i).Data());
    
    for(unsigned int j = 0; j < digits->size(); j++)
    {
      if(digits->at(j).hor == 1)
      {
        hrmsY.Fill(digits->at(j).z,digits->at(j).y*digits->at(j).y);
        hmeanY.Fill(digits->at(j).z,digits->at(j).y);
        hnY.Fill(digits->at(j).z);
      }
      else
      {
        hrmsX.Fill(digits->at(j).z,digits->at(j).x*digits->at(j).x);
        hmeanX.Fill(digits->at(j).z,digits->at(j).x);
        hnX.Fill(digits->at(j).z);
      }
    }
    
    for(unsigned int j = 0; j < binning.size()-1; j++)
    {
      mean = -1;
      rms = -1;
      n = hnX.GetBinContent(j+1);
      
      if(n  != 0)
      {
        mean = hmeanX.GetBinContent(j+1)/n;
        mean2 = hrmsX.GetBinContent(j+1)/n;
        var = mean2 - mean * mean;
        rms = sqrt(var);
      }
      if(n>1) n=2;
      
      hrmsX.SetBinContent(j+1,rms);
      hmeanX.SetBinContent(j+1,mean);
      hmultX.SetBinContent(j+1,n);
      
      mean = -1;
      rms = -1;
      n = hnY.GetBinContent(j+1);
      if(n  != 0)
      {
        mean = hmeanY.GetBinContent(j+1)/n;
        mean2 = hrmsY.GetBinContent(j+1)/n;
        var = mean2 - mean * mean;
        rms = sqrt(var);
      }
      if(n>1) n=2;
      
      hrmsY.SetBinContent(j+1,rms);
      hmeanY.SetBinContent(j+1,mean);
      hmultY.SetBinContent(j+1,n);
    }
    
    // filter
    // find wide regions 
    std::map<int,int> regionsX;
    std::map<int,int> regionsY;
    std::vector<int> regMtrX;
    std::vector<int> regMtrY;
    std::vector<int> reg0trX;
    std::vector<int> reg0trY;
    int last_vx = -1;
    int last_vy = -1;
    int last_bx = -1;
    int last_by = -1;
    
    for(unsigned int j = 0; j < binning.size()-1; j++)
    {
      if(hmultX.GetBinContent(j+1) == last_vx)
      {
        regionsX[last_bx]++;
      }
      else
      {
        regionsX[j+1] = 1;
        last_bx = j+1;
        last_vx = hmultX.GetBinContent(j+1);
      }
      
      if(hmultY.GetBinContent(j+1) == last_vy)
      {
        regionsY[last_by]++;
      }
      else
      {
        regionsY[j+1] = 1;
        last_by = j+1;
        last_vy = hmultY.GetBinContent(j+1);
      }
    }
    
    for (std::map<int,int>::iterator it=regionsX.begin(); it!=regionsX.end(); ++it)
    {
      if(hmultX.GetBinContent(it->first) == 2)
      {
        std::vector<int>::iterator ite;
        for (ite=regMtrX.begin(); ite!=regMtrX.end(); ++ite)
        {
          if(regionsX[*ite] < it->second)
            break;
        }
        regMtrX.insert(ite,it->first);
      }
      if(hmultX.GetBinContent(it->first) == 0)
      {
        std::vector<int>::iterator ite;
        for (ite=reg0trX.begin(); ite!=reg0trX.end(); ++ite)
        {
          if(regionsX[*ite] < it->second)
            break;
        }
        reg0trX.insert(ite,it->first);
      }
    }  
    
    for (std::map<int,int>::iterator it=regionsY.begin(); it!=regionsY.end(); ++it)
    {
      if(hmultY.GetBinContent(it->first) == 2)
      {
        std::vector<int>::iterator ite;
        for (ite=regMtrY.begin(); ite!=regMtrY.end(); ++ite)
        {
          if(regionsY[*ite] < it->second)
            break;
        }
        regMtrY.insert(ite,it->first);
      }
      if(hmultY.GetBinContent(it->first) == 0)
      {
        std::vector<int>::iterator ite;
        for (ite=reg0trY.begin(); ite!=reg0trY.end(); ++ite)
        {
          if(regionsY[*ite] < it->second)
            break;
        }
        reg0trY.insert(ite,it->first);
      }
    }
    
    const int minreg = 3;
    
    processMtrVct(regMtrX, reg0trX, regionsX, hmultX, minreg);
    process0trVct(reg0trX, regionsX, hmultX, minreg);
    
    processMtrVct(regMtrY, reg0trY, regionsY, hmultY, minreg);
    process0trVct(reg0trY, regionsY, hmultY, minreg);
    
    /*
    for (std::map<int,int>::iterator it=regionsX.begin(); it!=regionsX.end(); ++it)
    {
      if(it->second < minreg)
      {
        if(hmultX.GetBinContent(std::prev(it)->first) == hmultX.GetBinContent(std::next(it)->first))
        {
          std::prev(it)->second += it->second + std::next(it)->second;
          regionsX.erase(it, std::next(it));
        }
      }
    }
    
    for (std::map<int,int>::iterator it=regionsY.begin(); it!=regionsY.end(); ++it)
    {
      if(it->second < minreg)
      {
        if(hmultY.GetBinContent(std::prev(it)->first) == hmultY.GetBinContent(std::next(it)->first))
        {
          std::prev(it)->second += it->second + std::next(it)->second;
          regionsY.erase(it, std::next(it));
        }
      }
    }*/
    
    for (std::map<int,int>::iterator it=regionsX.begin(); it!=regionsX.end(); ++it)
    {
      for(int k = 0; k < it->second; k++)
      {
        h0trX.SetBinContent(it->first + k, hmultX.GetBinContent(it->first) == 0 ? 1 : 0);
        h1trX.SetBinContent(it->first + k, hmultX.GetBinContent(it->first) == 1 ? 1 : 0);
        hmtrX.SetBinContent(it->first + k, hmultX.GetBinContent(it->first) == 2 ? 1 : 0);
      }
    }
    
    for (std::map<int,int>::iterator it=regionsY.begin(); it!=regionsY.end(); ++it)
    {
      for(int k = 0; k < it->second; k++)
      {
        h0trY.SetBinContent(it->first + k, hmultY.GetBinContent(it->first) == 0 ? 1 : 0);
        h1trY.SetBinContent(it->first + k, hmultY.GetBinContent(it->first) == 1 ? 1 : 0);
        hmtrY.SetBinContent(it->first + k, hmultY.GetBinContent(it->first) == 2 ? 1 : 0);
      }
    }
    
    
    TLine lv1;
    lv1.SetLineStyle(2);
    lv1.SetLineColor(kBlack);
    lv1.SetLineWidth(2);
    lv1.SetX1(zv.GetVal());
    lv1.SetX2(zv.GetVal());
    
    TLine lv2;
    lv2.SetLineStyle(2);
    lv2.SetLineColor(kBlack);
    lv2.SetLineWidth(2);
    lv2.SetX1(zv.GetVal());
    lv2.SetX2(zv.GetVal());
    
    c.cd(1);
    h0trX.Scale(hrmsX.GetMaximum() == 0 ? 1 : hrmsX.GetMaximum());
    h1trX.Scale(hrmsX.GetMaximum() == 0 ? 1 : hrmsX.GetMaximum());
    hmtrX.Scale(hrmsX.GetMaximum() == 0 ? 1 : hrmsX.GetMaximum());
    h0trX.Draw("B");
    h1trX.Draw("BSAME");
    hmtrX.Draw("BSAME");
    hrmsX.Draw("HISTSAME");
    lv1.SetY1(hrmsX.GetMinimum() == -1 ? 0 : hrmsX.GetMinimum());
    lv1.SetY2(hrmsX.GetMaximum() == 0 ? 1 : hrmsX.GetMaximum());
    lv1.Draw();
    
    c.cd(2);
    h0trY.Scale(hrmsY.GetMaximum() == 0 ? 1 : hrmsY.GetMaximum());
    h1trY.Scale(hrmsY.GetMaximum() == 0 ? 1 : hrmsY.GetMaximum());
    hmtrY.Scale(hrmsY.GetMaximum() == 0 ? 1 : hrmsY.GetMaximum());
    h0trY.Draw("B");
    h1trY.Draw("BSAME");
    hmtrY.Draw("BSAME");
    hrmsY.Draw("HISTSAME");
    lv2.SetY1(hrmsY.GetMinimum() == -1 ? 0 : hrmsY.GetMinimum());
    lv2.SetY2(hrmsY.GetMaximum() == 0 ? 1 : hrmsY.GetMaximum());
    lv2.Draw();
    
    if(i == first)
    {
      c.SaveAs("rms.pdf(");
    }
    else
    {
      c.SaveAs("rms.pdf");
    }
    
    fd.Add(&hmeanX);
    fd.Add(&hnX);
    fd.Add(&hrmsX);
    
    fd.Add(&hmeanY);
    fd.Add(&hnY);
    fd.Add(&hrmsY);
    
    fd.Add(&xv);
    fd.Add(&yv);
    fd.Add(&zv);
    
    fout.cd();
    fd.Write();
    
    show(i,true,false,true);
    
    TCanvas* cev = (TCanvas*) gROOT->FindObject("cev");
    
    if(cev)
    {
      if(i == last)
      {
        cev->SaveAs("rms.pdf)");
      }
      else
      {
        cev->SaveAs("rms.pdf");
      }
    }
  }
}