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
              erase_element(prev2_el->first, regMtrX);
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
              erase_element(prev2_el->first, reg0trX);
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
              prev_el->second += this_el->second;
              regionsX.erase(this_el);
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

  bool display = true;
  bool save2root = true && display;
  bool save2pdf = true && display;

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
  TTree tv("tv","tv");
  
  double xvtx_true, yvtx_true, zvtx_true;
  int nint;
  double x_int[256];
  double y_int[256];
  double z_int[256];
    
  tv.Branch("xvtx",&xvtx_true,"xvtx/D");
  tv.Branch("yvtx",&yvtx_true,"yvtx/D");
  tv.Branch("zvtx",&zvtx_true,"zvtx/D");
  tv.Branch("nint",&nint,"nint/I"); 
  tv.Branch("x_int",x_int,"x_int[nint]/D"); 
  tv.Branch("y_int",y_int,"y_int[nint]/D"); 
  tv.Branch("z_int",z_int,"z_int[nint]/D"); 
  
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
  
  int first = 0;
  int last = 10;/*tDigit->GetEntries();*/
    
  vector<double> zvX;
  vector<double> zvY; 
  vector<double> zvtx;
  
  vector<TH1D> vharctg; 
  vector<TGraph>  vguv;
  vector<TH2D> vhHT; 
  
  TCanvas c;
  c.SaveAs("rms.pdf(");
  
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
    
    xvtx_true = xv.GetVal();
    yvtx_true = yv.GetVal();
    zvtx_true = zv.GetVal();
    
    TDirectoryFile* fd;
    
    if(save2root)
      fd = new TDirectoryFile(TString::Format("ev_%d",i).Data(),TString::Format("ev_%d",i).Data(),"",&fout);
    
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
    
    // eveluate rms, mean and multiplicity in Y and X as a function of the STT module
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
    
    
    
    // find regions of: 0 tracks, 1 track, multi tracks
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
    
    // order regions by width
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
    
    // merge regions less than threshold [in number of modules] (threshold to be optimized) 
    // to the adiacent ones starting from wider multi track regions
    // after multi tracks regions => wider 0 tracks regions
    // if regions of multi tracks and 0 tracks less than 
    // threshold are in the middle of 1 track region, the former
    // is merged to the latter
    const int minreg = 1;
    
    processMtrVct(regMtrX, reg0trX, regionsX, hmultX, minreg);
    process0trVct(reg0trX, regionsX, hmultX, minreg);
    
    processMtrVct(regMtrY, reg0trY, regionsY, hmultY, minreg);
    process0trVct(reg0trY, regionsY, hmultY, minreg);
    
    // find surface between:
    // - 0 tracks => 1 tracks : mono prong vertex
    // - 0 tracks => multi tracks: multi prong vertex
    // - 1 track => multi tracks: re-interaction or vertex with back-scattering
    zvX.clear();
    zvY.clear();
    
    for(std::map<int,int>::iterator it = regionsX.begin(); std::next(it) != regionsX.end(); ++it)
    {
      if(hmultX.GetBinContent(it->first) == 0 && hmultX.GetBinContent(std::next(it)->first) != 0)
      {
        zvX.push_back(hmultX.GetXaxis()->GetBinLowEdge(std::next(it)->first));
      }
      else if(hmultX.GetBinContent(it->first) == 1 && hmultX.GetBinContent(std::next(it)->first) == 2)
      {
        zvX.push_back(hmultX.GetXaxis()->GetBinLowEdge(std::next(it)->first));
      }
    }
    
    for(std::map<int,int>::iterator it = regionsY.begin(); std::next(it) != regionsY.end(); ++it)
    {
      if(hmultY.GetBinContent(it->first) == 0 && hmultY.GetBinContent(std::next(it)->first) != 0)
      {
        zvY.push_back(hmultY.GetXaxis()->GetBinLowEdge(std::next(it)->first));
      }
      else if(hmultY.GetBinContent(it->first) == 1 && hmultY.GetBinContent(std::next(it)->first) == 2)
      {
        zvY.push_back(hmultY.GetXaxis()->GetBinLowEdge(std::next(it)->first));
      }
    }
    
    // merge vertex between X and Y plan
    // if they are in between 200 mm => to be optimized
    std::map<double,int> vv;
    for(unsigned int j = 0; j < zvX.size(); j++)
    {
      vv[zvX[j]] = 1;
    }
    for(unsigned int j = 0; j < zvY.size(); j++)
    {
      vv[zvY[j]] = -1;
    }
    
    zvtx.clear();
    const double max_dist = 0.; //mm
    while(vv.size() != 0)
    {
      std::map<double,int>::iterator it = vv.begin();
      std::map<double,int>::iterator itt = it;
      bool found = false;
      while(++itt != vv.end())
      {
        if(abs(itt->first-it->first) < max_dist)
        {
          if(it->second != itt->second)
          {
            if(std::next(itt) != vv.end())
            {
              if(abs(itt->first-it->first) < abs(itt->first-std::next(it)->first) || itt->second == std::next(itt)->second)
              {
                found = true;
              }
            }
            else
            {
              found = true;
            }
            break;
          }
        }
      }
      if(found)
      {
        zvtx.push_back(0.5*(itt->first+it->first));
        vv.erase(itt);
        vv.erase(it);
      }
      else
      {
        zvtx.push_back(it->first);
        vv.erase(it);
      }
    }
    
    double u, v, d;
    double zc, yc;
    double zz, yy;
    
    vguv.clear();
    vharctg.clear();
    vhHT.clear();
    
    nint = zvtx.size();
    for(unsigned int k = 0; k < zvtx.size(); k++)
    {
      x_int[k] = hmeanX.GetBinContent(hmeanX.FindBin(zvtx[k]));
      y_int[k] = hmeanY.GetBinContent(hmeanY.FindBin(zvtx[k]));
      z_int[k] = zvtx[k];
      
      if(x_int[k] == -1. || y_int[k] == -1.)
      {
        /*
        std::cout << "event: " << i << " ; z_int[k]   = " << z_int[k] << std::endl;
        std::cout << "warning....: x_int[k]   = " << x_int[k] << " ; y_int[k]   = " << y_int[k] << std::endl;
        std::cout << "next   ....: x_int[k+1] = " << hmeanX.GetBinContent(hmeanX.FindBin(zvtx[k])+1) << " ; y_int[k+1] = " << hmeanY.GetBinContent(hmeanY.FindBin(zvtx[k])+1) << std::endl;
        */
        continue;
      }
      
      TH2D hHT(TString::Format("hHT_%d",k).Data(),TString::Format("zv: %.1f - reco: %.1f; zc (mm); yc (mm)",zv.GetVal(),z_int[k]).Data(),100,-2000,2000,100,-2000,2000);
      
      TH1D harctg(TString::Format("harctg_%d",k).Data(),TString::Format("zv: %.1f - reco: %.1f; #phi (rad)",zv.GetVal(),z_int[k]).Data(), 200, -TMath::Pi(), TMath::Pi());
      
      std::vector<double> vu;
      std::vector<double> vv;
      
      for(unsigned int j = 0; j < digits->size(); j++)
      {
        if(digits->at(j).hor == 1)
        {
          
          u = (digits->at(j).z - z_int[k]);
          v = (y_int[k] - digits->at(j).y);
          d = (u*u + v*v);
          
          u /=d;
          v /=d;
          
          vu.push_back(u);
          vv.push_back(v);
          
          double y = digits->at(j).y - y_int[k];
          double z = digits->at(j).z - z_int[k];
          double R = z*z + y*y; 
          
          for(int kk = 0; kk < hHT.GetNbinsX(); kk++)
          {
            zz = hHT.GetXaxis()->GetBinCenter(kk+1);
            yy = (-z*zz + 0.5*R)/y;
            hHT.Fill(zz, yy);
          }
          
          harctg.Fill(TMath::ATan2(v,u));
        }
      } 
      TGraph guv(vu.size(), vu.data(),vv.data());
      guv.SetName(TString::Format("huv_%d",k).Data());
      guv.SetTitle(TString::Format("zv: %.1f - reco: %.1f; u; v",zv.GetVal(),z_int[k]).Data()); 
      
      vguv.push_back(std::move(guv));
      vharctg.push_back(std::move(harctg));
      vhHT.push_back(std::move(hHT));
    }
    
    if(display)
    {
      
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
      
      // display of the event
      show(i,true,false,true);
      
      TCanvas* cev = (TCanvas*) gROOT->FindObject("cev");
      
      // histograms with regions
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
      
      TLine* lvx;
      TLine* lvy;
      
      c.Clear();
      c.Divide(1,2);
      
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
      
      c.cd(1);
      for(unsigned int j = 0; j < zvtx.size(); j++)
      {
        lvx = new TLine;
        lvx->SetLineColor(kViolet);
        lvx->SetLineWidth(2);
        lvx->SetLineStyle(8);
        lvx->SetX1(zvtx[j]);
        lvx->SetX2(zvtx[j]);
        lvx->SetY1(hrmsX.GetMinimum() == -1 ? 0 : hrmsX.GetMinimum());
        lvx->SetY2(hrmsX.GetMaximum() == 0 ? 1 : hrmsX.GetMaximum());
        lvx->Draw();
      }
      
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
      
      c.cd(2);
      for(unsigned int j = 0; j < zvtx.size(); j++)
      {
        lvy = new TLine;
        lvy->SetLineColor(kViolet);
        lvy->SetLineWidth(2);
        lvy->SetLineStyle(8);
        lvy->SetX1(zvtx[j]);
        lvy->SetX2(zvtx[j]);
        lvy->SetY1(hrmsY.GetMinimum() == -1 ? 0 : hrmsY.GetMinimum());
        lvy->SetY2(hrmsY.GetMaximum() == 0 ? 1 : hrmsY.GetMaximum());
        lvy->Draw();
      }
      
      // save to pdf
      if(save2pdf)
      {
        if(cev)
        {
          cev->SaveAs("rms.pdf");
        }
        
        c.SaveAs("rms.pdf");
      
        c.Clear();
        c.cd();
        
        for(unsigned int k = 0; k < vharctg.size(); k++)
        {
          vharctg[k].Draw();
          c.SaveAs("rms.pdf");
          vguv[k].Draw("ap*");
          c.SaveAs("rms.pdf");
          vhHT[k].Draw("colz");
          c.SaveAs("rms.pdf");
        }
      }
    }
    
    if(save2root)
    {
      fd->Add(&hmeanX);
      fd->Add(&hnX);
      fd->Add(&hrmsX);
      
      fd->Add(&hmeanY);
      fd->Add(&hnY);
      fd->Add(&hrmsY);
      
      fd->Add(&xv);
      fd->Add(&yv);
      fd->Add(&zv);
      
      for(unsigned int k = 0; k < vharctg.size(); k++)
      {
        fd->Add(&vharctg[k]);
        fd->Add(&vguv[k]);
        fd->Add(&vhHT[k]);
      }
      
      fout.cd();
      fd->Write();
    }
    
    tv.Fill();
  }
  c.SaveAs("rms.pdf)");
  fout.cd();
  tv.Write();
  fout.Close();
}