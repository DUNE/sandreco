void dispTrkFnd()
{
  gROOT->SetBatch(true);
  gSystem->Load("lib/libTrack3D.so");

  TFile f("vtx.root");
  TTree* tclX = (TTree*) f.Get("tclX");
  TTree* tclY = (TTree*) f.Get("tclY");
  
  std::vector<digit>* cluX = new std::vector<digit>;
  std::vector<digit>* cluY = new std::vector<digit>;
  TLorentzVector* PartPX = new TLorentzVector;
  TLorentzVector* PartPY = new TLorentzVector;
  std::vector<digit>* PartDX = new std::vector<digit>;
  std::vector<digit>* PartDY = new std::vector<digit>;
  int PartPDGX;
  int PartPDGY;
  double pHX;
  double pHY;
  int PartIDX;
  int PartIDY;
  int ix, iy;
  
  tclX->SetBranchAddress("evID", &ix);
  tclX->SetBranchAddress("cluX", &cluX);
  tclX->SetBranchAddress("PartIDX", &PartIDX);
  tclX->SetBranchAddress("PartPX", &PartPX);
  tclX->SetBranchAddress("PartDX", &PartDX);
  tclX->SetBranchAddress("PartPDGX", &PartPDGX);
  tclX->SetBranchAddress("pHX", &pHX);
  
  tclY->SetBranchAddress("evID", &iy);
  tclY->SetBranchAddress("cluY", &cluY);
  tclX->SetBranchAddress("PartIDY", &PartIDY);
  tclY->SetBranchAddress("PartPY", &PartPY);
  tclY->SetBranchAddress("PartDY", &PartDY);
  tclY->SetBranchAddress("PartPDGY", &PartPDGY);
  tclY->SetBranchAddress("pHY", &pHY);
  
  TCanvas c("c", "pattern recognition", 1200, 600);
  c.SaveAs("patternY.pdf(");
  
  std::vector<double> x;
  std::vector<double> y;
  TGraph* grMC;
  TGraph* grReco;
  
  tclY->GetEntry(0);
  
  int color = 0;
  int ievent = iy; 
  
  int nev = tclY->GetEntries();
  
  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
  
  c.Divide(2,1);
  c.cd(1)->DrawFrame(21500,-5000,26500,0,TString::Format("Event: %d [MC]",ievent).Data());
  c.cd(2)->DrawFrame(21500,-5000,26500,0,TString::Format("Event: %d [Reco]",ievent).Data());
  
  for(unsigned int i = 0; i < nev; i++)
  {
    tclY->GetEntry(i);
    
    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i) / nev * 100)
              << "%]" << std::flush;
    
    x.clear();
    y.clear();
    
    if(cluY->size() > 3)
    {
    
      for(unsigned int k = 0; k < cluY->size(); k++)
      {
        x.push_back(cluY->at(k).z);
        y.push_back(cluY->at(k).y);
      }
      grReco = new TGraph(x.size(), x.data(), y.data());
      
      for(unsigned int k = 0; k < PartDY->size(); k++)
      {
        x.push_back(PartDY->at(k).z);
        y.push_back(PartDY->at(k).y);
      }
      grMC = new TGraph(x.size(), x.data(), y.data());
      
      if(ievent == iy)
      {
        color++;
      }
      else
      {
        c.SaveAs("patternY.pdf");
        ievent = iy;
        color = 1;
        c.cd(1)->DrawFrame(21500,-5000,26500,0, TString::Format("Event: %d [MC]",ievent).Data());
        c.cd(2)->DrawFrame(21500,-5000,26500,0, TString::Format("Event: %d [Reco]",ievent).Data());
      }
      
    }
    c.cd(1);
    grMC->SetMarkerColor(color);
    grMC->Draw("*same");
    c.cd(2);
    grReco->SetMarkerColor(color);
    grReco->Draw("*same");
  }

  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  c.Clear();
  c.SaveAs("patternY.pdf)");
  
  tclX->GetEntry(0);
  
  color = 0;
  ievent = ix; 
  
  nev = tclX->GetEntries();
  
  std::cout << "Events: " << nev << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;
  
  c.SaveAs("patternX.pdf(");
  
  c.Divide(2,1);
  c.cd(1)->DrawFrame(21500,-2500,26500,2500,TString::Format("Event: %d [MC]",ievent).Data());
  c.cd(2)->DrawFrame(21500,-2500,26500,2500,TString::Format("Event: %d [Reco]",ievent).Data());
  
  for(unsigned int i = 0; i < nev; i++)
  {
    tclX->GetEntry(i);
    
    x.clear();
    y.clear();
    
    if(cluX->size() > 3)
    {
    
      for(unsigned int k = 0; k < cluX->size(); k++)
      {
        x.push_back(cluX->at(k).z);
        y.push_back(cluX->at(k).x);
      }
      grReco = new TGraph(x.size(), x.data(), y.data());
      
      for(unsigned int k = 0; k < PartDX->size(); k++)
      {
        x.push_back(PartDX->at(k).z);
        y.push_back(PartDX->at(k).x);
      }
      grMC = new TGraph(x.size(), x.data(), y.data());
      
      if(ievent == ix)
      {
        color++;
      }
      else
      {
        c.SaveAs("patternX.pdf");
        ievent = ix;
        color = 1;
        c.cd(1)->DrawFrame(21500,-2500,26500,2500,TString::Format("Event: %d [MC]",ievent).Data());
        c.cd(2)->DrawFrame(21500,-2500,26500,2500,TString::Format("Event: %d [Reco]",ievent).Data());
      }
    
    }
    c.cd(1);
    grMC->SetMarkerColor(color);
    grMC->Draw("*same");
    c.cd(2);
    grReco->SetMarkerColor(color);
    grReco->Draw("*same");
  }

  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;
  
  c.Clear();
  c.SaveAs("patternX.pdf)");
}