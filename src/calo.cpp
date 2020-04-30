void calo()
{
  gStyle->SetOptFit(1111);
  gROOT->SetBatch();
  
  RooMsgService::instance().setSilentMode(true);

  TFile f("../files/reco/mu_10GeV_1k.reco.root");
  TTree* t = (TTree*)f.Get("tDigit");
  TTree* tMC = (TTree*)f.Get("EDepSimEvents");
  TCanvas c;
  TH1D* h;
  TG4Event* ev = new TG4Event;
  std::vector<cell>* vc = new std::vector<cell>; 
  
  tMC->SetBranchAddress("Event",&ev);
  t->SetBranchAddress("cell",&vc);

  TFile fe250("../files/reco/e_250MeV_1k.reco.root");
  TFile fe500("../files/reco/e_500MeV_1k.reco.root");
  TFile fe750("../files/reco/e_750MeV_1k.reco.root");
  TFile fe1000("../files/reco/e_1GeV_1k.reco.root");
  TFile fe2000("../files/reco/e_2GeV_1k.reco.root");
  TFile fe5000("../files/reco/e_5GeV_1k.reco.root");

  TFile fg250("../files/reco/gamma_250MeV_1k.reco.root");
  TFile fg500("../files/reco/gamma_500MeV_1k.reco.root");
  TFile fg1000("../files/reco/gamma_1GeV_1k.reco.root");

  TFile fpi250("../files/reco/pi_250MeV_1k.reco.root");
  TFile fpi500("../files/reco/pi_500MeV_1k.reco.root");
  TFile fpi1000("../files/reco/pi_1GeV_1k.reco.root");

  TTree* te250 = (TTree*)fe250.Get("tDigit");
  TTree* te500 = (TTree*)fe500.Get("tDigit");
  TTree* te750 = (TTree*)fe750.Get("tDigit");
  TTree* te1000 = (TTree*)fe1000.Get("tDigit");
  TTree* te2000 = (TTree*)fe2000.Get("tDigit");
  TTree* te5000 = (TTree*)fe5000.Get("tDigit");

  TTree* tg250 = (TTree*)fg250.Get("tDigit");
  TTree* tg500 = (TTree*)fg500.Get("tDigit");
  TTree* tg1000 = (TTree*)fg1000.Get("tDigit");

  TTree* tpi250 = (TTree*)fpi250.Get("tDigit");
  TTree* tpi500 = (TTree*)fpi500.Get("tDigit");
  TTree* tpi1000 = (TTree*)fpi1000.Get("tDigit");

  RooRealVar de("de", "de", 5, 0, 20);
  RooFormulaVar npe_mean("npe_mean", "18.5*0.415*@0", RooArgList(de));
  RooRealVar npe("npe", "npe", 0, 200);

  t->Draw("cell.@pe_time1.size()>>h(200,0,200)", "cell.mod==18&&cell.cel==5",
          "E0");
  h = (TH1D*)gROOT->FindObject("h");

  RooDataHist histo_pmt1("histo_pmt1", "pmt1", npe, h);

  t->Draw("cell.@pe_time2.size()>>h(200,0,200)", "cell.mod==18&&cell.cel==5",
          "E0");
  h = (TH1D*)gROOT->FindObject("h");

  RooDataHist histo_pmt2("histo_pmt2", "pmt2", npe, h);

  RooRealVar de_mean("de_mean", "de mean", 5, 0, 20);
  RooRealVar de_sigma("de_sigma", "de sigma", 0.5, 0, 20);

  RooLandau landau("landau", "landau", de, de_mean, de_sigma);

  RooPoisson poisson("poisson", "poisson", npe, npe_mean);
  
  RooProdPdf prod("prod","landau x poisson",RooArgSet(landau,poisson));
  
  auto marg = prod.createProjection(RooArgSet(de));

  RooFitResult* r1 = marg->fitTo(histo_pmt1,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame1 = npe.frame();
  frame1->SetTitle("#p.e./cell for 10 GeV passing muon; number of p.e.");

  histo_pmt1.plotOn(frame1);
  marg->plotOn(frame1);
  Double_t vchi2 = frame1->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;

  frame1->Draw();
  c.SaveAs("calo.pdf(");
  
  r1->Print();

  RooFitResult* r2 = marg->fitTo(histo_pmt2,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame2 = npe.frame();
  frame2->SetTitle("#p.e./cell for 10 GeV passing muon; number of p.e.");

  histo_pmt2.plotOn(frame2);
  marg->plotOn(frame2);
  vchi2 = frame2->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;

  frame2->Draw();
  c.SaveAs("calo.pdf");
  
  r2->Print();
  
  RooRealVar time("time", "cell time", 9, 5, 15);
  RooRealVar tmean("tmean", "cell time mean", 5, 15);
  RooRealVar tsigma("tsigma", "cell time sigma", 0, 15);

  t->Draw(
      "0.5*(cell.tdc1+cell.tdc2-cell.l*5.85/"
      "1000.)>>h(100,7.5,11)",
      "cell.mod==18&&cell.cel==5&&cell.lay==0", "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("corrected cell time;t (ns)");
  RooDataHist histo_tlay0("histo_tlay0", "tlay0", time, h);

  t->Draw(
      "0.5*(cell.tdc1+cell.tdc2-cell.l*5.85/"
      "1000.)>>h(100,7.5,11)",
      "cell.mod==18&&cell.cel==5&&cell.lay==1", "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("corrected cell time;t (ns)");
  RooDataHist histo_tlay1("histo_tlay1", "tlay1", time, h);

  t->Draw(
      "0.5*(cell.tdc1+cell.tdc2-cell.l*5.85/"
      "1000.)>>h(100,7.5,11)",
      "cell.mod==18&&cell.cel==5&&cell.lay==2", "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("corrected cell time;t (ns)");
  RooDataHist histo_tlay2("histo_tlay2", "tlay2", time, h);

  t->Draw(
      "0.5*(cell.tdc1+cell.tdc2-cell.l*5.85/"
      "1000.)>>h(100,7.5,11)",
      "cell.mod==18&&cell.cel==5&&cell.lay==3", "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("corrected cell time;t (ns)");
  RooDataHist histo_tlay3("histo_tlay3", "tlay3", time, h);

  t->Draw(
      "0.5*(cell.tdc1+cell.tdc2-cell.l*5.85/"
      "1000.)>>h(100,7.5,11)",
      "cell.mod==18&&cell.cel==5&&cell.lay==4", "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("corrected cell time;t (ns)");
  RooDataHist histo_tlay4("histo_tlay4", "tlay4", time, h);
  
  RooGaussian gauss("gauss","gauss", time, tmean, tsigma);
  
  RooFitResult* rt0 = gauss.fitTo(histo_tlay0,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame_time0 = time.frame();
  frame_time0->SetTitle("cell (layer: 0) reco time for 10 GeV passing muon; t_{reco} (ns)");
  histo_tlay0.plotOn(frame_time0);
  gauss.plotOn(frame_time0);
  vchi2 = frame_time0->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_time0->Draw();
  c.SaveAs("calo.pdf");
  
  RooFitResult* rt1 = gauss.fitTo(histo_tlay1,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame_time1 = time.frame();
  frame_time1->SetTitle("cell (layer: 1) reco time for 10 GeV passing muon; t_{reco} (ns)");
  histo_tlay1.plotOn(frame_time1);
  gauss.plotOn(frame_time1);
  vchi2 = frame_time1->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_time1->Draw();
  c.SaveAs("calo.pdf");
  
  RooFitResult* rt2 = gauss.fitTo(histo_tlay2,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame_time2 = time.frame();
  frame_time2->SetTitle("cell (layer: 2) reco time for 10 GeV passing muon; t_{reco} (ns)");
  histo_tlay2.plotOn(frame_time2);
  gauss.plotOn(frame_time2);
  vchi2 = frame_time2->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_time2->Draw();
  c.SaveAs("calo.pdf");
  
  RooFitResult* rt3 = gauss.fitTo(histo_tlay3,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame_time3 = time.frame();
  frame_time3->SetTitle("cell (layer: 3) reco time for 10 GeV passing muon; t_{reco} (ns)");
  histo_tlay3.plotOn(frame_time3);
  gauss.plotOn(frame_time3);
  vchi2 = frame_time3->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_time3->Draw();
  c.SaveAs("calo.pdf");
  
  RooFitResult* rt4 = gauss.fitTo(histo_tlay4,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame_time4 = time.frame();
  frame_time4->SetTitle("cell (layer: 4) reco time for 10 GeV passing muon; t_{reco} (ns)");
  histo_tlay4.plotOn(frame_time4);
  gauss.plotOn(frame_time4);
  vchi2 = frame_time4->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_time4->Draw();
  c.SaveAs("calo.pdf");
  

  t->Draw(
      "0.5*(cell.tdc1+cell.tdc2-cell.l*5.85/"
      "1000.)>>h(100,7.5,11)",
      "cell.mod==18&&cell.cel==5&&cell.lay==4", "E0");

  t->Draw(
      "0.5*(cell.tdc1+cell.tdc2-cell.l*5.85/"
      "1000.)-0.044/0.299792458>>+h",
      "cell.mod==18&&cell.cel==5&&cell.lay==4", "E0");

  t->Draw(
      "0.5*(cell.tdc1+cell.tdc2-cell.l*5.85/"
      "1000.)-(0.044+0.044)/0.299792458>>+h",
      "cell.mod==18&&cell.cel==5&&cell.lay==4", "E0");

  t->Draw(
      "0.5*(cell.tdc1+cell.tdc2-cell.l*5.85/"
      "1000.)-(0.044+0.044+0.044)/0.299792458>>+h",
      "cell.mod==18&&cell.cel==5&&cell.lay==4", "E0");

  t->Draw(
      "0.5*(cell.tdc1+cell.tdc2-cell.l*5.85/"
      "1000.)-(0.044+0.044+0.044+0.049)/0.299792458>>+h",
      "cell.mod==18&&cell.cel==5&&cell.lay==4", "E0");
      
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("corrected cell time;t (ns)");
  
  RooDataHist histo_tall("histo_tall", "tall", time, h);
  
  RooFitResult* rtall = gauss.fitTo(histo_tall,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame_timeall = time.frame();
  histo_tall.plotOn(frame_timeall);
  gauss.plotOn(frame_timeall);
  vchi2 = frame_timeall->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_timeall->Draw();
  c.SaveAs("calo.pdf");
  
  rt0->Print();
  rt1->Print();
  rt2->Print();
  rt3->Print();
  rt4->Print();
  rtall->Print();
  
  /*
  TRandom3 rand(0);
  TH1D hSciPDF("hSciPDF","; p.e. time (ns)",100,10,50);
  double tm;
  
  for(int i = 0; i < 100000; i++)
  {
    tm = 1+2.022/0.299792458+2*5.85+3.08*TMath::Power((1./rand.Uniform()-1.),0.588) + rand.Gaus();
    hSciPDF.Fill(tm);
  }
  */
  
  RooRealVar partAndProgTime("partAndProgTime","particle and light propagation time", 18, 15, 30);
  RooRealVar peTime("peTime","pe time", 22, 10, 200);
  
  //RooRealVar scintTime("scintTime","scintillation time",3,0.,100.);
  RooGenericPdf scintillation("scintillation","scintillation","1./(0.588*3.08)*pow(@0/3.08,1./0.588-1)/pow(1+pow(@0/3.08,1/0.588),2)",RooArgList(peTime));
  
  //RooRealVar spreadTime("spreadTime","spreadTime",0.,-10.,10.);
  RooRealVar spread_m("spread_m","spread_m",0.,-10.,10.);
  RooRealVar spread_s("spread_s","spread_s",1.,0.,10);
  RooGaussian spread("spread","spread",peTime,partAndProgTime,spread_s);
  
  //RooFormulaVar peTime("peTime","@0+@1+@2",RooArgList(partTime,scintTime,spreadTime));
  RooFFTConvPdf pdf_peTime("pdf_peTime","scintillation (X) spread",peTime,scintillation,spread) ;   
  
  t->Draw("cell.pe_time1>>h(100,0,50)","cell.id==18005","E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle(";p.e. time (ns)");
  
  RooDataHist histo_petime("histo_petime", "petime", peTime, h);
  
  RooFitResult* rtpe = pdf_peTime.fitTo(histo_petime,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame_petime = peTime.frame();
  frame_petime->SetTitle("p.e. arrival time for 10 GeV passing muon; t_{arrival} (ns)");
  histo_petime.plotOn(frame_petime);
  pdf_peTime.plotOn(frame_petime);
  //scintillation.plotOn(frame_petime);
  //spread.plotOn(frame_petime);
  
  vchi2 = frame_petime->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_petime->Draw();
  
  c.SaveAs("calo.pdf");
  rtpe->Print();
  /*
  TH2D hCaloHitsL0("hCaloHitsL0","; z (mm); t (ns)",100,25800,26200,100,7,9);
  TH2D hCaloHitsL1("hCaloHitsL1","; z (mm); t (ns)",100,25800,26200,100,7,9);
  TH2D hCaloHitsL2("hCaloHitsL2","; z (mm); t (ns)",100,25800,26200,100,7,9);
  TH2D hCaloHitsL3("hCaloHitsL3","; z (mm); t (ns)",100,25800,26200,100,7,9);
  TH2D hCaloHitsL4("hCaloHitsL4","; z (mm); t (ns)",100,25800,26200,100,7,9);
  
  for(int i = 0; i < tMC->GetEntries(); i++)
  {
    tMC->GetEntry(i);
    t->GetEntry(i);
    
    for(unsigned int j = 0; j < vc->size(); j++)
    {
      if(vc->at(j).id == 18005)
      {
        for(unsigned int k = 0; k < vc->at(j).hindex1.size(); k++)
        {
          hCaloHitsL0.Fill(ev->SegmentDetectors["EMCalSci"][vc->at(j).hindex1.at(k)].GetStart().Z(),ev->SegmentDetectors["EMCalSci"][vc->at(j).hindex1.at(k)].GetStart().T());
        }
      }
      else if (vc->at(j).id == 18105)
      {
        for(unsigned int k = 0; k < vc->at(j).hindex1.size(); k++)
        {
          hCaloHitsL1.Fill(ev->SegmentDetectors["EMCalSci"][vc->at(j).hindex1.at(k)].GetStart().Z(),ev->SegmentDetectors["EMCalSci"][vc->at(j).hindex1.at(k)].GetStart().T());
        }
      }
      else if (vc->at(j).id == 18205)
      {
        for(unsigned int k = 0; k < vc->at(j).hindex1.size(); k++)
        {
          hCaloHitsL2.Fill(ev->SegmentDetectors["EMCalSci"][vc->at(j).hindex1.at(k)].GetStart().Z(),ev->SegmentDetectors["EMCalSci"][vc->at(j).hindex1.at(k)].GetStart().T());
        }
      }
      else if (vc->at(j).id == 18305)
      {
        for(unsigned int k = 0; k < vc->at(j).hindex1.size(); k++)
        {
          hCaloHitsL3.Fill(ev->SegmentDetectors["EMCalSci"][vc->at(j).hindex1.at(k)].GetStart().Z(),ev->SegmentDetectors["EMCalSci"][vc->at(j).hindex1.at(k)].GetStart().T());
        }
      }
      else if (vc->at(j).id == 18405)
      {
        for(unsigned int k = 0; k < vc->at(j).hindex1.size(); k++)
        {
          hCaloHitsL4.Fill(ev->SegmentDetectors["EMCalSci"][vc->at(j).hindex1.at(k)].GetStart().Z(),ev->SegmentDetectors["EMCalSci"][vc->at(j).hindex1.at(k)].GetStart().T());
        }
      }
    }
  }
  hCaloHitsL0.Draw("colz");
  c.SaveAs("calo.pdf");
  hCaloHitsL1.Draw("colz");
  c.SaveAs("calo.pdf");
  hCaloHitsL2.Draw("colz");
  c.SaveAs("calo.pdf");
  hCaloHitsL3.Draw("colz");
  c.SaveAs("calo.pdf");
  hCaloHitsL4.Draw("colz");
  c.SaveAs("calo.pdf");
  
  hCaloHitsL0.ProjectionY()->Draw();
  c.SaveAs("calo.pdf");
  hCaloHitsL1.ProjectionY()->Draw();
  c.SaveAs("calo.pdf");
  hCaloHitsL2.ProjectionY()->Draw();
  c.SaveAs("calo.pdf");
  hCaloHitsL3.ProjectionY()->Draw();
  c.SaveAs("calo.pdf");
  hCaloHitsL4.ProjectionY()->Draw();
  c.SaveAs("calo.pdf");*/
  
  TGraphErrors gr(5);
  gr.SetTitle("; electron energy (GeV); #Sigma(adc)");
  
  TGraphErrors grres(5);
  grres.SetTitle("; electron energy (GeV); #sigma_{#Sigma(adc)}/#Sigma(adc)");
  grres.SetMarkerStyle(34);
  
  RooRealVar adc_m("adc_m", "adc_m", 3000, 0, 60000);
  RooRealVar adc_s("adc_s", "adc_s", 3000, 0, 60000);
  
  te250->Draw("Sum$(cell.adc1+cell.adc2)>>h250(50,1000,5000)", "Sum$((cell.lay==4))==0&&cell.mod==18",
              "E0");
  h = (TH1D*)gROOT->FindObject("h250");
  h->SetTitle("sum ADC of 250 MeV e-;adc");
  
  RooRealVar adc250("adc250", "adc250", 3000, 0, 60000);
  RooGaussian gaussE250("gaussE250","gaussE250", adc250, adc_m, adc_s);
  RooDataHist histo250_dE("histo250_dE", "dEall", adc250, h);
  
  RooFitResult* r250dEall = gaussE250.fitTo(histo250_dE,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame250_dE = adc250.frame();
  frame250_dE->SetTitle("sum of adc for 250 MeV electron; #Sigma(adc)");
  histo250_dE.plotOn(frame250_dE);
  gaussE250.plotOn(frame250_dE);
  vchi2 = frame250_dE->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  gr.SetPoint(0, 0.25, adc_m.getVal());
  gr.SetPointError(0, 0., adc_s.getVal());
  grres.SetPoint(0, 0.25, adc_s.getVal()/adc_m.getVal());
  grres.SetPointError(0, 0., adc_s.getVal()/adc_m.getVal()*sqrt(pow(adc_s.getError()/adc_s.getVal(),2)+pow(adc_m.getError()/adc_m.getVal(),2)));
  
  frame250_dE->Draw();
  c.SaveAs("calo.pdf");
  
  te500->Draw("Sum$(cell.adc1+cell.adc2)>>h500(50,3000,8000)", "Sum$((cell.lay==4)*(cell.mod==18)*(cell.adc1+cell.adc2))<500&&cell.mod==18",
              "E0");
  h = (TH1D*)gROOT->FindObject("h500");
  h->SetTitle("sum ADC of 500 MeV e-;adc");
  
  RooRealVar adc500("adc500", "adc500", 3000, 0, 60000);
  RooGaussian gaussE500("gaussE500","gaussE500", adc500, adc_m, adc_s);
  RooDataHist histo500_dE("histo500_dE", "dEall", adc500, h);
  
  RooFitResult* r500dEall = gaussE500.fitTo(histo500_dE,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame500_dE = adc500.frame();
  histo500_dE.plotOn(frame500_dE);
  gaussE500.plotOn(frame500_dE);
  vchi2 = frame500_dE->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  gr.SetPoint(1, 0.5, adc_m.getVal());
  gr.SetPointError(1, 0., adc_s.getVal());
  grres.SetPoint(1, 0.5, adc_s.getVal()/adc_m.getVal());
  grres.SetPointError(1, 0., adc_s.getVal()/adc_m.getVal()*sqrt(pow(adc_s.getError()/adc_s.getVal(),2)+pow(adc_m.getError()/adc_m.getVal(),2)));
  
  frame500_dE->Draw();
  c.SaveAs("calo.pdf");
  
  te750->Draw("Sum$(cell.adc1+cell.adc2)>>h750(50,5000,11000)", "Sum$((cell.lay==4)*(cell.mod==18)*(cell.adc1+cell.adc2))<500&&cell.mod==18",
              "E0");
  h = (TH1D*)gROOT->FindObject("h750");
  h->SetTitle("sum ADC of 750 MeV e-;adc");
  
  RooRealVar adc750("adc750", "adc750", 3000, 0, 60000);
  RooGaussian gaussE750("gaussE750","gaussE750", adc750, adc_m, adc_s);
  RooDataHist histo750_dE("histo750_dE", "dEall", adc750, h);
  
  RooFitResult* r750dEall = gaussE750.fitTo(histo750_dE,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame750_dE = adc750.frame();
  histo750_dE.plotOn(frame750_dE);
  gaussE750.plotOn(frame750_dE);
  vchi2 = frame750_dE->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  gr.SetPoint(2, 0.75, adc_m.getVal());
  gr.SetPointError(2, 0., adc_s.getVal());
  grres.SetPoint(2, 0.75, adc_s.getVal()/adc_m.getVal());
  grres.SetPointError(2, 0., adc_s.getVal()/adc_m.getVal()*sqrt(pow(adc_s.getError()/adc_s.getVal(),2)+pow(adc_m.getError()/adc_m.getVal(),2)));
  
  frame750_dE->Draw();
  c.SaveAs("calo.pdf");
  
  te1000->Draw("Sum$(cell.adc1+cell.adc2)>>h1000(50,7000,14000)", "Sum$((cell.lay==4)*(cell.mod==18)*(cell.adc1+cell.adc2))<500&&cell.mod==18",
              "E0");
  h = (TH1D*)gROOT->FindObject("h1000");
  h->SetTitle("sum ADC of 1000 MeV e-;adc");
  
  RooRealVar adc1000("adc1000", "adc1000", 3000, 0, 60000);
  RooGaussian gaussE1000("gaussE1000","gaussE1000", adc1000, adc_m, adc_s);
  RooDataHist histo1000_dE("histo1000_dE", "dEall", adc1000, h);
  
  RooFitResult* r1000dEall = gaussE1000.fitTo(histo1000_dE,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame1000_dE = adc1000.frame();
  frame1000_dE->SetTitle("sum of adc for 1 GeV electron; #Sigma(adc)");
  histo1000_dE.plotOn(frame1000_dE);
  gaussE1000.plotOn(frame1000_dE);
  vchi2 = frame1000_dE->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  gr.SetPoint(3, 1., adc_m.getVal());
  gr.SetPointError(3, 0., adc_s.getVal());
  grres.SetPoint(3, 1., adc_s.getVal()/adc_m.getVal());
  grres.SetPointError(3, 0., adc_s.getVal()/adc_m.getVal()*sqrt(pow(adc_s.getError()/adc_s.getVal(),2)+pow(adc_m.getError()/adc_m.getVal(),2)));
  
  frame1000_dE->Draw();
  c.SaveAs("calo.pdf");
  
  te2000->Draw("Sum$(cell.adc1+cell.adc2)>>h2000(20,15000,25000)", "Sum$((cell.lay==4)*(cell.mod==18)*(cell.adc1+cell.adc2))<500&&cell.mod==18",
              "E0");
  h = (TH1D*)gROOT->FindObject("h2000");
  h->SetTitle("sum ADC of 2000 MeV e-;adc");
  
  RooRealVar adc2000("adc2000", "adc2000", 3000, 0, 60000);
  RooGaussian gaussE2000("gaussE2000","gaussE2000", adc2000, adc_m, adc_s);
  RooDataHist histo2000_dE("histo2000_dE", "dEall", adc2000, h);
  
  RooFitResult* r2000dEall = gaussE2000.fitTo(histo2000_dE,RooFit::Save(),RooFit::PrintEvalErrors(-1));
  
  RooPlot* frame2000_dE = adc2000.frame();
  histo2000_dE.plotOn(frame2000_dE);
  gaussE2000.plotOn(frame2000_dE);
  vchi2 = frame2000_dE->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  gr.SetPoint(4, 2., adc_m.getVal());
  gr.SetPointError(4, 0., adc_s.getVal());
  grres.SetPoint(4, 2., adc_s.getVal()/adc_m.getVal());
  grres.SetPointError(4, 0., adc_s.getVal()/adc_m.getVal()*sqrt(pow(adc_s.getError()/adc_s.getVal(),2)+pow(adc_m.getError()/adc_m.getVal(),2)));
  
  frame2000_dE->Draw();
  c.SaveAs("calo.pdf");
  
  r250dEall->Print();
  r500dEall->Print();
  r750dEall->Print();
  r1000dEall->Print();
  r2000dEall->Print();

  gr.Fit("pol1", "QS");
  gr.Draw("ap");
  
  c.Update();
  TPaveStats *st = (TPaveStats*)gr.FindObject("stats");
  st->SetX1NDC(0.2); //new x start position
  st->SetY1NDC(0.85); //new x end position
  st->SetX2NDC(0.4); //new x start position
  st->SetY2NDC(0.7); //new x end position
  gr.Draw("ap");

  c.SaveAs("calo.pdf");

  //TF1 fres("fres", "[0]+[1]/sqrt(x)+[2]/x",0,6);
  TF1 fres("fres", "[0]+[1]/sqrt(x)",0,6);
  fres.SetParameter(0, 0.0);
  fres.SetParameter(1, 0.057);
  //fres.SetParameter(2, 0.0);

  TF1 frestrue("frestrue", "[0]+[1]/sqrt(x)",0,6);
  frestrue.SetParameter(0, 0.0);
  frestrue.SetParameter(1, 0.057);
  frestrue.SetLineColor(kBlue);
  frestrue.SetLineStyle(2);

  grres.Fit("fres", "QS");
  grres.Draw("ap");
  
  c.Update();
  st = (TPaveStats*)grres.FindObject("stats");
  st->SetX1NDC(0.85); //new x start position
  st->SetY1NDC(0.85); //new x end position
  st->SetX2NDC(0.7); //new x start position
  st->SetY2NDC(0.7); //new x end position
  grres.Draw("ap");
  frestrue.Draw("same");

  c.SaveAs("calo.pdf)");
}
