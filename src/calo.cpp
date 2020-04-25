void calo()
{
  gStyle->SetOptFit(1111);
  gROOT->SetBatch();

  TFile f("../files/reco/mu_10GeV_1k.reco.root");
  TTree* t = (TTree*)f.Get("tDigit");
  TCanvas c;
  TH1D* h;

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

  RooFitResult* r1 = marg->fitTo(histo_pmt1,RooFit::Save());
  
  RooPlot* frame1 = npe.frame();

  histo_pmt1.plotOn(frame1);
  marg->plotOn(frame1);
  Double_t vchi2 = frame1->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;

  frame1->Draw();
  c.SaveAs("calo.pdf(");
  
  r1->Print();

  RooFitResult* r2 = marg->fitTo(histo_pmt2,RooFit::Save());
  
  RooPlot* frame2 = npe.frame();

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
  
  RooFitResult* rt0 = gauss.fitTo(histo_tlay0,RooFit::Save());
  
  RooPlot* frame_time0 = time.frame();
  histo_tlay0.plotOn(frame_time0);
  gauss.plotOn(frame_time0);
  vchi2 = frame_time0->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_time0->Draw();
  c.SaveAs("calo.pdf");
  
  rt0->Print();
  
  RooFitResult* rt1 = gauss.fitTo(histo_tlay1,RooFit::Save());
  
  RooPlot* frame_time1 = time.frame();
  histo_tlay1.plotOn(frame_time1);
  gauss.plotOn(frame_time1);
  vchi2 = frame_time1->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_time1->Draw();
  c.SaveAs("calo.pdf");
  
  rt1->Print();
  
  RooFitResult* rt2 = gauss.fitTo(histo_tlay2,RooFit::Save());
  
  RooPlot* frame_time2 = time.frame();
  histo_tlay2.plotOn(frame_time2);
  gauss.plotOn(frame_time2);
  vchi2 = frame_time2->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_time2->Draw();
  c.SaveAs("calo.pdf");
  
  rt2->Print();
  
  RooFitResult* rt3 = gauss.fitTo(histo_tlay3,RooFit::Save());
  
  RooPlot* frame_time3 = time.frame();
  histo_tlay3.plotOn(frame_time3);
  gauss.plotOn(frame_time3);
  vchi2 = frame_time3->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_time3->Draw();
  c.SaveAs("calo.pdf");
  
  rt3->Print();
  
  RooFitResult* rt4 = gauss.fitTo(histo_tlay4,RooFit::Save());
  
  RooPlot* frame_time4 = time.frame();
  histo_tlay4.plotOn(frame_time4);
  gauss.plotOn(frame_time4);
  vchi2 = frame_time4->chiSquare(2);
  std::cout << "chi2: " << vchi2 << std::endl;
  
  frame_time4->Draw();
  c.SaveAs("calo.pdf");
  
  rt4->Print();

  te250->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,1000,5000)", "cell.mod==18",
              "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 250 MeV e-;adc");

  TFitResultPtr r250 = h->Fit("gaus", "QS");
  h->Draw();

  c.SaveAs("calo.pdf");

  te500->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,3000,8000)", "cell.mod==18",
              "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 500 MeV e-;adc");

  TFitResultPtr r500 = h->Fit("gaus", "QS");
  h->Draw();

  c.SaveAs("calo.pdf");

  te750->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,5000,11000)", "cell.mod==18",
              "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 750 MeV e-;adc");

  TFitResultPtr r750 = h->Fit("gaus", "QS");
  h->Draw();

  c.SaveAs("calo.pdf");

  te1000->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,7000,14000)", "cell.mod==18",
               "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 1 GeV e-;adc");

  TFitResultPtr r1000 = h->Fit("gaus", "QS");
  h->Draw();

  c.SaveAs("calo.pdf");

  te2000->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,15000,25000)", "cell.mod==18",
               "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 2 GeV e-;adc");

  TFitResultPtr r2000 = h->Fit("gaus", "QS");
  h->Draw();

  c.SaveAs("calo.pdf");

  te5000->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,20000,60000)", "cell.mod==18",
               "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 5 GeV e-;adc");

  TFitResultPtr r5000 = h->Fit("gaus", "QS");
  h->Draw();

  c.SaveAs("calo.pdf");

  tg250->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,1000,5000)", "cell.mod==18",
              "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 250 MeV gamma;adc");

  TFitResultPtr rg250 = h->Fit("gaus", "QS");
  h->Draw();

  c.SaveAs("calo.pdf");

  tg500->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,3000,8000)", "cell.mod==18",
              "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 500 MeV gamma;adc");

  TFitResultPtr rg500 = h->Fit("gaus", "QS");
  h->Draw();

  c.SaveAs("calo.pdf");

  tg1000->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,7000,14000)", "cell.mod==18",
               "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 1 GeV gamma;adc");

  TFitResultPtr rg1000 = h->Fit("gaus", "QS");
  h->Draw();

  c.SaveAs("calo.pdf");

  TGraphErrors gr(6);
  gr.SetTitle("; E (MeV); adc");
  gr.SetPoint(0, 250, r250.Get()->GetParams()[1]);
  gr.SetPoint(1, 500, r500.Get()->GetParams()[1]);
  gr.SetPoint(2, 750, r750.Get()->GetParams()[1]);
  gr.SetPoint(3, 1000, r1000.Get()->GetParams()[1]);
  gr.SetPoint(4, 2000, r2000.Get()->GetParams()[1]);
  gr.SetPoint(5, 5000, r5000.Get()->GetParams()[1]);

  gr.SetPointError(0, 0., r250.Get()->GetParams()[2]);
  gr.SetPointError(1, 0., r500.Get()->GetParams()[2]);
  gr.SetPointError(2, 0., r750.Get()->GetParams()[2]);
  gr.SetPointError(3, 0., r1000.Get()->GetParams()[2]);
  gr.SetPointError(4, 0., r2000.Get()->GetParams()[2]);
  gr.SetPointError(5, 0., r5000.Get()->GetParams()[2]);

  TGraphErrors grg(3);
  grg.SetMarkerColor(kRed);
  grg.SetLineColor(kRed);
  grg.SetPoint(0, 250, rg250.Get()->GetParams()[1]);
  grg.SetPoint(1, 500, rg500.Get()->GetParams()[1]);
  grg.SetPoint(2, 1000, rg1000.Get()->GetParams()[1]);

  grg.SetPointError(0, 0., rg250.Get()->GetParams()[2]);
  grg.SetPointError(1, 0., rg500.Get()->GetParams()[2]);
  grg.SetPointError(2, 0., rg1000.Get()->GetParams()[2]);

  gr.Fit("pol1", "QS");
  gr.Draw("ap");
  grg.Draw("psame");

  c.SaveAs("calo.pdf");

  TGraph gr1(6);
  gr1.SetTitle("; E (GeV); #sigma_{E}/E");
  gr1.SetPoint(0, 0.25,
               r250.Get()->GetParams()[2] / r250.Get()->GetParams()[1]);
  gr1.SetPoint(1, 0.5, r500.Get()->GetParams()[2] / r500.Get()->GetParams()[1]);
  gr1.SetPoint(2, 0.75,
               r750.Get()->GetParams()[2] / r750.Get()->GetParams()[1]);
  gr1.SetPoint(3, 1.,
               r1000.Get()->GetParams()[2] / r1000.Get()->GetParams()[1]);
  gr1.SetPoint(4, 2.,
               r2000.Get()->GetParams()[2] / r2000.Get()->GetParams()[1]);
  gr1.SetPoint(5, 5.,
               r5000.Get()->GetParams()[2] / r5000.Get()->GetParams()[1]);
  gr1.SetMarkerStyle(34);

  TGraph gr1g(3);
  gr1g.SetPoint(0, 0.25,
                rg250.Get()->GetParams()[2] / rg250.Get()->GetParams()[1]);
  gr1g.SetPoint(1, 0.5,
                rg500.Get()->GetParams()[2] / rg500.Get()->GetParams()[1]);
  gr1g.SetPoint(3, 1.0,
                rg1000.Get()->GetParams()[2] / rg1000.Get()->GetParams()[1]);
  gr1g.SetMarkerStyle(34);
  gr1g.SetMarkerColor(kRed);

  TF1 fres("fres", "[0]+[1]/sqrt(x)");
  fres.SetParameter(0, 0.0);
  fres.SetParameter(1, 0.057);

  gr1.Fit("fres", "QS");
  gr1.Draw("ap");
  gr1g.Draw("psame");
  // fres.Draw("same");

  c.SaveAs("calo.pdf");

  tpi250->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,0,15000)", "cell.mod==18",
               "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 250 MeV pion;adc");

  TFitResultPtr rpi250 = h->Fit("gaus", "QS", "", 3000, 10000);
  h->Draw();

  c.SaveAs("calo.pdf");

  tpi500->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,0,20000)", "cell.mod==18",
               "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 500 MeV pion;adc");

  TFitResultPtr rpi500 = h->Fit("gaus", "QS", "", 4000, 30000);
  h->Draw();

  c.SaveAs("calo.pdf");

  tpi1000->Draw("Sum$(cell.adc1+cell.adc2)>>h(50,0,30000)", "cell.mod==18",
                "E0");
  h = (TH1D*)gROOT->FindObject("h");
  h->SetTitle("sum ADC of 1 GeV pion;adc");

  TFitResultPtr rpi1000 = h->Fit("gaus", "QS", "", 5000, 30000);
  h->Draw();

  c.SaveAs("calo.pdf");

  TGraphErrors grp(3);
  grp.SetTitle("; E (MeV); adc");
  grp.SetPoint(0, 250, rpi250.Get()->GetParams()[1]);
  grp.SetPoint(1, 500, rpi500.Get()->GetParams()[1]);
  grp.SetPoint(2, 1000, rpi1000.Get()->GetParams()[1]);

  grp.SetPointError(0, 0., rpi250.Get()->GetParams()[2]);
  grp.SetPointError(1, 0., rpi500.Get()->GetParams()[2]);
  grp.SetPointError(3, 0., rpi1000.Get()->GetParams()[2]);

  grp.Fit("pol1", "QS");
  grp.Draw("ap");

  c.SaveAs("calo.pdf)");
}
