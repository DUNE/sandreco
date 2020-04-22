
void ana()
{
  gROOT->SetBatch();
  
  gStyle->SetOptFit(1111);

  TFile f("numu_geoV12_100000.0.check.root");
  auto t = static_cast<TTree*>(f.Get("tCheck"));
  t->SetAlias("pt_true","sqrt(py_true*py_true+pz_true*pz_true)");
  t->SetAlias("pt_old","sqrt(py_old*py_old+pz_old*pz_old)");
  t->SetAlias("pt_new","sqrt(py_new*py_new+pz_new*pz_new)");
  t->SetAlias("dip_true","TMath::ATan2(px_true,sqrt(py_true*py_true+pz_true*pz_true))");
  t->SetAlias("dip_old","TMath::ATan2(px_old,sqrt(py_old*py_old+pz_old*pz_old))");
  t->SetAlias("dip_new","TMath::ATan2(px_new,sqrt(py_new*py_new+pz_new*pz_new))");
  
  TCanvas c;
  
  t->Draw("1-pt_true/pt_old>>h1_old(100,-1,1)","pt_true<500");
  auto h1_old = static_cast<TH1F*>(gDirectory->Get("h1_old"));
  auto r1_old = h1_old->Fit("gaus","QS","",-0.5,0.5);
  h1_old->Draw("E0");
  c.Print("ana.pdf(");
  
  t->Draw("1-pt_true/pt_old>>h2_old(100,-1,1)","pt_true>=500&&pt_true<1000");
  auto h2_old = static_cast<TH1F*>(gDirectory->Get("h2_old"));
  auto r2_old = h2_old->Fit("gaus","QS","",-0.5,0.5);
  h2_old->Draw("E0");
  
  t->Draw("1-pt_true/pt_old>>h3_old(100,-1,1)","pt_true>=1000&&pt_true<2000");
  auto h3_old = static_cast<TH1F*>(gDirectory->Get("h3_old"));
  auto r3_old = h3_old->Fit("gaus","QS","",-0.5,0.5);
  h3_old->Draw("E0");
  c.Print("ana.pdf");
  
  t->Draw("1-pt_true/pt_old>>h4_old(100,-1,1)","pt_true>=2000&&pt_true<3000");
  auto h4_old = static_cast<TH1F*>(gDirectory->Get("h4_old"));
  auto r4_old = h4_old->Fit("gaus","QS","",-0.5,0.5);
  h4_old->Draw("E0");
  c.Print("ana.pdf");
  
  t->Draw("1-pt_true/pt_old>>h5_old(100,-1,1)","pt_true>=3000&&pt_true<4000");
  auto h5_old = static_cast<TH1F*>(gDirectory->Get("h5_old"));
  auto r5_old = h5_old->Fit("gaus","QS","",-0.5,0.5);
  h5_old->Draw("E0");
  c.Print("ana.pdf");
  
  t->Draw("1-pt_true/pt_old>>h6_old(100,-1,1)","pt_true>=4000&&pt_true<5000");
  auto h6_old = static_cast<TH1F*>(gDirectory->Get("h6_old"));
  auto r6_old = h6_old->Fit("gaus","QS","",-0.5,0.5);
  h6_old->Draw("E0");
  c.Print("ana.pdf");
  
  t->Draw("1-pt_true/pt_old>>h7_old(100,-1,1)","pt_true>=5000&&pt_true<10000");
  auto h7_old = static_cast<TH1F*>(gDirectory->Get("h7_old"));
  auto r7_old = h7_old->Fit("gaus","QS","",-0.5,0.5);
  h7_old->Draw("E0");
  c.Print("ana.pdf");
  
  t->Draw("1-pt_true/pt_new>>h1_new(100,-1,1)","pt_true<500");
  auto h1_new = static_cast<TH1F*>(gDirectory->Get("h1_new"));
  auto r1_new = h1_new->Fit("gaus","QS","",-0.5,0.5);
  h1_new->Draw("E0");
  c.Print("ana.pdf(");
  
  t->Draw("1-pt_true/pt_new>>h2_new(100,-1,1)","pt_true>=500&&pt_true<1000");
  auto h2_new = static_cast<TH1F*>(gDirectory->Get("h2_new"));
  auto r2_new = h2_new->Fit("gaus","QS","",-0.5,0.5);
  h2_new->Draw("E0");
  
  t->Draw("1-pt_true/pt_new>>h3_new(100,-1,1)","pt_true>=1000&&pt_true<2000");
  auto h3_new = static_cast<TH1F*>(gDirectory->Get("h3_new"));
  auto r3_new = h3_new->Fit("gaus","QS","",-0.5,0.5);
  h3_new->Draw("E0");
  c.Print("ana.pdf");
  
  t->Draw("1-pt_true/pt_new>>h4_new(100,-1,1)","pt_true>=2000&&pt_true<3000");
  auto h4_new = static_cast<TH1F*>(gDirectory->Get("h4_new"));
  auto r4_new = h4_new->Fit("gaus","QS","",-0.5,0.5);
  h4_new->Draw("E0");
  c.Print("ana.pdf");
  
  t->Draw("1-pt_true/pt_new>>h5_new(100,-1,1)","pt_true>=3000&&pt_true<4000");
  auto h5_new = static_cast<TH1F*>(gDirectory->Get("h5_new"));
  auto r5_new = h5_new->Fit("gaus","QS","",-0.5,0.5);
  h5_new->Draw("E0");
  c.Print("ana.pdf");
  
  t->Draw("1-pt_true/pt_new>>h6_new(100,-1,1)","pt_true>=4000&&pt_true<5000");
  auto h6_new = static_cast<TH1F*>(gDirectory->Get("h6_new"));
  auto r6_new = h6_new->Fit("gaus","QS","",-0.5,0.5);
  h6_new->Draw("E0");
  c.Print("ana.pdf");
  
  t->Draw("1-pt_true/pt_new>>h7_new(100,-1,1)","pt_true>=5000&&pt_true<10000");
  auto h7_new = static_cast<TH1F*>(gDirectory->Get("h7_new"));
  auto r7_new = h7_new->Fit("gaus","QS","",-0.5,0.5);
  h7_new->Draw("E0");
  c.Print("ana.pdf");
  
  TGraphErrors gr_old(7);
  gr_old.SetPoint(0,0.25,r1_old->GetParams()[2]);
  gr_old.SetPoint(1,0.75,r2_old->GetParams()[2]);
  gr_old.SetPoint(2,1.5,r3_old->GetParams()[2]);
  gr_old.SetPoint(3,2.5,r4_old->GetParams()[2]);
  gr_old.SetPoint(4,3.5,r5_old->GetParams()[2]);
  gr_old.SetPoint(5,4.5,r6_old->GetParams()[2]);
  gr_old.SetPoint(6,7.5,r7_old->GetParams()[2]);
  gr_old.SetPointError(0,0.25,r1_old->GetErrors()[2]);
  gr_old.SetPointError(1,0.25,r2_old->GetErrors()[2]);
  gr_old.SetPointError(2,0.5,r3_old->GetErrors()[2]);
  gr_old.SetPointError(3,0.5,r4_old->GetErrors()[2]);
  gr_old.SetPointError(4,0.5,r5_old->GetErrors()[2]);
  gr_old.SetPointError(5,0.5,r6_old->GetErrors()[2]);
  gr_old.SetPointError(6,2.5,r7_old->GetErrors()[2]);
  gr_old.SetMarkerColor(kRed);
  gr_old.SetLineColor(kRed);
  
  //TF1 f_old("f_old","[0]+[1]/sqrt(x)",0.,10.);
  //auto r_old = gr_old.Fit("f_old","SQ","",0.,10.);
  
  TGraphErrors gr_new(7);
  gr_new.SetPoint(0,0.25,r1_new->GetParams()[2]);
  gr_new.SetPoint(1,0.75,r2_new->GetParams()[2]);
  gr_new.SetPoint(2,1.5,r3_new->GetParams()[2]);
  gr_new.SetPoint(3,2.5,r4_new->GetParams()[2]);
  gr_new.SetPoint(4,3.5,r5_new->GetParams()[2]);
  gr_new.SetPoint(5,4.5,r6_new->GetParams()[2]);
  gr_new.SetPoint(6,7.5,r7_new->GetParams()[2]);
  gr_new.SetPointError(0,0.25,r1_new->GetErrors()[2]);
  gr_new.SetPointError(1,0.25,r2_new->GetErrors()[2]);
  gr_new.SetPointError(2,0.5,r3_new->GetErrors()[2]);
  gr_new.SetPointError(3,0.5,r4_new->GetErrors()[2]);
  gr_new.SetPointError(4,0.5,r5_new->GetErrors()[2]);
  gr_new.SetPointError(5,0.5,r6_new->GetErrors()[2]);
  gr_new.SetPointError(6,2.5,r7_new->GetErrors()[2]);
  gr_new.SetMarkerColor(kBlue);
  gr_new.SetLineColor(kBlue);
  
  //TF1 f_new("f_new","[0]+[1]/sqrt(x)",0.,10.);
  //auto r_new = gr_new.Fit("f_new","SQ","",0.,10.);
  
  TH2D hcmp("hcmp",";P_{t} GeV; #sigma(1-P_{t}^{true}/P_{t}^{reco})",1,0.,10.,1,0.02,0.07);
  hcmp.GetYaxis()->SetTitleOffset(1.5);
  hcmp.SetStats(false);
  hcmp.Draw();
  gr_old.Draw("P SAME");
  //c.Modified();
  //c.Update();
  //auto ps_old = static_cast<TPaveStats*>(gr_old.GetListOfFunctions()->FindObject("stats"));
  //ps_old->SetX1NDC(0.20); 
  //ps_old->SetX2NDC(0.50);
  gr_new.Draw("P SAME");
  //c.Modified();
  //c.Update();
  //auto ps_new = static_cast<TPaveStats*>(gr_new.GetListOfFunctions()->FindObject("stats"));
  //ps_new->SetX1NDC(0.55); 
  //ps_new->SetX2NDC(0.85);
  
  TLegend leg(0.7,0.85,0.95,0.95);
  leg.AddEntry(&gr_old,"circular fit","pl");
  leg.AddEntry(&gr_new,"kalman filter (GENFIT)","pl");
  leg.Draw();
  c.Print("ana.pdf");
  
  t->SetLineColor(kRed);
  t->SetMarkerColor(kRed);
  t->Draw("dip_old-dip_true>>hdip_old(500,-0.02,0.02)","","E0");
  auto hdip_old = static_cast<TH1F*>(gDirectory->Get("hdip_old"));
  hdip_old->SetTitle(";#Delta#alpha (rad)");
  auto rdip_old = hdip_old->Fit("gaus","QS","",-0.001,0.001);
  hdip_old->Draw("E0");
  c.Print("ana.pdf");
  
  t->SetLineColor(kBlue);
  t->SetMarkerColor(kBlue);
  t->Draw("dip_new-dip_true>>hdip_new(500,-0.02,0.02)","","E0");
  auto hdip_new = static_cast<TH1F*>(gDirectory->Get("hdip_new"));
  hdip_new->SetTitle(";#Delta#alpha (rad)");
  auto rdip_new = hdip_new->Fit("gaus","QS","",-0.001,0.001);
  hdip_new->Draw("E0");
  c.Print("ana.pdf)");
}