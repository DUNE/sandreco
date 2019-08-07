// Copyright 2019 Matteo Tenti

#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TDatabasePDG.h>

TGeoManager* geo = 0;
TTree* t = 0;
TFile* f = 0;

bool initialized = false;
  
TCut isCC = "EvtCode.String().Contains(\"CC\") != 0";
TCut isDIS = "EvtCode.GetString().Contains(\"DIS\") != 0";
TCut isQE = "EvtCode.GetString().Contains(\"QE\") != 0";
TCut isRES = "EvtCode.GetString().Contains(\"RES\") != 0";
TCut isMEC = "EvtCode.GetString().Contains(\"MEC\") != 0";
TCut isHMULT = "EvtCode.GetString().Contains(\"hmult\") != 0";
TCut isIMD = "EvtCode.GetString().Contains(\"IMD\") != 0";
TCut isNUEEL = "EvtCode.GetString().Contains(\"NuEEL\") != 0";
TCut Htgt = "EvtCode.GetString().Contains(\"1000010010\") != 0";
TCut Ctgt = "EvtCode.GetString().Contains(\"1000060120\") != 0";
TCut Ntgt = "EvtCode.GetString().Contains(\"1000070140\") != 0";
TCut Otgt = "EvtCode.GetString().Contains(\"1000080160\") != 0";
TCut Artgt = "EvtCode.GetString().Contains(\"1000180400\") != 0";
TCut isLAr = "isLAr!=0";
TCut isStt = "isStt!=0";
TCut isMu = "particles.pdg==13";
TCut isPi = "abs(particles.pdg)==211";
TCut isPi0 = "particles.pdg==111";
TCut isP = "particles.pdg==2212";
TCut isN = "particles.pdg==2112";
TCut isGamma = "particles.pdg==22";
TCut isPrimary = "particles.primary==1";
TCut isGoodTrack = "particles.tr.@digits.size()>=3";

double p_loc[] = {0, 0, 0};
double p_mst[] = {0, 0, 0};

void init(const char* ifile)
{

  f = new TFile(ifile);
  TTree* tEvent = reinterpret_cast<TTree*>(f->Get("tEvent"));
  TTree* tReco = reinterpret_cast<TTree*>(f->Get("tReco"));
  TTree* tDigit = reinterpret_cast<TTree*>(f->Get("tDigit"));
  TTree* tEdep = reinterpret_cast<TTree*>(f->Get("EDepSimEvents"));
  TTree* tGenie = reinterpret_cast<TTree*>(f->Get("gRooTracker"));

  tEvent->AddFriend(tReco);
  tEvent->AddFriend(tDigit);
  tEvent->AddFriend(tEdep);
  tEvent->AddFriend(tGenie);

  t = tEvent;

  geo = reinterpret_cast<TGeoManager*>(f->Get("EDepSimGeometry"));
  
  if(f != 0 && t != 0 && geo != 0)
  {  
    geo->cd(
      "volWorld_PV/volDetEnclosure_PV_0/volKLOEFULLECALSENSITIVE_EXTTRK_NEWGAP_PV_0"
      "/KLOETrackingRegion_volume_PV_0/volKLOESTTFULLNEWCONFNEWGAPLAR_PV_0");
  
    geo->LocalToMaster(p_loc, p_mst);
    
    TGeoTubeSeg* lar = reinterpret_cast<TGeoTubeSeg*>(geo->GetVolume("LArTarget_volume_PV")->GetShape());
    double lar_rmin = lar->GetRmin();
    double lar_rmax = lar->GetRmax();
    double lar_dz = lar->GetDz();
    double phi1 = lar->GetPhi1();
    double phi2 = lar->GetPhi2();
    
    phi1 = phi1 > 180. ? phi1 - 360. : phi1;
    phi2 = phi2 > 180. ? phi2 - 360. : phi2;
    
    phi1 *= TMath::Pi()/180.;
    phi2 *= TMath::Pi()/180.;
    
    TGeoTube* stt = reinterpret_cast<TGeoTube*>(geo->GetVolume("StrawTubeFULL_for_STTPlane2FULL_0_0_lv_PV")->GetShape());
    double stt_dz = stt->GetDz();
    
    t->SetAlias("isLAr",TString::Format("r_kloe > %f && r_kloe < %f && abs(x - %f) < %f &&"
      "(TMath::ATan2(y - %f, z - %f) > %f || TMath::ATan2(y - %f, z - %f) < %f)", 
      lar_rmin, lar_rmax, p_mst[0], lar_dz, p_mst[1], p_mst[2], phi1, p_mst[1], p_mst[2], phi2).Data());
    
    t->SetAlias("isStt",TString::Format("r_kloe < %f && abs(x - %f) < %f", lar_rmin, p_mst[0], stt_dz).Data());
    
    t->SetAlias("pnu","sqrt(pxnu*pxnu+pynu*pynu+pznu*pznu)");
    t->SetAlias("r_kloe",TString::Format("sqrt((y - %f)*(y - %f)+(z - %f)*(z - %f))",p_mst[1],p_mst[1],p_mst[2],p_mst[2]));
    t->SetAlias("ppara_reco","particles.pxreco");
    t->SetAlias("ppara_true","particles.pxtrue");
    t->SetAlias("pperp_reco","sqrt(particles.pyreco*particles.pyreco+particles.pzreco*particles.pzreco)");
    t->SetAlias("pperp_true","sqrt(particles.pytrue*particles.pytrue+particles.pztrue*particles.pztrue)");
    t->SetAlias("dipangle_reco","TMath::ATan2(ppara_reco,pperp_reco)");
    t->SetAlias("dipangle_true","TMath::ATan2(ppara_true,pperp_true)");
    
    initialized = true;
  }
}

void make_plots(const char* ofile) 
{
  if(!initialized)
  {
    std::cout << "Tree not initialized" << std::endl;
    return;
  }

  gROOT->SetBatch();
  gStyle->SetOptFit(1111);
  
  TString pdffile = ofile;
  pdffile.ReplaceAll(".root",".pdf");

  TCanvas c;
  
  // X vertex
  t->SetLineColor(kRed);
  t->Draw("x>>hx(100,-1800,1800)",isStt,"");
  t->SetLineColor(kBlue);
  t->Draw("x",isLAr,"SAME");
  t->SetLineColor(kBlack);
  t->Draw("x",!isLAr&&!isStt,"SAME");
  c.Print(TString::Format("%s(",pdffile.Data()));
  
  // YZ vertex
  t->SetMarkerColor(kRed);
  t->Draw(TString::Format("y:z>>hyz(100,%f-2000.,%f+2000.,100,%f-2000.,%f+2000.)",p_mst[2],p_mst[2],p_mst[1],p_mst[1]),isStt,"");
  t->SetMarkerColor(kBlue);
  t->Draw("y:z",isLAr,"SAME");
  t->SetMarkerColor(kBlack);
  t->Draw("y:z",!isLAr&&!isStt,"SAME");
  c.Print(pdffile.Data());
  
  // nu energy
  c.SetLogy(true);
  t->Draw("Enu>>henu(1000,0,50000)","","");
  c.Print(pdffile.Data());
  
  // nu direction
  c.SetLogy(false);
  t->Draw("pxnu/pznu:pynu/pznu>>htxty(10,-0.11,-0.1,10,-0.01,0.01)","","COLZ TEXT");
  c.Print(pdffile.Data());
  
  // CC vs NC
  TH1D hCC("hCC","CC/NC",2,0.,2.);
  t->Draw("1>>+hCC",isCC);
  t->Draw("0>>+hCC",!isCC);
  hCC.GetXaxis()->CenterLabels();
  hCC.GetXaxis()->SetBinLabel(1,"NC");
  hCC.GetXaxis()->SetBinLabel(2,"CC");
  hCC.Draw("HIST TEXT");
  c.Print(pdffile.Data());
  
  // interaction type
  TH1D hintType("hintType","interaction type",7,0.,7.);
  t->Draw("0>>+hintType",isDIS);
  t->Draw("1>>+hintType",isQE);
  t->Draw("2>>+hintType",isRES);
  t->Draw("3>>+hintType",isMEC);
  t->Draw("4>>+hintType",isHMULT);
  t->Draw("5>>+hintType",isIMD);
  t->Draw("6>>+hintType",isNUEEL);
  hintType.GetXaxis()->CenterLabels();
  hintType.GetXaxis()->SetBinLabel(1,"DIS");
  hintType.GetXaxis()->SetBinLabel(2,"QE");
  hintType.GetXaxis()->SetBinLabel(3,"RES");
  hintType.GetXaxis()->SetBinLabel(4,"MEC");
  hintType.GetXaxis()->SetBinLabel(5,"HMULT");
  hintType.GetXaxis()->SetBinLabel(6,"IMD");
  hintType.GetXaxis()->SetBinLabel(7,"NUEEL");
  hintType.Draw("HIST TEXT");
  c.Print(pdffile.Data());
  
  // target nuclei
  TH1D htgt("htgt","target nuclei",5,0.,5.);
  t->Draw("0>>+htgt",Htgt);
  t->Draw("1>>+htgt",Ctgt);
  t->Draw("2>>+htgt",Ntgt);
  t->Draw("3>>+htgt",Otgt);
  t->Draw("4>>+htgt",Artgt);
  htgt.GetXaxis()->CenterLabels();
  htgt.GetXaxis()->SetBinLabel(1,"H");
  htgt.GetXaxis()->SetBinLabel(2,"C");
  htgt.GetXaxis()->SetBinLabel(3,"N");
  htgt.GetXaxis()->SetBinLabel(4,"O");
  htgt.GetXaxis()->SetBinLabel(5,"Ar");
  htgt.Draw("HIST TEXT");
  c.Print(pdffile.Data());
  
  // detector
  TH1D hdet("hdet","detector",3,0.,3.);
  t->Draw("0>>+hdet",isLAr);
  t->Draw("1>>+hdet",isStt);
  t->Draw("2>>+hdet",!isStt&&!isLAr);
  hdet.GetXaxis()->CenterLabels();
  hdet.GetXaxis()->SetBinLabel(1,"LAr");
  hdet.GetXaxis()->SetBinLabel(2,"Stt");
  hdet.GetXaxis()->SetBinLabel(3,"none");
  hdet.Draw("HIST TEXT");
  c.Print(pdffile.Data());
  
  // particles
  int npart = t->Draw("Primaries.PDGCode","","");
  double* pdg = t->GetV1();
  
  std::map<int,int> pdgpart;
  
  for(int i = 0; i < npart; i++)
  {
    pdgpart[static_cast<int>(pdg[i])]++;
  }
  
  int npdg = pdgpart.size();
  TH1D hpart("hpart","particles",npdg, 0., npdg);
  int ibin = 0;
  TDatabasePDG db;
  TParticlePDG* p;
  
  for(std::map<int,int>::iterator it=pdgpart.begin(); it!=pdgpart.end(); ++it)
  {
    hpart.SetBinContent(++ibin, it->second);
    p = db.GetParticle(it->first);
    if(p)
      hpart.GetXaxis()->SetBinLabel(ibin,p->GetName());
    else
      hpart.GetXaxis()->SetBinLabel(ibin,TString::Format("%d",it->first).Data());
  }
  hpart.Draw("HIST TEXT");
  c.Print(pdffile.Data());  
  
  // muon
  TH1D hmu_pperp_res("hmu_pperp_res","#sigma(P_{perp}) [#mu] (Stt); p_{t}/p_{r}-1",100,-1,1);
  t->Draw("pperp_true/pperp_reco-1>>hmu_pperp_res",isMu && isPrimary && isStt,"E0");
  hmu_pperp_res.Fit("gaus","Q","",-0.05,0.1);
  c.Print(pdffile.Data()); 
  
  TH2D hmu_pperp_res_vs_pperp("hmu_pperp_res_vs_pperp","#sigma(P_{perp}) vs P_{perp} [#mu] (Stt); p_{t}; p_{t}/p_{r}-1",100,0,10000,100,-1,1);
  t->Draw("pperp_true/pperp_reco-1:pperp_true>>hmu_pperp_res_vs_pperp",isMu && isPrimary && isStt,"COLZ");
  c.Print(pdffile.Data()); 
  
  TH2D hmu_pperp_res_vs_ndigits("hmu_pperp_res_vs_ndigits","#sigma(P_{perp}) vs n_digits [#mu] (Stt); n_digits; p_{t}/p_{r}-1",500,0,500,100,-1,1);
  t->Draw("pperp_true/pperp_reco-1:particles.tr.@digits.size()>>hmu_pperp_res_vs_ndigits",isMu && isPrimary && isStt,"COLZ");
  c.Print(pdffile.Data()); 
  
  /*
  TH1D hmu_dip_LAr("hmu_dip_LAr","dip angle  [#mu](LAr); #Delta#rho (rad)",50,-0.2,0.2);
  t->Draw("TMath::ATan2(ppara_reco,pperp_reco)-TMath::ATan2(ppara_true,pperp_true)>>hmu_dip_LAr",isMu && isPrimary && isLAr,"E0");
  hmu_dip_LAr.Fit("gaus","Q","",-0.02,0.02);
  c.Print(pdffile.Data());
  
  TH1D hmu_pperp_Stt("hmu_pperp_Stt","P_{perp} [#mu] (STT); p_t/p_r-1",100,-1,1);
  t->Draw("pperp_true/pperp_reco-1>>hmu_pperp_Stt",isMu && isPrimary && isStt,"E0");
  hmu_pperp_Stt.Fit("gaus","Q","",-0.05,0.1);
  c.Print(pdffile.Data()); 
  
  TH1D hmu_dip_Stt("hmu_dip_Stt","dip angle [#mu] (STT); #Delta#rho (rad)",50,-0.1,0.1);
  t->Draw("TMath::ATan2(ppara_reco,pperp_reco)-TMath::ATan2(ppara_true,pperp_true)>>hmu_dip_Stt",isMu && isPrimary && isStt,"E0");
  hmu_dip_Stt.Fit("gaus","Q","",-0.02,0.02);
  c.Print(pdffile.Data());
  
  TH1D hmu_charge("hmu_charge","charge mis-indentification [#mu] (STT)",2,-1.5,1.5);
  t->Draw("particles.charge*particles.charge_reco>>hmu_charge",isMu && isPrimary && isStt,"");
  hmu_charge.GetXaxis()->CenterLabels();
  hmu_charge.GetXaxis()->SetBinLabel(1,"WRONG");
  hmu_charge.GetXaxis()->SetBinLabel(2,"RIGHT");
  hmu_charge.Draw("HIST TEXT");
  c.Print(pdffile.Data());  
  
  // pi
  
  TH1D hpi_pperp_Stt("hpi_pperp_Stt","P_{perp} [#pi] (Stt); p_t/p_r-1",100,-1,1);
  t->Draw("pperp_true/pperp_reco-1>>hpi_pperp_Stt",isPi && isPrimary && isStt,"E0");
  hpi_pperp_Stt.Fit("gaus","Q","",-0.05,0.1);
  c.Print(pdffile.Data());  
  
  TH1D hpi_dip_Stt("hpi_dip_Stt","dip angle [#pi] (Stt); #Delta#rho (rad)",50,-0.2,0.2);
  t->Draw("TMath::ATan2(ppara_reco,pperp_reco)-TMath::ATan2(ppara_true,pperp_true)>>hpi_dip_Stt",isPi && isPrimary && isStt,"E0");
  hpi_dip_Stt.Fit("gaus","Q","",-0.01,0.01);
  c.Print(pdffile.Data());
  
  TH1D hpi_charge("hpi_charge","charge mis-indentification [#pi] (STT)",2,-1.5,1.5);
  t->Draw("particles.charge*particles.charge_reco>>hpi_charge",isPi && isPrimary && isStt,"");
  hpi_charge.GetXaxis()->CenterLabels();
  hpi_charge.GetXaxis()->SetBinLabel(1,"WRONG");
  hpi_charge.GetXaxis()->SetBinLabel(2,"RIGHT");
  hpi_charge.Draw("HIST TEXT");
  c.Print(pdffile.Data());
  
  // p
  
  TH1D hp_pperp_Stt("hp_pperp_Stt","P_{perp} [p] (Stt); p_t/p_r-1",50,-0.8,0.8);
  t->Draw("pperp_true/pperp_reco-1>>hp_pperp_Stt",isP && isPrimary && isStt,"E0");
  hp_pperp_Stt.Fit("gaus","Q","",-0.2,0.2);
  c.Print(pdffile.Data());  
  
  TH1D hp_dip_Stt("hp_dip_Stt","dip angle [p] (Stt); #Delta#rho (rad)",100,-0.4,0.4);
  t->Draw("TMath::ATan2(ppara_reco,pperp_reco)-TMath::ATan2(ppara_true,pperp_true)>>hp_dip_Stt",isP && isPrimary && isStt,"E0");
  hp_dip_Stt.Fit("gaus","Q","",-0.02,0.02);
  c.Print(pdffile.Data());
  
  TH1D hp_charge("hp_charge","charge mis-indentification [p] (STT)",2,-1.5,1.5);
  t->Draw("particles.charge*particles.charge_reco>>hp_charge",isP && isPrimary && isStt,"");
  hp_charge.GetXaxis()->CenterLabels();
  hp_charge.GetXaxis()->SetBinLabel(1,"WRONG");
  hp_charge.GetXaxis()->SetBinLabel(2,"RIGHT");
  hp_charge.Draw("HIST TEXT");
  c.Print(pdffile.Data());  
  
  t->Draw("pperp_reco:pperp_true >> hp(100,0,3000,100,0,3000)",isP && isPrimary && isStt,"colz");
  c.Print(pdffile.Data());
  */
  
  c.Print(TString::Format("%s)",pdffile.Data()));
}
