#include <TROOT.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TGraph.h>
#include <TCut.h>

#include <iostream>

void help()
{
  std::cout << "FastCheck <input> <output>" << std::endl;
  std::cout << "<input>  is input root file" << std::endl;
  std::cout << "<output> is pdf output file" << std::endl;
}

void print(TCanvas& c, TString fout, TTree* t, const char* var, TCut sel, const char* opt, const char* title, int firstlast = 0, long int nmax = 10000)
{
  t->Draw(var,sel,opt,nmax);
  TH1D* h1 = (TH1D*) gROOT->FindObject("htemp");
  h1->SetStats(false);
  h1->SetTitle(title);
  h1->Draw();
  if(firstlast == -1)
    fout += "(";
  else if(firstlast == 1)
    fout += ")";
  c.SaveAs(fout.Data());
}

int main(int argc, char* argv[])
{
  gROOT->SetBatch(true);
  
  TCanvas c;
  TGraph *g;
  TH1D* h1;
  TH2D* h2;
  int n;

  if(argc != 3)
  {
    help();
    exit(-1);
  }
  TFile fin(argv[1]);
  TString fout(argv[2]);
  
  TTree* tDigit = (TTree*) fin.Get("tDigit");
  TTree* tReco = (TTree*) fin.Get("tReco");
  TTree* tEvent = (TTree*) fin.Get("tEvent");
  
  TCut barrel_module = "cell.mod < 24";
  TCut endcap_module = "cell.mod > 24";
  TCut up_barrel_module = "cell.mod==0";
  TCut left_endcap_module = "cell.mod==30";
  TCut trackRecoOK = "track.ret_ln==0&&track.ret_cr==0";
  TCut EnuRecoOK = "Enureco>0.";
  TCut PrimaryPart = "particles.primary==1";
  TCut PartTrackRecoOK = "particles.tr.ret_ln==0&&particles.tr.ret_cr==0";
  
  if(tDigit)
  {  
    print(c, fout, tDigit, "cell.mod", "", "", "; module id", -1);
    print(c, fout, tDigit, "cell.lay", "", "", "; layer id");
    print(c, fout, tDigit, "cell.cel", barrel_module, "", "barrel modules; cell id");
    print(c, fout, tDigit, "cell.l", barrel_module, "", "barrel modules; cell length (mm)");
    print(c, fout, tDigit, "cell.cel", endcap_module, "", "endcap modules; cell id");
    print(c, fout, tDigit, "cell.l", endcap_module, "", "endcap modules; cell length (mm)");
    
    n = tDigit->Draw("cell.y:cell.x", "", "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1()); 
    g->SetTitle("; x cell position (mm); y cell position (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("cell.y:cell.z", "", "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("; z cell position (mm); y cell position (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("cell.y:cell.z", up_barrel_module, "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("module id == 0; z cell position (mm); y cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("TMath::ATan2(23910.000-cell.z,cell.y+2384.7300)/TMath::Pi()*180.:cell.mod", barrel_module, "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("; module id; angle from y axis (rad)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("cell.y:cell.lay", up_barrel_module, "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("module id == 0; lay id; y cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("cell.z:cell.cel", up_barrel_module, "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("module id == 0; cel id; z cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("cell.x:cell.z", left_endcap_module, "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("module id == 30; z cell position (mm); x cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("cell.x:cell.mod", endcap_module, "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("; module id; x cell position (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    c.SetLogy(true);
    print(c, fout, tDigit, "TMath::Log10(cell.adc1)", "", "", "; log_{10}(adc1)");
    print(c, fout, tDigit, "TMath::Log10(cell.adc2)", "", "", "; log_{10}(adc2)");
    print(c, fout, tDigit, "TMath::Log10(cell.tdc1)", "", "", "; log_{10}(tdc1/ns)");
    print(c, fout, tDigit, "TMath::Log10(cell.tdc2)", "", "", "; log_{10}(tdc2/ns)");
    print(c, fout, tDigit, "TMath::Log10(cell.pe_time1)", "", "", "; log_{10}(p.e. time1/ns)");
    print(c, fout, tDigit, "TMath::Log10(cell.pe_time2)", "", "", "; log_{10}(p.e. time2/ns)");
    c.SetLogy(false);
    
    n = tDigit->Draw("Stt.y:Stt.x", "", "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1()); 
    g->SetTitle("; x stt digit position (mm); y stt digit position (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("Stt.y:Stt.z", "", "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1()); 
    g->SetTitle("; z stt digit position (mm); y stt digit position (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    c.SetLogy(true);
    print(c, fout, tDigit, "TMath::Log10(Stt.t)", "", "", "; log_{10}(stt time/ns)");
    print(c, fout, tDigit, "TMath::Log10(Stt.de)", "", "", "; log_{10}(dE/MeV)");
    c.SetLogy(false);
    print(c, fout, tDigit, "Stt.hor", "", "", "; hor/ver");
  }
  
  if(tReco)
  {
    tReco->Draw("track.ret_ln:track.ret_cr", "", "colz"); 
    h2 = (TH2D*) gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle(";ret_cr;ret_ln");
    h2->Draw("colztext");
    c.SaveAs(fout.Data());
    
    c.SetLogy(true);
    print(c, fout, tReco, "TMath::Log10(track.r)", trackRecoOK, "", "; log_{10}(R/mm)");
    print(c, fout, tReco, "track.a", trackRecoOK, "", "; a (mm)");
    print(c, fout, tReco, "track.b", trackRecoOK, "", "; b");
    c.SetLogy(false);
    print(c, fout, tReco, "track.h", trackRecoOK, "", "; h");
    
    n = tReco->Draw("track.y0:track.x0", trackRecoOK, "goff"); 
    g = new TGraph(n,tReco->GetV2(),tReco->GetV1()); 
    g->SetTitle("; x0 (mm); y0 (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    n = tReco->Draw("track.y0:track.z0", trackRecoOK, "goff"); 
    g = new TGraph(n,tReco->GetV2(),tReco->GetV1()); 
    g->SetTitle("; z0 (mm); y0 (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    print(c, fout, tReco, "track.t0", trackRecoOK, "", "; t (ns)");
    print(c, fout, tReco, "TMath::Log10(track.chi2_cr)", trackRecoOK, "", "; log_{10}(#chi^{2}_{cr})");
    print(c, fout, tReco, "TMath::Log10(track.chi2_ln)", trackRecoOK, "", "; log_{10}(#chi^{2}_{ln})");
    
    n = tReco->Draw("cluster.y:cluster.x", "", "goff"); 
    g = new TGraph(n,tReco->GetV2(),tReco->GetV1()); 
    g->SetTitle("; x (mm); y (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    n = tReco->Draw("cluster.y:cluster.z", "", "goff"); 
    g = new TGraph(n,tReco->GetV2(),tReco->GetV1()); 
    g->SetTitle("; z (mm); y (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    print(c, fout, tReco, "TMath::Log10(cluster.t)", "", "", "; log_{10}(t/ns)");
    print(c, fout, tReco, "TMath::Log10(cluster.e)", "", "", "; log_{10}(E/MeV)");
    print(c, fout, tReco, "cluster.sx", "", "", "; sx");
    print(c, fout, tReco, "cluster.sy", "", "", "; sy");
    print(c, fout, tReco, "cluster.sz", "", "", "; sz");
    print(c, fout, tReco, "cluster.varx", "", "", "; varx");
    print(c, fout, tReco, "cluster.vary", "", "", "; vary");
    print(c, fout, tReco, "cluster.varz", "", "", "; varz");
  }
  
  if(tEvent)
  { 
    n = tEvent->Draw("y:x", "", "goff"); 
    g = new TGraph(n,tEvent->GetV2(),tEvent->GetV1()); 
    g->SetTitle("; x (mm); y (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    n = tEvent->Draw("y:z", "", "goff"); 
    g = new TGraph(n,tEvent->GetV2(),tEvent->GetV1()); 
    g->SetTitle("; z (mm); y (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    print(c, fout, tEvent, "t", "", "", "; t (ns)");
    print(c, fout, tEvent, "Enu", "", "", "neutrino; E (MeV)");
    print(c, fout, tEvent, "pxnu", "", "", "neutrino; px (MeV)");
    print(c, fout, tEvent, "pynu", "", "", "neutrino; py (MeV)");
    print(c, fout, tEvent, "pznu", "", "", "neutrino; pz (MeV)");
    print(c, fout, tEvent, "Enureco>0", "", "", "; reconstructed");
    print(c, fout, tEvent, "Enureco", EnuRecoOK, "", "neutrino; Ereco (MeV)");
    print(c, fout, tEvent, "pxnureco", EnuRecoOK, "", "neutrino; pxreco (MeV)");
    print(c, fout, tEvent, "pynureco", EnuRecoOK, "", "neutrino; pyreco (MeV)");
    print(c, fout, tEvent, "pznureco", EnuRecoOK, "", "neutrino; pzreco (MeV)");
    print(c, fout, tEvent, "particles.primary", "", "", "; primary");
    
    tEvent->Draw("Enureco:Enu>>h2(100,0,10000,100,0,10000)", EnuRecoOK, "colz"); 
    h2 = (TH2D*) gROOT->FindObject("h2");
    h2->SetStats(false);
    h2->SetTitle("; E (MeV); Ereco (MeV)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());
    
    c.SetLogy(true);    
    tEvent->Draw("particles.pdg>>hpdg(8000,-4000,4000)", "", ""); 
    h1 = (TH1D*) gROOT->FindObject("hpdg");
    h1->SetStats(false);
    h1->SetTitle("; pdg");
    h1->Draw();
    c.SaveAs(fout.Data());
    
    tEvent->Draw("particles.pdg>>hpdg_pri(8000,-4000,4000)", PrimaryPart, ""); 
    h1 = (TH1D*) gROOT->FindObject("hpdg_pri");
    h1->SetStats(false);
    h1->SetTitle("primaries ; pdg");
    h1->Draw();
    c.SaveAs(fout.Data());
    c.SetLogy(false);
    
    print(c, fout, tEvent, "particles.tr.ret_ln==0&&particles.tr.ret_cr==0", PrimaryPart, "", "primaries; recoOK==1");    
    
    tEvent->Draw("particles.charge_reco:particles.charge", PrimaryPart && PartTrackRecoOK, "colztext"); 
    h2 = (TH2D*) gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true charge; reco charge");
    h2->Draw("colztext");
    c.SaveAs(fout.Data());
    
    tEvent->Draw("particles.pxreco:particles.pxtrue", PrimaryPart && PartTrackRecoOK, "colz"); 
    h2 = (TH2D*) gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true px (MeV); reco px (MeV)");
    h2->Draw("colz");
    
    tEvent->Draw("particles.pyreco:particles.pytrue", PrimaryPart && PartTrackRecoOK, "colz"); 
    h2 = (TH2D*) gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true py (MeV); reco py (MeV)");
    h2->Draw("colz");
    
    tEvent->Draw("particles.pzreco:particles.pztrue", PrimaryPart && PartTrackRecoOK, "colz"); 
    h2 = (TH2D*) gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true pz (MeV); reco pz (MeV)");
    h2->Draw("colz");
    
    tEvent->Draw("particles.Ereco:particles.Etrue", PrimaryPart && PartTrackRecoOK, "colz"); 
    h2 = (TH2D*) gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true E (MeV); reco E (MeV)");
    h2->Draw("colz");
    
    tEvent->Draw("particles.xreco:particles.xtrue", PrimaryPart && PartTrackRecoOK, "colz"); 
    h2 = (TH2D*) gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true x (mm); reco x (mm)");
    h2->Draw("colz");
    
    tEvent->Draw("particles.yreco:particles.ytrue", PrimaryPart && PartTrackRecoOK, "colz"); 
    h2 = (TH2D*) gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true y (MeV); reco y (mm)");
    h2->Draw("colz");
    
    tEvent->Draw("particles.zreco:particles.ztrue", PrimaryPart && PartTrackRecoOK, "colz"); 
    h2 = (TH2D*) gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true z (mm); reco z (mm)");
    h2->Draw("colz");
    
    tEvent->Draw("particles.treco:particles.ttrue", PrimaryPart && PartTrackRecoOK, "colz"); 
    h2 = (TH2D*) gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true t (ns); reco t (ns)");
    h2->Draw("colz");
    
    
  }
  
  c.Clear();
  c.SaveAs(TString::Format("%s)",fout.Data()).Data());
  return 0;
}