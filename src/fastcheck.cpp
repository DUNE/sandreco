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

void print(TCanvas& c, TString fout, TTree* t, const char* var, TCut sel,
           const char* opt, const char* title, int firstlast = 0,
           long int nmax = 10000)
{
  t->Draw(var, sel, opt, nmax);
  TH1D* h1 = (TH1D*)gROOT->FindObject("htemp");
  h1->SetStats(false);
  h1->SetTitle(title);
  h1->Draw();
  if (firstlast == -1)
    fout += "(";
  else if (firstlast == 1)
    fout += ")";
  c.SaveAs(fout.Data());
}

int main(int argc, char* argv[])
{
  gROOT->SetBatch(true);

  TCanvas c;
  TGraph* g;
  TH1D* h1;
  TH2D* h2;
  int n;

  if (argc != 3) {
    help();
    exit(-1);
  }
  TFile fin(argv[1]);
  TString fout(argv[2]);

  TTree* tDigit = (TTree*)fin.Get("tDigit");
  TTree* tReco = (TTree*)fin.Get("tReco");
  TTree* tEvent = (TTree*)fin.Get("tEvent");

  TCut barrel_module = "cell.mod < 24";
  TCut endcap_module = "cell.mod > 24";
  TCut up_barrel_module = "cell.mod==0";
  TCut left_endcap_module = "cell.mod==30";
  TCut cellOK = "cell.adc1>0&&cell.adc2>0";
  TCut trackRecoOK = "track.ret_ln==0&&track.ret_cr==0";
  TCut EnuRecoOK = "Enureco>0.";
  TCut PrimaryPart = "particles.primary==1";
  TCut PartTrackRecoOK = "particles.tr.ret_ln==0&&particles.tr.ret_cr==0";

  std::cout << "tDigit histograms..." << std::flush;

  if (tDigit) {
    print(c, fout, tDigit, "cell.mod", "", "", "cells; module id", -1);
    print(c, fout, tDigit, "cell.lay", "", "", "cells; layer id");
    print(c, fout, tDigit, "cell.cel", barrel_module, "",
          "barrel modules; cell id");
    print(c, fout, tDigit, "cell.l", barrel_module, "",
          "barrel modules; cell length (mm)");
    print(c, fout, tDigit, "cell.cel", endcap_module, "",
          "endcap modules; cell id");
    print(c, fout, tDigit, "cell.l", endcap_module, "",
          "endcap modules; cell length (mm)");

    n = tDigit->Draw("cell.y:cell.x", "", "goff");
    g = new TGraph(n, tDigit->GetV2(), tDigit->GetV1());
    g->SetTitle("cells; x cell position (mm); y cell position (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    n = tDigit->Draw("cell.y:cell.z", "", "goff");
    g = new TGraph(n, tDigit->GetV2(), tDigit->GetV1());
    g->SetTitle("cells; z cell position (mm); y cell position (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    n = tDigit->Draw("cell.y:cell.z", up_barrel_module, "goff");
    g = new TGraph(n, tDigit->GetV2(), tDigit->GetV1());
    g->SetTitle("module id == 0; z cell position (mm); y cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());

    n = tDigit->Draw(
        "TMath::ATan2(23910.000-cell.z,cell.y+2384.7300)/"
        "TMath::Pi()*180.:cell.mod",
        barrel_module, "goff");
    g = new TGraph(n, tDigit->GetV2(), tDigit->GetV1());
    g->SetTitle("barrel modules; module id; angle from y axis (rad)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    n = tDigit->Draw("cell.y:cell.lay", up_barrel_module, "goff");
    g = new TGraph(n, tDigit->GetV2(), tDigit->GetV1());
    g->SetTitle("module id == 0; lay id; y cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());

    n = tDigit->Draw("cell.z:cell.cel", up_barrel_module, "goff");
    g = new TGraph(n, tDigit->GetV2(), tDigit->GetV1());
    g->SetTitle("module id == 0; cel id; z cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());

    n = tDigit->Draw("cell.x:cell.z", left_endcap_module, "goff");
    g = new TGraph(n, tDigit->GetV2(), tDigit->GetV1());
    g->SetTitle("module id == 30; z cell position (mm); x cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());

    n = tDigit->Draw("cell.x:cell.mod", endcap_module, "goff");
    g = new TGraph(n, tDigit->GetV2(), tDigit->GetV1());
    g->SetTitle("endcap module; module id; x cell position (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    c.SetLogy(true);
    print(c, fout, tDigit, "TMath::Log10(cell.adc1)", "", "",
          "cells; log_{10}(adc1)");
    print(c, fout, tDigit, "TMath::Log10(cell.adc2)", "", "",
          "cells; log_{10}(adc2)");
    c.SetLogy(false);

    c.SetLogz(true);
    tDigit->Draw("adc2:adc1>>htemp(200,0,2000,200,0,2000)", "", "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("cells; adc1; adc2");
    h2->Draw("colz");
    c.SaveAs(fout.Data());
    c.SetLogz(false);

    tDigit->Draw("(adc2>0):(adc1>0)>>htemp(2,-0.5,1.5,2,-0.5,1.5)", "",
                 "colztext");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("cells; adc1>0; adc2>0");
    h2->Draw("colztext");
    c.SaveAs(fout.Data());

    c.SetLogy(true);
    print(c, fout, tDigit, "TMath::Log10(cell.tdc1)", "cell.adc1>0", "",
          "OK cells; log_{10}(tdc1/ns)");
    print(c, fout, tDigit, "TMath::Log10(cell.tdc2)", "cell.adc2>0", "",
          "OK cells; log_{10}(tdc2/ns)");
    print(c, fout, tDigit, "TMath::Log10(abs(cell.tdc1-cell.tdc2))", "", "",
          "OK cells; log_{10}(|tdc1-tdc2|/ns)");
    print(c, fout, tDigit, "TMath::Log10(cell.tdc1+cell.tdc2)", "", "",
          "OK cells; log_{10}((tdc1+tdc2)/ns)");
    c.SetLogy(false);

    c.SetLogz(true);
    tDigit->Draw("adc2:adc1>>htemp(100,0,5000,100,0,5000)", cellOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("cells; adc1; adc2");
    h2->Draw("colz");
    c.SaveAs(fout.Data());
    c.SetLogz(false);

    c.SetLogz(true);
    tDigit->Draw("tdc2:tdc1>>htemp(100,0,100,100,0,100)", cellOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("cells; tdc1 (ns); tdc2 (ns)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());
    c.SetLogz(false);

    c.SetLogy(true);
    print(c, fout, tDigit, "TMath::Log10(cell.@pe_time1.size())", "", "",
          "cells; log_{10}(#p.e.)");
    print(c, fout, tDigit, "TMath::Log10(cell.@pe_time2.size())", "", "",
          "cells; log_{10}(#p.e.)");
    print(c, fout, tDigit, "TMath::Log10(cell.pe_time1)", "", "",
          "cells; log_{10}(p.e. time1/ns)");
    print(c, fout, tDigit, "TMath::Log10(cell.pe_time2)", "", "",
          "cells; log_{10}(p.e. time2/ns)");
    c.SetLogy(false);

    n = tDigit->Draw("Stt.y:Stt.x", "", "goff", 1000);
    g = new TGraph(n, tDigit->GetV2(), tDigit->GetV1());
    g->SetTitle("stt; x stt digit position (mm); y stt digit position (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    n = tDigit->Draw("Stt.y:Stt.z", "", "goff", 1000);
    g = new TGraph(n, tDigit->GetV2(), tDigit->GetV1());
    g->SetTitle("stt; z stt digit position (mm); y stt digit position (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    n = tDigit->Draw("Stt.x:Stt.z", "", "goff", 1000);
    g = new TGraph(n, tDigit->GetV2(), tDigit->GetV1());
    g->SetTitle("stt; z stt digit position (mm); x stt digit position (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    c.SetLogy(true);
    print(c, fout, tDigit, "TMath::Log10(Stt.t)", "", "",
          "stt; log_{10}(stt time/ns)");
    print(c, fout, tDigit, "TMath::Log10(Stt.de)", "", "",
          "stt; log_{10}(dE/MeV)");
    c.SetLogy(false);
    print(c, fout, tDigit, "Stt.hor", "", "", "stt; hor/ver");
  }

  std::cout << " done" << std::endl;

  std::cout << "tReco  histograms..." << std::flush;

  if (tReco) {
    tReco->Draw("track.ret_ln:track.ret_cr", "", "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("tracks;ret_cr;ret_ln");
    h2->Draw("colztext");
    c.SaveAs(fout.Data());

    c.SetLogy(true);
    print(c, fout, tReco, "TMath::Log10(track.r)", trackRecoOK, "",
          "tracks; log_{10}(R/mm)");
    c.SetLogy(false);

    n = tReco->Draw("track.yc:track.zc", trackRecoOK, "goff");
    g = new TGraph(n, tReco->GetV2(), tReco->GetV1());
    g->SetTitle("tracks; zc (mm); yc (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    c.SetLogy(true);
    print(c, fout, tReco, "track.a", trackRecoOK, "", "tracks; a (mm)");
    print(c, fout, tReco, "track.b", trackRecoOK, "", "tracks; b");
    c.SetLogy(false);
    print(c, fout, tReco, "track.h", trackRecoOK, "", "tracks; h");

    n = tReco->Draw("track.y0:track.x0", trackRecoOK, "goff");
    g = new TGraph(n, tReco->GetV2(), tReco->GetV1());
    c.DrawFrame(-2500, -5000, 2500, 0, "tracks; x0 (mm); y0 (mm)");
    g->SetTitle("tracks; x0 (mm); y0 (mm)");
    g->Draw("psame");
    c.SaveAs(fout.Data());

    n = tReco->Draw("track.y0:track.z0", trackRecoOK, "goff");
    g = new TGraph(n, tReco->GetV2(), tReco->GetV1());
    c.DrawFrame(21500, -5000, 26500, 0, "tracks; z0 (mm); y0 (mm)");
    g->SetTitle("tracks; z0 (mm); y0 (mm)");
    g->Draw("psame");
    c.SaveAs(fout.Data());

    print(c, fout, tReco, "TMath::Log10(track.t0)", trackRecoOK, "",
          "tracks; log_{10}(t/ns)");
    print(c, fout, tReco, "TMath::Log10(track.chi2_cr)", trackRecoOK, "",
          "tracks; log_{10}(#chi^{2}_{cr})");
    print(c, fout, tReco, "TMath::Log10(track.chi2_ln)", trackRecoOK, "",
          "tracks; log_{10}(#chi^{2}_{ln})");

    tReco->Draw("TMath::Log10(track.chi2_ln):TMath::Log10(track.chi2_cr)",
                trackRecoOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("tracks; log_{10}(#chi^{2}_{cr}); log_{10}(#chi^{2}_{ln})");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    n = tReco->Draw("cluster.y:cluster.x", "", "goff");
    g = new TGraph(n, tReco->GetV2(), tReco->GetV1());
    g->SetTitle("clusters; x (mm); y (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    n = tReco->Draw("cluster.y:cluster.z", "", "goff");
    g = new TGraph(n, tReco->GetV2(), tReco->GetV1());
    g->SetTitle("clusters; z (mm); y (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    print(c, fout, tReco, "TMath::Log10(cluster.t)", "", "",
          "clusters; log_{10}(t/ns)");
    print(c, fout, tReco, "TMath::Log10(cluster.e)", "", "",
          "clusters; log_{10}(E/MeV)");
    print(c, fout, tReco, "cluster.sx", "", "", "clusters; sx");
    print(c, fout, tReco, "cluster.sy", "", "", "clusters; sy");
    print(c, fout, tReco, "cluster.sz", "", "", "clusters; sz");

    tReco->Draw("cluster.sy:cluster.sx", "", "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("cluster; sx; sy");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    tReco->Draw("cluster.sy:cluster.sz", "", "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("cluster; sz; sy");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    tReco->Draw("cluster.sx:cluster.sz", "", "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("cluster; sz; sx");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    c.SetLogy(true);
    print(c, fout, tReco, "TMath::Log10(cluster.varx)", "cluster.varx>0", "",
          "clusters; log_{10}(varx)");
    print(c, fout, tReco, "TMath::Log10(cluster.vary)", "cluster.vary>0", "",
          "clusters; log_{10}(vary)");
    print(c, fout, tReco, "TMath::Log10(cluster.varz)", "cluster.varz>0", "",
          "clusters; log_{10}(varz)");
    c.SetLogy(false);

    tReco->Draw(
        "TMath::Log10(cluster.vary):TMath::Log10(cluster.varx)>>htemp(100,-20,"
        "10,100,-20,10)",
        "", "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("cluster;log_{10}(vary);log_{10}(varx)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    tReco->Draw(
        "TMath::Log10(cluster.vary):TMath::Log10(cluster.varz)>>htemp(100,-20,"
        "10,100,-20,10)",
        "", "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("cluster;log_{10}(vary);log_{10}(varz)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    tReco->Draw(
        "TMath::Log10(cluster.varx):TMath::Log10(cluster.varz)>>htemp(100,-20,"
        "10,100,-20,10)",
        "", "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("cluster;log_{10}(varx);log_{10}(varz)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());
  }

  std::cout << " done" << std::endl;

  std::cout << "tEvent histograms..." << std::flush;

  if (tEvent) {
    n = tEvent->Draw("y:x", "", "goff");
    g = new TGraph(n, tEvent->GetV2(), tEvent->GetV1());
    g->SetTitle("events; x (mm); y (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    n = tEvent->Draw("y:z", "", "goff");
    g = new TGraph(n, tEvent->GetV2(), tEvent->GetV1());
    g->SetTitle("events; z (mm); y (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());

    print(c, fout, tEvent, "t", "", "", "events; t (ns)");
    print(c, fout, tEvent, "Enu", "", "", "neutrino; E (MeV)");
    print(c, fout, tEvent, "pxnu", "", "", "neutrino; px (MeV)");
    print(c, fout, tEvent, "pynu", "", "", "neutrino; py (MeV)");
    print(c, fout, tEvent, "pznu", "", "", "neutrino; pz (MeV)");
    print(c, fout, tEvent, "Enureco>0>>htemp(2,-0.5,1.5)", "", "",
          "; reconstructed");

    c.SetLogy(true);
    print(c, fout, tEvent, "Enureco>>htemp(100,0,20000)", EnuRecoOK, "",
          "neutrino; Ereco (MeV)");
    print(c, fout, tEvent, "pxnureco>>htemp(100,-5000,5000)", EnuRecoOK, "",
          "neutrino; pxreco (MeV)");
    print(c, fout, tEvent, "pynureco>>htemp(100,-10000,5000)", EnuRecoOK, "",
          "neutrino; pyreco (MeV)");
    print(c, fout, tEvent, "pznureco>>htemp(100,-5000,20000)", EnuRecoOK, "",
          "neutrino; pzreco (MeV)");
    c.SetLogy(false);

    print(c, fout, tEvent, "particles.primary>>htemp(2,-0.5,1.5)", "", "",
          "particles; primary");

    tEvent->Draw("Enureco:Enu>>h2(100,0,10000,100,0,10000)", EnuRecoOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("h2");
    h2->SetStats(false);
    h2->SetTitle("neutrino; E (MeV); Ereco (MeV)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    c.SetLogy(true);
    tEvent->Draw("particles.pdg>>hpdg(8000,-4000,4000)", "", "");
    h1 = (TH1D*)gROOT->FindObject("hpdg");
    h1->SetStats(false);
    h1->SetTitle("particles; pdg");
    h1->Draw();
    c.SaveAs(fout.Data());

    tEvent->Draw("particles.pdg>>hpdg_pri(8000,-4000,4000)", PrimaryPart, "");
    h1 = (TH1D*)gROOT->FindObject("hpdg_pri");
    h1->SetStats(false);
    h1->SetTitle("primaries ; pdg");
    h1->Draw();
    c.SaveAs(fout.Data());
    c.SetLogy(false);

    print(c, fout, tEvent,
          "(particles.tr.ret_ln==0&&particles.tr.ret_cr==0)>>htemp(2,-0.5,1.5)",
          PrimaryPart, "", "primaries; recoOK==1");

    tEvent->Draw(
        "particles.charge_reco:particles.charge>>htemp(3,-1.5,1.5,3,-1.5,1.5)",
        PrimaryPart && PartTrackRecoOK, "colztext");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("reco primaries; true charge; reco charge");
    h2->Draw("colztext");
    c.SaveAs(fout.Data());

    tEvent->Draw(
        "particles.pxreco:particles.pxtrue>>htemp(200,-2000,2000,200,-2000,"
        "2000)",
        PrimaryPart && PartTrackRecoOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("reco primaries; true px (MeV); reco px (MeV)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    tEvent->Draw(
        "particles.pyreco:particles.pytrue>>htemp(200,-2000,2000,200,-2000,"
        "2000)",
        PrimaryPart && PartTrackRecoOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("reco primaries; true py (MeV); reco py (MeV)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    tEvent->Draw(
        "particles.pzreco:particles.pztrue>>htemp(200,-1000,15000,200,-1000,"
        "15000)",
        PrimaryPart && PartTrackRecoOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("reco primaries; true pz (MeV); reco pz (MeV)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    tEvent->Draw(
        "particles.Ereco:particles.Etrue>>htemp(200,0,15000,200,0,15000)",
        PrimaryPart && PartTrackRecoOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("reco primaries; true E (MeV); reco E (MeV)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    tEvent->Draw(
        "particles.xreco:particles.xtrue>>htemp(200,-2000,2000,200,-2000,2000)",
        PrimaryPart && PartTrackRecoOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true x (mm); reco x (mm)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    tEvent->Draw(
        "particles.yreco:particles.ytrue>>htemp(200,-4500,-200,200,-4500,-200)",
        PrimaryPart && PartTrackRecoOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true y (MeV); reco y (mm)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    tEvent->Draw(
        "particles.zreco:particles.ztrue>>htemp(200,22000,26000,200,22000,"
        "26000)",
        PrimaryPart && PartTrackRecoOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true z (mm); reco z (mm)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());

    tEvent->Draw("particles.treco:particles.ttrue",
                 PrimaryPart && PartTrackRecoOK, "colz");
    h2 = (TH2D*)gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle("primaries; true t (ns); reco t (ns)");
    h2->Draw("colz");
    c.SaveAs(fout.Data());
  }

  std::cout << " done" << std::endl;

  c.Clear();
  c.SaveAs(TString::Format("%s)", fout.Data()).Data());
  return 0;
}
