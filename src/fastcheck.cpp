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
  
  if(tDigit)
  {
    print(c, fout, tDigit, "cell.mod", "", "", "; module id", -1);
    print(c, fout, tDigit, "cell.lay", "", "", "; layer id");
    print(c, fout, tDigit, "cell.cel", "cell.mod < 24", "", "; cell id");
    print(c, fout, tDigit, "cell.l", "cell.mod < 24", "", "; cell length (mm)");
    print(c, fout, tDigit, "cell.cel", "cell.mod > 24", "", "; cell id");
    print(c, fout, tDigit, "cell.l", "cell.mod > 24", "", "; cell length (mm)");
    
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
    
    n = tDigit->Draw("cell.y:cell.z", "cell.mod==0", "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("module id == 0; z cell position (mm); y cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("TMath::ATan2(23910.000-cell.z,cell.y+2384.7300)/TMath::Pi()*180.:cell.mod", "cell.mod < 24", "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("; module id; angle from y axis (rad)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("cell.y:cell.lay", "cell.mod==0", "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("module id == 0; lay id; y cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("cell.z:cell.cel", "cell.mod==0", "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("module id == 0; cel id; z cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("cell.x:cell.z", "cell.mod==30", "goff"); 
    g = new TGraph(n,tDigit->GetV2(),tDigit->GetV1());
    g->SetTitle("module id == 30; z cell position (mm); y cell position (mm)");
    g->Draw("ap*");
    c.SaveAs(fout.Data());
    
    n = tDigit->Draw("cell.x:cell.mod", "cell.mod >= 24", "goff"); 
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
    n = tReco->Draw("track.ret_ln:track.ret_cr", "", "colz"); 
    h2 = (TH2D*) gROOT->FindObject("htemp");
    h2->SetStats(false);
    h2->SetTitle(";ret_cr;ret_ln");
    h2->Draw("colztext");
    c.SaveAs(fout.Data());
    
    c.SetLogy(true);
    print(c, fout, tReco, "TMath::Log10(track.r)", "track.ret_ln==0&&track.ret_cr==0", "", "; log_{10}(R/mm)");
    print(c, fout, tReco, "track.a", "track.ret_ln==0&&track.ret_cr==0", "", "; a (mm)");
    print(c, fout, tReco, "track.b", "track.ret_ln==0&&track.ret_cr==0", "", "; b");
    c.SetLogy(false);
    print(c, fout, tReco, "track.h", "track.ret_ln==0&&track.ret_cr==0", "", "; h");
    
    n = tReco->Draw("track.y0:track.x0", "", "goff"); 
    g = new TGraph(n,tReco->GetV2(),tReco->GetV1()); 
    g->SetTitle("; x0 (mm); y0 (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    n = tReco->Draw("track.y0:track.z0", "", "goff"); 
    g = new TGraph(n,tReco->GetV2(),tReco->GetV1()); 
    g->SetTitle("; z0 (mm); y0 (mm)");
    g->Draw("ap");
    c.SaveAs(fout.Data());
    
    print(c, fout, tReco, "track.t0", "track.ret_ln==0&&track.ret_cr==0", "", "; t (ns)");
    print(c, fout, tReco, "TMath::Log10(track.chi2_cr)", "track.ret_ln==0&&track.ret_cr==0", "", "; log_{10}(#chi^{2}_{cr})");
    print(c, fout, tReco, "TMath::Log10(track.chi2_ln)", "track.ret_ln==0&&track.ret_cr==0", "", "; log_{10}(#chi^{2}_{ln})");
    
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
    
    print(c, fout, tReco, "TMath::Log10(cluster.t)", "track.ret_ln==0&&track.ret_cr==0", "", "; log_{10}(t/ns)");
    print(c, fout, tReco, "TMath::Log10(cluster.e)", "track.ret_ln==0&&track.ret_cr==0", "", "; log_{10}(E/MeV)");
    print(c, fout, tReco, "cluster.sx", "track.ret_ln==0&&track.ret_cr==0", "", "; sx");
    print(c, fout, tReco, "cluster.sy", "track.ret_ln==0&&track.ret_cr==0", "", "; sy");
    print(c, fout, tReco, "cluster.sz", "track.ret_ln==0&&track.ret_cr==0", "", "; sz");
    print(c, fout, tReco, "TMath::Log10(cluster.varx)", "track.ret_ln==0&&track.ret_cr==0", "", "; log_{10}(varx)");
    print(c, fout, tReco, "TMath::Log10(cluster.vary)", "track.ret_ln==0&&track.ret_cr==0", "", "; log_{10}(vary)");
    print(c, fout, tReco, "cluster.varz", "track.ret_ln==0&&track.ret_cr==0", "", "; varz");
  }
  
  if(tEvent)
  {
  }
  
  c.Clear();
  c.SaveAs(TString::Format("%s)",fout.Data()).Data());
  return 0;
}