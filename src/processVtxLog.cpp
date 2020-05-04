void processVtxLog()
{
  gSystem->Exec("echo \"minreg/I:epsilon/D:wideMultReg/I:gefficiency/D:efficiency:rmsX:rmsY:rmsZ:rms3D\" > tmp");
  gSystem->Exec("cat findvtx.log | grep -v \"Events\" | grep -v \"minreg\" | grep -v \"sampl\" | grep -v \"44.4                   5\" >> tmp");

  gROOT->SetBatch(true);

  TTree t;
  t.ReadFile("tmp");
  
  int n = t.Draw("rmsX:epsilon","minreg==1&&wideMultReg==0","goff");
  TGraph gRMSX(n,t.GetV2(),t.GetV1());
  
  t.Draw("rmsY:epsilon","minreg==1&&wideMultReg==0","goff");
  TGraph gRMSY(n,t.GetV2(),t.GetV1());
  gRMSY.SetMarkerColor(kRed);
  gRMSY.SetLineColor(kRed);
  
  t.Draw("rmsZ:epsilon","minreg==1&&wideMultReg==0","goff");
  TGraph gRMSZ(n,t.GetV2(),t.GetV1());
  gRMSZ.SetMarkerColor(kBlue);
  gRMSZ.SetLineColor(kBlue);
  
  t.Draw("100+efficiency*300:epsilon","minreg==1&&wideMultReg==0","goff");
  TGraph gEff(n,t.GetV2(),t.GetV1());
  gEff.SetMarkerColor(kOrange);
  gEff.SetLineColor(kOrange);
  
  t.Draw("100+gefficiency*300:epsilon","minreg==1&&wideMultReg==0","goff");
  TGraph gGEff(n,t.GetV2(),t.GetV1());
  gGEff.SetMarkerColor(kGreen);
  gGEff.SetLineColor(kGreen);
  
  TCanvas c;
  c.DrawFrame(0,100,2.2,400,";#epsilon;RMS (mm)");
  gRMSX.Draw("*lsame");
  gRMSY.Draw("*lsame");
  gRMSZ.Draw("*lsame");
  gEff.Draw("*lsame");
  gGEff.Draw("*lsame");
  
  TLine l(0.5,100,0.5,400);
  l.SetLineStyle(2);
  l.Draw();
  
  TLegend leg(0.64,0.69,0.89,0.89);
  leg.AddEntry(&gRMSX,"RMSX","lp");
  leg.AddEntry(&gRMSY,"RMSY","lp");
  leg.AddEntry(&gRMSZ,"RMSZ","lp");
  leg.AddEntry(&gEff,"mtr vtx","lp");
  leg.AddEntry(&gGEff,"good mtr vtx","lp");
  leg.Draw();
  
  TGaxis ax(2.2,100,2.2,400,0,1,510,"+L");
  ax.SetTitle("efficiency");
  ax.SetTextColor(kOrange);
  ax.SetLabelColor(kOrange);
  ax.SetLineColor(kOrange);
  ax.Draw();
  
  c.SaveAs("optPar.pdf(");
  
  n = t.Draw("rmsX:minreg","epsilon==0.5&&wideMultReg==0","goff");
  TGraph gRMSX_mr(n,t.GetV2(),t.GetV1());
  
  t.Draw("rmsY:minreg","epsilon==0.5&&wideMultReg==0","goff");
  TGraph gRMSY_mr(n,t.GetV2(),t.GetV1());
  gRMSY_mr.SetMarkerColor(kRed);
  gRMSY_mr.SetLineColor(kRed);
  
  t.Draw("rmsZ:minreg","epsilon==0.5&&wideMultReg==0","goff");
  TGraph gRMSZ_mr(n,t.GetV2(),t.GetV1());
  gRMSZ_mr.SetMarkerColor(kBlue);
  gRMSZ_mr.SetLineColor(kBlue);
  
  t.Draw("50+efficiency*200:minreg","epsilon==0.5&&wideMultReg==0","goff");
  TGraph gEff_mr(n,t.GetV2(),t.GetV1());
  gEff_mr.SetMarkerColor(kOrange);
  gEff_mr.SetLineColor(kOrange);
  
  t.Draw("50+gefficiency*200:minreg","epsilon==0.5&&wideMultReg==0","goff");
  TGraph gGEff_mr(n,t.GetV2(),t.GetV1());
  gGEff_mr.SetMarkerColor(kGreen);
  gGEff_mr.SetLineColor(kGreen);
  
  t.Draw("50+gefficiency*efficiency*200:minreg","epsilon==0.5&&wideMultReg==0","goff");
  TGraph gGEffprod_mr(n,t.GetV2(),t.GetV1());
  gGEffprod_mr.SetLineStyle(2);
  
  c.DrawFrame(0.5,50,10.5,250,";min_{tr};RMS");
  gRMSX_mr.Draw("*lsame");
  gRMSY_mr.Draw("*lsame");
  gRMSZ_mr.Draw("*lsame");
  gEff_mr.Draw("*lsame");
  gGEff_mr.Draw("*lsame");
  gGEffprod_mr.Draw("*lsame");
  
  ax.SetX1(10.5);
  ax.SetX2(10.5);
  ax.SetY1(50);
  ax.SetY2(250);
  ax.Draw();
  leg.Draw();
  l.SetX1(2);
  l.SetX2(2);
  l.SetY1(50);
  l.SetY2(250);
  l.Draw();
    
  c.SaveAs("optPar.pdf");
  
  n = t.Draw("rmsX:minreg","epsilon==0.5&&wideMultReg==1","goff");
  TGraph gRMSX_mr_wdm(n,t.GetV2(),t.GetV1());
  
  t.Draw("rmsY:minreg","epsilon==0.5&&wideMultReg==1","goff");
  TGraph gRMSY_mr_wdm(n,t.GetV2(),t.GetV1());
  gRMSY_mr_wdm.SetMarkerColor(kRed);
  gRMSY_mr_wdm.SetLineColor(kRed);
  
  t.Draw("rmsZ:minreg","epsilon==0.5&&wideMultReg==1","goff");
  TGraph gRMSZ_mr_wdm(n,t.GetV2(),t.GetV1());
  gRMSZ_mr_wdm.SetMarkerColor(kBlue);
  gRMSZ_mr_wdm.SetLineColor(kBlue);
  
  t.Draw("50+efficiency*200:minreg","epsilon==0.5&&wideMultReg==1","goff");
  TGraph gEff_mr_wdm(n,t.GetV2(),t.GetV1());
  gEff_mr_wdm.SetMarkerColor(kOrange);
  gEff_mr_wdm.SetLineColor(kOrange);
  
  t.Draw("50+gefficiency*200:minreg","epsilon==0.5&&wideMultReg==1","goff");
  TGraph gGEff_mr_wdm(n,t.GetV2(),t.GetV1());
  gGEff_mr_wdm.SetMarkerColor(kGreen);
  gGEff_mr_wdm.SetLineColor(kGreen);
  
  t.Draw("50+gefficiency*efficiency*200:minreg","epsilon==0.5&&wideMultReg==1","goff");
  TGraph gGEffprod_mr_wdm(n,t.GetV2(),t.GetV1());
  gGEffprod_mr_wdm.SetLineStyle(2);
  
  c.DrawFrame(0.5,50,10.5,250,";min_{tr};RMS");
  gRMSX_mr_wdm.Draw("*lsame");
  gRMSY_mr_wdm.Draw("*lsame");
  gRMSZ_mr_wdm.Draw("*lsame");
  gEff_mr_wdm.Draw("*lsame");
  gGEff_mr_wdm.Draw("*lsame");
  gGEffprod_mr_wdm.Draw("*lsame");
  
  ax.Draw();
  leg.Draw();
    
  c.SaveAs("optPar.pdf)");
  
  gSystem->Exec("rm -f tmp");
  gROOT->SetBatch(false);
}
