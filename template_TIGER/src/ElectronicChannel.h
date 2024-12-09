#ifndef ElectronicChannel_h
#define ElectronicChannel_h
#include "Common.h"
using namespace std; 

namespace TIGER{
  class ElectronicChannel{
  public:
    //Constructor
    //ElectronicChannel():above_thr(false){
    //    h_time_raw     = new TH1D();
    //    h_time_apv     = new TH1D();//
    //    h_time_tiger_T = new TH1D();
    //    h_time_tiger_E = new TH1D();
    //    h_time_cap     = new TH1D();//
    //    h_time_int     = new TH1D();//?
    //    h_time_res     = new TH1D();//
    //    h_time_tot     = new TH1D();
    //    h_time_cur     = new TH1D();//
    //};
    ElectronicChannel():                                                                                  
        above_thr(false){                        
        h_time_raw     = new TH1D("  raw","Current       raw",n_ns, 0,  n_ns);         //raw signal induced by electrons
        h_time_tiger_T = new TH1D("tig_T","Voltage[mV] tig_T",n_ns, 0,  n_ns);         //tiger Tbranch signal
        h_time_tiger_E = new TH1D("tig_E","Voltage[mV] tig_E",n_ns, 0,  n_ns);         //tiger Ebranch signal
        //h_time_int     = new TH1D("  int","Charge[fC]    int",n_ns, 0, n_ns);          //integrator circuit
        h_time_tot     = new TH1D("  tot","Charge[fC]    tot",n_ns, 0, n_ns);          //integrated charge over time
        h_time_cur     = new TH1D("  cur","Current       cur",n_ns, 0, n_ns);          //current after the resistive
      
        h_time_raw->GetXaxis()->SetTitle("Time [ns]");
        h_time_tiger_T->GetXaxis()->SetTitle("Time [ns]");
        h_time_tiger_E->GetXaxis()->SetTitle("Time [ns]");
        //h_time_int->GetXaxis()->SetTitle("Time [ns]");
        h_time_tot->GetXaxis()->SetTitle("Time [ns]");
        h_time_cur->GetXaxis()->SetTitle("Time [ns]");
      
        h_time_raw->GetYaxis()->SetTitle("Current [1 fC / 1 ns]");
        h_time_tiger_T->GetYaxis()->SetTitle("Voltage [1 mV / 1 ns]");
        h_time_tiger_E->GetYaxis()->SetTitle("");
        //h_time_int->GetYaxis()->SetTitle("Charge [fC]");
        h_time_tot->GetYaxis()->SetTitle("Charge [fC]");
        h_time_cur->GetYaxis()->SetTitle("Current [1 fC / 1 ns]");
      
    };
    //Destructor
    ~ElectronicChannel() {
      h_time_raw->~TH1D();
      h_time_tiger_T->~TH1D();
      h_time_tiger_E->~TH1D();
      //h_time_int->~TH1D();
      h_time_tot->~TH1D();
      h_time_cur->~TH1D();
      /*
      delete h_time_raw;
      delete h_time_apv;
      delete h_time_tiger_T;
      delete h_time_tiger_E;
      delete h_time_cap;
      delete h_time_int;
      delete h_time_res;
      delete h_time_tot;
      delete h_time_cur;
      */
    };
    //Function
    int       Get_ChannelID     ()          {return channelID;};
    double    Get_Charge        ()          {return charge;};
    double    Get_Time          ()          {return time;};
    double    Get_dTime         ()          {return dtime;};
    double    Get_t_thr_E       ()          {return t_thr_E;};
    double    Get_t_Q_E         ()          {return t_Q_E;};
    TH1D*     Get_Histo_raw     ()          {return h_time_raw;};
    TH1D*     Get_Histo_tiger_E ()          {return h_time_tiger_E;};
    TH1D*     Get_Histo_tiger_T ()          {return h_time_tiger_T;};
    //TH1D*     Get_Histo_int     ()          {return h_time_int;};
    TH1D*     Get_Histo_tot     ()          {return h_time_tot;};
    TH1D*     Get_Histo_cur     ()          {return h_time_cur;};
    bool      Get_AboveThr      ()          {return above_thr;};
    
    void      Set_Charge        (double io) {charge=io;};
    void      Set_Time          (double io) {time=io;};
    void      Set_dTime         (double io) {dtime=io;};
    void      Set_t_thr_E       (double io) {t_thr_E=io;};
    void      Set_t_Q_E         (double io) {t_Q_E=io;};
    void      Set_AboveThr      (bool   io) {above_thr=io;};
    
    void      Reset             () { Reset_Time(); Set_Charge(0); Set_Time(0); Set_dTime(0); };
    void      Reset_Time        () {
      h_time_raw->Reset();
      //h_time_int->Reset();
      h_time_cur->Reset();
      h_time_tiger_E->Reset();
      h_time_tiger_T->Reset();
    };
    void      Fill_Time(double io, int charge) { h_time_raw->Fill(io, -charge); }; // charge in fC

    void      Print_Time_TIGER    (int ch_id) {
      TCanvas *c = new TCanvas("ccc","ccc",600,600);
      c->Divide(2,3);
      c->cd(1); h_time_raw->Draw();
      c->cd(2); h_time_tot->Draw();
      c->cd(3); h_time_cur->Draw();
      //c->cd(8); h_time_int->Draw();
      c->cd(5);
      h_time_tiger_T->Draw();
      float V_thr_T = thrT_TIGER;
      float t_thr_T = Get_Time();
      TLine *l_t_thr_T = new TLine(t_thr_T,h_time_tiger_T->GetMinimum()-10,t_thr_T,V_thr_T);
      TLine *l_V_thr_T = new TLine(0,V_thr_T,t_thr_T,V_thr_T);
      l_V_thr_T->SetLineColor(kRed);
      l_t_thr_T->SetLineColor(kRed);
      l_V_thr_T->Draw("same");
      l_t_thr_T->Draw("same");
      c->cd(6);
      h_time_tiger_E->Draw();
      float V_thr_E = thrE_TIGER;
      float t_thr_E = Get_Time();
      float V_Q     = Get_Charge()*gain_TIGER;
      TLine *l_t_thr_E = new TLine(Get_t_thr_E(),Get_Histo_tiger_E()->GetMinimum()-10,Get_t_thr_E(),V_thr_E);
      TLine *l_V_thr_E = new TLine(0,V_thr_E,Get_t_thr_E(),V_thr_E);
      TLine *l_t_Q     = new TLine(Get_t_Q_E(),Get_Histo_tiger_E()->GetMinimum()-10,Get_t_Q_E(),V_Q);
      TLine *l_V_Q     = new TLine(0,V_Q,Get_t_Q_E(),V_Q);
      l_V_thr_E->SetLineColor(kRed);
      l_t_thr_E->SetLineColor(kRed);
      l_t_Q->SetLineColor(kRed);
      l_V_Q->SetLineColor(kRed);
      l_V_thr_E->Draw("same");
      l_t_thr_E->Draw("same");
      l_t_Q->Draw("same");
      l_V_Q->Draw("same");
      c->cd(4);
      TPaveText *text_tiger = new TPaveText(0.2,0.2,0.8,0.8);
      text_tiger->AddText(Form("Tiger T thr = %.0f mV",thrT_TIGER));
      text_tiger->AddText(Form("Tiger E thr = %.0f mV",thrE_TIGER));
      text_tiger->AddText(Form("Time measured = %.2f ns",Get_Time()));
      text_tiger->AddText(Form("Charge measured = %.2f fC",Get_Charge()));
      text_tiger->AddText(Form("Voltage at Qmeas = %.2f mV",Get_Charge()*gain_TIGER));
      text_tiger->AddText(Form("Time at Qmeas = %.2f ns",Get_t_Q_E()));
      text_tiger->AddText(Form("Time at thr_E = %.2f ns",Get_t_thr_E()));
      text_tiger->AddText(Form("Max charge collected = %.2f fC",Get_Histo_tot()->GetMaximum()));
      text_tiger->Draw();
      c->SaveAs(Form("data/Channel_histos_%i.pdf", ch_id));

      delete c;
    };

    //Variables
    TH1D*     h_time_raw;
    TH1D*     h_time_tiger_E;
    TH1D*     h_time_tiger_T;
    //TH1D*     h_time_int;
    TH1D*     h_time_tot;
    TH1D*     h_time_cur;
    //constant
    const int nbin = 1000;
  private:
    //Variable
    bool      PrintNTuple;
    int       channelID;
    double    charge; //fC
    double    time;   //ns
    double    dtime;  //ns
    bool      above_thr;
    double    t_thr_E;
    double    t_Q_E;
  };
}
#endif
