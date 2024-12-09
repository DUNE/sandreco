#ifndef UTILS_H
#define UTILS_H

#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TMultiGraph.h"
//#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TROOT.h"
using namespace std;

//create a class ?

//I want to store the wf vector here, possibly in the same way it is done in sandReco
namespace WF_DATA {
    class DATA {
    public:
        //Construction
        DATA() {
            ReadWFfile();
        };
        //Destruction
        ~DATA() {};

        //Function 
        void SetTimeArr(int io) {
            time.clear();
            for (int t = 0; t < io; ++t) {
                time.push_back(t);
            }
        };

        void ReadWFfile() {
            TString wf_path = "/mnt/e/data/drift_chamber/wf_sample_data";
            TString fname = "/map_drift_0_.root";

            cout << "Opening " << wf_path + fname << endl;
            TFile* infile = new TFile(wf_path + fname);
            TTree* tree = (TTree*)infile->Get("out_tree");
            //save to class all_wf
            vector<double>* wf_tmpPtr = &wf_tmp;
            tree->SetBranchAddress("ind_wf", &wf_tmpPtr);

            nentries = tree->GetEntries();
            for (Int_t i = 0; i < nentries; i++) {
                tree->GetEntry(i);
                SetTimeArr(wf_tmp.size());
                time_arr.push_back(time);
                all_wf.push_back(wf_tmp);
            }
            cout << "entries: " << nentries << endl;
        };

        vector< vector<double> > Get_TimeWF(int io) {
            //read wf n°io from wf vector and returns it
            time_wf_arr.clear();
            time_wf_arr.push_back(time_arr[io]); //maybe i can create a wf data type
            time_wf_arr.push_back(all_wf[io]);
            return time_wf_arr;
        };

        vector<double> Get_WF(int io) {
            //read wf n°io from wf vector and returns it
            return all_wf[io];
        };

        int Get_Nentries() { return nentries; };
    private:
        int nentries = 0;
        vector<double>              wf_tmp;
        vector<double>              time;
        vector< vector<double> >    all_wf;
        vector< vector<double> >    time_arr;
        vector< vector<double> >    time_wf_arr;
        
    };//class
};//namespace
#endif
