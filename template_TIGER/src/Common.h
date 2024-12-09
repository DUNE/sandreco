#ifndef Common_h
#define Common_h
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TMath.h"
#include "TROOT.h"
using namespace std;
// XYZ are refered to the detector frame where the anode is on the XY plane and the electrinc drift field is along the Z direction
namespace TIGER {

	//Simulation
	const bool   ToT = false;
	const bool   SH = true;
	const bool   NO_Noise = false;
	const bool   NO_Saturation = false;
	//Readout - Induction TIGER
	#define n_ns 1625//1500//875//675 //500
	const double jitter_TIGER = 0;
	const double timestep_TIGER = 6.25; //ns
	const double noise_TIGER = 10; //mV
	const double thrT_TIGER = 10; //mV
	const double thrE_TIGER = 20; //mV
	const double gain_TIGER = 12.7; //mV/fC
	const int    integration_time_TIGER = 8;
	const double saturation_TIGER = 800;
	const bool   enable_tfine_tiger = true;
}
#endif
