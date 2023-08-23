// #include <TFile.h>
// #include <TTree.h>
#include <TApplication.h>
// #include <TCanvas.h>
// #include <TMarker.h>
// #include <TLine.h>
// #include <TG4Event.h>

// #include "STTCluster.h"
// #include "STTKFTrack.h"
// #include "STTUtils.h"
#include "SANDEventDisplay.h"

#include <iostream>

void usage();

using namespace std;

int main(int argc, char* argv[])
{
  int eventNumber = 0;
  TString simFileName = "";
  TString digitFileName = "";

  // read command line arguments

  for (int iParam = 1; iParam < argc; iParam++) {

    if (argv[iParam][0] == '-') {
      if (!strcmp(argv[iParam], "--sim")) {
        iParam++;
        if (iParam < argc && argv[iParam][0] != '-') {
          simFileName = argv[iParam];
          continue;
        } else {
          cerr << "simulated file name missing" << endl;
          exit(1);
        }
      } else if (!strcmp(argv[iParam], "--digit")) {
        iParam++;
        if (iParam < argc && argv[iParam][0] != '-') {
          digitFileName = argv[iParam];
          continue;
        } else {
          cerr << "digit file name missing" << endl;
          exit(1);
        }
      }

      switch (argv[iParam][1]) {
        case 'e':
          iParam++;
          if (iParam < argc && argv[iParam][0] != '-')
            eventNumber = ((TString)argv[iParam]).Atoll();
          else {
            cerr << "event number is missing" << endl;
            exit(0);
          }
          break;
        case 'h':
          usage();
          exit(0);
        default:
          cerr << "unknown command " << argv[iParam] << endl;
          usage();
          exit(0);
      }
    }
  }

  // TFile* fIn = 0;
  // TFile* fMc = 0;
  // int ientry = -1;

  // TTree* tIn = 0;
  // TTree* tMc = 0;
  // TGeoManager* geo = 0;

  // try {
  //   ientry = std::atoi(argv[1]);
  //   fIn = TFile::Open(argv[2]);
  //   fMc = TFile::Open(argv[3]);

  //   try {
  //     tIn = static_cast<TTree*>(fIn->Get("tKF"));
  //   }
  //   catch(...)
  //   {
  //     std::cout << "<ERROR> tKF Tree not found in " << fIn->GetName() <<
  //     std::endl; return 1;
  //   }

  //   try {
  //     geo = static_cast<TGeoManager*>(fMc->Get("EDepSimGeometry"));
  //     tMc = static_cast<TTree*>(fMc->Get("EDepSimEvents"));
  //   }
  //   catch(...)
  //   {
  //     std::cout << "<ERROR> EDepSimGeometry TGeoManager not found in " <<
  //     fMc->GetName() << std::endl; return 1;
  //   }
  // }
  // catch(...)
  // {
  //   std::cout << "<ERROR> Usage: " << argv[0] << " <entry> <input file> <edep
  //   file>" << std::endl; return 1;
  // }

  TApplication application("SAND event viewer", &argc, argv);
  // SANDEventDisplay::Init(geo);

  // // auto c = SANDEventDisplay::Create2DDisplaySideBySide();
  // // c->SetTitle(TString::Format("Entry: %d",ientry).Data());

  // std::vector<STTKFTrack>* fKFTracks = 0;
  // std::vector<STTCluster>* fCluster = 0;

  // TG4Event *fEvent = 0;

  // // tIn->SetBranchAddress("tracks",&fKFTracks);
  // // tIn->SetBranchAddress("clusters",&fCluster);

  // tMc->SetBranchAddress("Event", &fEvent);

  // // SANDEventDisplay eventDisplay(fEvent, gClient->GetRoot(), 1600, 1200);
  // tMc->GetEntry(ientry);

  SANDEventDisplay eventDisplay(gClient->GetRoot(), 1600, 800);
  //   SANDEventDisplay eventDisplay(gClient->GetRoot(),
  //   gClient->GetDisplayWidth(), gClient->GetDisplayHeight());
  eventDisplay.SetEventNumber(eventNumber);
  if (simFileName != "") eventDisplay.SetSimData(simFileName);
  if (digitFileName != "") eventDisplay.SetDigitData(digitFileName);

  eventDisplay.Run();

  // cout << gClient->GetDisplayWidth() << " " << gClient->GetDisplayHeight() <<
  // endl;

  // tIn->GetEntry(ientry);

  // eventDisplay.DrawEvent(ientry);

  // EColor color[] = {EColor::kBlack, EColor::kRed, EColor::kBlue};
  // auto index = 0;
  //
  // cout << fKFTracks->size() << endl;
  //
  // for(auto& tr: *fKFTracks)
  // {
  //   if(tr.GetSteps().size() < 10)
  //     continue;
  //
  //   auto trackColor = color[index++ % 3];
  //
  //   for(unsigned int istep = 0; istep < tr.GetSteps().size()-1; ++istep)
  //   {
  //     auto clusterID = tr.GetStep(istep).GetClusterIDForThisState();
  //     auto thisCluster = STTUtils::GetClusterPointer(clusterID, *fCluster);
  //
  //     if(thisCluster != 0)
  //     {
  //       // cout << index << " " << trackColor << endl;
  //
  //       auto& cluster = *thisCluster;
  //       auto m = SANDEventDisplay::GetMarkerFromCluster(cluster, trackColor);
  //
  //       c->cd(cluster.GetOrientation() == STTPlane::EOrientation::kHorizontal
  //       ? 2 : 1); m->Draw();
  //     }
  //     else
  //     {
  //       std::cout << "<ERROR>: cluster with ID: " << clusterID << " not
  //       found" << std::endl;
  //     }
  //
  //     auto lines = SANDEventDisplay::GetFilterToPredictionLines(tr, istep,
  //     trackColor);
  //
  //     c->cd(1)->cd();
  //     lines[0]->Draw();
  //     c->cd(2)->cd();
  //     lines[1]->Draw();
  //   }

  // }

  //  c->Update();

  application.Run();
  return 0;
}

//----------------------------------------------------------------------------
void usage()
{
  cout << "SAND EVENT VIEWER PACKAGE" << endl;
  cout << "usage:\n";
  cout << "   eventDisplay [-h] [-e <eventNumber>] [--sim <fileName> | --digit "
          "<fileName>]\n";
  cout << "command are:\n";
  cout << "  -e     <eventNumber> : Event number\n";
  cout << " --sim   <fileName>    : file with simulated information\n";
  cout << " --digit <fileName>    : file with digits\n";
  cout << "  -h                   : Print the help" << endl;
}
