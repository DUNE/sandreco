#include "SANDEventDisplay.h"
#include "SANDDisplayUtils.h"
// #include "STTStrawTubeTracker.h"
// #include "STTUtils.h"
#include "utils.h"
#include <iostream>

#include <TG4Event.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TPolyMarker.h>
#include <TClonesArray.h>
#include <TPolyLine.h>
#include <TDatabasePDG.h>
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGToolBar.h>
#include <TGMenu.h>
#include <TGDockableFrame.h>
#include <TGNumberEntry.h>
#include <TGButtonGroup.h>
#include <TGLabel.h>
#include <TColor.h>
#include <TMath.h>

using namespace std;

ClassImp(SANDEventDisplay)

enum ETestCommandIdentifiers {
   M_FILE_OPEN,
   M_SAVE_PICTURE,
   M_SAVE_PICTURE_PNG,
   M_SAVE_PICTURE_PDF,
   M_PREVIOUS_EVENT,
   M_EXECUTE_EVENT,
   M_NEXT_EVENT,
   M_ZOOM_IN,
   M_ZOOM_OUT,
   M_ZOOM_EVENT,
   M_ZOOM_VERTEX,
   M_ZOOM_DETECTOR,
   M_FILE_EXIT
};

const char *EVXpmNames[] = {
  "open.xpm",
  "save-picture.xpm",
  "",
  "previous.xpm",
  "reload.xpm",
  "next.xpm",
  "",
  "zoom-in.xpm",
  "zoom-out.xpm",
  "zoom-event.xpm",
  "zoom-vertex.xpm",
  "zoom-detector.xpm",
  "",
  "exit.xpm",
  0
};

ToolBarData_t EVTbData[] = {
  { "", "Open Root event file",     kFALSE, M_FILE_OPEN,         NULL },
  { "", "Save picture",             kFALSE, M_SAVE_PICTURE,      NULL },
  { "",              0,             0,      -1,                  NULL },
  { "", "Previous event",           kFALSE, M_PREVIOUS_EVENT,    NULL },
  { "", "Execute event",            kFALSE, M_EXECUTE_EVENT,     NULL },
  { "", "Next event",               kFALSE, M_NEXT_EVENT,        NULL },
  { "",              0,             0,      -1,                  NULL },
  { "", "Zoom In",                  kFALSE, M_ZOOM_IN,           NULL },
  { "", "Zoom Out",                 kFALSE, M_ZOOM_OUT,          NULL },
  { "", "Zoom to event",            kFALSE, M_ZOOM_EVENT,        NULL },
  { "", "Zoom to vertex",           kFALSE, M_ZOOM_VERTEX,       NULL },
  { "", "Zoom to detector",         kFALSE, M_ZOOM_DETECTOR,     NULL },
  { "",              0,             0,      -1,                  NULL },
  { "", "Exit Application",         kFALSE, M_FILE_EXIT,         NULL },
  { NULL,            NULL,          0,      0,                   NULL }
};



SANDEventDisplay::SANDEventDisplay(const TGWindow *p, int w, int h) : TGMainFrame(p, w, h)
{
  // initialization

  fColNum = 20;
  fSetHitsMode = 2;
  fGeomInitialized = false;
  fDrawHits = kDigitHits;

  fEvent = NULL;
  fEventNumber = 0;
  fFileSimData = NULL;
  fFileDigitData = NULL;
  fTreeSimData = NULL;
  fTreeDigitData = NULL;
  fTubeDigitVect = NULL;
  fCellDigitVect = NULL;
  fPDGcode = TDatabasePDG::Instance();

  fPadZY = NULL;
  fPadZX = NULL;

  // fDisplayCanvas = new TCanvas("SAND viewer", "SAND viewer", 2000, 1000);

  // auto sandCenter = STTUtils::GetSANDInnerVolumeCenterPosition();
  // auto sandRadius = STTUtils::GetSANDInnerVolumeRadius()/10;
  // auto sandLenght = STTUtils::GetSANDInnerVolumeLength()/10;
  //
  // for (int i = 0; i < 3; ++i) sandCenter[i] /= 10;

  // fHistoHitsZYviewTrue->SetStats(false);
  // fHistoHitsZXviewTrue->SetStats(false);
  //
  // fHistoHitsZYviewTrue->GetXaxis()->SetTitle("z position, cm");
  // fHistoHitsZYviewTrue->GetYaxis()->SetTitle("y position, cm");
  // fHistoHitsZXviewTrue->GetXaxis()->SetTitle("z position, cm");
  // fHistoHitsZXviewTrue->GetYaxis()->SetTitle("x position, cm");
  //
  // fHistoHitsZYviewTrue->SetTitleOffset(1.5,"Y");
  // fHistoHitsZXviewTrue->SetTitleOffset(1.5,"Y");

  fTracksArrayZYTrue = new TClonesArray("TPolyLine");
  fTracksArrayZXTrue = new TClonesArray("TPolyLine");
  fVerticesZYTrue = new TPolyMarker();
  fVerticesZXTrue = new TPolyMarker();

  fVerticesZYTrue->SetMarkerStyle(29);
  fVerticesZYTrue->SetMarkerColor(6);
  fVerticesZYTrue->SetMarkerSize(2);

  fVerticesZXTrue->SetMarkerStyle(29);
  fVerticesZXTrue->SetMarkerColor(6);
  fVerticesZXTrue->SetMarkerSize(2);

  DefineColors();
  DrawButtons();

  SetWindowName("SAND Event Viewer");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
}


SANDEventDisplay::~SANDEventDisplay()
{
  SafeDelete(fTracksArrayZYTrue);
  SafeDelete(fTracksArrayZXTrue);
  SafeDelete(fVerticesZYTrue);
  SafeDelete(fVerticesZXTrue);
  SafeDelete(fFileSimData);
  SafeDelete(fFileDigitData);
}


//---------------------------------------------------------------------------
void SANDEventDisplay::Run()
{
  if (fTreeSimData) fTreeSimData->GetEntry(fEventNumber);
  if (fTreeDigitData) fTreeDigitData->GetEntry(fEventNumber);
  DrawEvent();
}

//----------------------------------------------------------------------------
void SANDEventDisplay::DrawEvent()
{
  InitObjects(); // clear all drawing objects
  DrawDetector();
  FillEventHits();
  FillDigitHits();
  DrawTracks();
  fEntryEventId->SetIntNumber(fEventNumber);
  fDisplayCanvas->cd(); fDisplayCanvas->Update();
}


//----------------------------------------------------------------------------
void SANDEventDisplay::FillEventHits()
{
  // event primary vertices

  for (auto &vtx : fEvent->Primaries) {
    fVerticesZYTrue->SetNextPoint(vtx.GetPosition().Z()/10, vtx.GetPosition().Y()/10);
    fVerticesZXTrue->SetNextPoint(vtx.GetPosition().Z()/10, vtx.GetPosition().X()/10);
  }

  // event hits

  EVHits_t hit;

  for (auto &segment : fEvent->SegmentDetectors) {
    for (auto &seg : segment.second) {

      if (seg.GetEnergyDeposit() > 0.00025) {

        // searching particle for this hit

        int pdgCode = 0;
        for (auto &contrib : seg.Contrib) {
          for (auto &trj : fEvent->Trajectories) {
            if (trj.GetTrackId() == contrib) {pdgCode = trj.GetPDGCode(); break;}
          }
          if (pdgCode) break;
        }

        hit.particle = pdgCode;
        hit.z = seg.GetStart().Z()/10;
        hit.e = seg.GetEnergyDeposit();
        hit.x = seg.GetStart().Y()/10;
        fEventHitsZY.push_back(hit);
        hit.x = seg.GetStart().X()/10;
        fEventHitsZX.push_back(hit);
      }
    }
  }

  // event trajectories

  for (auto &trj : fEvent->Trajectories) {
    int itrack = fTracksArrayZYTrue->GetEntries();

    TPolyLine *trackZY = (TPolyLine*)fTracksArrayZYTrue->ConstructedAt(itrack);
    TPolyLine *trackZX = (TPolyLine*)fTracksArrayZXTrue->ConstructedAt(itrack);

    trackZY->SetPolyLine(0);
    trackZX->SetPolyLine(0);
    int pdgCode = trj.GetPDGCode();

    if (fPDGcode->GetParticle(pdgCode)) {

      double charge = fPDGcode->GetParticle(pdgCode)->Charge();

      if (abs(pdgCode) == 11) {
        trackZY->SetLineColor(2); trackZX->SetLineColor(2);
      }
      else if (abs(pdgCode) == 13) {
        trackZY->SetLineColor(4); trackZX->SetLineColor(4);
      }
      else {
        trackZY->SetLineColor(3); trackZX->SetLineColor(3);
      }

      if (!charge) {
        trackZY->SetLineStyle(2); trackZX->SetLineStyle(2);
        trackZY->SetLineColor(9); trackZX->SetLineColor(9);
        if (pdgCode == 2112) {
          trackZY->SetLineColor(1); trackZX->SetLineColor(1);
        }
      }
      else {trackZY->SetLineStyle(1); trackZX->SetLineStyle(1);}

      if (charge) {
        for (auto &trk : trj.Points) {
          trackZY->SetNextPoint(trk.GetPosition().Z()/10, trk.GetPosition().Y()/10);
          trackZX->SetNextPoint(trk.GetPosition().Z()/10, trk.GetPosition().X()/10);
        }
      }
    }

    if (!trackZY->Size()) fTracksArrayZYTrue->Remove(trackZY);
    if (!trackZX->Size()) fTracksArrayZXTrue->Remove(trackZX);
  }
}


//----------------------------------------------------------------------------
void SANDEventDisplay::FillDigitHits()
{
  EVHits_t hit;

  // STT hits

  for (auto &tube : *fTubeDigitVect) {
    hit.particle = 0;
    hit.z = tube.z/10;
    hit.e = tube.de;
    if (tube.hor) {
      hit.x = tube.y/10;
      fTubeDigitHitsZY.push_back(hit);
    }
    else {
      hit.x = tube.x/10;
      fTubeDigitHitsZX.push_back(hit);
    }
  }

  // ECAL hits

  for (auto &cell : *fCellDigitVect) {
    hit.particle = 0;
    hit.z = cell.z/10;

    int pe = 0;

    for (auto &ps : cell.ps1) {
      for (auto &pe1 : ps.photo_el) ++pe;
    }

    for (auto &ps : cell.ps2) {
      for (auto &pe1 : ps.photo_el) ++pe;
    }

    hit.e = pe;

    if (cell.det == 2) {
      hit.x = cell.y/10;
      fCellDigitHitsZY.push_back(hit);
    }
    else if (cell.det == 1) {
      hit.x = cell.x/10;
      fCellDigitHitsZX.push_back(hit);
    }
  }

}


//----------------------------------------------------------------------------
void SANDEventDisplay::DrawTracks()
{
  // fHistoHitsZYviewTrue->SetTitle(Form("%i", entry));

  TEllipse *ellipse;
  TBox *box;

  // define track colors

  double de = (0.01-0.000250)/20.;

  fPadZY->cd();

  if (fDrawHits == kSimHits) {
    for (auto &hit : fEventHitsZY) {
      int nc = int((hit.e-0.000250)/de);
      hit.color = hit.e > 0.01 ? fPalette[fColNum-1] : fPalette[nc];
      ellipse = new TEllipse(hit.z, hit.x, 1, 1, 0, 360, 0);
      SANDDisplayUtils::DrawEllipse(ellipse, hit.color, 1001);
    }
  }
  else if (fDrawHits == kDigitHits) {
    for (auto &hit : fTubeDigitHitsZY){
      ellipse = new TEllipse(hit.z, hit.x, 0.5, 0.5, 0, 360, 0);
      SANDDisplayUtils::DrawEllipse(ellipse, 1, 1001);
    }
    for (auto &hit : fCellDigitHitsZY){
      ellipse = new TEllipse(hit.z, hit.x, 2.5, 2.5, 0, 360, 0);
      SANDDisplayUtils::DrawEllipse(ellipse, 1, 1001);
    }
  }

  fPadZX->cd();

  if (fDrawHits == kSimHits) {
    for (auto &hit : fEventHitsZX) {
      int nc = int((hit.e-0.000250)/de);
      hit.color = hit.e > 0.01 ? fPalette[fColNum-1] : fPalette[nc];
      ellipse = new TEllipse(hit.z, hit.x, 1, 1, 0, 360, 0);
      SANDDisplayUtils::DrawEllipse(ellipse, hit.color, 1001);
    }
  }
  else if (fDrawHits == kDigitHits) {
    for (auto &hit : fTubeDigitHitsZX){
      ellipse = new TEllipse(hit.z, hit.x, 0.5, 0.5, 0, 360, 0);
      SANDDisplayUtils::DrawEllipse(ellipse, 1, 1001);
    }
    for (auto &hit : fCellDigitHitsZX){
      ellipse = new TEllipse(hit.z, hit.x, 2.5, 2.5, 0, 360, 0);
      SANDDisplayUtils::DrawEllipse(ellipse, 1, 1001);
    }
  }

  fPadZY->cd();

  for (int i = 0; i < fTracksArrayZYTrue->GetEntries(); ++i)
  ((TPolyLine*)fTracksArrayZYTrue->At(i))->Draw();

  fVerticesZYTrue->Draw();

  fPadZX->cd();

  for (int i = 0; i < fTracksArrayZXTrue->GetEntries(); ++i)
  ((TPolyLine*)fTracksArrayZXTrue->At(i))->Draw();

  fVerticesZXTrue->Draw();

  //  canvas2->Print(Form("pictures/detector_view/detector_view_stt_%d.pdf", entry));
}


//----------------------------------------------------------------------------
void SANDEventDisplay::InitObjects()
{

  fVerticesZYTrue->SetPolyMarker(-1);
  fVerticesZXTrue->SetPolyMarker(-1);
  fTracksArrayZYTrue->Clear("C");
  fTracksArrayZXTrue->Clear("C");
  fEventHitsZY.clear();
  fEventHitsZX.clear();
  fTubeDigitHitsZY.clear();
  fTubeDigitHitsZX.clear();
  fCellDigitHitsZY.clear();
  fCellDigitHitsZX.clear();

  SafeDelete(fPadZY); SafeDelete(fPadZX);
  fDisplayCanvas->Clear(); fDisplayCanvas->cd();

  // SafeDelete(fPadZY); SafeDelete(fPadZX);

  fPadZY = new TPad("ZY plane", "ZY plane", 0., 0, 0.5, 1);
  fPadZX = new TPad("ZX plane", "ZX plane", 0.5, 0, 1, 1);

  fPadZY->Draw(); fPadZX->Draw();
}


//----------------------------------------------------------------------------
void SANDEventDisplay::DrawDetector()
{
  // set all sizes in cm

  // auto sandCenter = STTUtils::GetSANDInnerVolumeCenterPosition();
  // auto sandRadius = STTUtils::GetSANDInnerVolumeRadius()/10;
  // auto sandLenght = STTUtils::GetSANDInnerVolumeLength()/10;

  double sandCenter[3] = {0, -2384.73/10, 23910/10};
  double sandRadius = 200;
  double sandLenght = 338;

  // for (int i = 0; i < 3; ++i) sandCenter[i] /= 10;

  fPadZY->cd();
  fPadZY->DrawFrame(sandCenter[2] - 1.1 * sandRadius, sandCenter[1] - 1.1 * sandRadius, sandCenter[2] + 1.1 * sandRadius, sandCenter[1] + 1.1 * sandRadius);

  TEllipse *el1 = new TEllipse(sandCenter[2], sandCenter[1], sandRadius);
  el1->SetLineColor(2);
  el1->Draw();

  fPadZX->cd();
  fPadZX->DrawFrame(sandCenter[2] - 1.1 * sandRadius, sandCenter[0] - 1.1 * 0.5 * sandLenght, sandCenter[2] + 1.1 * sandRadius, sandCenter[0] + 1.1 * 0.5 * sandLenght);

  TBox *box1 = new TBox(sandCenter[2] - sandRadius, sandCenter[0] - 0.5*sandLenght, sandCenter[2] + sandRadius, sandCenter[0] + 0.5*sandLenght);
  box1->SetFillStyle(0);
  box1->SetLineColor(2);
  box1->Draw();

  // TString name = "Plane ZY";
  // TH2F *histo = (TH2F*)gROOT->FindObject(name);
  // SafeDelete(histo);
  //
  // histo = new TH2F(name, "", 100, sandCenter[2] - 1.1*sandRadius, sandCenter[0] - 1.1*0.5*sandLenght,
  // 100, sandCenter[2] + 1.1*sandRadius, sandCenter[0] + 1.1*0.5*sandLenght);

  // histo = new TH2F(name, "", 100, sandCenter[2] - 1.1*sandRadius, sandCenter[0] - 1.1*0.5*sandLenght,
  // 100, sandCenter[2] + 1.1*sandRadius, sandCenter[0] + 1.1*0.5*sandLenght);
  // histo->GetXaxis()->SetTitle("Z (cm)");
  // histo->GetYaxis()->SetTitle("Y (cm)");

  // SetHistoStyle(histo);
  //
  // histo->Draw();

  // fPadZX->cd();
  //
  // name = "Plane ZX";
  // histo = (TH2F*)gROOT->FindObject(name);
  // SafeDelete(histo);
  //
  // histo = new TH2F(name, "", 100, sandCenter[2] - 1.1*sandRadius, sandCenter[1] - 1.1*sandRadius,
  // 100, sandCenter[2] + 1.1*sandRadius, sandCenter[1] + 1.1*sandRadius);
  // histo->GetXaxis()->SetTitle("Z (cm)");
  // histo->GetYaxis()->SetTitle("X (cm)");
  //
  // SetHistoStyle(histo);
  //
  // histo->Draw();



  // auto c = new TCanvas("myDisplay","", 2000, 1000);
  // c->Divide(2,1);
  // c->cd(1)->DrawFrame(sandCenter[2]/10 - 1.1 * sandRadius, sandCenter[0]/10 - 1.1 * 0.5 * sandLenght, sandCenter[2]/10 + 1.1 * sandRadius, sandCenter[0]/10 + 1.1 * 0.5 * sandLenght);
  //
  // TBox *box1 = new TBox(sandCenter[2]/10 - sandRadius, sandCenter[0]/10 - 0.5*sandLenght, sandCenter[2]/10 + sandRadius, sandCenter[0]/10 + 0.5*sandLenght);
  // box1->SetFillStyle(0);
  // box1->SetLineColor(2);
  // box1->Draw();
  //
  // c->cd(2)->DrawFrame(sandCenter[2]/10 - 1.1 * sandRadius, sandCenter[1]/10 - 1.1 * sandRadius, sandCenter[2]/10 + 1.1 * sandRadius, sandCenter[1]/10 + 1.1 * sandRadius);
  //
  // TEllipse *el1 = new TEllipse(sandCenter[2]/10, sandCenter[1]/10, sandRadius);
  // el1->SetLineColor(2);
  // el1->Draw();
  //
  // return c;
}

// TMarker* SANDEventDisplay::GetMarkerFromCluster(const STTCluster& cluster, EColor color){
//
//   auto trk = cluster.GetRecoParameters().front().trk;
//   auto z = cluster.GetZ();
//   auto x = trk.m * z + trk.q;
//   auto m = new TMarker(z,x,5);
//   m->SetMarkerColor(color);
//
//   return m;
// }
//
// TLine** SANDEventDisplay::GetFilterToPredictionLines(const STTKFTrack& tr, int stepIndex, EColor color)
// {
//   auto step1 = tr.GetStep(stepIndex);
//   auto step2 = tr.GetStep(stepIndex+1);
//   auto state1 = step1.GetStage(STTKFTrackStep::STTKFTrackStateStage::kFiltering).GetStateVector();
//   auto state2 = step2.GetStage(STTKFTrackStep::STTKFTrackStateStage::kPrediction).GetStateVector();
//   auto x1 = state1.X();
//   auto y1 = state1.Y();
//   auto z1 = STTStrawTubeTracker::GetZPlane(step1.GetPlaneID());
//   auto x2 = state2.X();
//   auto y2 = state2.Y();
//   auto z2 = STTStrawTubeTracker::GetZPlane(step2.GetPlaneID());
//
//   auto lines = new TLine*[2];
//   lines[0] = new TLine(z1, x1, z2, x2);
//   lines[1] = new TLine(z1, y1, z2, y2);
//
//   lines[0]->SetLineColor(color);
//   lines[1]->SetLineColor(color);
//
//   return lines;
// }


//----------------------------------------------------------------------------
void SANDEventDisplay::DrawButtons()
{
  // layout hints

  fLayoutExpX      = new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2);
  fLayoutExpY      = new TGLayoutHints(kLHintsExpandY, 2, 2, 2, 2);
  fLayoutLeftExpY  = new TGLayoutHints(kLHintsLeft|kLHintsExpandY, 2, 4, 0, 0);
  fLayoutRightExpY = new TGLayoutHints(kLHintsRight|kLHintsExpandY,4, 2, 0, 0);
  fLayoutLeftExpX  = new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 2, 4, 0, 0);
  fLayoutRightExpX = new TGLayoutHints(kLHintsRight|kLHintsExpandX,4, 2, 0, 0);
  fLayoutExpXExpY  = new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,2,2,2,2);

  SetCleanup(kDeepCleanup);

  // Correct close window

  Connect("CloseWindow()", "SANDEventDisplay", this, "ExitApplication()");

  // Create main frame

  // fMainFrame = new TGVerticalFrame(this, gClient->GetDisplayWidth(), gClient->GetDisplayHeight());
  fMainFrame = new TGVerticalFrame(this, 1600, 800);
  fMainFrame->SetCleanup(kDeepCleanup);

  // Create a menu bar

  AddMenuBar(fMainFrame);

  TGHorizontalFrame *workframe = new TGHorizontalFrame(fMainFrame);

  // fTab = new TGTab(workframe);

  // Create a Tab with run event buttons

  // AddRunEventFrame(fTab->AddTab("Run Events"));
  AddRunEventFrame(workframe);

  // Create a Tab for drawing options

  // AddOptionsFrame(fTab->AddTab("Options"));

  // workframe->AddFrame(fTab, fLayoutExpY);

  AddCanvasFrame(workframe);               // draw canvases

  // AddFrame(workframe, fLayoutLeftExpY);

  // Status bar

  // Int_t parts[] = {20, 20, 20, 20, 20};
  // fStatusBar = new TGStatusBar(fMainFrame);
  // fStatusBar->SetParts(parts, 5);
  // fStatusBar->Draw3DCorner(kFALSE);
  //
  fMainFrame->AddFrame(workframe, fLayoutExpXExpY);
  // fMainFrame->AddFrame(fStatusBar, new TGLayoutHints(kLHintsExpandX,0,0,1,0));

  AddFrame(fMainFrame, fLayoutExpXExpY);

}


//----------------------------------------------------------------------------
void SANDEventDisplay::AddMenuBar(TGVerticalFrame *workframe)
{
  // create menu bar

  TGLayoutHints *LayoutMenuBar = new TGLayoutHints(kLHintsTop|kLHintsExpandX);
  TGLayoutHints *LayoutMenuBarItem = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);

  // Section File

  TGPopupMenu *MenuFile = new TGPopupMenu(gClient->GetRoot());
  // MenuFile->AddEntry("Open...",           M_FILE_OPEN);
  MenuFile->AddEntry("Save picture .png", M_SAVE_PICTURE_PNG);
  MenuFile->AddEntry("Save picture .pdf", M_SAVE_PICTURE_PDF);
  MenuFile->AddSeparator();
  MenuFile->AddEntry("Exit",              M_FILE_EXIT);

  MenuFile->DisableEntry(M_SAVE_PICTURE_PNG);
  MenuFile->DisableEntry(M_SAVE_PICTURE_PDF);
  MenuFile->DisableEntry(M_FILE_EXIT);

  // MenuFile->Associate(this);

  // menu bar

  TGDockableFrame *MenuDock = new TGDockableFrame(workframe);
  MenuDock->SetWindowName("Event Viewer Menu");

  TGMenuBar *MenuBar = new TGMenuBar(MenuDock, 1, 1, kHorizontalFrame);
  MenuBar->AddPopup("&File", MenuFile, LayoutMenuBarItem);

  MenuDock->AddFrame(MenuBar, LayoutMenuBar);
  workframe->AddFrame(MenuDock, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 1, 0));

  // Toolbar

  Int_t spacing = 8;
  TGToolBar *toolBar = new TGToolBar(workframe, 60, 20,
				     kHorizontalFrame | kRaisedFrame);
  for (Int_t i = 0; EVXpmNames[i]; i++) {
    TString iconname = (TString)getenv("SandReco_DIR") +
      "/icons/" + EVXpmNames[i];
    // TString iconname = TString("icons/") + TString(EVXpmNames[i]);
    EVTbData[i].fPixmap = iconname.Data();
    if (strlen(EVXpmNames[i]) == 0) {spacing = 8; continue;}
    toolBar->AddButton(this, &EVTbData[i], spacing);
    spacing = 0;
  }

  workframe->AddFrame(toolBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));
}


//----------------------------------------------------------------------------
void SANDEventDisplay::AddRunEventFrame(TGHorizontalFrame *workframe)
{
  // Create a vertical frame widget with buttons

  TGVerticalFrame *buttonframe = new TGVerticalFrame(workframe);

  // Filter event button

  // TGGroupFrame *GroupShowHits = new TGGroupFrame(buttonframe,"Show hits");
  // TGComboBox *fComboShowHits = new TGComboBox(GroupShowHits, 100);
  // fComboShowHits->Resize(50, 20);
  // fComboShowHits->Associate(this);
  // GroupShowHits->AddFrame(fComboShowHits, fLayoutExpX);
  // buttonframe->AddFrame(GroupShowHits, fLayoutExpX);

  // Tracks visualization

  // TGVButtonGroup *GroupTrkDir = new TGVButtonGroup(buttonframe,
	// 					   "Draw directions of");
  // fListTrkDir = new TGListBox(GroupTrkDir, 125);
  // fListTrkDir->SetMultipleSelections();
  // fListTrkDir->Resize(50, 70);
  // fListTrkDir->Associate(this);
  // GroupTrkDir->AddFrame(fListTrkDir, fLayoutExpX);
  //
  // fCheckTrkDir = new TGCheckButton(GroupTrkDir, "show OpRec tracking", 123);
  // fCheckTrkDir->Associate(this);
  //
  // if (fSetDrawMuonSeg) fCheckTrkDir->SetState(kButtonDown, kFALSE);
  //
  // buttonframe->AddFrame(GroupTrkDir, fLayoutExpX);

  TGHorizontalFrame *FrameEventId = new TGHorizontalFrame(buttonframe);
  TGLabel *LblEventId = new TGLabel(FrameEventId, "Event Number");
  fEntryEventId       = new TGNumberEntry(FrameEventId, fEventNumber, 7, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMin, 0);
  FrameEventId->AddFrame(LblEventId,    fLayoutLeftExpY);
  FrameEventId->AddFrame(fEntryEventId, fLayoutRightExpY);
  buttonframe->AddFrame(FrameEventId, fLayoutExpX);
  fEntryEventId->GetNumberEntry()->Connect("TextChanged(char*)", "SANDEventDisplay", this, "SetEventId()");

  // Energy mode

  // TGVButtonGroup *GroupEnergyMode = new TGVButtonGroup(buttonframe, "Hits drawing");
  TGGroupFrame *GroupEnergyMode = new TGGroupFrame(buttonframe, "Hits drawing", kVerticalFrame|kRaisedFrame);
  fRadioEnergyMode[0] = new TGRadioButton(GroupEnergyMode, new TGHotString("color energy mode"), 130);
  fRadioEnergyMode[1] = new TGRadioButton(GroupEnergyMode, new TGHotString("bar energy mode"), 131);
  fRadioEnergyMode[2] = new TGRadioButton(GroupEnergyMode, new TGHotString("particles mode (MC only)"), 132);

  for (int i = 0; i < 3; ++i) {
    GroupEnergyMode->AddFrame(fRadioEnergyMode[i], new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2));
    fRadioEnergyMode[i]->Connect("Clicked()", "SANDEventDisplay", this, "SetHitsMode()");
  }

  fRadioEnergyMode[fSetHitsMode]->SetState(kButtonDown);
  for (int i = 0; i < 2; ++i) fRadioEnergyMode[i]->SetState(kButtonDisabled);

  buttonframe->AddFrame(GroupEnergyMode, fLayoutExpX);

  // Hits type

  // TGVButtonGroup *GroupHitsType = new TGVButtonGroup(buttonframe, "Hits type");
  TGGroupFrame *GroupHitsType = new TGGroupFrame(buttonframe, new TGString("HitsType"), kVerticalFrame|kRaisedFrame);
  fRadioHitsType[0] = new TGRadioButton(GroupHitsType, new TGHotString("simulated hits"), 140);
  fRadioHitsType[1] = new TGRadioButton(GroupHitsType, new TGHotString("digitized hits"), 141);
  fRadioHitsType[2] = new TGRadioButton(GroupHitsType, new TGHotString("reconstructed hits"), 142);

  for (int i = 0; i < 3; ++i) {
    GroupHitsType->AddFrame(fRadioHitsType[i], new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2));
    fRadioHitsType[i]->Connect("Clicked()", "SANDEventDisplay", this, "SetHitsType()");
  }

  for (int i = 0; i < 3; ++i) fRadioHitsType[i]->SetState(kButtonDisabled);

  buttonframe->AddFrame(GroupHitsType, fLayoutExpX);


  // fEntryEventId->Connect("ValueSet(long)", "SANDEventDisplay", this, "SetEventId()");

  // TGHorizontalFrame *FrameEventNumber = new TGHorizontalFrame(buttonframe);
  // TGLabel *LblEventNumber = new TGLabel(FrameEventNumber, "Event number");
  // fEntryEventNumber       = new TGTextEntry(FrameEventNumber, new TGTextBuffer(100));
  // fEntryEventNumber->Resize(95, fEntryEventNumber->GetDefaultHeight());
  // fEntryEventNumber->SetAlignment(kTextRight);
  // fEntryEventNumber->Connect("ProcessedEvent(Event_t*)","EventViewer",this,"SelectText(Event_t*)");
  // fEntryEventNumber->Connect("ReturnPressed()","EventViewer",this,"SetEventNumber()");
  //
  // FrameEventNumber->AddFrame(LblEventNumber,    fLayoutLeftExpY);
  // FrameEventNumber->AddFrame(fEntryEventNumber, fLayoutRightExpY);
  // buttonframe->AddFrame(FrameEventNumber, fLayoutExpX);

  // Using selected brick

  // fCheckSelectedBrickMode = new TGCheckButton(buttonframe,
	// 				      "Zoom to selected brick", 104);
  // fCheckSelectedBrickMode->Associate(this);
  // if (fSetZoomToSelectedBrick)
  //   fCheckSelectedBrickMode->SetState(kButtonDown, kTRUE);
  // buttonframe->AddFrame(fCheckSelectedBrickMode, fLayoutExpX);
  //
  // fCheckDrawVetoTracks = new TGCheckButton(buttonframe,
	// 				   "Draw veto tracks", 174);
  // fCheckDrawVetoTracks->Associate(this);
  // fCheckDrawVetoTracks->SetState(kButtonDown, kTRUE);
  // buttonframe->AddFrame(fCheckDrawVetoTracks, fLayoutExpX);
  //
  // fCheckRemoveSelectedTracks = new TGCheckButton(buttonframe,
	// 					 "Remove selected tracks",175);
  // fCheckRemoveSelectedTracks->Associate(this);
  // fCheckRemoveSelectedTracks->SetState(kButtonUp, kTRUE);
  // buttonframe->AddFrame(fCheckRemoveSelectedTracks, fLayoutExpX);

 // Scanned tracks

  // TGVButtonGroup    *GroupScannedTrk = new TGVButtonGroup(buttonframe,"Scanned tracks");
  // TGHorizontalFrame *FrameSelectTrk  = new TGHorizontalFrame(GroupScannedTrk);
  //
  // TGTextButton *selectbutton[3];
  // selectbutton[0] = new TGTextButton(FrameSelectTrk, " All ",      171);
  // selectbutton[1] = new TGTextButton(FrameSelectTrk, " Prev Trk ", 172);
  // selectbutton[2] = new TGTextButton(FrameSelectTrk, " Next Trk ", 173);
  //
  // for (Int_t i = 0; i < 3; i++) {
  //   selectbutton[i]->Associate(this);
  //   FrameSelectTrk->AddFrame(selectbutton[i], fLayoutExpX);
  // }
  //
  // GroupScannedTrk->AddFrame(FrameSelectTrk, fLayoutExpX);
  // buttonframe->AddFrame(GroupScannedTrk, fLayoutExpX);

  // add navigate view buttons

  AddNavigateButtons(buttonframe);

  // create frame for fitting hits

  //  fTTCSConnection->FittingHitsFrame(buttonframe);

  workframe->AddFrame(buttonframe, new TGLayoutHints(kLHintsLeft, 2, 4, 0, 0));
  // workframe->AddFrame(buttonframe, fLayoutExpX);
}


//----------------------------------------------------------------------------
void SANDEventDisplay::AddNavigateButtons(TGVerticalFrame *workframe)
{
  // View

  TGPictureButton *pictButton[4];
  TGHorizontalFrame *viewNavigate1 = new TGHorizontalFrame(workframe);
  TGHorizontalFrame *viewNavigate2 = new TGHorizontalFrame(workframe);
  TString icpath = "icons/";

  pictButton[0] = new TGPictureButton(viewNavigate1, gClient->GetPicture(icpath+"up.xpm"));
  pictButton[1] = new TGPictureButton(viewNavigate2, gClient->GetPicture(icpath+"left.xpm"));
  pictButton[2] = new TGPictureButton(viewNavigate2, gClient->GetPicture(icpath+"down.xpm"));
  pictButton[3] = new TGPictureButton(viewNavigate2, gClient->GetPicture(icpath+"right.xpm"));

  // pictButton[0]->Connect("Clicked()","SANDEventDisplay",this,"ShiftDetector(=2)");
  // pictButton[1]->Connect("Clicked()","SANDEventDisplay",this,"ShiftDetector(=-1)");
  // pictButton[3]->Connect("Clicked()","SANDEventDisplay",this,"ShiftDetector(=1)");
  // pictButton[2]->Connect("Clicked()","SANDEventDisplay",this,"ShiftDetector(=-2)");

  viewNavigate1->AddFrame(pictButton[0], fLayoutExpY);
  viewNavigate2->AddFrame(pictButton[1], fLayoutExpY);
  viewNavigate2->AddFrame(pictButton[2], fLayoutExpY);
  viewNavigate2->AddFrame(pictButton[3], fLayoutExpY);

  // navigating axis

//   TGVButtonGroup *GroupShiftAxis = new TGVButtonGroup(workframe,
// 						      "Navigating axis");
  // TGHorizontalFrame *FrameShiftAxis  = new TGHorizontalFrame(workframe);
  // TGLabel *LblShiftAxis = new TGLabel(FrameShiftAxis, "Axis:");
  // TGHButtonGroup *GroupShiftAxis = new TGHButtonGroup(FrameShiftAxis);
  //
  // TGRadioButton *RadioShiftAxis[3];
  // RadioShiftAxis[0] = new TGRadioButton(GroupShiftAxis,"X   ", 140);
  // RadioShiftAxis[1] = new TGRadioButton(GroupShiftAxis,"Y   ", 141);
  // RadioShiftAxis[2] = new TGRadioButton(GroupShiftAxis,"XY  ", 142);
  // RadioShiftAxis[fSetShiftAxis]->SetState(kButtonDown);
  //
  // for (Int_t i = 0; i < 3; i++) RadioShiftAxis[i]->Associate(this);

  // FrameShiftAxis->AddFrame(LblShiftAxis, fLayoutLeftExpY);
  // FrameShiftAxis->AddFrame(GroupShiftAxis, fLayoutExpX);

  workframe->AddFrame(viewNavigate1, new TGLayoutHints(kLHintsCenterX | kLHintsExpandY, 2, 2, 2, 2));
  workframe->AddFrame(viewNavigate2, new TGLayoutHints(kLHintsCenterX | kLHintsExpandY, 2, 2, 2, 2));
  // workframe->AddFrame(FrameShiftAxis, fLayoutExpX);
}


//----------------------------------------------------------------------------
void SANDEventDisplay::AddCanvasFrame(TGHorizontalFrame *workframe)
{
  // TGCompositeFrame *tf;

  // Making optimal window size

  int displayWidth = 1600; // gClient->GetDisplayWidth() - 251;
  int displayHeight = 800; // gClient->GetDisplayHeight() - 205;

  // Create frame and canvases for detector view and information

  // TGVerticalFrame *CanvasFrame = new TGVerticalFrame(workframe, displayWidth, displayHeight);

  // fTabCanvas = new TGTab(CanvasFrame);

  // Create a Tab with detector view

  // tf = fTabCanvas->AddTab("Detector");
  TRootEmbeddedCanvas *DisplayDetector = new TRootEmbeddedCanvas("DisplayDetector", workframe, displayWidth, displayHeight);

  // DisplayDetector->GetCanvas()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "SANDEventDisplay", this, "CanvasInfo(Int_t, Int_t, Int_t, TObject*)");
  // tf->AddFrame(fDisplayDetector, fLayoutExpXExpY);
  // AddFrame(fDisplayDetector);
  // CanvasFrame->AddFrame(fDisplayDetector, fLayoutExpXExpY);

  fDisplayCanvas = DisplayDetector->GetCanvas();

  // fTTCSConnection->SetCanvas(fDisplayCanvas);

  // create CS frame

  // fTTCSConnection->CreateCSCanvas(fTabCanvas, displayWidth, displayHeight);

  // Create a Tab with event information

  // tf = fTabCanvas->AddTab("Info");
  // fTextInfo = new TGTextView(tf, displayWidth, displayHeight);
  // const TGFont *font = gClient->GetFont("-*-courier-bold-r-*-*-14-140-75-75-*-*-*-*");
  // fTextInfo->SetFont(font->GetFontStruct());
  // tf->AddFrame(fTextInfo, fLayoutExpXExpY);

  // CanvasFrame->AddFrame(fTabCanvas, fLayoutExpXExpY);
  workframe->AddFrame(DisplayDetector, fLayoutExpXExpY);
}


//----------------------------------------------------------------------------
bool SANDEventDisplay::ProcessMessage(Long_t msg, Long_t param1, Long_t)
{
  // Handle messages send to the TDisplayMainFrame object.
  // E.g. all menu button messages.

  TString picture_format;
  // const char *filetypes[] = {"ROOT files", "*.root", 0, 0};

  // cout << param1 << " " << GET_MSG(msg) << " " << GET_SUBMSG(msg) << " " << kC_COMMAND << " " << kCM_RADIOBUTTON << endl;
  // cout << M_PREVIOUS_EVENT << " " << kCM_MENU << endl;


  switch (GET_MSG(msg)) {
    case kC_COMMAND:
    switch (GET_SUBMSG(msg)) {
      //   case kCM_COMBOBOX:
      //   if (param1 == 100 && fSM) DrawDetector(kFALSE);
      //   break;
      //   case kCM_CHECKBUTTON:
      //   if (param1 == 101) fSetNavigateMode=fCheckNavigateMode->IsDown();
      //   else if (param1 == 103) {
      //     fSetBrickFinderInfo = fCheckBrickFinderInfo->IsDown();
      //     DrawDetector();
      //   }
      //   else if (param1 == 104) SetViewToBrick();
      //   else if (param1 == 105) {
      //     fSetPresentationMode = fCheckPresentationMode->IsDown();
      //     DrawDetector(kFALSE);
      //   }
      //   else if (param1 >= 120 && param1 <= 123 && fSM) DrawDetector(kFALSE);
      //   else if (param1 == 174) {
      //     fTTCSConnection->SetDrawVetoTracks(fCheckDrawVetoTracks->IsDown());
      //     DrawDetector(kFALSE);
      //   }
      //   else if (param1 == 175) {
      //     fTTCSConnection->SetRemoveSelectedTracks(fCheckRemoveSelectedTracks->IsDown());
      //     DrawDetector(kFALSE);
      //   }
      //   break;
//      case kCM_RADIOBUTTON:
      //   if (param1 == 100) fSetZoomTo = 0;
      //   else if (param1 == 101) fSetZoomTo = 1;
      //   else if (param1 == 102) fSetZoomTo = 2;
      // if (param1 == 130) fSetHitsMode = 0;
      // else if (param1 == 131) fSetHitsMode = 1;
      // else if (param1 == 132) fSetHitsMode = 2;
      // else if (param1 == 140) fDrawHits = kSimHits;
      // else if (param1 == 141) fDrawHits = kDigitHits;
      // else if (param1 == 142) fDrawHits = kRecHits;
      //   else if (param1 == 140) fSetShiftAxis = 0;
      //   else if (param1 == 141) fSetShiftAxis = 1;
      //   else if (param1 == 142) fSetShiftAxis = 2;
      // if (param1 >= 130 && param1 <= 142) Run();
      // break;
      //   case kCM_LISTBOX:
      //   if (param1 == 125 && fSM) DrawDetector(kFALSE);
      //   break;
      //   case kCM_BUTTON:
      //   if (param1 == 171) fTTCSConnection->SetStep(0);
      //   else if (param1 == 172) fTTCSConnection->SetStep(-1);
      //   else if (param1 == 173) fTTCSConnection->SetStep(1);
      //   if (param1 >= 171 && param1 <= 173) DrawDetector(kFALSE);
      // break; // We do not needed "break" in this case!
      case kCM_BUTTON:
      switch (param1) {
        // case M_FILE_OPEN: {
        //   static TString dir(fInputDataDirName);
        //   TGFileInfo fi;
        //   fi.fFileTypes = filetypes;
        //   fi.fIniDir    = StrDup(dir.Data());
        //
        //   new TGFileDialog(gClient->GetRoot(), fMainFrame, kFDOpen, &fi);
        //
        //   dir = fi.fIniDir;
        //
        //   // open new file
        //
        //   if (fi.fFilename) {
        //     fEventManager->ClearEventList();
        //     if (fEventManager->InitTreeFile((TString)fi.fFilename)) {
        //       NextEvent();
        //       UpdateFileInfo();
        //     }
        //   }
        // }
        // break;
        // case M_SAVE_PICTURE:
        // SaveDisplayPicture(".");
        // break;
        // case M_SAVE_PICTURE_PNG:
        // picture_format = fSetPictureFormat;
        // fSetPictureFormat = "png";
        // SaveDisplayPicture(".");
        // fSetPictureFormat = picture_format;
        // break;
        // case M_SAVE_PICTURE_PDF:
        // picture_format = fSetPictureFormat;
        // fSetPictureFormat = "pdf";
        // SaveDisplayPicture(".");
        // fSetPictureFormat = picture_format;
        // break;
        case M_PREVIOUS_EVENT:
        PreviousEvent();
        break;
        // case M_EXECUTE_EVENT:
        // ProcessEvent();
        // break;
        case M_NEXT_EVENT:
        NextEvent();
        break;
        // case M_ZOOM_IN:
        // ZoomDetector(1.07);
        // break;
        // case M_ZOOM_OUT:
        // ZoomDetector(0.93);
        // break;
        // case M_ZOOM_EVENT:
        // FitZoomToEvent();
        // break;
        // case M_ZOOM_VERTEX:
        // FitZoomToVertex();
        // break;
        // case M_ZOOM_DETECTOR:
        // ZoomDetector();
        // break;
        case M_FILE_EXIT:
        ExitApplication();
        break;
        default:
        break;
      }
      default:
      break;
    }
    break;
    default:
    break;
  }

  return true;
}


//----------------------------------------------------------------------------
void SANDEventDisplay::NextEvent()
{
  fEventNumber = TMath::Min(fEventNumber+1, fTreeSimData->GetEntries()-1);
  Run();
}


//----------------------------------------------------------------------------
void SANDEventDisplay::PreviousEvent()
{
  if (fEventNumber > 0) --fEventNumber;
  Run();
}


//---------------------------------------------------------------------------
void SANDEventDisplay::SetHitsMode()
{
  // Handle radio buttons.

  TGButton *btn = (TGButton *) gTQSender;
  int id = btn->WidgetId();

  if (id == 130) fSetHitsMode = 0;
  else if (id == 131) fSetHitsMode = 1;
  else if (id == 132) fSetHitsMode = 2;

  if (id >= 130 && id <= 132) {
    for (int i = 0; i < 3; i++)
    if (fRadioEnergyMode[i]->WidgetId() != id && fRadioEnergyMode[i]->GetState() != kButtonDisabled) fRadioEnergyMode[i]->SetState(kButtonUp);
  }

  Run();
}


//---------------------------------------------------------------------------
void SANDEventDisplay::SetHitsType()
{
  // Handle radio buttons.

  TGButton *btn = (TGButton *) gTQSender;
  int id = btn->WidgetId();

  if (id == 140) fDrawHits = kSimHits;
  else if (id == 141) fDrawHits = kDigitHits;
  else if (id == 142) fDrawHits = kRecHits;

  if (id >= 140 && id <= 142) {
    for (int i = 0; i < 3; i++)
    if (fRadioHitsType[i]->WidgetId() != id && fRadioHitsType[i]->GetState() != kButtonDisabled) fRadioHitsType[i]->SetState(kButtonUp);
  }

  Run();
}


//----------------------------------------------------------------------------
void SANDEventDisplay::SetEventId()
{
  long long eventId = atoi(fEntryEventId->GetNumberEntry()->GetText());
  fEventNumber = TMath::Range(0, fTreeSimData->GetEntries()-1, eventId);
  fEntryEventId->GetNumberEntry()->SetIntNumber(fEventNumber);
  Run();
}


//----------------------------------------------------------------------------
void SANDEventDisplay::SetSimData(TString fileName)
{
  TGeoManager* geo = 0;

  try {
    fFileSimData = new TFile(fileName, "READ");
    try {
      if (!fGeomInitialized) geo = static_cast<TGeoManager*>(fFileSimData->Get("EDepSimGeometry"));
      fTreeSimData = static_cast<TTree*>(fFileSimData->Get("EDepSimEvents"));
    }
    catch(...)
    {
      cout << "<ERROR> MC info not found in " << fFileSimData->GetName() << endl;
      return;
    }
  }
  catch(...)
  {
    cout << "<ERROR> MC file can not be opened " << fileName << endl;
    return;
  }

  if (!fGeomInitialized) {
    sand_reco::init(geo);
    fGeomInitialized = true;
  }

  fRadioHitsType[kSimHits]->SetState(kButtonEngaged);
  fRadioHitsType[fDrawHits]->SetState(kButtonDown);

  fTreeSimData->SetBranchAddress("Event", &fEvent);
}


//----------------------------------------------------------------------------
void SANDEventDisplay::SetDigitData(TString fileName)
{
  TGeoManager* geo = 0;

  try {
    fFileDigitData = new TFile(fileName, "READ");
    try {
      if (!fGeomInitialized) geo = static_cast<TGeoManager*>(fFileSimData->Get("EDepSimGeometry"));
      fTreeDigitData = static_cast<TTree*>(fFileDigitData->Get("tDigit"));
    }
    catch(...)
    {
      cout << "<ERROR> Digit info not found in " << fFileDigitData->GetName() << endl;
      return;
    }
  }
  catch(...)
  {
    cout << "<ERROR> Digit file can not be opened " << fileName << endl;
    return;
  }

  if (!fGeomInitialized) {
    sand_reco::init(geo);
    fGeomInitialized = true;
  }

  fRadioHitsType[kDigitHits]->SetState(kButtonEngaged);
  fRadioHitsType[fDrawHits]->SetState(kButtonDown);

  fTubeDigitVect = new std::vector<dg_tube>;
  fCellDigitVect = new std::vector<dg_cell>;

  fTreeDigitData->SetBranchAddress("dg_tube", &fTubeDigitVect);
  fTreeDigitData->SetBranchAddress("dg_cell", &fCellDigitVect);
}


//----------------------------------------------------------------------------
void SANDEventDisplay::DefineColors()
{
  for (int iColor = 0; iColor < fColNum; iColor++) {

    float cl3   =  1 - 2.*iColor/fColNum;
    float cl1   = -1 + 2.*iColor/fColNum;
    float cl2_1 = 2.*(iColor-1)/fColNum;
    float cl2_2 = 1 - 2.*(iColor-8)/fColNum;
    float cl2   = cl2_1;

    if (cl3<0) {cl3 = 0; cl2 = cl2_2;}
    if (cl1<0) {cl1 = 0; cl2 = cl2_1;}
    if (cl2<0) {cl2 = 0;}

    if (!gROOT->GetColor(300+iColor))
      fColor = new TColor(300+iColor, cl1, cl2, cl3, "");
    else {
      fColor = gROOT->GetColor(300+iColor);
      fColor->SetRGB(cl1, cl2, cl3);
    }
    fPalette[iColor] = 300 + iColor;
  }
}


//----------------------------------------------------------------------------
void SANDEventDisplay::SetHistoStyle(TH2F *histo)
{
  histo->GetYaxis()->CenterTitle(true);

  histo->GetXaxis()->SetLabelFont(63);
  histo->GetXaxis()->SetLabelSize(23);
  histo->GetYaxis()->SetLabelFont(63);
  histo->GetYaxis()->SetLabelSize(23);

  histo->GetXaxis()->SetTitleFont(63);
  histo->GetXaxis()->SetTitleSize(23);
  histo->GetYaxis()->SetTitleFont(63);
  histo->GetYaxis()->SetTitleSize(27);

  histo->SetTitleOffset(1.2,"X");
  histo->SetTitleOffset(1.2,"Y");

  histo->GetXaxis()->SetNdivisions(507, true);
  histo->GetYaxis()->SetNdivisions(505, true);

  histo->SetStats(false);
}


//----------------------------------------------------------------------------
void SANDEventDisplay::ExitApplication()
{
  gApplication->Terminate();
}
