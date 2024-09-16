#ifndef SANDEVENTDISPLAY_H
#define SANDEVENTDISPLAY_H

#include <TCanvas.h>
#include <TGeoManager.h>
#include <TH2F.h>
#include <TLine.h>
#include <TMarker.h>
// #include <TG4Event.h>
#include <TApplication.h>
#include <TGFrame.h>
#include <TGTab.h>
#include <TRootEmbeddedCanvas.h>
#include <vector>

// #include "STTCluster.h"
// #include "STTKFTrack.h"

// #include "utils.h"
#include "struct.h"

class TClonesArray;
class TPolyMarker;
class TDatabasePDG;
class TCanvas;
class TGNumberEntry;
class TGRadioButton;
class TFile;
class TTree;
class TColor;
class TG4Event;

enum DetectorType_t { kECAL, kGRAIN, kSTT };

enum HitsType_t { kSimHits, kDigitHits, kRecHits };

typedef struct {
  DetectorType_t detector;  // detector type
  int particle;             // particle code
  double x;                 // X or Y coordinates of the digits
  double z;                 // Z coordinates of the digits
  double e;                 // energy deposited in the digits
  double driftdist;         // drift distance
  Color_t color;            // Drawing color in Event Viewer
} EVHits_t;

class SANDEventDisplay : public TGMainFrame
{
 public:
  SANDEventDisplay(const TGWindow *p, int w, int h);
  ~SANDEventDisplay();

  // static void Init(TGeoManager* geo) {
  //   sand_reco::init(geo);
  //   STTStrawTubeTracker::Init();
  // };
  // static TCanvas* Create2DDisplaySideBySide();
  // static TMarker* GetMarkerFromCluster(const STTCluster& cluster, EColor
  // color); static TLine** GetFilterToPredictionLines(const STTKFTrack& tr, int
  // stepIndex, EColor color);

  void SetEventId();
  void SetEventNumber(long long eventNumber) { fEventNumber = eventNumber; }
  void SetSimData(TString fileName);
  void SetDigitData(TString fileName);
  void Run();
  void NextEvent();
  void PreviousEvent();
  void SetHitsMode();
  void SetHitsType();
  bool ProcessMessage(Long_t msg, Long_t param1, Long_t);
  // void   AddNavigateButtons(TGVerticalFrame *workframe);
  void ExitApplication();

 private:
  int fColNum;
  int fPalette[100];
  int fSetHitsMode;
  long long fEventNumber;
  bool fGeomInitialized;
  HitsType_t fDrawHits;
  TG4Event *fEvent;
  TDatabasePDG *fPDGcode;
  TString fSimFileName;
  TString fDiditFileName;
  TFile *fFileSimData;
  TFile *fFileDigitData;
  TTree *fTreeSimData;
  TTree *fTreeDigitData;
  TColor *fColor;

  std::vector <dg_tube> *fTubeDigitVect;
  std::vector <dg_wire> *fWireDigitVect;
  std::vector <dg_cell> *fCellDigitVect;

  TClonesArray *fTracksArrayZYTrue;
  TClonesArray *fTracksArrayZXTrue;
  TPolyMarker *fVerticesZYTrue;
  TPolyMarker *fVerticesZXTrue;
  TCanvas *fDisplayCanvas;
  TPad *fPadZY;
  TPad *fPadZX;

  TGVerticalFrame *fMainFrame;
  TGTab *fTab;
  TGLayoutHints *fLayoutExpX;
  TGLayoutHints *fLayoutExpY;
  TGLayoutHints *fLayoutLeftExpY;
  TGLayoutHints *fLayoutRightExpY;
  TGLayoutHints *fLayoutLeftExpX;
  TGLayoutHints *fLayoutRightExpX;
  TGLayoutHints *fLayoutExpXExpY;
  TGNumberEntry *fEntryEventId;
  TGRadioButton *fRadioEnergyMode[3];
  TGRadioButton *fRadioHitsType[3];

  // TRootEmbeddedCanvas *fDisplayDetector;

  std::vector <EVHits_t>    fEventHitsZY;
  std::vector <EVHits_t>    fEventHitsZX;
  std::vector <EVHits_t>    fTubeDigitHitsZY;
  std::vector <EVHits_t>    fTubeDigitHitsZX;
  std::vector <EVHits_t>    fWireDigitHitsZY;
  std::vector <EVHits_t>    fWireDigitHitsZX;
  std::vector <EVHits_t>    fCellDigitHitsZY;
  std::vector <EVHits_t>    fCellDigitHitsZX;

  void   InitObjects(); // init drawing objects
  void   FillEventHits();
  void   FillDigitHits();
  void   DrawEvent();
  void   DrawDetector();
  void   DrawTracks();
  void   SetHistoStyle(TH2F *histo);
  void   DrawButtons();
  void   DefineColors();
  void   AddMenuBar(TGVerticalFrame *workframe);
  void   AddRunEventFrame(TGHorizontalFrame *workframe);
  void   AddCanvasFrame(TGHorizontalFrame *workframe);
  void   AddNavigateButtons(TGVerticalFrame *workframe);
  // void   DrawEllipse(TEllipse *ellipse, Color_t color, Style_t style);
  // void   DrawBox(TBox *box, Color_t color, Style_t style);

  ClassDef(SANDEventDisplay, 0)
};

#endif
