#include "STTKFTrack.h"
#include "SANDTrackerUtils.h"

void STTKFTrackStep::SetStage(STTKFTrackStateStage stage, STTKFState state) {
  switch (stage)
  {
    case STTKFTrackStateStage::kPrediction:
      fPrediction = state;
      break;
    case STTKFTrackStateStage::kFiltering:
      fFiltered = state;
      break;
    case STTKFTrackStateStage::kSmoothing:
      fSmoothed = state;
      break;
  }
}

const STTKFState& STTKFTrackStep::GetStage(STTKFTrackStateStage stage) const {
  switch (stage)
  {
    case STTKFTrackStateStage::kPrediction:
      return fPrediction;
    case STTKFTrackStateStage::kFiltering:
      return fFiltered;
    case STTKFTrackStateStage::kSmoothing:
      return fSmoothed;
    default:
      // STTTRACKRECO_LOG("ERROR", "Unknown track state stage. Code should never reach this part");
      return fPrediction;
  }
}