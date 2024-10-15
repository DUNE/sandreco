#include "SANDKFTrack.h"
#include "SANDTrackerUtils.h"

void SANDKFTrackStep::SetStage(SANDKFTrackStateStage stage, SANDKFState state) {
  switch (stage)
  {
    case SANDKFTrackStateStage::kPrediction:
      fPrediction = state;
      break;
    case SANDKFTrackStateStage::kFiltering:
      fFiltered = state;
      break;
    case SANDKFTrackStateStage::kSmoothing:
      fSmoothed = state;
      break;
  }
}

const SANDKFState& SANDKFTrackStep::GetStage(SANDKFTrackStateStage stage) const {
  switch (stage)
  {
    case SANDKFTrackStateStage::kPrediction:
      return fPrediction;
    case SANDKFTrackStateStage::kFiltering:
      return fFiltered;
    case SANDKFTrackStateStage::kSmoothing:
      return fSmoothed;
    default:
      // SANDTRACKRECO_LOG("ERROR", "Unknown track state stage. Code should never reach this part");
      return fPrediction;
  }
}