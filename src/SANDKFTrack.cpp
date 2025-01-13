#include "SANDKFTrack.h"
#include "SANDTrackerUtils.h"

namespace sand_reco
{

namespace kf
{
void TrackStep::setStage(TrackStateStage stage, State state) {
  switch (stage)
  {
    case TrackStateStage::kPrediction:
      prediction_ = state;
      break;
    case TrackStateStage::kFiltering:
      filtered_ = state;
      break;
    case TrackStateStage::kSmoothing:
      smoothed_ = state;
      break;
  }
}

const State& TrackStep::getStage(TrackStateStage stage) const {
  switch (stage)
  {
    case TrackStateStage::kPrediction:
      return prediction_;
    case TrackStateStage::kFiltering:
      return filtered_;
    case TrackStateStage::kSmoothing:
      return smoothed_;
    default:
      // SANDTRACKRECO_LOG("ERROR", "Unknown track state stage. Code should never reach this part");
      return prediction_;
  }
}
}
}