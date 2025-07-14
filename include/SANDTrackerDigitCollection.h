#ifndef SANDTrackerDIGITCOLLECTION_H
#define SANDTrackerDIGITCOLLECTION_H

#include "struct.h"
#include "SANDGeoManager.h"

#include <TTreeReader.h>

namespace sand_reco
{
namespace tracker
{
// digit id -> dg_wire.did
class DigitID : public SingleElStruct<long>
{
 public:
  DigitID(long id) : SingleElStruct<long>(id){};
  DigitID() : SingleElStruct<long>(){};
};

// digit index -> index inside Digit vector
class DigitIndex : public SingleElStruct<long>
{
 public:
  DigitIndex(long id) : SingleElStruct<long>(id){};
  DigitIndex() : SingleElStruct<long>(){};
};

// Digit
using Digit = dg_wire;

// Digit map: key: digit id; value: index in gTreeReaderDigit
using DigitMap = std::map<DigitID, DigitIndex>;

/**********************************************
 * Class to read and access digits (Digit)
 * from the input tree through a TTreeReaderValue
 ***********************************************/
class DigitCollection
{
 private:
  // vector of Digits
  static std::vector<Digit> sand_fg_tracker_digits_;

  // digit map -> key: digit id; value: index in gTreeReaderDigit
  static DigitMap fg_map_digit_;

 public:
  DigitCollection(){};
  ~DigitCollection(){};

  // fill digit map
  static void fillMap(const std::vector<Digit>* digits)
  {
    sand_fg_tracker_digits_ = *digits;
    for (auto i = 0u; i < sand_fg_tracker_digits_.size(); i++) {
      fg_map_digit_[DigitID(static_cast<long>(
          sand_fg_tracker_digits_.at(i).did))] = DigitIndex(i);
    }
  };

  // get digit vector
  static const std::vector<Digit> &getDigits()
  {
    return sand_fg_tracker_digits_;
  };

  // get i-th digit
  static const Digit &getDigit(const DigitID &id)
  {
    // std::cout << "DIGIT COLLECTION: " << id() << " " << fg_map_digit_[id]() << std::endl;
    return sand_fg_tracker_digits_.at(fg_map_digit_[id]());
  };
};
} // namespace tracker
} // namespace sand_reco
#endif
