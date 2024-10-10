#ifndef SANDTrackerDIGITCOLLECTION_H
#define SANDTrackerDIGITCOLLECTION_H

#include "struct.h"
#include "SANDGeoManager.h"

#include <TTreeReader.h>

// digit id -> dg_wire.did
class SANDTrackerDigitID : public SingleElStruct<unsigned long>
{
 public:
  SANDTrackerDigitID(unsigned long id) : SingleElStruct<unsigned long>(id){};
  SANDTrackerDigitID() : SingleElStruct<unsigned long>(){};
};

// digit index -> index inside SANDTrackerDigit vector
class SANDTrackerDigitIndex : public SingleElStruct<unsigned long>
{
 public:
  SANDTrackerDigitIndex(unsigned long id) : SingleElStruct<unsigned long>(id){};
  SANDTrackerDigitIndex() : SingleElStruct<unsigned long>(){};
};

// SANDTrackerDigit
using SANDTrackerDigit = dg_wire;

// Digit map: key: digit id; value: index in gTreeReaderDigit
using SANDTrackerDigitMap = std::map<SANDTrackerDigitID, SANDTrackerDigitIndex>;

/**********************************************
 * Class to read and access digits (SANDTrackerDigit)
 * from the input tree through a TTreeReaderValue
 ***********************************************/
class SANDTrackerDigitCollection
{
 private:
  // vector of SANDTrackerDigits
  static std::vector<SANDTrackerDigit> SANDfgTrackerDigits;

  // digit map -> key: digit id; value: index in gTreeReaderDigit
  static SANDTrackerDigitMap fgMapDigit;

 public:
  SANDTrackerDigitCollection(){};
  ~SANDTrackerDigitCollection(){};

  // fill digit map
  static void FillMap(const std::vector<SANDTrackerDigit>* digits)
  {
    SANDfgTrackerDigits = *digits;
    for (auto i = 0u; i < SANDfgTrackerDigits.size(); i++) {
      fgMapDigit[SANDTrackerDigitID(static_cast<unsigned long>(
          SANDfgTrackerDigits.at(i).did))] = SANDTrackerDigitIndex(i);
    }
  };

  // get digit vector
  static const std::vector<SANDTrackerDigit> &GetDigits()
  {
    return SANDfgTrackerDigits;
  };

  // get i-th digit
  static const SANDTrackerDigit &GetDigit(const SANDTrackerDigitID &id)
  {
    // std::cout << "DIGIT COLLECTION: " << id() << " " << fgMapDigit[id]() << std::endl;
    return SANDfgTrackerDigits.at(fgMapDigit[id]());
  };
};

#endif
