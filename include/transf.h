#include <TLorentzVector.h>

#ifndef TRANSF_H
#define TRANSF_H

TLorentzVector LocalToGlobalCoordinates(TLorentzVector pos);
TLorentzVector GlobalToLocalCoordinates(TLorentzVector pos);

#endif