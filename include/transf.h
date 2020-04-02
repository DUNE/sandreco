//for transforming coordinates 
TLorentzVector LocalToGlobalCoordinates(TLorentzVector pos){
        TLorentzVector finalpos;
        finalpos.SetXYZT(pos.X(),
                        pos.Y() -2384.73,
                        pos.Z() +23910,
                        pos.T());
        return finalpos;
}

TLorentzVector GlobalToLocalCoordinates(TLorentzVector pos){
        TLorentzVector finalpos;
        finalpos.SetXYZT(pos.X(),
                        pos.Y() +2384.73,
                        pos.Z() -23910,
                        pos.T());
        return finalpos;
}



