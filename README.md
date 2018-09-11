============================================
Digitization
- Create digits of STT and cells of calorimeter
============================================
mtenti@neut05:KLOEcal (master)$ root -l
root [0] .L digitization/digitization.cpp+
void Digitize(const char* finname, const char* foutname)
root [1] Digitize("muon-LArTarget_volume-geometry_v12.1E7-CC/muon-LArTarget_volume-geometry_v12.1E7-CC.*.root","output.digits.root")

============================================
Reconstruction
- Track find and fit of STT track
- Clustering of calorimeter cells
============================================
root [0] .L reconstruction/reconstruction.cpp+
void Reconstruct(const char* fDigit, const char* fTrueMC, const char* fOut)
root [1] Reconstruct("output.digits.root","muon-LArTarget_volume-geometry_v12.1E7-CC/muon-LArTarget_volume-geometry_v12.1E7-CC.*.root","output.reco.root")

============================================
Analysis
- Evaluate parameters of particle
============================================
mtenti@neut05:KLOEcal (master)$ root -l
root [0] .L analysis/analysis.cpp+
void analysis(const char* fReco, const char* fTrueMC, const char* fOut)
root [1] analysis("output.reco.root","muon-LArTarget_volume-geometrroot [2] analysis(LArTarget_volume-geometry_v12.1E7-CC.*.root","output.analysis.root")


