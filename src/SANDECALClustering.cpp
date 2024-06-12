#include "SANDClustering.h"

int clustering(std::string const& input)
{
    // TH1F* ClusterEnergy = new TH1F("ClusterEnergy", "Total Energy reconstructed
    // in the ECal;Energy [MeV];Entries/bins", 100, 0., 1500.);
    //gROOT->ProcessLine("gInterpreter->AddIncludePath("-I/opt/exp_software/neutrino/EDEPSIM/include/EDepSim")");
  // gSystem->Load("libStruct.so");
    const char* finname = input.c_str();
    TFile f(finname, "READ");
    TTree* t = (TTree*)f.Get("tDigit");
    int nEvents = t->GetEntries();
    std::vector<dg_cell>* cell = new std::vector<dg_cell>;
    std::vector<cluster> f_clust, og_clust;
    TString output = input;
    output.ReplaceAll(".root", ".Clusters.root");
    TFile fout(output, "RECREATE");
    TTree tout("tCluster", "Clustering");
    tout.Branch("cluster", "std::vector<cluster>", &f_clust);
    t->SetBranchAddress("dg_cell", &cell);
    //infofile.open("info.txt");
    // info_new.open("muon_6marzo.txt");
    //nEvents = 10; //DC
    for (int i = 0; i < nEvents; i++) {
        //infofile<< "-----------" << std::endl;
        // infofile<< "Entry " << i << std::endl;
        t->GetEntry(i);
        std::vector<cluster> clust = Clusterize(std::move(cell));
        double CluEn = 0;
        // for (auto const& clu_info : clust) {
        //     Clust_info(clu_info);
            
        //     // if(clu_info.sy >0.2){
        //     //     info_new<< "Entry with cluster.sy> 0.2 : " << i << std::endl;  
        //     // }
        // }
        f_clust = clust;
        tout.Fill();
        clust.clear();
        f_clust.clear();
        og_clust.clear();
    }
    fout.cd();
    tout.Write();
    fout.Close();
    delete cell;
    return 0;
}

bool initializeFiles(int argc, char* argv[], std::string& digitFileName) {
    // Process command-line arguments to build the file names
    for (int i = 1; i < argc; i += 2) {
        std::string flag = argv[i];
        std::string fileName = argv[i + 1];

        if (flag == "-d") {
            digitFileName = fileName;
        } else {
            std::cerr << "Error: Unknown flag: " << flag << std::endl;
            return false;
        }
    }

    // Check if the required files are present
    if (digitFileName.empty()) {
        std::cerr << "Error: Missing required flags. Please use '-d' with corresponding file name." << std::endl;
        return false;
    }

    // Validate file names
    if (!endsWith(digitFileName, ".digit.root")) {
        std::cerr << "Error: Invalid arguments. Please use '-d' before the .digit.root file" << std::endl;
        return false;
    }

    return true; // Files initialized successfully
}

int main(int argc, char* argv[]){
  
  if (argc < 2) {
        std::cerr << "Usage: Clustering" << " -d <file.digit.root>" << std::endl;
        return 1; // Exit with an error code
    }
    std::string digitFileName;
    if (!initializeFiles(argc, argv, digitFileName)) {
        return 1; // Exit with an error code
    }
    
    clustering(digitFileName);
  return 0;
}
