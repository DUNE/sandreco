#include "SANDClustering.h"

int clustering(std::string const& input, std::string const& edep_input)
{
  
  const char* finname = input.c_str();
  TFile f(finname, "READ");
  TTree* t = (TTree*)f.Get("tDigit");
  if (f.IsZombie()){
    std::cout << "Error in opening file\n";
    exit(1);
  }
  
  if (t == nullptr ){
    std::cout << "Error in retrieving objects from root file: tDigit" << std::endl; 
    exit(-1);
  }

  const char* edep_name = edep_input.c_str();
  TFile f_edep(edep_name, "READ");
  TGeoManager* geo = (TGeoManager*)f_edep.Get("EDepSimGeometry");
  SANDGeoManager sand_geo;
  sand_geo.init(geo);

  int nEvents = t->GetEntries();
  std::vector<dg_cell>* cell = new std::vector<dg_cell>;
  std::vector<cluster> f_clust, og_clust;
  TString output = input;
  output.ReplaceAll(".root", ".Clusters.root");
  TFile fout(output, "RECREATE");
  TTree tout("tCluster", "Clustering");
  tout.Branch("cluster", "std::vector<cluster>", &f_clust);
  t->SetBranchAddress("dg_cell", &cell);
  
  std::cout << "Events: " << nEvents << " [";
  std::cout << std::setw(3) << int(0) << "%]" << std::flush;

  for (int i = 0; i < nEvents; i++) {
    t->GetEntry(i);

    std::cout << "\b\b\b\b\b" << std::setw(3) << int(double(i) / nEvents * 100)
              << "%]" << std::flush;

    std::vector<cluster> clust = clusterize(&sand_geo, *cell);
    
    f_clust = clust;
    tout.Fill();
    clust.clear();
    f_clust.clear();
    og_clust.clear();
  }
  std::cout << "\b\b\b\b\b" << std::setw(3) << 100 << "%]" << std::flush;
  std::cout << std::endl;

  fout.cd();
  tout.Write();
  fout.Close();
  delete cell;
  return 0;
}

bool initializeFiles(int argc, char* argv[], std::string& digitFileName, std::string& edepFileName)
{

  for (int i = 1; i < argc; i += 2) {
    std::string flag = argv[i];
    std::string fileName = argv[i + 1];

    if (flag == "-d") {
      digitFileName = fileName;
    } else if (flag == "-e") {
      edepFileName = fileName;
    } else {
      std::cerr << "Error: Unknown flag: " << flag << std::endl;
      return false;
    }
  }

  if (digitFileName.empty()) {
    std::cerr << "Error: Missing required flags. Please use '-d' with "
                 "corresponding file name."
              << std::endl;
    return false;
  }
  if (edepFileName.empty()) {
    std::cerr << "Error: Missing required flags. Please use '-e' with "
                 "corresponding file name."
              << std::endl;
    return false;
  }

  if (!endsWith(digitFileName, ".digit.root")) {
    std::cerr << "Error: Invalid arguments. Please use '-d' before the "
                 ".digit.root file"
              << std::endl;
    return false;
  }
  if (!endsWith(edepFileName, ".edep.root")) {
    std::cerr << "Error: Invalid arguments. Please use '-e' before the "
                 ".edep.root file"
              << std::endl;
    return false;
  }

  return true;
}

int main(int argc, char* argv[])
{

  if (argc < 2) {
    std::cerr << "Usage: SANDECALClustering -d <file.digit.root> -e <file.edep.root>" << std::endl;
    return 1;
  }
  std::string digitFileName;
  std::string ecalFileName;
  if (!initializeFiles(argc, argv, digitFileName, ecalFileName)) {
    return 1;
  }

  clustering(digitFileName, ecalFileName);
  return 0;
}
