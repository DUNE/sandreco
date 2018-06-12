double rmass = 0.; // kg

bool initok = false;

TGeoManager* geo = 0;

void init(const char* file)
{
  geo = TGeoManager::Import(file);
  initok = true;
}


void searchRad(TGeoVolume* vol, const char* volname)
{
  TString name(vol->GetName());
  
  if(name.Contains(volname))
  {
    rmass += vol->Weight();
    std::cout << "Node: " << vol->GetName() << std::endl;
  } 
  
  for(int i = 0; i < vol->GetNdaughters(); i++)
  {
    searchRad(vol->GetNode(i)->GetVolume(), volname);
  }
}

void radmass(const char* volname)
{
  if(!initok || geo == 0)
  {
    std::cout << "Geometry not initializated" << std::endl;
    return;
  }
  
  rmass = 0.;
  
  searchRad(geo->GetTopVolume(), volname);
  
  std::cout << "Total radiator mass: " << rmass << " kg" << std::endl;
  
}
