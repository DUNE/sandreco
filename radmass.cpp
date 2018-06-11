double rmass = 0.; // kg

void searchRad(TGeoVolume* vol)
{
  TString name(vol->GetName());
  
  if(name.Contains("Radiator"))
  {
    rmass += vol->Weight();
    std::cout << "Node: " << vol->GetName() << std::endl;
  } 
  
  for(int i = 0; i < vol->GetNdaughters(); i++)
  {
    searchRad(vol->GetNode(i)->GetVolume());
  }
}

void radmass(const char* file)
{
  TGeoManager* geo = TGeoManager::Import(file);
  
  rmass = 0.;
  
  searchRad(geo->GetTopVolume());
  
  std::cout << "Total radiator mass: " << rmass << " kg" << std::endl;
  
}
