

// frontST0_vol_hor_vol_PV_0
// sttmod0_ST_hor_vol_PV_0
// '(sttmod|frontST)([0-9]+)_(ST|vol)_(hor|ver)_vol_PV_0'

// ST front :
// 'frontST([0-9]+)_vol_(hor|ver)_ST_stGas_([a-zA-Z]{2})19_vol_PV_([0-9]+)'
// ST module:
// 'sttmod([0-9]+)_ST_(hor|ver)_ST_stGas_([a-zA-Z]{2})19_vol_PV_([0-9]+)'
// both:
// '(sttmod|frontST)([0-9]+)_(ST|vol)_(hor|ver)_ST_stGas_([a-zA-Z]{2})19_vol_PV_([0-9]+)'

// slab
// 'sttmod([0-9]+)_slab_vol_PV_0'

TPRegexp rST(
    "(sttmod|frontST)([0-9]+)_(ST|vol)_(hor|ver)_ST_stGas_([a-zA-Z]{2})19_vol_"
    "PV_([0-9]+)");
TPRegexp rSTplane("(sttmod|frontST)([0-9]+)_(ST|vol)_(hor|ver)_vol_PV_0");
TPRegexp rSlab("sttmod([0-9]+)_slab_vol_PV_0");

void showElemPos(TGeoNode* nod, TGeoHMatrix m, TString path);

bool isSlab(TString name)
{
  return name.Contains(rSlab);
}

bool isST(TString name)
{
  return name.Contains(rST);
}

bool isSTPlane(TString name)
{
  return name.Contains(rSTplane);
}

int getSTId(TString name)
{
  int id = -999;

  TObjArray* obj = name.Tokenize("_");

  if (obj->GetEntries() != 9) {
    std::cout << "Error: tokenizing " << name.Data() << std::endl;
  } else {
    TString sid = ((TObjString*)obj->At(8))->GetString();

    id = sid.Atoi();
  }
  delete obj;

  return id;
}

int getPlaneID(TString name)
{
  int mod = 0;
  int type = 0;

  TObjArray* obj = name.Tokenize("_");

  if (obj->GetEntries() != 6) {
    std::cout << "Error: tokenizing " << name.Data() << std::endl;
  } else {
    TString stype = ((TObjString*)obj->At(2))->GetString();
    TString smod = ((TObjString*)obj->At(0))->GetString();

    TObjArray* oarr;

    if (smod.Contains("frontST")) {
      oarr = smod.Tokenize("frontST");
      mod = 90 + ((TObjString*)oarr->At(0))->GetString().Atoi();
    } else if (smod.Contains("sttmod")) {
      oarr = smod.Tokenize("sttmod");
      mod = ((TObjString*)oarr->At(0))->GetString().Atoi();
    } else
      std::cout << "Error evaluating module id for: " << name.Data()
                << std::endl;

    if (stype.Contains("ver"))
      type = 1;
    else if (stype.Contains("hor"))
      type = 2;
    else
      std::cout << "Error evaluating type for: " << name.Data() << std::endl;

    delete oarr;
  }

  delete obj;

  return type + 10 * mod;
}

void getSTinfo(TGeoNode* nod, TGeoHMatrix mat, int pid,
               std::map<double, int>& stX, std::map<int, TVector2>& stPos)
{
  int ic = pid - int(double(pid) / 10.) * 10 - 1;

  if (ic != 0 && ic != 1)
    std::cout << "Error: ic expected 0 or 1 -> " << ic << std::endl;

  for (int i = 0; i < nod->GetNdaughters(); i++) {
    TString name = nod->GetDaughter(i)->GetName();

    if (!isST(name))
      std::cout << "Error: expected ST but not -> " << name.Data() << std::endl;

    TGeoMatrix* thismat = nod->GetDaughter(i)->GetMatrix();
    TGeoHMatrix mymat = mat * (*thismat);

    int id = getSTId(name);

    TVector2 v;
    v.SetX(mymat.GetTranslation()[2]);
    v.SetY(mymat.GetTranslation()[ic]);

    stX[v.Y()] = id;
    stPos[id] = v;
  }
}

void getSlabInfo(TGeoNode* nod, TGeoHMatrix mat, int& pid, double& z)
{
  TString name = nod->GetName();

  TObjArray* obj = name.Tokenize("_");

  if (obj->GetEntries() != 5) {
    std::cout << "Error: tokenizing " << name.Data() << std::endl;
  } else {
    TString smod = ((TObjString*)obj->At(0))->GetString();

    TObjArray* oarr = smod.Tokenize("sttmod");
    pid = ((TObjString*)oarr->At(0))->GetString().Atoi();

    z = mat.GetTranslation()[2];

    delete oarr;
  }

  delete obj;
}

void getSTPlaneinfo(TGeoNode* nod, TGeoHMatrix mat,
                    std::map<int, std::map<double, int> >& stX,
                    std::map<int, std::map<int, TVector2> >& stPos,
                    std::map<int, double>& slabPos)
{
  TString name = nod->GetName();
  TGeoMatrix* thismat = nod->GetMatrix();
  TGeoHMatrix mymat = mat * (*thismat);

  int pid = 0;
  double x = 0;
  double z = 0;

  if (isSTPlane(name)) {
    pid = getPlaneID(name);

    std::map<double, int> mstX;
    std::map<int, TVector2> mstPos;

    getSTinfo(nod, mymat, pid, mstX, mstPos);

    stX[pid] = mstX;
    stPos[pid] = mstPos;
  } else if (isSlab(name)) {
    getSlabInfo(nod, mymat, pid, z);
    slabPos[pid] = z;
  } else {

    for (int i = 0; i < nod->GetNdaughters(); i++) {
      getSTPlaneinfo(nod->GetDaughter(i), mymat, stX, stPos, slabPos);
    }
  }
}

void showElemPos(TGeoNode* nod, TGeoHMatrix mat, TString path)
{
  TString name = nod->GetName();
  TString mypath = path + "/" + name;
  TGeoMatrix* thismat = nod->GetMatrix();
  TGeoHMatrix mymat = mat * (*thismat);

  // std::cout << nod->GetName() << ": " << mypath.Data() << std::endl;
  // mymat.Print();

  int mod = 0;
  int type = 0;
  int id = 0;
  double x = 0;
  double z = 0;

  if ((name.Contains("stGas_Ar19_vol_PV_") ||
       name.Contains("stGas_Xe19_vol_PV_")) &&
      !name.Contains("coat") && !name.Contains("wire") &&
      !name.Contains("air") && !name.Contains("mylar")) {
    TObjArray* obj = name.Tokenize("_");

    if (obj->GetEntries() != 9) {
      std::cout << "Error: tokenizing " << name.Data() << std::endl;
    } else {
      TString stype = ((TObjString*)obj->At(2))->GetString();
      TString smod = ((TObjString*)obj->At(0))->GetString();
      TString sid = ((TObjString*)obj->At(8))->GetString();

      id = sid.Atoi();

      TObjArray* oarr;

      z = mymat.GetTranslation()[2];

      if (smod.Contains("frontST")) {
        oarr = smod.Tokenize("frontST");
        mod = 90 + ((TObjString*)oarr->At(0))->GetString().Atoi();
      } else if (smod.Contains("sttmod")) {
        oarr = smod.Tokenize("sttmod");
        mod = ((TObjString*)oarr->At(0))->GetString().Atoi();
      } else
        std::cout << "Error evaluating module id for: " << name.Data()
                  << std::endl;

      if (stype.Contains("ver")) {
        x = mymat.GetTranslation()[0];
        type = 1;
      } else if (stype.Contains("hor")) {
        x = mymat.GetTranslation()[1];
        type = 2;
      } else
        std::cout << "Error evaluating type for: " << name.Data() << std::endl;

      delete oarr;
    }
    delete obj;

    std::cout << name.Data() << " -> mod: " << mod << " type: " << type
              << " id: " << id << " x(y): " << x << " z: " << z << std::endl;
  }

  for (int i = 0; i < nod->GetNdaughters(); i++) {
    showElemPos(nod->GetDaughter(i), mymat, mypath);
  }
}

void showElemPos()
{
  TFile f("../files/reco/numu_internal_10k.0.reco.root");
  TGeoManager* gGeoManager = (TGeoManager*)f.Get("EDepSimGeometry");

  TGeoHMatrix mat = *gGeoIdentity;

  std::map<int, std::map<double, int> > stX;
  std::map<int, std::map<int, TVector2> > stPos;
  std::map<int, double> slabPos;

  getSTPlaneinfo(gGeoManager->GetTopVolume()->GetNode(0), mat, stX, stPos,
                 slabPos);

  slabPos[85] = slabPos[84];
  slabPos[86] = slabPos[84];
  slabPos[87] = slabPos[84];
  slabPos[88] = slabPos[84];
  slabPos[89] = slabPos[84];
  slabPos[90] = 22000.;
  slabPos[91] = 22000.;
  /*
  for(std::map<int, std::map<int, TVector2> >::iterator it = stPos.begin(); it
  != stPos.end(); ++it)
  {
    for(std::map<int, TVector2>::iterator iter = it->second.begin(); iter !=
  it->second.end(); ++iter)
    {
      std::cout << it->first << " " << iter->first << " " << iter->second.X() <<
  " " << iter->second.Y() << std::endl;
    }
  }*/

  for (std::map<int, double>::iterator it = slabPos.begin();
       it != slabPos.end(); ++it) {
    std::cout << it->first << " " << it->second << std::endl;
  }

  // showElemPos(gGeoManager->GetTopVolume()->GetNode(0), mat, path);
}