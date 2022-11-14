#include "SANDECALCellInfo.h"
#include "SANDSTTTubeInfo.h"

#include <TGeoManager.h>
#include <TVector3.h>

#include <map>

#ifndef SANDGEOMANAGER_H
#define SANDGEOMANAGER_H

namespace sand_geometry {
    namespace stt {
        const char* const path_internal_volume =
            "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
            "MagIntVol_volume_PV_0/sand_inner_volume_PV_0";
        const char* const name_internal_volume = "sand_inner_volume_PV";

        const char* const rST_string =
            "(horizontalST_(Ar|Xe)|STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_vv_ST)_PV_([0-9]+"
            ")(/|)";
        // "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_ST_stGas_(Xe|Ar)19_vol_PV_(["
        // "0-9]+)";
        const char* const r2ST_string =
            "STT_([0-9]+)_(Trk|C3H6|C)Mod_(ST_|)(hh|vv)_2straw_PV_([0-9]+)(/|)";
        const char* const rSTplane_string =
            "STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_(hh|vv)_PV_([0-9]+)(/|)";
        // "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_vol_PV_0";
        const char* const rSTmod_string =
            "STT_([0-9]+)_(Trk|C3H6|C)Mod_PV_([0-9]+)(/|)";
    }

    namespace ecal {
        const char* const path_barrel_template =
            "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
            "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_lv_PV_%d";
        const char* const path_endcapL_template =
            "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
            "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_0";
        const char* const path_endcapR_template =
            "volWorld_PV_1/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
            "MagIntVol_volume_PV_0/kloe_calo_volume_PV_0/ECAL_end_lv_PV_1";
        const char* const barrel_module_name = "ECAL_lv_PV";
        const char* const endcap_module_name = "ECAL_end_lv_PV";
        
        const int nLay = 5;
        const int nCel = 12;
        const int nMod = 24;

        const int nCel_ec = 90;

        // thickness of the layers in mm
        const double dzlay[nLay] = {44., 44., 44., 44., 54.};

        //endcap module id
        const int endcap_module_ids[2] = {30, 40};
    }
}

class SANDGeoManager : public TObject {
    private:
        enum class SANDGeoType {kFromEdepSim, kFromFluka};
        std::map<int, SANDECALCellInfo> cellmap_;
        std::map<int, SANDSTTTubeInfo> sttmap_;
        std::vector<double> get_levels_z(double half_module_height);
        int encode_ecal_barrel_cell_local_id(int layer, int cell);
        int encode_ecal_endcap_cell_local_id(int layer, int cell);
        std::pair<int, int> decode_ecal_barrel_cell_local_id(int id);
        std::pair<int, int> decode_ecal_endcap_cell_local_id(int id);
        std::map<int, TVector3> get_ecal_barrel_cell_center_local_position(const std::vector<double>& zlevels, double m, double q);
        std::map<int, TVector3> get_ecal_endcap_cell_center_local_position(const std::vector<double>& zlevels, double rmin, double rmax);
        void set_ecal_cell_center_position(TGeoManager* const geo, SANDGeoType geo_type);
    
    public:
        SANDGeoManager(): cellmap_(), sttmap_() {};
        ~SANDGeoManager() {cellmap_.clear(); sttmap_.clear();};
        void init(TGeoManager* const geo, SANDGeoType geo_type);


    ClassDef(SANDGeoManager, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDGeoManager + ;
#endif

#endif