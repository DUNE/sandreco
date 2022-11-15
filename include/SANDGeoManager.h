#include "SANDECALCellInfo.h"
#include "SANDSTTTubeInfo.h"

#include <TGeoManager.h>
#include <TVector3.h>
#include <TPRegexp.h>

#include <map>

#ifndef SANDGEOMANAGER_H
#define SANDGEOMANAGER_H

namespace sand_geometry {
    namespace stt {
        const char* const path_internal_volume =
            "volWorld_PV/rockBox_lv_PV_0/volDetEnclosure_PV_0/volSAND_PV_0/"
            "MagIntVol_volume_PV_0/sand_inner_volume_PV_0";
        const char* const name_internal_volume = "sand_inner_volume_PV";

        const char* const stt_single_tube_regex_string =
            "(horizontalST_(Ar|Xe)|STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_vv_ST)_PV_([0-9]+"
            ")(/|)";
        // "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_ST_stGas_(Xe|Ar)19_vol_PV_(["
        // "0-9]+)";
        const char* const stt_two_tubes_regex_string =
            "STT_([0-9]+)_(Trk|C3H6|C)Mod_(ST_|)(hh|vv)_2straw_PV_([0-9]+)(/|)";
        const char* const stt_plane_regex_string =
            "STT_([0-9]+)_(Trk|C3H6|C)Mod(_ST|)_(hh|vv)_PV_([0-9]+)(/|)";
        // "_(C3H6|C|Tr)Mod_([0-9]+)_(ST_|)(hor|ver|hor2)_vol_PV_0";
        const char* const stt_module_regex_string =
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
        
        const int number_of_layers = 5;
        const int number_of_cells_per_barrel_layer = 12;
        const int number_of_barrel_modules = 24;
        const int number_of_cells_per_endcap_layer = 90;

        // thickness of the layers in mm
        const double layer_thickness[number_of_layers] = {44., 44., 44., 44., 54.};

        //endcap module id
        const int endcap_module_ids[2] = {30, 40};

        namespace fluka {
            // ecal dimension for fluka
            const double barrel_module_xmin = 262.55;
            const double barrel_module_xmax = 292.85;
            const double barrel_module_thickness = 115.0;
            const double barrel_module_length = 4300;

            const double endcap_rmax = 2000.0;  // ad essere precisi nella realtà è 1980
            const double endcap_thickness = 115.0;
        }
    }
}

class SANDGeoManager : public TObject {
    public:
        enum class SANDGeoType {kFromEdepSim, kFromFluka};

    private:
        TGeoManager* geo_;
        SANDGeoType geo_type_;
        std::map<int, SANDECALCellInfo> cellmap_;
        std::map<int, SANDSTTTubeInfo> sttmap_;
        TPRegexp stt_single_tube_regex_;
        TPRegexp stt_two_tubes_regex_;
        TPRegexp stt_plane_regex_;
        TPRegexp stt_module_regex_;
        std::map<int, std::map<double, int> > stt_tube_tranverse_position_map_;

        // ECAL
        std::vector<double> get_levels_z(double half_module_height);
        int encode_ecal_barrel_cell_local_id(int layer, int cell);
        int encode_ecal_endcap_cell_local_id(int layer, int cell);
        std::pair<int, int> decode_ecal_barrel_cell_local_id(int id);
        std::pair<int, int> decode_ecal_endcap_cell_local_id(int id);
        std::map<int, TVector3> get_ecal_barrel_cell_center_local_position(const std::vector<double>& zlevels, double m, double q);
        std::map<int, TVector3> get_ecal_endcap_cell_center_local_position(const std::vector<double>& zlevels, double rmin, double rmax);
        void set_ecal_info();

        // STT
        int encode_stt_tube_id(int stt_plane_global_id, int stt_tube_local_id);
        void decode_stt_tube_id(int stt_tube_global_id, int& stt_plane_global_id, int& stt_tube_local_id);
        int encode_stt_plane_id(int stt_module_id, int stt_plane_local_id, int stt_plane_type);
        void decode_stt_plane_id(int stt_plane_global_id, int& stt_module_id, int& stt_plane_local_id, int& stt_plane_type);
        bool is_stt_tube(const TString& volume_name);
        bool is_stt_plane(const TString& volume_name);
        int get_stt_plane_id(const TString& volume_path);
        void set_stt_tube_info(const TGeoNode* const node, const TGeoHMatrix& matrix, int stt_plane_id);
        void set_stt_info(const TGeoHMatrix& matrix);
        void set_stt_info();
    
    public:
        SANDGeoManager(): cellmap_(), sttmap_(), 
            stt_single_tube_regex_(sand_geometry::stt::stt_single_tube_regex_string),
            stt_two_tubes_regex_(sand_geometry::stt::stt_two_tubes_regex_string),
            stt_plane_regex_(sand_geometry::stt::stt_plane_regex_string),
            stt_module_regex_(sand_geometry::stt::stt_module_regex_string),
            stt_tube_tranverse_position_map_() {};
        ~SANDGeoManager() {cellmap_.clear(); 
                           sttmap_.clear();
                           stt_tube_tranverse_position_map_.clear();};
        void init(TGeoManager* const geo, SANDGeoType geo_type);
        const SANDECALCellInfo& get_ecal_cell_info(int ecal_cell_id) {return cellmap_.at(ecal_cell_id);}
        const SANDSTTTubeInfo& get_stt_tube_info(int stt_tube_id) {return sttmap_.at(stt_tube_id);}
        int get_ecal_cell_id(double x, double y, double z);
        int get_stt_tube_id(double x, double y, double z);

    ClassDef(SANDGeoManager, 1);
};

#ifdef __MAKECINT__
#pragma link C++ class SANDGeoManager + ;
#endif

#endif