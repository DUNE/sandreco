import uproot4 as upr
import numpy as np
import pandas as pd
from dg_wire import dg_wire

from Helix import Helix

class Reader:
    def __init__(self, arg_file_name, arg_tree_name):
        self.file_name = arg_file_name
        self.tree_name = arg_tree_name
        self.file_upr = upr.open(arg_file_name)
        self.tree_upr = self.file_upr[arg_tree_name]
        self.dataframe_helix = self.tree_upr.arrays([
                                'edep_file_input',
                                'event_index',
                                'reco_object/pt_true',
                                'reco_object/pt_reco',
                                'reco_object/true_helix/true_helix.R_',
                                'reco_object/reco_helix/reco_helix.R_',
                                'reco_object/true_helix/true_helix.dip_',
                                'reco_object/reco_helix/reco_helix.dip_',
                                'reco_object/true_helix/true_helix.Phi0_',
                                'reco_object/reco_helix/reco_helix.Phi0_',
                                'reco_object/true_helix/true_helix.x0_',
                                'reco_object/reco_helix/reco_helix.x0_',
                                ], library='pd').rename(columns={
                                'edep_file_input'                         : "edep_file",
                                'event_index'                             : "event_index",
                                'reco_object/pt_true'                     : "pt_true",
                                'reco_object/pt_reco'                     : "pt_reco",
                                'reco_object/true_helix/true_helix.R_'    : 'R_true',
                                'reco_object/reco_helix/reco_helix.R_'    : 'R_reco',
                                'reco_object/true_helix/true_helix.dip_'  : 'dip_true',
                                'reco_object/reco_helix/reco_helix.dip_'  : 'dip_reco',
                                'reco_object/reco_helix/reco_helix.Phi0_' : 'Phi0_reco',
                                'reco_object/true_helix/true_helix.Phi0_' : 'Phi0_true',
                                'reco_object/reco_helix/reco_helix.x0_fX' : 'x0_reco',
                                'reco_object/true_helix/true_helix.x0_fX' : 'x0_true',
                                'reco_object/reco_helix/reco_helix.x0_fY' : 'y0_reco',
                                'reco_object/true_helix/true_helix.x0_fY' : 'y0_true',
                                'reco_object/reco_helix/reco_helix.x0_fZ' : 'z0_reco',
                                'reco_object/true_helix/true_helix.x0_fZ' : 'z0_true',
                                })
        self.dataframe_wires = self.tree_upr.arrays([ 
                                'edep_file_input',
                                'event_index',
                                'reco_object/fired_wires/fired_wires.did',
                                'reco_object/fired_wires/fired_wires.x',
                                'reco_object/fired_wires/fired_wires.y',
                                'reco_object/fired_wires/fired_wires.z',
                                'reco_object/fired_wires/fired_wires.de',
                                'reco_object/fired_wires/fired_wires.adc',
                                'reco_object/fired_wires/fired_wires.tdc',
                                'reco_object/fired_wires/fired_wires.hor',
                                'reco_object/fired_wires/fired_wires.wire_length',
                                'reco_object/fired_wires/fired_wires.t_hit',
                                'reco_object/fired_wires/fired_wires.drift_time',
                                'reco_object/fired_wires/fired_wires.signal_time',
                                'reco_object/fired_wires/fired_wires.t_hit_measured',
                                'reco_object/fired_wires/fired_wires.signal_time_measured',
                                'reco_object/fired_wires/fired_wires.drift_time_measured',
                                'reco_object/fired_wires/fired_wires.t_hit'],library='pd').rename(columns={
                                'edep_file_input'                                          : "edep_file",
                                'event_index'                                              : "event_index",
                                'reco_object/fired_wires/fired_wires.did'                  : "did",
                                'reco_object/fired_wires/fired_wires.x'                    : "x",
                                'reco_object/fired_wires/fired_wires.y'                    : "y",
                                'reco_object/fired_wires/fired_wires.z'                    : 'z',
                                'reco_object/fired_wires/fired_wires.de'                   : 'de',
                                'reco_object/fired_wires/fired_wires.adc'                  : 'adc',
                                'reco_object/fired_wires/fired_wires.tdc'                  : 'tdc',
                                'reco_object/fired_wires/fired_wires.t_hit'                : 't_hit',
                                'reco_object/fired_wires/fired_wires.drift_time'           : 'drift_time',
                                'reco_object/fired_wires/fired_wires.signal_time'          : 'signal_time',
                                'reco_object/fired_wires/fired_wires.signal_time_measured' : 'signal_time_measured',
                                'reco_object/fired_wires/fired_wires.t_hit_measured'       : 't_hit_measured',
                                'reco_object/fired_wires/fired_wires.drift_time_measured'  : 'drift_time_measured',
                                'reco_object/fired_wires/fired_wires.t_hit'                : 't_hit',
                                'reco_object/fired_wires/fired_wires.hor'                  : 'hor',
                                'reco_object/fired_wires/fired_wires.wire_length'          : 'wire_length'
                                })
    
    def get_true_helix(self, ev_index):
        row =  self.dataframe_helix[self.dataframe_helix.event_index==ev_index]
        return Helix(row.R_true[0], row.dip_true[0], row.Phi0_true[0], 1, (row.x0_true[0], row.y0_true[0], row.z0_true[0]))
    
    def get_reco_helix(self, ev_index):
        row =  self.dataframe_helix[self.dataframe_helix.event_index==ev_index]
        return Helix(row.R_reco[0], row.dip_reco[0], row.Phi0_reco[0], 1, (row.x0_reco[0], row.y0_reco[0], row.z0_reco[0]))

    def get_track_segmants_dataframe(self, plane = "XZ"):
        df = self.tree_upr.arrays([f'reco_object/track_segments_{plane}/track_segments_{plane}.dx_',
                                   f'reco_object/track_segments_{plane}/track_segments_{plane}.dy_',
                                   f'reco_object/track_segments_{plane}/track_segments_{plane}.ax_',
                                   f'reco_object/track_segments_{plane}/track_segments_{plane}.ay_',
                                   f'reco_object/track_segments_{plane}/track_segments_{plane}.m_',
                                   f'reco_object/track_segments_{plane}/track_segments_{plane}.q_',
                                   ], library = 'pd').rename(columns={
                                   f'reco_object/track_segments_{plane}/track_segments_{plane}.dx_':"dx",
                                   f'reco_object/track_segments_{plane}/track_segments_{plane}.dy_':"dy",
                                   f'reco_object/track_segments_{plane}/track_segments_{plane}.ax_':"ax",
                                   f'reco_object/track_segments_{plane}/track_segments_{plane}.ay_':"ay",
                                   f'reco_object/track_segments_{plane}/track_segments_{plane}.m_':"m",
                                   f'reco_object/track_segments_{plane}/track_segments_{plane}.q_':"q",
                                    })
        return df
    
    def get_wire_info(self, ev_index, return_type = 'pd_series'):
        if return_type == 'dg_wire':
            return [dg_wire("", r.did, r.x, r.y, r.z, 0, 0, 0, r.hor, r.wire_length, 0) 
                    for _, r in self.dataframe_wires[self.dataframe_wires.event_index==ev_index].iterrows()]
        else:
           return self.dataframe_wires[self.dataframe_wires.event_index==ev_index]

    
