import uproot4 as upr
import numpy as np
import pandas as pd

from dg_wire import dg_wire
from Helix import Helix
from Circle import Circle

class Reader:
    drift_velocity = 0.05
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
        self.dataframe_fitted_values_xz = self.tree_upr.arrays([
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.name',
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.id',
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.initial_guess',
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.value',
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.error',
                                ], library='pd').rename(columns={
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.name'          : "name",
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.id'            : "id",
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.initial_guess' : "initial_guess",
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.value'         : "value",
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.error'         : "error"
                                })
        self.dataframe_fitted_values_zy = self.tree_upr.arrays([
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.name',
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.id',
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.initial_guess',
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.value',
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.error',
                                ], library='pd').rename(columns={
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.name'          : "name",
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.id'            : "id",
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.initial_guess' : "initial_guess",
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.value'         : "value",
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.error'         : "error"
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
    
    def _define_drift_circle_(self, wire):
        if wire.hor == 1:
            return Circle(wire.z, wire.y, wire.drift_time_measured * self.drift_velocity)
        elif wire.hor == 0:
            return Circle(wire.x, wire.z, wire.drift_time_measured * self.drift_velocity)
        else:
            raise ValueError(f"cannot define a wire circle for {wire.hor}")

    def get_wire_info(self, ev_index):
        df_wire_ev = self.dataframe_wires[self.dataframe_wires.event_index==ev_index]
        df_wire_ev['drift_circle'] = [self._define_drift_circle_(wire) for _, wire in df_wire_ev.iterrows()]
        return df_wire_ev

    
