import uproot4 as upr
import numpy as np
import pandas as pd
import ROOT
import time
import glob

from dg_wire import dg_wire
from Helix import Helix
from Circle import Circle


class tDigitReader:
    def __init__(self, arg_file_name, arg_tree_name):
        self.file_name = arg_file_name
        self.tree_name = arg_tree_name
        self.file_upr = upr.open(arg_file_name)
        self.tree_upr = self.file_upr[arg_tree_name]
        self.dataframe_wires = self.tree_upr.arrays([ 
                                'edep_file_input',
                                'edep_event_index',
                                'fired_wires/fired_wires.did',
                                'fired_wires/fired_wires.x',
                                'fired_wires/fired_wires.y',
                                'fired_wires/fired_wires.z',
                                'fired_wires/fired_wires.de',
                                'fired_wires/fired_wires.adc',
                                'fired_wires/fired_wires.tdc',
                                'fired_wires/fired_wires.hor',
                                # 'fired_wires/fired_wires.hindex',
                                'fired_wires/fired_wires.wire_length',
                                'fired_wires/fired_wires.t_hit',
                                'fired_wires/fired_wires.drift_time',
                                'fired_wires/fired_wires.signal_time',
                                'fired_wires/fired_wires.t_hit_measured',
                                'fired_wires/fired_wires.signal_time_measured',
                                'fired_wires/fired_wires.drift_time_measured',
                                'fired_wires/fired_wires.missing_coordinate',
                                ],library='pd').rename(columns={
                                'fired_wires/fired_wires.did'                  : "did",
                                'fired_wires/fired_wires.x'                    : "x",
                                'fired_wires/fired_wires.y'                    : "y",
                                'fired_wires/fired_wires.z'                    : 'z',
                                'fired_wires/fired_wires.de'                   : 'de',
                                'fired_wires/fired_wires.adc'                  : 'adc',
                                'fired_wires/fired_wires.tdc'                  : 'tdc',
                                'fired_wires/fired_wires.t_hit'                : 't_hit',
                                # 'fired_wires/fired_wires.hindex'               : 'hindex',
                                'fired_wires/fired_wires.drift_time'           : 'drift_time',
                                'fired_wires/fired_wires.signal_time'          : 'signal_time',
                                'fired_wires/fired_wires.signal_time_measured' : 'signal_time_measured',
                                'fired_wires/fired_wires.t_hit_measured'       : 't_hit_measured',
                                'fired_wires/fired_wires.drift_time_measured'  : 'drift_time_measured',
                                'fired_wires/fired_wires.missing_coordinate'   : 'missing_coordinate',
                                'fired_wires/fired_wires.hor'                  : 'hor',
                                'fired_wires/fired_wires.wire_length'          : 'wire_length'
                                })