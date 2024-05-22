import uproot4 as upr
import numpy as np
import pandas as pd
import ROOT
import time
import glob

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
                                'pt_true',
                                'pt_reco',
                                ], library='pd').rename(columns={
                                'edep_file_input'                         : "edep_file",
                                'event_index'                             : "event_index",
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
                                'edep_file_input',
                                'event_index',
                                'reco_object/fit_infos_xz/fit_infos_xz.MinValue',
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.name',
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.id',
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.initial_guess',
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.value',
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.error',
                                ], library='pd').rename(columns={
                                'edep_file_input'                                                                                      : "edep_file",
                                'event_index'                                                                                          : "event_index",
                                'reco_object/fit_infos_xz/fit_infos_xz.MinValue'                                                       : 'MinValue',
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.name'          : "name",
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.id'            : "id",
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.initial_guess' : "initial_guess",
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.value'         : "value",
                                'reco_object/fit_infos_xz/fit_infos_xz.fitted_parameters/fit_infos_xz.fitted_parameters.error'         : "error"
                                })
        self.dataframe_fitted_values_zy = self.tree_upr.arrays([
                                'edep_file_input',
                                'event_index',
                                'reco_object/fit_infos_zy/fit_infos_zy.MinValue',
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.name',
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.id',
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.initial_guess',
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.value',
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.error',
                                ], library='pd').rename(columns={
                                'edep_file_input'                                                                                      : "edep_file",
                                'event_index'                                                                                          : "event_index",
                                'reco_object/fit_infos_zy/fit_infos_zy.MinValue'                                                       : "MinValue",
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.name'          : "name",
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.id'            : "id",
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.initial_guess' : "initial_guess",
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.value'         : "value",
                                'reco_object/fit_infos_zy/fit_infos_zy.fitted_parameters/fit_infos_zy.fitted_parameters.error'         : "error"
                                })
    
    def get_true_helix(self, ev_index: int):
        row =  self.dataframe_helix[self.dataframe_helix.event_index==ev_index]
        return Helix(row.R_true.values[0], row.dip_true.values[0], row.Phi0_true.values[0], 1, (row.x0_true.values[0], row.y0_true.values[0], row.z0_true.values[0]))
    
    def get_reco_helix(self, ev_index: int):
        row =  self.dataframe_helix[self.dataframe_helix.event_index==ev_index]
        return Helix(row.R_reco.values[0], row.dip_reco.values[0], row.Phi0_reco.values[0], 1, (row.x0_reco.values[0], row.y0_reco.values[0], row.z0_reco.values[0]))
    
    def get_fit_xz(self, ev_index: int):
        row = self.dataframe_fitted_values_xz[self.dataframe_fitted_values_xz.event_index==ev_index]
        return row
    
    def get_fit_zy(self, ev_index: int):
        row = self.dataframe_fitted_values_zy[self.dataframe_fitted_values_zy.event_index==ev_index]
        return row
    
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

    # def get_wire_info(self, ev_index):
    #     df_wire_ev = self.dataframe_wires[self.dataframe_wires.event_index==ev_index]
    #     df_wire_ev['drift_circle'] = [self._define_drift_circle_(wire) for _, wire in df_wire_ev.iterrows()]
    #     return df_wire_ev
    
    def get_wire_info(self, ev_index):
        self.dataframe_wires.loc[self.dataframe_wires.event_index==ev_index, 'drift_circle'] = \
        [self._define_drift_circle_(wire) for _, wire in self.dataframe_wires[self.dataframe_wires.event_index==ev_index].iterrows()]
        return self.dataframe_wires[self.dataframe_wires.event_index==ev_index]

class Reader_production:

    def __init__(self, arg_production = ""):
        self.production = arg_production
        self.files = self._sort_files_(arg_production)
        self.NaN_events = []
    
    def _sort_files_(self, production: str) -> list:
        if not isinstance(production, str):
            raise ValueError("Input has to be a string")
        
        if "*" not in production:
            raise ValueError("Input has to have * in the name")
        
        list_of_files = glob.glob(production)
        if len(list_of_files) == 0:
            raise ValueError("List of file production is empty")
        
        print("nof of files in the production {}".format(len(list_of_files)))
        
        return sorted(list_of_files, key=lambda x: int(x.split(".")[-3]) if x.split(".")[-3].isdigit() else float('inf'))
    
    def _try_dataframe_helix_(self, file: str) -> pd.DataFrame:
        try:
            reader = Reader(file, "tReco")
            return reader.dataframe_helix
        except:
            print("skipping file {}".format(file))
    
    def _try_dataframe_wires_(self, file: str) -> pd.DataFrame:
        try:
            reader = Reader(file, "tReco")
            return reader.dataframe_wires
        except:
            print("skipping file {}".format(file))
    
    def _try_dataframe_fit_xz_(self, file: str) -> pd.DataFrame:
        try:
            reader = Reader(file, "tReco")
            return reader.dataframe_fitted_values_xz
        except:
            print("skipping file {}".format(file))
    
    def _try_dataframe_fit_zy_(self, file: str) -> pd.DataFrame:
        try:
            reader = Reader(file, "tReco")
            return reader.dataframe_fitted_values_zy
        except:
            print("skipping file {}".format(file))
    
    def get_dataframe_production(self, up_2_index = 999) -> pd.DataFrame:
        dfs = []
        for file in self.files[:up_2_index]:
            print("processing file {}".format(file))
            dfs.append(self._try_dataframe_helix_(file))
        return pd.concat(dfs, ignore_index=True)
    
    def get_dataframe_wire(self, up_2_index = 999) -> pd.DataFrame:
        dfs = []
        for file in self.files[:up_2_index]:
            print("processing file {}".format(file))
            dfs.append(self._try_dataframe_wires_(file))
        return pd.concat(dfs, ignore_index=True)
    
    def get_dataframe_fit_xz(self, up_2_index = 999) -> pd.DataFrame:
        dfs = []
        for file in self.files[:up_2_index]:
            print("processing file {}".format(file))
            dfs.append(self._try_dataframe_fit_xz_(file))
        return pd.concat(dfs, ignore_index=True)
    
    def get_dataframe_fit_zy(self, up_2_index = 999) -> pd.DataFrame:
        dfs = []
        for file in self.files[:up_2_index]:
            print("processing file {}".format(file))
            dfs.append(self._try_dataframe_fit_zy_(file))
        return pd.concat(dfs, ignore_index=True)
    
    def get_nan_rows(self, df: pd.DataFrame) -> tuple:
        nan_rows = df[df.isnull().any(axis=1)]
        rows_to_exclude = [tuple(row) for row in nan_rows[['edep_file', 'event_index']].values]
        return pd.DataFrame(rows_to_exclude, columns=['edep_file', 'event_index']) 



class ROOT_tools:

    def __init__(self, arg_file = ""):
        self.file = arg_file
    
    def _check_TH1D_(self, h: ROOT.TH1D):
        if not isinstance(h, ROOT.TH1D) or h.GetEntries() == 0:
            raise ValueError("input hist has to be a non empty ROOT.TH1D")
        else:
            return True

    def FillTH1D(self, iterable: list, histogram_name: str, title: str, nbins: int , x_min: int, x_max: int) -> ROOT.TH1D:
        # fill TH1D from python-like vector
        # prevent memory leak
        if len(iterable)==0: raise ValueError("input vector is empty")
        histogram_name = histogram_name + " time:" + str(int(time.time()))
        histogram = ROOT.TH1D(histogram_name, title, nbins, x_min, x_max)
        for entry in iterable: histogram.Fill(entry)
        return histogram

    def FitTH1D_w_gauss(self, hist: ROOT.TH1D, gauss_range: tuple, fit_range: tuple) -> tuple:
        # fit TH1D with gauss function and return hist, mean and sigma of the fit
        self._check_TH1D_(hist)
        gaussian_func = ROOT.TF1("gaussian_func", "gaus", gauss_range[0], gauss_range[1])
        gaussian_func.SetRange(fit_range[0], fit_range[1])
        hist.Fit(gaussian_func, "R")
        mean = gaussian_func.GetParameter(1)
        sigma = gaussian_func.GetParameter(2)
        return (hist, mean, sigma)
    
    def FitTH1D_w_chi2(self, hist: ROOT.TH1D, chi2_range: tuple, fit_range: tuple) -> tuple:
        # fit TH1D with chi2 function and return hist, mean and sigma of the fit
        self._check_TH1D_(hist)
        chi2_func = ROOT.TF1("chi2_func", "[0]*ROOT::Math::chisquared_pdf(x, [1])", chi2_range[0], chi2_range[1])
        # Initial guess for the parameters: scale and degrees of freedom
        chi2_func.SetParameters(1, 1)
        # Set the range for the fit
        chi2_func.SetRange(fit_range[0], fit_range[1])
        # Fit the histogram with the chi-squared function
        hist.Fit(chi2_func, "R")
        # Extract the parameters: scale and degrees of freedom
        scale = chi2_func.GetParameter(0)
        dof = chi2_func.GetParameter(1)
        return (hist, scale, dof)
    
    def PlotTH1D(self, hist: ROOT.TH1D, canvas_name: str, canvas_dimensions = (500,800)) -> None:
        self._check_TH1D_(hist)
        canvas_name = "{} time: {}".format(canvas_name, str(int(time.time()))[-4:])
        ROOT.gStyle.SetOptFit(1011)
        canvas = ROOT.TCanvas(canvas_name, "Canvas", canvas_dimensions[0], canvas_dimensions[1])
        hist.Draw()
        canvas.Draw()
    