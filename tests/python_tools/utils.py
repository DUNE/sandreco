import ROOT as r
import uproot4 as upr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import re
import uncertainties as unc
from uncertainties import ufloat
import time
import glob

from Helix import Helix
from Line2D import Line2D
from Circle import Circle

# FUNCTIONS TO GENERATE TRACK POINTS IN THE SPACE _________________________________________________

def GetFittingEfficiency(df, x_variable, bins, y_variable='TMinuitFinalStatus'):
    eff = []
    bins = np.array(bins)
    x_centers = (bins[1:] + bins[:-1])/2
    x_widths = (bins[1:] - bins[:-1])/2
    for x, w in zip(x_centers, x_widths):
        eff.append(df[(df[x_variable]>=x-w) & 
                      (df[x_variable]<=x+w)][y_variable].mean())
    eff = np.array(eff)
    y = {"x":x_centers, "efficiency":eff, "widths":x_widths}
    return y

def GetChiSquareFit(df, x_variable, bins, y_variable = 'MinValue'):
    min = []
    bins = np.array(bins)
    x_centers = (bins[1:] + bins[:-1])/2
    x_widths = (bins[1:] - bins[:-1])/2
    for x, w in zip(x_centers, x_widths):
        min.append(df[(df[x_variable]>=x-w) & 
                      (df[x_variable]<=x+w)][y_variable].mean())
    min = np.array(min)
    y = {"x":x_centers, y_variable:min, "widths":x_widths}
    return y

# FUCTIONS FOR HISTOGRAMS _________________________________________________________________________

def GaussFitTH1D(hist, fit_range):
    # fit TH1D with gauss function
    # __________________________________
    gaussian_func = r.TF1("gaussian_func", "gaus", -1, 1)
    gaussian_func.SetRange(fit_range[0], fit_range[1])
    hist.Fit(gaussian_func, "R")
    mean = gaussian_func.GetParameter(1)
    sigma = gaussian_func.GetParameter(2)
    return (mean, sigma)


def FillTH1D(vector, histogram_name, title, nbins, x_min, x_max):
    # fill TH1D from python-like vector
    # __________________________________
    # prevent potential memory leaks if the existing histogram isn't properly deleted.
    histogram_name = histogram_name+"_" + str(int(time.time()))
    histogram = r.TH1D(histogram_name, title, nbins, x_min, x_max)
    for entry in vector: histogram.Fill(entry)
    return histogram

def GetResolution(dataframe, variable, nbins = 100, hist_range = (-0.5,0.5), fit_range = (-0.1,0.1)):
    # Get resolution of variable using a hist
    # defined in hist_range and fit the hist
    # in fit_range
    # __________________________________
    variable_true = variable + "_true"
    variable_reco = variable + "_reco"
    if(variable=='dip'):
        residuals = (dataframe[variable_true] - dataframe[variable_reco])
    else:
        residuals = (1/dataframe[variable_true] - 1/dataframe[variable_reco]) / (1/dataframe[variable_reco])
        # residuals = (dataframe[variable_true] - dataframe[variable_reco]) / (dataframe[variable_reco])

    # fill TH1D
    hist = FillTH1D(residuals, variable+"_residuals", variable+"_residuals", nbins, hist_range[0], hist_range[1])

    # fit histogram
    # hist.Fit("gaus","S","",fit_range[0],fit_range[1])
    # hist.Draw()
    mean, sigma = GaussFitTH1D(hist, fit_range)

    print(variable+" Mean: {:.6f}, Sigma: {:.6f}".format(mean, sigma))
    # PlotTH1D(hist)

    return (mean, sigma, hist)

def PlotTH1D(hist, canvas_name = "canvas_"):
    canvas_name = canvas_name + str(int(time.time()))  # Generate a unique canvas name using a timestamp
    canvas = r.TCanvas(canvas_name, "Canvas", 800, 500)
    hist.Draw()
    canvas.Draw()

def GetExpRes_NOSMEAR(pt, dip_angle, n_digits, track_length = 1.5):
    # track_length = 1.5 m taken as average track length
    B = 0.6
    sigma = 200E-6
    pt = pt*1e-3
    return pt*np.sqrt(np.sqrt(720/(n_digits+4))*((sigma*pt)/(0.3*B*track_length**2*np.cos(dip_angle)**2)))

def GetDataframeInRange(df, variable, central_value, bin_width):
    # return sliced df whose variable "variable" is in range [variable_true - width, variable_true + width]
    return df[abs(df[variable] - central_value) <= bin_width/2]

def GetVariableValueInRange(df, variable_input, bins_input, function, variable_output):
    # get variable_output function value for each bins_input of the variable_input
    # return (input_variable_bin_center, output_variable_function_value)
    bin_center = (bins_input[1:]+bins_input[:-1])/2
    bin_width  = (bins_input[1:]-bins_input[:-1])
    var_in,var_out = [],[]
    for central_value, width in zip(bin_center, bin_width):
        df_inrange = GetDataframeInRange(df, variable_input, central_value, width)
        # mean_value_in_range = df_inrange[variable_output].apply(function)
        try:
            mean_value_in_range = function(df_inrange[variable_output])
            print(f"input_variable {variable_input}, range [{central_value - width/2},{central_value + width/2}], output_variable function value : {mean_value_in_range}")
            var_in.append(central_value)
            var_out.append(mean_value_in_range)
        except RuntimeWarning as rw:
            print("Caught a RuntimeWarning:", rw)
        except Exception as e:
            print("Caught an exception:", e)
            
    return (var_in, var_out)

# FUNTIONS TO READ RECONSTRUCTION TREE OUTPUT

def DF_Helix(file):
    
    file = upr.open(file)

    df = file['tReco'].arrays([ 
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
                        ],
                      library='pd').rename(columns={
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

    df['pt_res'] = (df['pt_true']-df['pt_reco'])/df['pt_true']

    df['dip_res'] = df['dip_true'] - df['dip_reco']

# TMinuit: 0 means fit success, 4 fit failure
    # df['TMinuitFinalStatus'].replace({0:1, 4:0}, inplace=True)

    event_chi = file['tReco'].arrays(['reco_object/fit_infos/fit_infos.MinValue'],library='np')['reco_object/fit_infos/fit_infos.MinValue']
    
    chi_v = []
    for i in range(len(event_chi)):
      chi_v.append(event_chi[i][-2:].mean())
    
    df['MinValue'] = chi_v

    return df


def GetTrackSegments(file, plane = "XZ"):

    df = upr.open(file)['tReco'].arrays([f'reco_object/track_segments_{plane}/track_segments_{plane}.dx_',
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

def GetFiredWires(file):

    file = upr.open(file)

    fired_wires = file['tReco'].arrays([ 
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
     'reco_object/fired_wires/fired_wires.hor'                  : 'hor'})
    
    fired_wires['true_drift_distance'] = fired_wires['drift_time'] * 0.05
    fired_wires['measured_drift_distance'] = fired_wires['drift_time_measured'] * 0.05

    return fired_wires

def Try_DF_TrueHelix(file):
    try: 
        return DF_Helix(file)
    except:
        print("skipping: ",file)

def Try_GetFiredWires(file):
    try: 
        return GetFiredWires(file)
    except:
        print("skipping: ",file)

def ConcatFiredWires(dfs):
    # concatenate fired wires dataframes
    # to get a concatenated df with 
    # coninuos indiex over the event number
    nof_events_previous = 0
    new_dfs = []
    for i, df in enumerate(dfs):
        try:
            df_reset = df.reset_index()
            # update event index making it continuos
            df_reset['entry'] += nof_events_previous
            new_dfs.append(df_reset)
            nof_events_previous = df_reset['entry'].max()
        except:
            print('skipping : ',i)
    return pd.concat(new_dfs).set_index(['entry','subentry'])

sand_center = [0,-2384.73,23910]
sand_radius = 2000
theta = np.linspace(0, 2*np.pi, 100000)
z_sand = sand_center[2] + sand_radius * np.cos(theta)
y_sand = sand_center[1] + sand_radius * np.sin(theta)

def extract_number_from_filename(filename):
    # Definisci il pattern regex per estrarre il numero dal nome del file
    pattern = r'\.(\d+)\.'

    # Cerca il pattern nel nome del file
    match = re.search(pattern, filename)

    if match:
        # Se trovi il match, restituisci il numero come intero
        return int(match.group(1))
    else:
        # Se non trovi nessun match, restituisci None o gestisci l'errore come preferisci
        return None

# FUNCTIONS FOR PLOTTING __________________________________________________________________________

def PlotImpactParameterZY(df, file_name, event_number,
              fired_wires,
              all_wires = None,
              hits = None,
              ZYcircleFit = None,
              ZYcircleFirst_guess = None,
              track_segments_YZ = None,
              y_range = [-4500, -300], 
              z_range = [21500,26000],
              plot_helix = True,
              plot_all_wires = False,
              plot_impact_par=False):
    
    h_points_true = GetHelixPoints2(df[(df.edep_file==file_name)&(df.event_index==event_number)].R_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].Phi0_true.values, 
                                    1, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].x0_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].y0_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].z0_true.values,
                                    np.linspace(-1e4,0,int(1e3)))
    
    h_points_reco = GetHelixPoints2(df[(df.edep_file==file_name)&(df.event_index==event_number)].R_reco.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_reco.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].Phi0_reco.values, 
                                    1, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].x0_reco.values,
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].y0_reco.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].z0_reco.values,
                                    np.linspace(-1e4,0,int(1e3)))
    
    fig,ax = plt.subplots(figsize=(20,9))
    
    # true and reco helix
    if plot_helix:
        ax.plot(h_points_true[2], h_points_true[1], linewidth=1, label='true helix', color='blue')
        ax.plot(h_points_reco[2], h_points_reco[1], linewidth=1, label='reco helix', color='green')
    
    if(ZYcircleFit is not None):
        ax.plot(ZYcircleFit[0], ZYcircleFit[1], linewidth=1, label='fit of TDCs of horizontal wires', color='magenta')
    
    if(ZYcircleFirst_guess is not None):
        ax.plot(ZYcircleFirst_guess[0], ZYcircleFirst_guess[1], linewidth=1, label='first guess', color='green')
    
    # fired digits
    wires_h = fired_wires[(fired_wires.edep_file==file_name)&(fired_wires.event_index==event_number)&(fired_wires.hor==True)]

    ax.scatter(x=wires_h.z, y=wires_h.y, s=15, marker="x",label='fired wires',color='orange')

    for i,w in wires_h.iterrows():
        # z_true,y_true = GetCirclePoints(w.z, w.y, w.true_drift_distance)
        # ax.plot(z_true, y_true, 'r-')
        z_measured,y_measured = GetCirclePoints(w.z, w.y, w.measured_drift_distance)
        if (i==0) :
            ax.plot(z_measured, y_measured, 'b-',label='drift circles', linestyle='dashed')
        else :
            ax.plot(z_measured, y_measured, 'b-', linestyle='dashed')
    
    # track segments
    if(track_segments_YZ is not None):
        for index, row in track_segments_YZ.iterrows():
            ax.plot([row['ax'], row['ay']], [row['ax'] + row['dx'], row['ay'] + row['dy']], 'g-')
    # for m,q in zip([0.331407, 0.380008],[-9119.57, -10255.4]):
    #     yy = [(m*z_i + q) for z_i in z_range]
    #     ax.plot(z_range, yy, 'g-')
    
    ax.scatter(x=all_wires[all_wires.orientation==0].z, 
               y=all_wires[all_wires.orientation==0].y,
               s=1, marker="x", label='wires', color='blue',alpha=0.5)
    ax.set_xlim(z_range)
    ax.set_ylim(y_range)

    ax.set_xlabel("Z [mm]")
    ax.set_ylabel("Y [mm]")

    plt.legend()
    # fig.suptitle(
    #              f'file number : {extract_number_from_filename(file_name)}, '+
    #              f'Event number : {event_number}, '+
    #              f'pt true : {df[(df.edep_file==file_name)&(df.event_index==event_number)].pt_true[0]/1e3:.3f} GeV, '+
    #              f'pt reco : {df[(df.edep_file==file_name)&(df.event_index==event_number)].pt_reco[0]/1e3:.3f} GeV, '+
    #              f'dip true : {df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_true[0]:.3f}'+
    #              f'dip reco : {df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_reco[0]:.3f}',fontsize=16)

def PlotImpactParameterXZ(df, file_name, event_number,
              fired_wires,
              all_wires = None,
              hits = None,
              XZFit = None,
              XZFirst_guess = None,
              track_segments_ZX = None,
              x_range = [-1700,1700], 
              z_range = [21500,26000],
              plot_helix = True,
              plot_all_wires = False,
              plot_impact_par=False):
    
    h_points_true = GetHelixPoints2(df[(df.edep_file==file_name)&(df.event_index==event_number)].R_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].Phi0_true.values, 
                                    1, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].x0_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].y0_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].z0_true.values,
                                    np.linspace(-1e4,0,int(1e3)))

    
    h_points_reco = GetHelixPoints2(df[(df.edep_file==file_name)&(df.event_index==event_number)].R_reco.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_reco.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].Phi0_reco.values, 
                                    1, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].x0_reco.values,
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].y0_reco.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].z0_reco.values,
                                    np.linspace(-1e4,0,int(1e3)))
    
    fig,ax = plt.subplots(figsize=(20,9))
    # true and reco helix
    if plot_helix:
        ax.plot(h_points_true[0], h_points_true[2], linewidth=1, label='true helix', color='blue')
        ax.plot(h_points_reco[0], h_points_reco[2], linewidth=1, label='reco helix', color='green')
    
    if(XZFit is not None):
        ax.plot(XZFit[0], XZFit[1], linewidth=1, label=' fit of TDCs of vertical wires', color='magenta')
    
    if(XZFirst_guess is not None):
        ax.plot(XZFirst_guess[0], XZFirst_guess[2], linewidth=1, label='first guess', color='green')
    
    # fired digits
    wires_v = fired_wires[(fired_wires.edep_file==file_name)&(fired_wires.event_index==event_number)&(fired_wires.hor==False)]

    ax.scatter(x=wires_v.x, y=wires_v.z, s=15, marker="x",label='fired wires',color='orange')

    for i,w in wires_v.iterrows():
        # x_true,z_true = GetCirclePoints(w.x, w.z, w.true_drift_distance)
        # ax.plot(x_true, z_true, 'r-')
        x_measured,z_measured = GetCirclePoints(w.x, w.z, w.measured_drift_distance)
        if(i==0) : 
            ax.plot(x_measured, z_measured, 'b-', label='drift circles', linestyle='dashed')
        else:
            ax.plot(x_measured, z_measured, 'b-', linestyle='dashed')
    
    # track segments
    if(track_segments_ZX is not None):
        for index, row in track_segments_ZX.iterrows():
            ax.plot([row['ax'], row['ay']], [row['ax'] + row['dx'], row['ay'] + row['dy']], 'g-')
    
    ax.scatter(x=all_wires[all_wires.orientation==1].x, 
               y=all_wires[all_wires.orientation==1].z,
               s=1, marker="x", label='wires', color='blue',alpha=0.5)
    ax.set_xlim(x_range)
    ax.set_ylim(z_range)

    ax.set_xlabel("X [mm]")
    ax.set_ylabel("Z [mm]")

    plt.legend()
    # fig.suptitle(
    #              f'file number : {extract_number_from_filename(file_name)}, '+
    #              f'Event number : {event_number}, '+
    #              f'pt true : {df[(df.edep_file==file_name)&(df.event_index==event_number)].pt_true[0]/1e3:.3f} GeV, '+
    #              f'pt reco : {df[(df.edep_file==file_name)&(df.event_index==event_number)].pt_reco[0]/1e3:.3f} GeV, '+
    #              f'dip true : {df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_true[0]:.3f}'+
    #              f'dip reco : {df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_reco[0]:.3f}',fontsize=16)

def PlotHelix(df, file_name, event_number,
              fired_wires, 
              pdf = None, 
              all_wires = None,
              hits = None,
              XZsinFit = None,
              ZYcircleFit = None,
              x_range = [-1700,1700], 
              y_range = [-4500, -300], 
              z_range = [21500,26000],
              plot_helix = True,
              plot_all_wires = False,
              plot_impact_par=False):

    h_points_true = GetHelixPoints2(df[(df.edep_file==file_name)&(df.event_index==event_number)].R_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].Phi0_true.values, 
                                    1, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].x0_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].y0_true.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].z0_true.values,
                                    np.linspace(-1e4,0,int(1e3)))

    h_points_reco = GetHelixPoints2(df[(df.edep_file==file_name)&(df.event_index==event_number)].R_reco.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_reco.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].Phi0_reco.values, 
                                    1, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].x0_reco.values,
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].y0_reco.values, 
                                    df[(df.edep_file==file_name)&(df.event_index==event_number)].z0_reco.values,
                                    np.linspace(-1e4,0,int(1e3)))
    
    
    fig,ax = plt.subplots(1,2,figsize=(18,9))

    # true and reco helix
    if plot_helix:
        ax[0].plot(h_points_true[2], h_points_true[1], linewidth=1, label='true helix', color='blue')
        ax[0].plot(h_points_reco[2], h_points_reco[1], linewidth=1, label='reco helix', color='green')

    if plot_helix:
        ax[1].plot(h_points_true[0], h_points_true[2], linewidth=1, label='true helix', color='blue')
        ax[1].plot(h_points_reco[0], h_points_reco[2], linewidth=1, label='reco helix', color='green')
    
    if(XZsinFit is not None):
        ax[1].plot(XZsinFit[0], XZsinFit[1], linewidth=1, label='sine fit of wire vertical coordinates', color='magenta')

    if(ZYcircleFit is not None):
        ax[0].plot(ZYcircleFit[0], ZYcircleFit[1], linewidth=1, label='circle fit of horizontal wire coordinates', color='magenta')
    
    # sand zy
    ax[0].plot(z_sand, y_sand, linestyle='-',color='red')
    
    # sand zx
    ax[1].vlines(-1650, sand_center[2] - sand_radius, sand_center[2] + sand_radius, color='red', label='sand')
    ax[1].vlines(+1650, sand_center[2] - sand_radius, sand_center[2] + sand_radius, color='red')
    ax[1].hlines(sand_center[2] - sand_radius,-1650, 1650, color='red')
    ax[1].hlines(sand_center[2] + sand_radius,-1650, 1650, color='red')
    
    # fired digits
    wires_h = fired_wires[(fired_wires.edep_file==file_name)&(fired_wires.event_index==event_number)&(fired_wires.hor==True)]
    wires_v = fired_wires[(fired_wires.edep_file==file_name)&(fired_wires.event_index==event_number)&(fired_wires.hor==False)]
    
    ax[0].scatter(x=wires_h.z, y=wires_h.y, s=15, marker="x",label='fired wires',color='orange')
    ax[1].scatter(x=wires_v.x, y=wires_v.z, s=15, marker="x",label='fired wires',color='orange')

    if(plot_impact_par):
        for i,w in wires_h.iterrows():
            z_true,y_true = GetCirclePoints(w.z, w.y, w.true_drift_distance)
            z_measured,y_measured = GetCirclePoints(w.z, w.y, w.measured_drift_distance)
            z_estimated,y_estimated = GetCirclePoints(w.z, w.y, w.impact_par_estimated)
            ax[0].plot(z_true, y_true, 'g-')
            # ax[0].plot(z_measured, y_measured, 'b-')
            # ax[0].plot(z_estimated, y_estimated, 'r-')

        for i,w in wires_v.iterrows():
            z_true,y_true = GetCirclePoints(w.x, w.z, w.true_drift_distance)
            z_measured,y_measured = GetCirclePoints(w.x, w.z, w.measured_drift_distance)
            z_estimated,y_estimated = GetCirclePoints(w.x, w.z, w.impact_par_estimated)
            ax[1].plot(z_true, y_true, 'g-')
            # ax[1].plot(z_measured, y_measured, 'b-')
            # ax[1].plot(z_estimated, y_estimated, 'r-')
    
    # all wires
    if (all_wires is not None):
        ax[0].scatter(x=all_wires[all_wires.orientation==0].z, 
                      y=all_wires[all_wires.orientation==0].y,
                      s=1, marker="x", label='wires', color='blue',alpha=0.5)
        # plot vertical wires in the zoom
        if plot_all_wires==True:
            for index, vertical_w in all_wires[all_wires.orientation==1].iterrows():
                if(vertical_w.z > z_range[0] and vertical_w.z < z_range[1]):
                    ax[0].vlines(vertical_w.z, vertical_w.y - vertical_w.length/2, vertical_w.y + vertical_w.length/2, 
                                color='blue',linewidth=0.005)

        ax[1].scatter(x=all_wires[all_wires.orientation==1].x, 
                      y=all_wires[all_wires.orientation==1].z,
                      s=1, marker="x", label='wires', color='blue',alpha=0.5)
        # plot horizontal wires in the zoom
        if plot_all_wires==True:    
            for index, horizontal_w in all_wires[all_wires.orientation==0].iterrows():
                if(horizontal_w.z > z_range[0] and horizontal_w.z < z_range[1]):
                    ax[1].hlines(horizontal_w.z, horizontal_w.x - horizontal_w.length/2, horizontal_w.x + horizontal_w.length/2, 
                                 color='blue',linewidth=0.005)  

    # hits
    if (hits is not None):
        for index, row in hits.iterrows():
            ax[0].plot([row.start_z, row.stop_z], [row.start_y, row.stop_y], 
                       marker='o', linestyle='-',color='green',markersize=1)

            ax[1].plot([row.start_x, row.stop_x], [row.start_z, row.stop_z], 
                       marker='o', linestyle='-',color='green',markersize=1)
            
            # # closest hit point 2 wire
            # ax[0].scatter(x=[row.hpoint_z],y=[row.hpoint_y], color='blue', marker='v', s=5)
            # ax[1].scatter(x=[row.hpoint_x],y=[row.hpoint_z], color='blue', marker='v', s=5)
    
    ax[0].set_xlim(z_range)
    ax[0].set_ylim(y_range)
    ax[1].set_xlim(x_range)
    ax[1].set_ylim(z_range)

    ax[0].set_xlabel("Z [mm]")
    ax[0].set_ylabel("Y [mm]")
    ax[1].set_xlabel("X [mm]")
    ax[1].set_ylabel("Z [mm]")

    plt.legend()
    fig.suptitle(
                 f'file number : {extract_number_from_filename(file_name)}, '+
                 f'Event number : {event_number}, '+
                 f'pt true : {df[(df.edep_file==file_name)&(df.event_index==event_number)].pt_true.values[0]/1e3:.3f} GeV, '+
                 f'pt reco : {df[(df.edep_file==file_name)&(df.event_index==event_number)].pt_reco.values[0]/1e3:.3f} GeV, '+
                 f'dip true : {df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_true.values[0]:.3f} '+
                 f'dip reco : {df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_reco.values[0]:.3f} ',fontsize=16)
    if(pdf): 
        pdf.savefig()
