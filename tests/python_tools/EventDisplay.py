import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes._axes import Axes
import re

from dg_wire import dg_wire
from Circle import Circle
from Helix import Helix

import numpy as np

class EventDisplay:
    sand_center = [0, -2384.73, 23910]
    sand_radius = 2000
    default_x_range = (-1700, 1700)
    default_y_range = (-4500, -300)
    default_z_range = (21500, 26000)
    theta = np.linspace(0, 2*np.pi, 100000)
    z_sand = sand_center[2] + sand_radius * np.cos(theta)
    y_sand = sand_center[1] + sand_radius * np.sin(theta)
    
    def __init__(self, arg_file_name : str,  # string
                       arg_event_idx : int,  # int
                       fig : Figure , # figure
                       ax_0 : Axes,  # axis for plotting view xz
                       ax_1 : Axes,  # axis for plotting view yz
                       arg_fired_wires_info = None,
                ):
        '''
        arg_file_name -> string: name of the edepsim file
        arg_event_idx -> int: event number
        ax_0 -> AxesSubplot: to be filled with xz plot
        ax_1 -> AxesSubplot: to be filled with zy plot
        arg_fired_wires_info -> dataframe with columns dg_wire attributes: vector of dg_wire for the fired wires
        '''
        self.file_name = arg_file_name
        self.event_idx = arg_event_idx
        self.fig = fig
        self.ax_0 = ax_0
        self.ax_1 = ax_1
        self.fired_wires_info = self._handle_wires_input_(arg_fired_wires_info)
        # self.hor_fired_wires_info = []
        # self.ver_fired_wires_info = []
    
    def _handle_wires_input_(self, wires):
        if isinstance(wires, pd.DataFrame):
            required_columns = ['edep_file', 'event_index', 'did', 'x', 'y', 'z', 'de', 'adc', 'tdc',
                                'hor', 'wire_length', 't_hit', 'drift_time', 'signal_time',
                                't_hit_measured', 'signal_time_measured', 'drift_time_measured', 't_hit', 'drift_circle']
            for column in wires.columns:
                if column not in required_columns:
                    raise ValueError(f"wires dataframe has no a required column : {column}")
            self.hor_fired_wires_info = wires[wires.hor == True]
            self.ver_fired_wires_info = wires[wires.hor == False]
        else:
            raise ValueError("wires info are expected to be pandas dataframe")

    def _handle_shape_(self, points):
        """
        Manipulates the shape of the points array according to the following criteria:
        - If the shape is (3, n) or (2, n), leaves it unchanged.
        - If the shape is (n, 3) or (n, 2), transpose it so it becomes (3, n) or (2, n).
        - Otherwise, raises an exception.

        Parameters:
        points (numpy.ndarray): Array of points.

        Returns:
        numpy.ndarray or None: Manipulated array of points, or None if input is None.

        Raises:
        ValueError: If the shape of the points array is not among (3, n), (2, n), (n, 3), or (n, 2).
        """
        if points is None:
            return None

        shape_ = np.shape(points)
        if len(shape_) == 2:
            if shape_[0] in [2, 3] and shape_[1] > 1:
                return points
            elif shape_[1] in [2, 3] and shape_[0] > 1:
                return np.transpose(points)

        raise ValueError("The shape of the points array is invalid. It must be (3, n), (2, n), (n, 3), or (n, 2).")


    def _is_wire_in_zoom_(self, wire, 
                      coordinate_to_check='x', 
                      zoom_x_range = default_x_range, 
                      zoom_y_range = default_y_range, 
                      zoom_z_range = default_z_range):
        """
        Check if the wire is within the specified zoom box along the specified axis.

        Args:
            wire (dg_wire): The wire to check.
            coordinate_to_check (str): The axis of the coordinates to check ('x', 'y', or 'z').
            zoom_x_range (tuple): The zoom range along the x-axis.
            zoom_y_range (tuple): The zoom range along the y-axis.
            zoom_z_range (tuple): The zoom range along the z-axis.

        Returns:
            bool: True if the wire is within the zoom box, False otherwise.
        """

        if coordinate_to_check not in ['x', 'y', 'z']:
            raise ValueError("coordinate_to_check has to be one of 'x', 'y' o 'z'")

        range_middle = {
            'x': (zoom_x_range[1] + zoom_x_range[0]) / 2.0,
            'y': (zoom_y_range[1] + zoom_y_range[0]) / 2.0,
            'z': (zoom_z_range[1] + zoom_z_range[0]) / 2.0
        }[coordinate_to_check]

        range_half_length = {
            'x': (zoom_x_range[1] - zoom_x_range[0]) / 2.0,
            'y': (zoom_y_range[1] - zoom_y_range[0]) / 2.0,
            'z': (zoom_z_range[1] - zoom_z_range[0]) / 2.0
        }[coordinate_to_check]

        wire_coordinate = {
            'x': wire.x,
            'y': wire.y,
            'z': wire.z
        }[coordinate_to_check]

        return (abs(wire_coordinate - range_middle) <= range_half_length)
    
    def get_wires_in_zoom(self,
                          df_wires,
                          coordinate_1 = "z",
                          coordinate_2 = None,
                          zoom_x_range = default_x_range,
                          zoom_y_range = default_y_range,
                          zoom_z_range = default_z_range):
        df_wires_ = df_wires.loc[df_wires.apply(self._is_wire_in_zoom_, 
                                    axis = 1, 
                                    args=(coordinate_1, zoom_x_range, zoom_y_range, zoom_z_range,))]
        if coordinate_2 == None:
            return df_wires_
        else:
            df_wires__ = df_wires_.loc[df_wires_.apply(self._is_wire_in_zoom_, 
                                    axis = 1, 
                                    args=(coordinate_2, zoom_x_range, zoom_y_range, zoom_z_range,))]
            return df_wires__

    def plot_wires(self, 
                   zoom_x_range = default_x_range,
                   zoom_y_range = default_y_range,
                   zoom_z_range = default_z_range,
                   color = 'blue',
                   alpha = 0.5,
                   plot_drift_circles = False,
                   plot_opposit_view = False):
        
        hor_wires_in_zoom = self.get_wires_in_zoom(self.hor_fired_wires_info, "z", "y", 
                                                   zoom_x_range, 
                                                   zoom_y_range, 
                                                   zoom_z_range)
        ver_wires_in_zoom = self.get_wires_in_zoom(self.ver_fired_wires_info, "x", "z", 
                                                   zoom_x_range, 
                                                   zoom_y_range,    
                                                   zoom_z_range)

        self.ax_0.scatter(x=ver_wires_in_zoom['x'], y=ver_wires_in_zoom['z'], 
                          s=15, marker="x", label='Vertical Wires', color=color, alpha=alpha)
        self.ax_1.scatter(x=hor_wires_in_zoom['z'], y=hor_wires_in_zoom['y'], 
                          s=15, marker="x", label='Horizontal Wires', color=color, alpha=alpha)
        
        if plot_drift_circles:
            self.plot_drift_circles(hor_wires_in_zoom, ver_wires_in_zoom, color=color, alpha=alpha)

        if plot_opposit_view:
            hor_wires_in_zoom_z = self.get_wires_in_zoom(self.hor_fired_wires_info, "z", 
                                                         self.default_x_range, 
                                                         self.default_y_range, 
                                                         zoom_z_range)
            ver_wires_in_zoom_z = self.get_wires_in_zoom(self.ver_fired_wires_info, "z", 
                                                         self.default_x_range, 
                                                         self.default_y_range, 
                                                         zoom_z_range)
            self.ax_0.hlines(hor_wires_in_zoom_z['z'], zoom_x_range[0], zoom_x_range[1], color=color, linewidth=0.005)
            self.ax_1.vlines(ver_wires_in_zoom_z['z'], zoom_y_range[0], zoom_y_range[1], color=color, linewidth=0.005)

    def plot_sand(self):
        self.ax_0.vlines(-1650, self.sand_center[2] - self.sand_radius, self.sand_center[2] + self.sand_radius, color='red', label='sand')
        self.ax_0.vlines(+1650, self.sand_center[2] - self.sand_radius, self.sand_center[2] + self.sand_radius, color='red')
        self.ax_0.hlines(self.sand_center[2] - self.sand_radius, -1650, 1650, color='red')
        self.ax_0.hlines(self.sand_center[2] + self.sand_radius, -1650, 1650, color='red')
        self.ax_1.plot(self.z_sand, self.y_sand, linestyle='-', color='red')
    
    def plot_wires_from_csv_table(self, df, color = 'blue', alpha = 0.5):
        if not isinstance(df, pd.DataFrame):
            raise ValueError("Input must be a Pandas DataFrame.")
        required_columns = {"x", "y", "z", "orientation"}
        # hecks if the input is a Pandas DataFrame and contains columns with names "x", "y", "z", and "orientation".
        if not required_columns.issubset(df.columns):
            missing_columns = required_columns - set(df.columns)
            raise ValueError(f"DataFrame must contain columns: {missing_columns}")
        
        self.ax_0.scatter(x=df[df.orientation==1]['x'], y=df[df.orientation==1]['z'], s=1, label='Horizontal Wires', color=color, alpha=alpha)
        self.ax_1.scatter(x=df[df.orientation==0]['z'], y=df[df.orientation==0]['y'], s=1, label='Vertical Wires', color=color, alpha=alpha)

    def plot_helix_points(self, helix_points, label = 'helix', color = 'blue'):
        helix_points = self._handle_shape_(helix_points)
        self.ax_0.plot(helix_points[0], helix_points[2], linewidth = 1, label = label, color = color)
        self.ax_1.plot(helix_points[2], helix_points[1], linewidth = 1, label = label, color = color)
    
    def plot_helix(self, helix : Helix, label = 'helix', color = 'blue'):
        self.plot_helix_points(helix.get_helix_points(), label, color)

    def plot_drift_circles(self,
                   hor_wires_in_zoom,
                   ver_wires_in_zoom,
                   color = 'blue',
                   label = 'drift circles',
                   alpha = 1):
        # first get wires in the zoom
        for i, wire in hor_wires_in_zoom.iterrows():
            circle_points = self._handle_shape_(wire['drift_circle'].get_points())
            if(i==0):
                self.ax_1.plot(circle_points[0], circle_points[1], 
                               color=color, alpha = alpha, label = label, linestyle = 'dashed')
            else:    
                self.ax_1.plot(circle_points[0], circle_points[1], 
                               color=color, alpha = alpha, linestyle = 'dashed')
        
        for i, wire in ver_wires_in_zoom.iterrows():
            circle_points = self._handle_shape_(wire['drift_circle'].get_points())
            if(i==0):
                self.ax_0.plot(circle_points[0], circle_points[1], 
                               color=color, alpha = alpha, label = label, linestyle = 'dashed')
            else:
                self.ax_0.plot(circle_points[0], circle_points[1], linestyle = 'dashed',
                               color=color, alpha = alpha)
    
    def plot(self, points, ax = 0, label = '', linestyle = 'dashed'):
        self._handle_shape_(points)
        if ax == 0:
            self.ax_0.plot(points[0], points[1], 'b-', label = label, linestyle = linestyle)
        elif ax == 1:
            self.ax_1.plot(points[0], points[1], 'b-', label = label, linestyle = linestyle)
        else:
            raise ValueError(f"Cannot plot line on specified axis: {ax}")
    
    def set_sand_range(self):
        self.ax_0.set_xlim(self.default_x_range)
        self.ax_0.set_ylim(self.default_z_range)
        self.ax_1.set_xlim(self.default_z_range)
        self.ax_1.set_ylim(self.default_y_range)
    
    def set_ax_ranges(self, zoom_x_range, zoom_y_range, zoom_z_range):
        self.ax_0.set_xlim(zoom_x_range)
        self.ax_0.set_ylim(zoom_z_range)
        self.ax_1.set_xlim(zoom_z_range)
        self.ax_1.set_ylim(zoom_y_range)
    
    def set_ax_labels(self):
        self.ax_0.set_xlabel("X [mm]")
        self.ax_0.set_ylabel(r"Z [mm] $\nu$ beam axis")
        self.ax_1.set_xlabel(r"Z [mm] $\nu$ beam axis")
        self.ax_1.set_ylabel("Y [mm]")
    
    def extract_number_from_filename(self, filename):
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
    
    def set_title(self):
        self.fig.suptitle(
             f'file number : {self.extract_number_from_filename(self.file_name)}, '+
             f'Event number : {self.event_idx}, '
            #  +
            #  f'pt true : {df[(df.edep_file==file_name)&(df.event_index==event_number)].pt_true[0]/1e3:.3f} GeV, '+
            #  f'pt reco : {df[(df.edep_file==file_name)&(df.event_index==event_number)].pt_reco[0]/1e3:.3f} GeV, '+
            #  f'dip true : {df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_true[0]:.3f}'+
            #  f'dip reco : {df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_reco[0]:.3f}'
             ,fontsize=16)

    def set_figure(self, 
                   zoom_x_range = default_x_range, 
                   zoom_y_range = default_y_range, 
                   zoom_z_range = default_z_range):
        self.set_title()
        self.set_ax_labels()
        self.set_ax_ranges(zoom_x_range, zoom_y_range, zoom_z_range)
        plt.legend()