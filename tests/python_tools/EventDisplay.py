import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes._axes import Axes
import re

from dg_wire import dg_wire
from Line2D import Line2D
from Circle import Circle
from Helix import Helix

class EventDisplay:
    sand_center = [0, -2384.73, 23910]
    sand_radius = 2000
    default_x_range = (-1700, 1700)
    default_y_range = (-4500, -300)
    default_z_range = (21500, 26000)
    zoom_x_range = default_x_range
    zoom_y_range = default_y_range
    zoom_z_range = default_z_range
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
            required_columns = ['edep_file_input', 'edep_event_index', 'did', 'x', 'y', 'z', 'de', 'adc', 'tdc',
                                'hor', 'wire_length', 't_hit', 'drift_time', 'signal_time',
                                't_hit_measured', 'signal_time_measured', 'drift_time_measured', 't_hit', 
                                # 'drift_circle'
                                ]

            # wires_columns = [
            #     (column.split('.')[-1] if '.' in column else column) 
            #     for column in wires.columns
            # ]


            for column in required_columns:
                if column not in wires.columns:
                    raise ValueError(f"wires dataframe has no a required column : {column}")
            self.hor_fired_wires_info = wires[(wires.hor == True)&(wires.edep_event_index == self.event_idx)]
            self.ver_fired_wires_info = wires[(wires.hor == False)&(wires.edep_event_index == self.event_idx)]
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
    
    @classmethod
    def set_zoom(cls, zoom_x_range = None, zoom_y_range = None, zoom_z_range = None):
        if zoom_x_range is not None:
            cls.zoom_x_range = zoom_x_range
        if zoom_y_range is not None:
            cls.zoom_y_range = zoom_y_range
        if zoom_z_range is not None:
            cls.zoom_z_range = zoom_z_range
    
    @classmethod
    def reset_zoom(cls):
        cls.zoom_x_range = cls.default_x_range
        cls.zoom_y_range = cls.default_y_range
        cls.zoom_z_range = cls.default_z_range
    
    def get_wires_in_zoom(self,
                          df_wires,
                          coordinate_1 = "z",
                          coordinate_2 = "x",
                          zoom_x_range = default_x_range,
                          zoom_y_range = default_y_range,
                          zoom_z_range = default_z_range):
        df_wires_ = df_wires.loc[df_wires.apply(self._is_wire_in_zoom_, 
                                                axis = 1, 
                                                args=(coordinate_1, zoom_x_range, zoom_y_range, zoom_z_range,))]
        if coordinate_1 == coordinate_2:
            return df_wires_
        else:
            df_wires__ = df_wires_.loc[df_wires_.apply(self._is_wire_in_zoom_, 
                                                       axis = 1, 
                                                       args=(coordinate_2, zoom_x_range, zoom_y_range, zoom_z_range,))]
            return df_wires__

    def plot_wires(self, 
                   color = 'blue',
                   alpha = 0.5,
                   label = "wires",
                   marker_size = 40,
                   marker_line_width = 2,
                   plot_drift_circles = False,
                   plot_opposit_view = False,
                   add_assumed_signal_time_axis = False):

        hor_wires_in_zoom = self.get_wires_in_zoom(self.hor_fired_wires_info, "z", "y", 
                                                   self.zoom_x_range, 
                                                   self.zoom_y_range,
                                                   self.zoom_z_range)
        ver_wires_in_zoom = self.get_wires_in_zoom(self.ver_fired_wires_info, "x", "z", 
                                                   self.zoom_x_range, 
                                                   self.zoom_y_range,    
                                                   self.zoom_z_range)

        self.ax_0.scatter(x=ver_wires_in_zoom['x'], y=ver_wires_in_zoom['z'], 
                          s=marker_size, marker="x", label=label, color=color, alpha=alpha, linewidth=marker_line_width)
        self.ax_1.scatter(x=hor_wires_in_zoom['z'], y=hor_wires_in_zoom['y'], 
                          s=marker_size, marker="x", label=label, color=color, alpha=alpha, linewidth=marker_line_width)
        
        if plot_drift_circles:
            self.plot_drift_circles(hor_wires_in_zoom, ver_wires_in_zoom, color=color, alpha=alpha)

        if plot_opposit_view:
            hor_wires_in_zoom_z = self.get_wires_in_zoom(self.hor_fired_wires_info, "z", "z",
                                                         self.default_x_range, 
                                                         self.default_y_range,
                                                         self.zoom_z_range)
            ver_wires_in_zoom_z = self.get_wires_in_zoom(self.ver_fired_wires_info, "z", "z",
                                                         self.default_x_range, 
                                                         self.default_y_range,
                                                         self.zoom_z_range)
            self.ax_0.hlines(hor_wires_in_zoom_z['z'], self.zoom_x_range[0], self.zoom_x_range[1], color="blue", linewidth=1)
            self.ax_1.vlines(ver_wires_in_zoom_z['z'], self.zoom_y_range[0], self.zoom_y_range[1], color="blue", linewidth=1)
        
        if add_assumed_signal_time_axis:    
            # add signal time axes
            ax_0_second_y = self.ax_0.twinx()
            ax_0_second_y.set_ylim(self.zoom_z_range)
            ax_0_second_y.set_yticks(hor_wires_in_zoom_z['z'])
            assumed_signal_time_hor = [f'{value:.2f} ns' for value in hor_wires_in_zoom_z['signal_time_measured']]
            ax_0_second_y.set_yticklabels(assumed_signal_time_hor)
            # ax_0_second_y.set_ylabel('Assumed signal time for hor wires')
            ax_0_second_y.tick_params(axis='y', labelcolor='b', labelrotation=-30, labelsize = 20)

            ax_1_second_x = self.ax_1.twiny()
            ax_1_second_x.set_xlim(self.zoom_z_range)
            ax_1_second_x.set_xticks(ver_wires_in_zoom_z['z'])
            assumed_signal_time_ver = [f'{value:.2f} ns' for value in ver_wires_in_zoom_z['signal_time_measured']]
            ax_1_second_x.set_xticklabels(assumed_signal_time_ver)
            # ax_0_second_y.set_ylabel('Assumed signal time for ver wires')
            ax_1_second_x.tick_params(axis='x', labelcolor='b', labelrotation=-30, labelsize = 20)

    def plot_sand(self):
        self.ax_0.vlines(-1650, self.sand_center[2] - self.sand_radius, self.sand_center[2] + self.sand_radius, color='red', label='SAND')
        self.ax_0.vlines(+1650, self.sand_center[2] - self.sand_radius, self.sand_center[2] + self.sand_radius, color='red')
        self.ax_0.hlines(self.sand_center[2] - self.sand_radius, -1650, 1650, color='red')
        self.ax_0.hlines(self.sand_center[2] + self.sand_radius, -1650, 1650, color='red')
        self.ax_1.plot(self.z_sand, self.y_sand, linestyle='-', color='red', label = 'SAND')

    def plot_wires_from_csv_table(self, df, color = 'blue', alpha = 0.5, size = 10):
        if not isinstance(df, pd.DataFrame):
            raise ValueError("Input must be a Pandas DataFrame.")
        required_columns = {"x", "y", "z", "orientation"}
        # hecks if the input is a Pandas DataFrame and contains columns with names "x", "y", "z", and "orientation".
        if not required_columns.issubset(df.columns):
            missing_columns = required_columns - set(df.columns)
            raise ValueError(f"DataFrame must contain columns: {missing_columns}")
        
        self.ax_0.scatter(x=df[df.orientation==1]['x'], y=df[df.orientation==1]['z'], s=size, label='Wires', color=color, alpha=alpha)
        self.ax_1.scatter(x=df[df.orientation==0]['z'], y=df[df.orientation==0]['y'], s=size, label='Wires', color=color, alpha=alpha)

    def plot_helix_points(self, helix_points, label, color, linewidth):
        helix_points = self._handle_shape_(helix_points)
        self.ax_0.plot(helix_points[0], helix_points[2], linewidth = linewidth, label = label, color = color)
        self.ax_1.plot(helix_points[2], helix_points[1], linewidth = linewidth, label = label, color = color)
    
    def plot_helix(self, helix : Helix, label = 'helix', color = 'blue', linewidth = 3):
        self.plot_helix_points(helix.get_helix_points(), label, color, linewidth)

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
                               color = color, alpha = alpha, label = label, linestyle = 'dashed')
            else:    
                self.ax_1.plot(circle_points[0], circle_points[1], 
                               color = color, alpha = alpha, linestyle = 'dashed')
        
        for i, wire in ver_wires_in_zoom.iterrows():
            circle_points = self._handle_shape_(wire['drift_circle'].get_points())
            if(i==0):
                self.ax_0.plot(circle_points[0], circle_points[1], 
                               color = color, alpha = alpha, label = label, linestyle = 'dashed')
            else:
                self.ax_0.plot(circle_points[0], circle_points[1], linestyle = 'dashed',
                               color = color, alpha = alpha)
    
    def plot(self, points, ax, color, label, linestyle, linewidth):
        self._handle_shape_(points)
        if ax == 0:
            self.ax_0.plot(points[0], points[1], color = color, label = label, linestyle = linestyle, linewidth = linewidth)
        elif ax == 1:
            self.ax_1.plot(points[0], points[1], color = color, label = label, linestyle = linestyle, linewidth = linewidth)
        else:
            raise ValueError(f"Cannot plot line on specified axis: {ax}"    )
    
    def plot_line(self, line : Line2D,  ax = 0, color = 'blue', label = '', linestyle = 'dashed', linewidth = 3):
        self.plot(points = line.get_points(np.linspace(self.zoom_x_range[0], self.zoom_x_range[1])),
                  ax = ax, color = color, label = label, linestyle = linestyle, linewidth = linewidth)
    
    def plot_circle(self, circle : Circle,  ax = 1, color = 'blue', label = '', linestyle = 'dashed', linewidth = 3):
        self.plot(points = circle.get_points(n_points=5000), 
                  ax = ax, color = color, label = label, linestyle = linestyle, linewidth = linewidth)
    
    def set_sand_range(self):
        self.ax_0.set_xlim(self.default_x_range)
        self.ax_0.set_ylim(self.default_z_range)
        self.ax_1.set_xlim(self.default_z_range)
        self.ax_1.set_ylim(self.default_y_range)
    
    def set_ax_ranges(self):
        self.ax_0.set_xlim(self.zoom_x_range)
        self.ax_0.set_ylim(self.zoom_z_range)
        self.ax_1.set_xlim(self.zoom_z_range)
        self.ax_1.set_ylim(self.zoom_y_range)
    
    def set_ax_labels(self, 
                      labelsize = 30,
                      numbersize = 30):
        self.ax_0.set_xlabel("X [mm]", fontsize=labelsize)
        self.ax_0.set_ylabel(r"Z [mm] $\nu$ beam axis", fontsize=labelsize)
        self.ax_0.tick_params(labelsize=numbersize)
        self.ax_1.set_xlabel(r"Z [mm] $\nu$ beam axis", fontsize=labelsize)
        self.ax_1.set_ylabel("Y [mm]", fontsize=labelsize)
        self.ax_1.tick_params(labelsize=numbersize)
    
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
    
    def set_title(self, titlesize = 30):
        self.fig.suptitle(
             f'file number : {self.extract_number_from_filename(self.file_name)}, '+
             f'Event number : {self.event_idx} '
            #  +
            #  f'pt true : {df[(df.edep_file==file_name)&(df.event_index==event_number)].pt_true[0]/1e3:.3f} GeV, '+
            #  f'pt reco : {df[(df.edep_file==file_name)&(df.event_index==event_number)].pt_reco[0]/1e3:.3f} GeV, '+
            #  f'dip true : {df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_true[0]:.3f}'+
            #  f'dip reco : {df[(df.edep_file==file_name)&(df.event_index==event_number)].dip_reco[0]:.3f}'
             ,fontsize=titlesize)

    def set_figure(self, titlesize = 30, labelsize = 30, numbersize = 30, legend_fontsize = 30):
        self.set_title(titlesize)
        self.set_ax_labels(labelsize, numbersize)
        self.set_ax_ranges()
        plt.legend(fontsize=legend_fontsize)