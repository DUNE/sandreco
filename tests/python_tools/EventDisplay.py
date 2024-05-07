import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from dg_wire import dg_wire
from Circle import Circle

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
    
    def __init__(self, arg_file_name,  # string
                       arg_event_idx,  # int
                       ax_0,  # axis for plotting view xz
                       ax_1,  # axis for plotting view yz
                       arg_true_helix_points = None,  # (x, y, z) with x, y, z vectors
                       arg_reco_helix_points = None,  # (x, y, z) with x, y, z vectors
                       arg_xz_points = None,  # (x, z)
                       arg_zy_points = None,  # (z, y)
                       arg_fired_wires_info = None,  # vector of dg_wire
                ):
        '''
        arg_file_name -> string: name of the edepsim file
        arg_event_idx -> int: event number
        ax_0 -> AxesSubplot: to be filled with xz plot
        ax_1 -> AxesSubplot: to be filled with zy plot
        arg_true_helix_points -> object (3, n): 3 vector of points x,y,z representig the true helix
        arg_reco_helix_points -> object (3, n): 3 vector of points x,y,z representig the true helix
        arg_xz_points -> object (2, n): 2 vector of points x,z representig the projection of the reco helix
        arg_zy_points -> object (2, n): 2 vector of points z,y representig the projection of the reco helix
        arg_fired_wires_info -> dataframe with columns dg_wire attributes: vector of dg_wire for the fired wires
        '''
        self.file_name = arg_file_name
        self.event_idx = arg_event_idx
        self.ax_0 = ax_0
        self.ax_1 = ax_1
        self.true_helix_points = self._handle_shape(arg_true_helix_points)
        self.reco_helix_points = self._handle_shape(arg_reco_helix_points)
        self.xz_points = self._handle_shape(arg_xz_points)
        self.zy_points = self._handle_shape(arg_zy_points)
        self.fired_wires_info = self._handle_wires_type(arg_fired_wires_info)
        
        self.hor_fired_wires_info = [wire for  wire in self.fired_wires_info if wire.hor == True]
        self.ver_fired_wires_info = [wire for  wire in self.fired_wires_info if wire.hor == False]
        self.hor_drift_circles = [Circle(w.z, w.y, w.drift_time_measured * 0.05) for w in self.hor_fired_wires_info]
        self.ver_drift_circles = [Circle(w.x, w.z, w.drift_time_measured * 0.05) for w in self.ver_fired_wires_info]
        print(type(self.ax_0))
    
    def _handle_wires_type(self, wires):
        if isinstance(wires, pd.DataFrame):
            # If the input is a Pandas DataFrame, create a list of dg_wire objects from the DataFrame rows
            return [dg_wire("", r.did, r.x, r.y, r.z, 0, 0, 0, r.hor, r.wire_length, 0) for _, r in wires.iterrows()]
        elif isinstance(wires, list) and all(isinstance(w, dg_wire) for w in wires):
            # If the input is a list of dg_wire objects, return the list directly
            return wires
        else:
            # If the input is neither a Pandas DataFrame nor a list of dg_wire objects, raise an exception
            raise ValueError(f"Cannot manage input wire type {type(wires)}")



    def _handle_shape(self, points):
        '''
        all inputs are expected to have shape (3, n) or (2, n)
        where 3 or 2 is the coordinate x, y, z and n is the
        number of entries
        '''
        if points is not None:
            # Check if the points have shape (3, n)
            shape_ = np.shape(points)
            nof_entries_ = len(points)
            if shape_[1] == nof_entries_:
                return points
            # If not, assume shape is (n, 3) and reshape
            elif shape_[0] == nof_entries_:
                return np.transpose(points)
        return None

    def _is_wire_in_zoom(self, wire, 
                      coordinate_to_check='x', 
                      zoom_x_range=default_x_range, 
                      zoom_y_range=default_y_range, 
                      zoom_z_range=default_z_range):
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
            raise ValueError("coordinate_to_check deve essere 'x', 'y' o 'z'")

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


    def plot_sand(self):
        print(type(self.ax_0))
        self.ax_0.vlines(-1650, self.sand_center[2] - self.sand_radius, self.sand_center[2] + self.sand_radius, color='red', label='sand')
        self.ax_0.vlines(+1650, self.sand_center[2] - self.sand_radius, self.sand_center[2] + self.sand_radius, color='red')
        self.ax_0.hlines(self.sand_center[2] - self.sand_radius, -1650, 1650, color='red')
        self.ax_0.hlines(self.sand_center[2] + self.sand_radius, -1650, 1650, color='red')

        self.ax_1.plot(self.z_sand, self.y_sand, linestyle='-', color='red')

    def plot_wires(self,
                    wires, 
                    zoom_x_range = default_x_range, 
                    zoom_y_range = default_y_range, 
                    zoom_z_range = default_z_range,
                    color = 'blue',
                    alpha = 0.5,
                    plot_opposit_view = False):
        wires = self._handle_wires_type(wires)
        ver_wires_x, ver_wires_y, ver_wires_z = [], [], []
        hor_wires_x, hor_wires_y, hor_wires_z = [], [], []
        for wire in wires:
            if wire.hor:  # True for horizontal wires, False for vertical wires
                if (self._is_wire_in_zoom(wire, coordinate_to_check="z", zoom_z_range=zoom_z_range) and 
                    self._is_wire_in_zoom(wire, coordinate_to_check="y", zoom_y_range=zoom_y_range)):
                    hor_wires_x.append(wire.x)
                    hor_wires_y.append(wire.y)
                    hor_wires_z.append(wire.z)
            else:
                if (self._is_wire_in_zoom(wire, coordinate_to_check="x", zoom_x_range=zoom_x_range) and 
                    self._is_wire_in_zoom(wire, coordinate_to_check="z", zoom_z_range=zoom_z_range)):
                    ver_wires_x.append(wire.x)
                    ver_wires_y.append(wire.y)
                    ver_wires_z.append(wire.z)
                
        # Plotting all wires
        self.ax_0.scatter(x=ver_wires_x, y=ver_wires_z, s=15, marker="x", label='Vertical Wires', color=color, alpha=alpha)
        self.ax_1.scatter(x=hor_wires_z, y=hor_wires_y, s=15, marker="x", label='Horizontal Wires', color=color, alpha=alpha)
        if(plot_opposit_view): 
            self._plot_wire_opposit_view(wires, zoom_x_range, zoom_y_range, zoom_z_range, color = color)

    def _plot_wire_opposit_view(self,
                                wires,
                                zoom_x_range=default_x_range, 
                                zoom_y_range=default_y_range, 
                                zoom_z_range=default_z_range,
                                color='blue'):
        hor_wires, ver_wires = [], []
        for wire in wires:
            if self._is_wire_in_zoom(wire, coordinate_to_check="z", zoom_z_range=zoom_z_range):
                if wire.hor:
                    hor_wires.append(wire.z)
                else:
                    ver_wires.append(wire.z)
        if hor_wires:
            self.ax_0.hlines(hor_wires, zoom_x_range[0], zoom_x_range[1], color=color, linewidth=0.005)
        if ver_wires:
            self.ax_1.vlines(ver_wires, zoom_y_range[0], zoom_y_range[1], color=color, linewidth=0.005)
    
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

    def plot_true_helix(self, color = 'blue'):
        self.ax_0.plot(self.true_helix_points[0], self.true_helix_points[2], linewidth = 1, label = 'true helix', color = color)
        self.ax_1.plot(self.true_helix_points[2], self.true_helix_points[1], linewidth = 1, label = 'true helix', color = color)
    
    def plot_reco_helix(self, color = 'green'):
        self.ax_0.plot(self.reco_helix_points[0], self.reco_helix_points[2], linewidth = 1, label = 'true helix', color = color)
        self.ax_1.plot(self.reco_helix_points[2], self.reco_helix_points[1], linewidth = 1, label = 'true helix', color = color)
    
    def plot_hits(self, hits, plane = 'zy'):
        for index, row in hits.iterrows():
            self.ax_1.plot([row.start_z, row.stop_z], [row.start_y, row.stop_y], 
                         marker = 'o', 
                         linestyle = '-', 
                         color = 'green', 
                         markersize = 1)
        for index, row in hits.iterrows():
            self.ax_0.plot([row.start_x, row.stop_x], [row.start_z, row.stop_z], 
                         marker = 'o', 
                         linestyle = '-', 
                         color = 'green', 
                         markersize = 1)

    def plot_drift_circles(self, plane = 'zy'):
        for circle in self.hor_drift_circles:
            z_measured, y_measured = circle.get_points(nof_points = 20)
            self.ax_1.plot(z_measured, y_measured, 
                           'b-', 
                           linestyle = 'dashed')
    
        for circle in self.ver_drift_circles:
            z_measured, y_measured = circle.get_points(nof_points = 20)
            self.ax_0.plot(z_measured, y_measured, 
                           'b-',
                           linestyle = 'dashed')
        
