import math
import numpy as np
import Line2D as Line2D

class Circle:
    def __init__(self, arg_center_x, arg_center_y, arg_R):
        self.center_x = arg_center_x
        self.center_y = arg_center_y
        self.R = arg_R
    
    def x_l(self, angle):
        return self.center_x + self.R * math.cos(angle)
    
    def y_l(self, angle):
        return self.center_y + self.R * math.sin(angle)
    
    def dx_derivative(self, angle):
        return -1. * self.R * math.sin(angle)
    
    def dy_derivative(self, angle):
        return self.R * math.cos(angle)
    
    def get_point_at(self, angle):
        return (self.x_l(angle), self.y_l(angle))
    
    def get_points(self, n_points=100):
        angles = np.linspace(0, 2*np.pi, n_points)
        x = [self.x_l(i) for i in angles]
        y = [self.y_l(i) for i in angles]
        return (x,y)
    
    def get_derivative_at(self, angle):
        return (self.dx_derivative(angle), self.dy_derivative(angle))
    
    def get_derivative_from_point(self, arg_x, arg_y):
        angle = self.get_angle_from_point(arg_x, arg_y)
        return self.get_derivative_at(angle)
    
    def get_angle_from_point(self, arg_x, arg_y):
        if arg_y - self.center_y > 0:
            return math.atan2(arg_y - self.center_y, arg_x - self.center_x)
        else:
            return math.atan2(arg_y - self.center_y, arg_x - self.center_x) + 2 * math.pi
    
    def get_upper_semi_circle_points(self, n_points = 5000):
        angles = np.linspace(0, np.pi, n_points)
        x = [self.x_l(i) for i in angles]
        y = [self.y_l(i) for i in angles]
        return (x,y)

    def get_lower_semi_circle_points(self, n_points = 5000):
        angles = np.linspace(np.pi, 2 * np.pi, n_points)
        x = [self.x_l(i) for i in angles]
        y = [self.y_l(i) for i in angles]
        return (x,y)
    
    def distance_to_point(self, point):
        return abs(((self.center_x - point[0]) ** 2 + (self.center_y - point[1]) ** 2) ** 0.5 - self.R)
    
    def center(self):
        return (self.center_x, self.center_y)
    
    def print_circle_info(self):
        print("circle center (x,y): (", self.center_x, ", ", self.center_y, "), R: ", self.R)
    
    def center_x(self):
        return self.center_x
    
    def center_y(self):
        return self.center_y
    
    def R(self):
        return self.R
