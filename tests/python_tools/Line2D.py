import math

class Line2D:
    def __init__(self, arg_dx, arg_dy, arg_ax, arg_ay):
        self.dx = arg_dx
        self.dy = arg_dy
        self.ax = arg_ax
        self.ay = arg_ay
        self.m = arg_dy / arg_dx
        self.q = arg_ay - self.m * arg_ax
        self.direction = (arg_dx, arg_dy)
        self.p0 = (arg_ax, arg_ay)
    
    @classmethod
    def from_mq(cls, arg_m, arg_q, arg_ax=0.):
        dy_ = arg_m
        dx_ = 1.
        ax_ = arg_ax
        ay_ = arg_q + arg_m * arg_ax
        direction_ = (dx_, dy_)
        p0_ = (ax_, ay_)
        return cls(dx_, dy_, ax_, ay_)
    
    def distance_to_point(self, p1):
        diff = (p1[0] - self.p0[0], p1[1] - self.p0[1])
        projection_length = diff[0] * self.direction[0] + diff[1] * self.direction[1]
        projection = (projection_length * self.direction[0], projection_length * self.direction[1])
        diff_projection = (diff[0] - projection[0], diff[1] - projection[1])
        return math.sqrt(diff_projection[0] ** 2 + diff_projection[1] ** 2)
    
    def get_x_from_y(self, arg_y):
        return (arg_y - self.q) / self.m
    
    def dx(self):
        return self.dx
    
    def dy(self):
        return self.dy
    
    def ax(self):
        return self.ax
    
    def ay(self):
        return self.ay
    
    def m(self):
        return self.m
    
    def q(self):
        return self.q
    
    def direction(self):
        return self.direction
    
    def p0(self):
        return self.p0
    
    def get_points(self, x):
        return (x, self.m * x + self.q)
