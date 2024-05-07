import math
import particle

class Helix:
    def __init__(self, arg_R, arg_dip, arg_Phi0, arg_h, arg_x0, low_lim=-999., up_lim=-999.):
        self.R = arg_R
        self.dip = arg_dip
        self.Phi0 = arg_Phi0
        self.h = arg_h
        self.x0 = arg_x0
        self.low_lim = low_lim
        self.up_lim = up_lim
    
    @classmethod
    def from_trajectory(cls, momentum, x0, pdg):
        pt = math.sqrt(momentum.Z() * momentum.Z() + momentum.Y() * momentum.Y())
        pl = momentum.X()
        charge = particle.Particle.from_pdgid(13).charge
        
        dip = math.atan2(pl, pt)
        R = pt / (0.3 * 0.6)  # [m] = [GeV]/[T] or [mm] = [MeV]/[T]
        # h = 1 if charge < 0 else -1
        h = -1 if charge < 0 else 1
        Phi0 = math.atan2(momentum.Y(), momentum.Z()) + h * math.pi * 0.5
        
        return cls(R, dip, Phi0, h, x0)
    
    def x_h(self, s):
        return self.x0[0] + s * math.sin(self.dip)
    
    def y_h(self, s):
        return self.x0[1] + self.R * (math.sin(self.Phi0 - self.h * s * math.cos(self.dip) / self.R) - math.sin(self.Phi0))
    
    def z_h(self, s):
        return self.x0[2] + self.R * (math.cos(self.Phi0 - self.h * s * math.cos(self.dip) / self.R) - math.cos(self.Phi0))
    
    def dx_over_ds(self):
        return -math.sin(self.dip)
    
    def dy_over_ds(self, s):
        return math.cos(self.Phi0 + self.h * s * math.cos(self.dip) / self.R) * self.h * math.cos(self.dip)
    
    def dz_over_ds(self, s):
        return -math.sin(self.Phi0 + self.h * s * math.cos(self.dip) / self.R) * self.h * math.cos(self.dip)
    
    def get_point_at(self, s):
        return (self.x_h(s), self.y_h(s), self.z_h(s))
    
    def get_tangent_vector(self, s):
        return (self.dx_over_ds(), self.dy_over_ds(s), self.dz_over_ds(s))
    
    def get_phi_from_z(self, z):
        return math.acos((z - self.x0[2]) / self.R + math.cos(self.Phi0)) - self.Phi0
    
    def get_s_from_phi(self, Phi):
        return Phi * self.R / self.h / math.cos(self.dip)
    
    def get_helix_points(self, s_min=0, s_max=2e3, step=10.):
        if s_min >= s_max:
            print("s_min cannot be larger than s_max")
            raise ValueError
        points = []
        while s_min < s_max:
            points.append(self.get_point_at(s_min))
            s_min += step
        return points
    
    def set_helix_param(self, p):
        self.R = p[0]
        self.dip = p[1]
        self.Phi0 = p[2]
        self.h = p[3]
        self.x0 = (p[4], p[5], p[6])
    
    def set_R(self, arg):
        self.R = arg
    
    def set_dip(self, arg):
        self.dip = arg
    
    def set_Phi0(self, arg):
        self.Phi0 = arg
    
    def set_x0(self, arg):
        self.x0 = arg
    
    def set_low_lim(self, arg_low):
        self.low_lim = arg_low
    
    def set_up_lim(self, arg_up):
        self.up_lim = arg_up
    
    def set_helix_range_from_digit(self, digit):
        z_min = digit.z - 8
        z_max = digit.z + 8
        Phi_min = self.get_phi_from_z(z_max)
        Phi_max = self.get_phi_from_z(z_min)
        self.set_low_lim(self.get_s_from_phi(Phi_min))
        self.set_up_lim(self.get_s_from_phi(Phi_max))
    
    def print_helix_pars(self):
        print("R   ->", self.R)
        print("dip ->", self.dip)
        print("phi ->", self.Phi0)
        print("h   ->", self.h)
        print("x0, y0, z0 ->", self.x0[0], self.x0[1], self.x0[2])
        print("zc, yc     ->", self.x0[2] - self.R * math.cos(self.Phi0), self.x0[1] - self.R * math.sin(self.Phi0))
    
    def R(self):
        return self.R
    
    def dip(self):
        return self.dip
    
    def Phi0(self):
        return self.Phi0
    
    def h(self):
        return self.h
    
    def x0(self):
        return self.x0
    
    def low_lim(self):
        return self.low_lim
    
    def up_lim(self):
        return self.up_lim
    
    def center(self):
        return (self.x0[2] - self.R * math.cos(self.Phi0), self.x0[1] - self.R * math.sin(self.Phi0))