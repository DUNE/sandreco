class dg_wire:
    def __init__(self, det, did, x, y, z, t0, de, adc, hor, wire_length, hindex):
        self.det = det
        self.did = did
        self.x = x
        self.y = y
        self.z = z
        self.t0 = t0
        self.de = de
        self.adc = adc
        self.tdc = 1e9
        self.hor = hor
        self.wire_length = wire_length
        self.hindex = hindex
        
        # True quantities
        self.t_hit = 1e9
        self.signal_time = 1e9
        self.drift_time = 1e9
        
        # Measured quantities
        self.t_hit_measured = 1e9  # via global trigger
        self.signal_time_measured = 1e9  # exploit different wire orientation
        self.drift_time_measured = 1e9  # tdc - signal_time_measured - t_hit_measured